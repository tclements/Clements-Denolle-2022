using Distributed
addprocs()

@everywhere begin

    using Arrow 
    using DataFrames
    using Dates 
    using Glob 
    using SCEDCCorr
    using SeisDvv
    using SeisNoise
    using Serialization

    function robust_pmap(f::Function, args...; num_retries::Int=3)

        # This function returns true if we should retry and fails if we should error out.
        function retry_check(delay_state, err)
            # Below each condition to retry is listed along with an explanation about why
            # retrying should/might work.
            should_retry = (
                # Worker death is normally stocastic, if not then doesn't matter how many
                # retries as it will rapidly kill all workers
                err isa ProcessExitedException ||
                # If we are in the middle of fetching data and the process is killed we could
                # get an ArgumentError saying that the stream was closed or unusable.
                # So same as above.
                err isa ArgumentError && occursin("stream is closed or unusable", err.msg) ||
                # In general IO errors can be transient and related to network blips
                err isa Base.IOError
            )
            should_retry && info(LOGGER, "Retrying computation that failed due to $err")
            return should_retry
        end
    
        return pmap(f, CachingPool(workers()), args...;
            retry_check=retry_check,
            retry_delays=ExponentialBackOff(n=num_retries)
        )
    end

    function do_stretch(
        C::CorrData, 
        freqmin::Real,
        freqmax::Real, 
        minlag::Real, 
        maxlag::Real, 
        dvmin::Real, 
        dvmax::Real;
        ntrial::Int=101,
        ntrial2::Int=21,
        remove_max::Bool=false,
        interval::Period=Day(1),
    )
    
        remove_nan!(C)
    
        # preprocess
        if length(C.t) == 0
            return nothing
        end
        clean_up!(C,freqmin,freqmax)
        if remove_max
            abs_max!(C)
        end

        # smooth correlation 
        smooth!(C,interval=interval)
    
        # get reference
        lags = -C.maxlag:1/C.fs:C.maxlag
        times = Dates.unix2datetime.(C.t)
        minpad = max(5 / C.fs, dvmax * minlag)
        maxpad = max(5 / C.fs, dvmax * maxlag)
        indpos = findall((lags .>= minlag - minpad) .& (lags .<= maxlag + maxpad))
        indneg = findall((lags .>= -maxlag - maxpad) .& (lags .<= -minlag + minpad))
        window = findall((lags[indpos] .>= minlag) .& (lags[indpos] .<= maxlag))
        N = length(indpos)
        A = C.corr
        ref = SeisNoise.stack(C,allstack=true)
        if remove_max
            abs_max!(ref)
        end
        r = ref.corr
    
        dvvpos, dvvneg = zeros(size(A,2)), zeros(size(A,2))
        ccpos, ccneg = zeros(size(A,2)), zeros(size(A,2))
        cdppos, cdpneg = zeros(size(A,2)), zeros(size(A,2))
        errpos, errneg = zeros(size(A,2)), zeros(size(A,2))
        allCpos, allCneg = zeros(ntrial,size(A,2)), zeros(ntrial,size(A,2))
    
        # inital stretching 
        for ii in 1:size(A,2)
            dvvpos[ii],ccpos[ii],cdppos[ii], ϵ, errpos[ii], allCpos[:,ii] = stretching(
                r[indpos],
                A[indpos,ii],
                lags[indpos],
                window,
                freqmin,
                freqmax,
                dvmin=dvmin,
                dvmax=dvmax,
                ntrial=ntrial,
            )
            dvvneg[ii],ccneg[ii],cdpneg[ii], ϵ, errneg[ii], allCneg[:,ii] = stretching(
                r[indneg],
                A[indneg,ii],
                lags[indneg],
                window,
                freqmin,
                freqmax,
                dvmin=dvmin,
                dvmax=dvmax,
                ntrial=ntrial,
            )
        end
    
        # restretch reference using initial dv/v - use weights based on CC 
        ccpos[ccpos .< 0.0] .= 0.0
        ccneg[ccneg .< 0.0] .= 0.0 
    
        # load into DataFrame, then write to disk
        df = DataFrame(
            DATE = Date.(times), 
            DVVPOS = dvvpos, 
            CCPOS = ccpos,
            CDPPOS = cdppos, 
            ERRPOS = errpos,
            DVVNEG = dvvneg, 
            CCNEG = ccneg,
            CDPNEG = cdpneg,
            ERRNEG = errneg,
        )
        return df
    end

    function stretch_freqs(
        path::String, 
        DVVDIR::String,
        freqmin::AbstractArray,
        freqmax::AbstractArray, 
        tmin::AbstractArray, 
        tmax::AbstractArray, 
        dvmin::Real, 
        dvmax::Real;
        ntrial::Int=101,
        interval::Period=Day(1),
    )
        C = deserialize(path)
        if all(isnan.(C.corr))
            return nothing
        end
        for ii = 1:length(freqmin)
            println("Computing dv/v for $(path) $(freqmin[ii])-$(freqmax[ii]) Hz $(now())")
            df = do_stretch(
                deepcopy(C),
                freqmin[ii],
                freqmax[ii],
                tmin[ii],
                tmax[ii],
                dvmin,
                dvmax,
                ntrial=ntrial,
                interval=interval,
            )
            DVVFREQ = joinpath(DVVDIR,"$(freqmin[ii])-$(freqmax[ii])")
            mkpath(DVVFREQ)
            dvvpath = joinpath(DVVFREQ,"$(basename(path)).arrow")
            Arrow.write(dvvpath,df)
        end
        return nothing 
    end

    # dv/v parameters 
    freqmin = 2. .^ (0:3)
    freqmax = 2 .* freqmin 
    tmin = round.(0 ./ freqmin, digits=2) 
    tmax = round.(20 ./ freqmin, digits=2)
    dv = 0.05
    dvmin, dvmax = -dv, dv
    dvvdays = 90 
    interval = Day(dvvdays)

    # directories 
    ROOTDIR = "/media/FOUR/data/"
    SCDIR = joinpath(ROOTDIR,"ONECORR")
    NCDIR = joinpath(ROOTDIR,"NCCORR")
    DVVDIR = joinpath(ROOTDIR,"DVV-$dvvdays-DAY")
    if !isdir(DVVDIR)
        mkpath(DVVDIR)
    end

    corrs = [glob("*",SCDIR);glob("*",NCDIR)]
end
robust_pmap(x -> stretch_freqs(
    x,
    DVVDIR,
    freqmin,
    freqmax,
    tmin,
    tmax,
    dvmin,
    dvmax,
    ntrial=101,
    interval=interval,
    ),
    corrs,
)
