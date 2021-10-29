using Distributed
addprocs()

@everywhere begin 
    using Arrow
    using DataFrames
    using Dates
    using Glob
    using Serialization
    using SCEDCCorr
    using SeisIO
    using SeisNoise
    using SeisDvv
    using Statistics

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

    function do_mwcs(C::CorrData, 
        freqmin::Real,
        freqmax::Real, 
        minlag::Real, 
        tmin::Real,
        window_length::Real,
        window_step::Real,
        smoothing_half_win::Int,
        fs::Real,
    )

        # filter data 
        remove_nan!(C)
        clean_up!(C,freqmin,freqmax)

        # check for nonunique days and sort 
        ind = indexin(unique(C.t), C.t)
        C.t = C.t[ind]
        C.corr = C.corr[:,ind]
        ind = sortperm(C.t)
        C.t = C.t[ind]
        C.corr = C.corr[:,ind]

        # get correct lags 
        lags = -C.maxlag:1/C.fs:C.maxlag
        ind = findall(abs.(lags) .<= tmin)
        A = deepcopy(C.corr)
        T = eltype(A)
        Ncol,Nrows = size(A)
        ref = mean(A,dims=2)
        dttright = zeros(T,Nrows)
        errright = similar(dttright)
        dttleft = similar(dttright)
        errleft = similar(dttright)
        dttboth = similar(dttright)
        errboth = similar(dttright)

        for ii = 1:Nrows
            # left dtt first 
            cur = A[:,ii]
            time_axis, dt, err, mcoh = mwcs(
                ref[ind],
                cur[ind],
                freqmin,
                freqmax,
                fs,
                -tmin,
                window_length,
                window_step,
                smoothing_half_win,
            )

            # left dv/v 
            m, em, a, ea, m0, em0 = mwcs_dvv(
                time_axis,
                dt,
                err,
                mcoh,
                "static",
                0.,
                0.,
                minlag,
                tmin-minlag,
                "left",
            )
            dttleft[ii] = m0
            errleft[ii] = em0

            # right dv/v 
            m, em, a, ea, m0, em0 = mwcs_dvv(
                time_axis,
                dt,
                err,
                mcoh,
                "static",
                0.,
                0.,
                minlag,
                tmin-minlag,
                "right",
            )
            dttright[ii] = m0
            errright[ii] = em0

            # both dv/v 
            m, em, a, ea, m0, em0 = mwcs_dvv(
                time_axis,
                dt,
                err,
                mcoh,
                "static",
                0.,
                0.,
                minlag,
                tmin-minlag,
                "both",
            )
            dttboth[ii] = m0
            errboth[ii] = em0
        end

        # load into DataFrame, then write to disk
        df = DataFrame(
            DATE = Date.(u2d.(C.t)), 
            DVVPOS = dttright .* -100, 
            ERRPOS = errright,
            DVVNEG = dttleft .* -100, 
            ERRNEG = errleft,
            DVVBOTH = dttboth .* -100, 
            ERRBOTH = errboth,
        )

        return df
    end

    function stretch_freqs(
        path::String, 
        DVVDIR::String,
        freqmin::AbstractArray,
        freqmax::AbstractArray, 
        minlag::AbstractArray, 
        tmin::AbstractArray,
        window_length::AbstractArray,
        window_step::AbstractArray,
        smoothing_half_win::AbstractArray,
        fs::Real,

    )
        C = deserialize(path)
        for ii = 1:length(freqmin)
            println("Computing dv/v for $(path) $(freqmin[ii])-$(freqmax[ii]) Hz $(now())")
            df = do_mwcs(
                deepcopy(C),
                freqmin[ii],
                freqmax[ii],
                minlag[ii],
                tmin[ii],
                window_length[ii],window_step[ii],
                smoothing_half_win[ii],
                fs,
            )
            DVVFREQ = joinpath(DVVDIR,"$(freqmin[ii])-$(freqmax[ii])")
            mkpath(DVVFREQ)
            dvvpath = joinpath(DVVFREQ,"$(basename(path)).arrow")
            Arrow.write(dvvpath,df)
        end
        return nothing 
    end

    ROOTDIR = "/media/FOUR/data/"
    CORRDIR = joinpath(ROOTDIR,"NCCORR")
    DVVDIR = joinpath(ROOTDIR,"NC-MWCS-DVV")
    if !isdir(DVVDIR)
        mkpath(DVVDIR)
    end

    # dv/v parameters 
    fs = 40. 
    freqmin = [2.]
    freqmax = [4.]
    minlag = [0.]
    tmin = [5.]
    window_length = [2.75]
    window_step = [0.25]
    smoothing_half_win = [10]

    corrs = glob("*",CORRDIR)
end
robust_pmap(x -> stretch_freqs(
    x,
    DVVDIR,
    freqmin,
    freqmax,
    minlag,
    tmin,
    window_length,
    window_step,
    smoothing_half_win,
    fs,
    ),
    corrs,
)
