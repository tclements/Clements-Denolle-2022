using Distributed
addprocs()

@everywhere begin 
    using Arrow
    using AWS
    using AWSS3
    using DataFrames
    using Dates
    using Serialization
    using SCEDCCorr
    using SeisNoise
    using SeisDvv

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
        remove_max::Bool=false
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
    
        # get reference
        lags = -C.maxlag:1/C.fs:C.maxlag
        times = Dates.unix2datetime.(C.t)
        minpad = max(5 / C.fs, dvmax * minlag)
        maxpad = max(5 / C.fs, dvmax * maxlag)
        indpos = findall((lags .>= minlag - minpad) .& (lags .<= maxlag + maxpad))
        indneg = findall((lags .>= -maxlag - maxpad) .& (lags .<= -minlag + minpad))
        indboth = findall(abs.(lags) .<= maxlag + maxpad)
        window = findall((lags[indpos] .>= minlag) .& (lags[indpos] .<= maxlag))
        bothwindow = findall(abs.(lags[indboth]) .<= maxlag)
        N = length(indpos)
        A = C.corr
        ref = SeisNoise.stack(C,allstack=true)
        if remove_max
            abs_max!(ref)
        end
        r = ref.corr
    
        dvvpos, dvvneg, dvvboth = zeros(size(A,2)), zeros(size(A,2)), zeros(size(A,2))
        ccpos, ccneg, ccboth = zeros(size(A,2)), zeros(size(A,2)), zeros(size(A,2))
        cdppos, cdpneg, cdpboth = zeros(size(A,2)), zeros(size(A,2)), zeros(size(A,2))
        errpos, errneg, errboth = zeros(size(A,2)), zeros(size(A,2)), zeros(size(A,2))
        allCpos, allCneg, allCboth = zeros(ntrial,size(A,2)), zeros(ntrial,size(A,2)), zeros(ntrial,size(A,2))
    
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
            dvvboth[ii],ccboth[ii],cdpboth[ii], ϵ, errboth[ii], allCboth[:,ii] = stretching(
                r[indboth],
                A[indboth,ii],
                lags[indboth],
                bothwindow,
                freqmin,
                freqmax,
                dvmin=dvmin,
                dvmax=dvmax,
                ntrial=ntrial,
            )
        end
    
        # set negative CC to zero  
        ccpos[ccpos .< 0.0] .= 0.0
        ccneg[ccneg .< 0.0] .= 0.0 
        ccboth[ccboth .< 0.0] .= 0.0 

        cdppos[cdppos .< 0.0] .= 0.0
        cdpneg[cdpneg .< 0.0] .= 0.0 
        cdpboth[cdpboth .< 0.0] .= 0.0 
    
        # load into DataFrame, then write to disk
        df = DataFrame(
            DATE = Date.(times), 
            DVVPOS = dvvpos, 
            CCPOS = ccpos,
            CDPPOS = cdppos, 
            DVVNEG = dvvneg, 
            CCNEG = ccneg,
            CDPNEG = cdpneg,
            DVVBOTH = dvvboth, 
            CCBOTH = ccboth,
            CDPBOTH = cdpboth,
        )
        return df
    end

    function stretch_freqs(
        aws::AWSConfig,
        bucket::String,
        path::String, 
        freqmin::AbstractArray,
        freqmax::AbstractArray, 
        tmin::AbstractArray, 
        tmax::AbstractArray, 
        dvmin::Real, 
        dvmax::Real;
        ntrial::Int=101,
    )
        C = s3_get_autocorr(aws,bucket,path)
        for ii = 1:length(freqmin)
            println("Computing dv/v for $(path) $(freqmin[ii])-$(freqmax[ii]) Hz $(now())")
            df = do_stretch(deepcopy(C),freqmin[ii],freqmax[ii],tmin[ii],tmax[ii],dvmin,dvmax,ntrial=ntrial)
            io = IOBuffer()
            Arrow.write(io,df)
            dvvpath = "NC-DVV/$(freqmin[ii])-$(freqmax[ii])/$(basename(path)).arrow"
            s3_put(aws,bucket,dvvpath,io.data)
        end
        return nothing 
    end

    # dv/v parameters 
    freqmin = 2. .^ (0:3)
    freqmax = 2 .* freqmin 
    tmin = round.(4 ./ freqmin, digits=2) 
    tmax = round.(16 ./ freqmin, digits=2)
    dv = 0.1
    dvmin, dvmax = -dv, dv  

    # get all autocorrelations from s3 
    aws = AWSConfig(region="us-west-2")
    upload_bucket = "bh-auto-corr"
    autolist = collect(s3_list_keys(aws,upload_bucket,"NCCORR/"))
end

pmap(
    x->stretch_freqs(
        aws,
        upload_bucket,
        x,
        freqmin,
        freqmax,
        tmin,
        tmax,
        dvmin,
        dvmax,
    ),
    autolist,
)