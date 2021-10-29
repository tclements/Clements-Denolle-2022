using Distributed 
addprocs()

@everywhere begin 
    using Dates
    using Glob 
    using JLD2
    using SeisIO
    using SeisNoise 
    using Statistics

    function seischannel2psd(
        file::String,
        fs::Real,
        RESP::SeisChannel,
        cc_len::Real,
        cc_step::Real,
        freqmin::Real,
        freqmax::Real;
        responsefreq::Real=0.1,
    )
        println("Reading file $file, $(now())")
        S = read_data("mseed", file)
        # if can't merge, return NaNs 
        try 
            merge!(S)
            sort!(S)
            taper!(S,t_max=100.)
            ungap!(S,m=false)
            demean!(S)
            detrend!(S) 
        catch e
            for ii in 1:S.n
                S.x[ii] .= NaN
            end
            return S
        end

        # find and replace nans with zeros 
        nans = Dict()
        for ii in 1:S.n
            nans[ii] = findall(isnan.(S.x[ii]))
            S.x[ii][nans[ii]] .= 0 
        end

        # highpass filter 
        highpass!(S,responsefreq,zerophase=true,corners=2)

        # remove instrument response
        for ii in 1:S.n
            S.loc[ii] = RESP.loc
            S.gain[ii] = RESP.gain
            S.resp[ii] = RESP.resp
        end
        remove_resp!(S)

        for ii = 1:S.n
            S.x[ii] ./= S.gain[ii]
        end
        
        # replace gaps with zeros
        for ii in 1:S.n
            S.x[ii][nans[ii]] .= 0 
        end 

        resample!(S,fs=fs) 

        # convert to RawData
        R = RawData(S[1],cc_len,cc_step) 
        detrend!(R)
        taper!(R,max_percentage=0.1)
        highpass!(R,freqmin,zerophase=true)
        F = rfft(R)

        nfft = size(R.x,1)
        nF = size(F.fft,1)

        # get psd for each sample 
        psd = abs2.(F.fft) .* 2 ./ nfft ./ fs .* 1.142857
        return median(psd,dims=2)[:]
    end

    function julday2Date(s::String)
        year = parse(Int,s[1:4])
        days = parse(Int,s[6:end])
        return Date(year) + Day(days - 1)
    end

    # parameters used for spectrum 
    fs = 40.0 
    cc_len = 1200.0 
    cc_step = 300.0 
    freqmin = 0.5 
    freqmax = 9.5 
    responsefreq = 0.4

    # gather files 
    files = glob("*/*/*","/media/FOUR/data/continuous_waveforms/")
    files = [f for f in files if occursin("LJR__BHZ",f)]

    # filter based on file sizes
    filesizes = [filesize(f) for f in files]
    med = median(filesizes) 
    ind = findall(filesizes .> med / 10)
    files = files[ind]
    N = length(files)

    # load instrument response for CI.LJR..BHZ 
    RESP = FDSNsta("CI.LJR..BHZ",s="2002-01-01",t="2021-05-01",src="SCEDC")
    ind = findfirst(RESP.id .== "CI.LJR..BHZ")
    RESP = RESP[ind]
end

tpsd = Array{Date}(undef,N)
# get median spectrum for each day 
out = pmap( 
    x -> seischannel2psd(
        x,
	    fs,
	    RESP,
	    cc_len,
	    cc_step,
	    freqmin,
	    freqmax,
        responsefreq=responsefreq,
    ),
    files,
)
psd = zeros(eltype(out[1]),size(out[1],1),size(out,1))
for ii in 1:N
    psd[:,ii] .= out[ii]
    tpsd[ii] = julday2Date(basename(dirname(files[ii])))
end

# get frequencies of interest 
fpsd = SeisNoise.FFTW.rfftfreq(Int(cc_len * fs),fs)
find = findall((fpsd .>= freqmin) .& (fpsd .<= freqmax))
psd = psd[find,:]
fpsd = fpsd[find]

# remove NaNs 
nanind = unique([c[2] for c in findall(isnan.(psd))])
keepind = findall(.!in(nanind),1:size(psd,2))
psd = psd[:,keepind]
tpsd = tpsd[keepind]

# remove high amplitudes 
ampind = unique([c[2] for c in findall(psd .> 1e5)])
keepind = findall(.!in(ampind),1:size(psd,2))
psd = psd[:,keepind]
tpsd = tpsd[keepind]

# average over 1/128 octaves 
fbins = 2 .^ (range(-1,stop=3.25,length=Int(4.25*128)))
Nbins = length(fbins) - 1 
smoothPSD = zeros(eltype(psd), Nbins, size(psd,2))
fsmooth = zeros(Nbins)
for ii in 1:Nbins
    fbinind = findall((fpsd .>= fbins[ii]) .& (fpsd .<= fbins[ii + 1]))
    smoothPSD[ii,:] = mean(psd[fbinind,:],dims=1)
    fsmooth[ii] = mean(fbins[ii:ii+1])
end
psd = smoothPSD 
fpsd = fsmooth 

# save to JLD2 
@save "/media/FOUR/data/LJR-psd.jld2" psd tpsd fpsd