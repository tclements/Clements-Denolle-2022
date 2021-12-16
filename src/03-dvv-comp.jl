using Arrow
using DataFrames
using Dates 
using Glob 
using Statistics 

ROOTDIR = "/media/FOUR/data/"
DVVDIR = joinpath(ROOTDIR,"DVV-90-DAY")
COMPDVV = joinpath(ROOTDIR,"DVV-90-DAY-COMP")
if !isdir(COMPDVV)
    mkpath(COMPDVV)
end

# load seismic data
freqdirs = glob("*",DVVDIR)

# loop through all stations
for ii in 1:length(freqdirs)
    freqdir = freqdirs[ii]
    freqmin, freqmax = parse.(Float64,(split(basename(freqdir),"-")))
    freqarrow = sort(glob("*",freqdir))
    freqstas = unique([join(split(basename(f),".")[1:2],".") for f in freqarrow])
    stas = [sort(freqarrow[ii:ii+2]) for ii = 1:3:length(freqarrow)-1]

    OUTDIR = joinpath(COMPDVV,"$freqmin-$freqmax")
    if !isdir(OUTDIR)
        mkpath(OUTDIR)
    end

    for jj = 1:length(stas)
        sta = stas[jj]
        freqsta = join(split(basename(sta[1]),".")[1:2],".")
        println("$freqsta $freqmin-$freqmax $(now())")

        # get dvv
        EN = DataFrame(Arrow.columntable(Arrow.Table(sta[1])))
        EZ = DataFrame(Arrow.columntable(Arrow.Table(sta[2])))
        NZ = DataFrame(Arrow.columntable(Arrow.Table(sta[3])))

        # get unqiue rows 
        unique!(EN,:DATE)
        unique!(EZ,:DATE)
        unique!(NZ,:DATE)

        # remove high error days 
        # EN = EN[(EN[!,:ERRPOS] .< maxerr) .& (EN[!,:ERRNEG] .< maxerr),:]
        # EZ = EZ[(EZ[!,:ERRPOS] .< maxerr) .& (EZ[!,:ERRNEG] .< maxerr),:]
        # NZ = NZ[(NZ[!,:ERRPOS] .< maxerr) .& (NZ[!,:ERRNEG] .< maxerr),:]

        dates = intersect(EN[!,:DATE],EZ[!,:DATE],NZ[!,:DATE])
        if length(dates) == 0 
            continue
        end
        EN = EN[[d in dates for d in EN[!,:DATE]],:]
        EZ = EZ[[d in dates for d in EZ[!,:DATE]],:]
        NZ = NZ[[d in dates for d in NZ[!,:DATE]],:]

        # compute CC from Hobiger, 2014
        CC = EN[!,:CCPOS] .^ 2 .+ EZ[!,:CCPOS] .^ 2 .+ NZ[!,:CCPOS] .^ 2 .+
            EN[!,:CCNEG] .^ 2 .+ EZ[!,:CCNEG] .^ 2 .+ NZ[!,:CCNEG] .^ 2

        # compute dv/v from Hobiger, 2014
        DVV = (EN[!,:CCPOS] .^ 2 .* EN[!,:DVVPOS]) .+ (EZ[!,:CCPOS] .^ 2 .*
            EZ[!,:DVVPOS]) .+ (NZ[!,:CCPOS] .^ 2 .* NZ[!,:DVVPOS]) .+
            (EN[!,:CCNEG] .^ 2 .* EN[!,:DVVNEG]) .+ (EZ[!,:CCNEG] .^ 2 .*
                    EZ[!,:DVVNEG]) .+ (NZ[!,:CCNEG] .^ 2 .* NZ[!,:DVVNEG])
        DVV ./= CC
        DVV[isnan.(DVV)] .= 0.0 
        CC = (EN[!,:CCPOS] .^ 3 .+ EZ[!,:CCPOS] .^ 3 .+ NZ[!,:CCPOS] .^ 3 .+
            EN[!,:CCNEG] .^ 3 .+ EZ[!,:CCNEG] .^ 3 .+ NZ[!,:CCNEG] .^ 3) ./ CC
        CC[isnan.(CC)] .= 0.0

        # write to dataframe 
        df = DataFrame("DATE"=>dates,"DVV"=>DVV,"CC"=>CC)
        savepath = joinpath(OUTDIR,"$freqsta.arrow")
        Arrow.write(savepath,df)

    end 
end
