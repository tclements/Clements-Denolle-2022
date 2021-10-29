using Arrow
using CSV
using DataFrames 
using Dates
using DSP
using Glob 
using LsqFit
using NetCDF
using Plots
using Plots.Measures
using Statistics

import GR
GR.inline("png")

function CDM(A::AbstractArray,k::Int)
    Amean = [cumsum(A[1:k-1]) ./ (1:k-1) ;rolling_mean(A,k)]
    return cumsum(A .- Amean)
end

function rolling_mean(A::AbstractArray,k::Int)
    B = cumsum([zero(eltype(A)); A])
    return (B[k+1:end] - B[1:end-k]) ./ k
end

# create lowpass butterworth filter 
responsetype = Lowpass(0.01,fs=1.)
designmethod = Butterworth(4)

# load precip data 
filename = "/media/FOUR/data/ppt.nc"
lon = ncread(filename,"lon")
lat = ncread(filename,"lat")
tppt = ncread(filename,"t")
tppt = Date.(Dates.unix2datetime.(tppt))
ppt = ncread(filename,"ppt")
tpptday = (tppt .- tppt[1]) ./ Day(1)

# load tmean data 
filename = "/media/FOUR/data/tmean.nc"
tmean = ncread(filename,"tmean")
ttmean = ncread(filename,"t")
ttmean = Date.(Dates.unix2datetime.(ttmean))
ttmeanday = (ttmean .- ttmean[1]) ./ Day(1)

# get CI station locations 
SCdf = DataFrame(CSV.File("/home/timclements/CALI/CIstations.csv"))
NCdf = DataFrame(CSV.File("/home/timclements/CALI/NCstations.csv"))
CAdf = vcat(NCdf, SCdf)

# number of days 
days = 90 

# load dv/v data 
DVVlist = glob("*","/media/FOUR/data/DVV-$days-DAY-COMP/2.0-4.0")
FIGDIR = "/media/FOUR/data/FIGURES-FIT-DVV/$days-DAY"
if !isdir(FIGDIR)
    mkpath(FIGDIR)
end

# df to hold parameters 
fitdf = DataFrame()

for ii = 1:length(DVVlist)
    file = DVVlist[ii]
    df = Arrow.Table(file) |> Arrow.columntable |> DataFrame
    sort!(df,:DATE)
    netsta = replace(basename(file),".arrow"=>"")
    net, sta = split(netsta,'.')

    if size(df,1) < 365
        continue
    end

    # get precip nearest station lat, lon 
    staind = findfirst(CAdf[!,:Station] .== sta)
    if isnothing(staind)
        continue
    end
    stalat = CAdf[staind,:Latitude]
    stalon = CAdf[staind,:Longitude]
    lonind = argmin(abs.(lon .- stalon))
    latind = argmin(abs.(lat .- stalat))
    precip = ppt[latind,lonind,:]
    temp = tmean[latind,lonind,:]
    temp .-= mean(temp)

    # skip stations that are outside PRISM for now 
    if iszero(temp)
        continue
    end
    tempfilt = filtfilt(digitalfilter(responsetype, designmethod), temp)

    # find best fitting temperature 
    maxlag = isodd(length(df[!,:DVV])) ? size(df[!,:DVV],1) : size(df[!,:DVV],1) - 1
    lags = -maxlag:maxlag
    lagind = findall((lags .<= 360) .& (lags .>= 0))
    tempind = findall(in(df[!,:DATE]),ttmean)
    tempdvvcorr = xcorr(tempfilt[tempind],df[!,:DVV])[lagind]
    tempcor = cor(df[!,:DVV],tempfilt[tempind])
    tempshift = argmax(tempdvvcorr) 

    # find common times between DV/V and shifted temperature 
    shiftind = findall(in(df[!,:DATE]),ttmean .- Day(tempshift))
    dvvshiftind = findall(in(ttmean .- Day(tempshift)),df[:,:DATE])
    tempcor = cor(df[dvvshiftind,:DVV],tempfilt[shiftind])

    # find best fitting precip model against dv/v 
    days = 180 : 14 * 365
    N = length(days)
    cors = zeros(N)
    pptminind = findfirst(tppt .== df[1,:DATE])
    pptmaxind = findfirst(tppt .== df[end,:DATE])
    # find best CDM date range 
    for jj in 1:N
        CDMii = CDM(precip[pptminind - days[jj]:end],days[jj])
        pptind = findall(in(df[!,:DATE]),tppt[pptminind - days[jj]:end])
        cors[jj] = cor(df[!,:DVV],CDMii[pptind])
    end
    CDMbest = CDM(precip[pptminind - days[argmin(cors)]:end],days[argmin(cors)]) ./ 1000
    CDMbest .-= CDMbest[1]
    pptfit = tppt[pptminind - days[argmin(cors)]:end]
    precipind = findall(in(df[!,:DATE]),pptfit)

    # fit precip + temperature model to dv/v
    @. multimodel(x, p) = p[1] + p[2] * x[:,1] + p[3] * x[:,2]
    dvvfit = curve_fit(
        multimodel,
        [CDMbest[precipind][dvvshiftind] tempfilt[shiftind]],
        df[dvvshiftind,:DVV],
        [0.,-1.,1.],
    )
    fitcor = cor(
        df[dvvshiftind,:DVV],
        tempfilt[shiftind] .* 
        coef(dvvfit)[3] .+ 
        CDMbest[precipind][dvvshiftind] .* 
        coef(dvvfit)[2] .+ 
        coef(dvvfit)[1]
    )

    println("$(netsta): GWL $(round(minimum(cors),digits=2)), TEMP $(round(tempcor,digits=2)), FIT $(round(fitcor,digits=2))")
    # plot predicted dv/v 
    ylims = (quantile(df[!,:DVV],0.05) * 2,quantile(df[!,:DVV],0.95) * 2)
    scatter(
        Date.(df[!,:DATE]),
        df[!,:DVV],
        alpha=df[!,:CC] ./ 10,
        label="",
        seriescolor=:Reds_9,
        marker_z=df[!,:CC],
        ylims=ylims,
        colorbar=false,
        ylabel="dv/v [%]",
        right_margin = 3cm,
        left_margin = 1.5cm,
        yaxis=:flip,
        yguidefontcolor=:darkred,
        ytickfont=font(10,:darkred),
        ytick_direction=:out,
        yforeground_color_axis=:darkred,
        dpi=500,
        size=(800,400),
        xlims=(df[1,:DATE],df[end,:DATE]),
    )
    plot!(
        pptfit[precipind][dvvshiftind],
        tempfilt[shiftind] .* coef(dvvfit)[3] .+ CDMbest[precipind][dvvshiftind] .* coef(dvvfit)[2] .+ coef(dvvfit)[1] ,
        c=:dodgerblue,
        xlims=(df[1,:DATE],df[end,:DATE]),
        linewidth=2.5,
        alpha=0.85,
        label="dv/v Fit",
    )
    savefig(joinpath(FIGDIR,"$(netsta)-PREDICT-DVV.png"))

    append!(
        fitdf,
        DataFrame(
            Dict(
                :NETSTA=>netsta,
                :PRECIPCOR=>minimum(cors),
                :TEMPCOR=>tempcor,
                :FITCOR=>fitcor,
                :BIAS=>coef(dvvfit)[1],
                :PRECIPSCALE=>coef(dvvfit)[2],
                :TEMPSCALE=>coef(dvvfit)[3],
                :TEMPSHIFT=>-tempshift,
                :PRECIPSHIFT=>days[argmin(cors)],
                :LAT=>stalat,
                :LON=>stalon,
            )
        )
    )
end

# write to arrow 
Arrow.write("/media/FOUR/data/dvvfit-$days-day.arrow",fitdf)