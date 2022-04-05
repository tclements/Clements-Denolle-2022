using Arrow 
using CSV
using Dates
using DataFrames 
using DelimitedFiles
using DSP 
using Glob 
using GLM 
using GMT
using Interpolations 
using NetCDF
using Statistics
using StatsModels 

function DVVslope(df::DataFrame,ind::AbstractArray)
    t = (df[ind,:DATE] .- df[ind[1],:DATE]) ./ Day(1) .+ 1

    # create linear model before 
    data = DataFrame(X=t ./ t[end],Y=df[ind,:DVV])
    ols = lm(@formula(Y ~ X), data, wts=df[ind,:CC])

    # get relevant data 
    inter, slope = coef(ols)
    rtwo = r2(ols)
    interSTD, slopeSTD = stderror(ols)
    return inter, slope, rtwo
end

"""
    smooth_withfiltfilt(A::AbstractArray; window_len::Int=7, window::Symbol=:bartlett)
Apply DSP.filtfilt() to smooth the waveform.
# Arguments
- "A::AbstractArray": periodogram of spectrum
- "window_len::Int=7": window length
- "window::Symbol=:bartlett": window type: (see https://juliadsp.github.io/DSP.jl/stable/windows/)
"""
function smooth_withfiltfilt(A::AbstractArray; window_len::Int=11, window::Symbol=:rect)
    w = getfield(DSP.Windows, window)(window_len)
    A = DSP.filtfilt(w ./ sum(w), A)
    return A
end

#load station locations 
NCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/NCstations.csv")))
SCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/CIstations.csv")))
CAdf = vcat(NCdf, SCdf)

# lat/lon for LA map 
LAminlon = -119.5
LAmaxlon = -116.5 
LAminlat = 32.75 
LAmaxlat = 34.75 

# load GRACE data 
# from http://www2.csr.utexas.edu/grace/RL06_mascons.html
filename = joinpath(@__DIR__,"../data/CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02.nc")
tlwe = Date(2002,1,1) .+ Day.(floor.(ncread(filename,"time")))
lon = ncread(filename,"lon")
lat = ncread(filename,"lat")
# lwe ./= 100 # convert to m 
ind2004 = findfirst(tlwe .== Date(2004,10,16))
ind2005a = findfirst(tlwe .== Date(2005,5,16))
ind2005b = findfirst(tlwe .== Date(2005,10,16))
ind2011 = findfirst(tlwe .== Date(2011,9,16))
ind2016 = findfirst(tlwe .== Date(2016,8,21))
ind2020 = findfirst(tlwe .== Date(2020,10,16))
G2004 = gmtread(filename,layer=ind2004,varname="lwe_thickness")
G2005a = gmtread(filename,layer=ind2005a,varname="lwe_thickness")
G2005b = gmtread(filename,layer=ind2005b,varname="lwe_thickness")
G2011 = gmtread(filename,layer=ind2011,varname="lwe_thickness")
G2016 = gmtread(filename,layer=ind2016,varname="lwe_thickness")
G2020 = gmtread(filename,layer=ind2020,varname="lwe_thickness")

# load fitted values 
fitdf = Arrow.Table(joinpath(@__DIR__,"../data/hydro-model-90-day.arrow")) |> Arrow.columntable |> DataFrame

# load temperature values 
filename = joinpath(@__DIR__,"../data/tmean.nc")
pptlon = ncread(filename,"lon")
pptlat = ncread(filename,"lat")
tmean = ncread(filename,"tmean")
ttmean = ncread(filename,"t")
ttmean = Date.(Dates.unix2datetime.(ttmean))
ttmeanday = (ttmean .- ttmean[1]) ./ Day(1)

# plot maps 2-4 Hz 
freqmin = 2.0
freqmax = 4.0

# Oct 2004 - May 2005 
maxgap = 15 # maximum data gap [days]
mindays = 100 # minimum number of days needed for analysis 
gray = makecpt(color=150, range=(-10000,10000), no_bg=:true);
arrowfiles = glob("*",joinpath(@__DIR__,"../data/DVV-90-DAY-COMP/$freqmin-$freqmax/"))
ΔDVVdf = DataFrame()

for jj in 1:length(arrowfiles)
    # read dv/v for each station 
    DVV = Arrow.Table(arrowfiles[jj]) |> Arrow.columntable |> DataFrame
    ind = findall(Date(2005,1,1) .< DVV[:,:DATE] .< Date(2005,5,1))
    # only keep stations with > 2 years of data after 2012 
    if length(ind) < mindays
        continue 
    end

    # check for gaps in the data 
    tDVV = (DVV[ind,:DATE] .- DVV[ind[1],:DATE]) ./ Day(1) .+ 1
    if maximum(diff(tDVV)) > maxgap
        continue 
    end

    # get station name and location 
    netsta = replace(basename(arrowfiles[jj]),".arrow"=>"")
    net, sta = split(netsta,".")
    staind = findfirst(CAdf[!,:Station] .== sta)
    if isnothing(staind)
        println("No lat/lon for $netsta")
        continue 
    end
    stalat = CAdf[staind,:Latitude]
    stalon = CAdf[staind,:Longitude]

    println("Calculating Δdv/v $netsta $(now())")

    # remove temperature effects
    fitind = findfirst(fitdf[:,:NETSTA] .== netsta)
    lonind = argmin(abs.(pptlon .- fitdf[fitind,:LON]))
    latind = argmin(abs.(pptlat .- fitdf[fitind,:LAT]))
    temp = tmean[latind,lonind,:]
    temp .-= mean(temp)
    tempfilt = smooth_withfiltfilt(temp,window_len=45)
    tempind = findall(DVV[1,:DATE] .<= ttmean .<= DVV[end,:DATE])
    # ttemp = ttmean[tempind]

    # filter temperature data 
    dvvT = tempfilt[tempind .- round(Int,fitdf[jj,:E5])] 
    dvvT .-= mean(dvvT)
    dvvT ./= std(dvvT)
    dvvT .*= fitdf[fitind,:E4]

    # remove temperature from dv/v 
    dvvtempind = findall(in(ttmean[tempind]),DVV[ind,:DATE])
    DVV[ind,:DVV] .-= dvvT[dvvtempind]  

    # create linear model 
    data = DataFrame(X=tDVV ./ tDVV[end],Y=DVV[ind,:DVV])
    ols = lm(@formula(Y ~ X), data, wts=DVV[ind,:CC] .^ 3)

    # get relevant data 
    inter, slope = coef(ols)
    rtwo = r2(ols)
    interSTD, slopeSTD = stderror(ols)
    mindate = DVV[ind[1],:DATE]
    maxdate = DVV[ind[end],:DATE]

    # create dataframe
    df = DataFrame(
        MINDATE=mindate,
        MAXDATE=maxdate,
        NOBS=size(data,1),
        INTER=inter,
        SLOPE=slope,
        R2=rtwo,
        INTERSTD=interSTD,
        SLOPESTD=slopeSTD,
        NETSTA=netsta,
        NET=net,
        STA=sta,
        LON=stalon,
        LAT=stalat,
    )
    append!(ΔDVVdf,df)
end

# plot with GMT 
minlon = minimum(ΔDVVdf[:,:LON]) - 0.5
maxlon = maximum(ΔDVVdf[:,:LON]) + 0.5
minlat = minimum(ΔDVVdf[:,:LAT]) - 0.5
maxlat = maximum(ΔDVVdf[:,:LAT]) + 0.5

# plot Δdv/v 
colorlimit=1.0
CDVV = GMT.makecpt(T=(-colorlimit,colorlimit,colorlimit / 101), cmap=:vik, reverse=true);
CGRACE = GMT.makecpt(T=(-30,30,30 / 101), cmap=:vik);
GMT.grdimage(
    "@earth_relief_15s",
    cmap=gray, 
    J=:guess,
    shade=true,
    region=(minlon,maxlon,minlat,maxlat),
    coast=true,
    colorbar=false,
    )
GMT.coast!(
    N="a",
    ocean=:white,
)
GMT.grdimage!(
    -(G2005a .- G2004),
    cmap=CGRACE,
    alpha=50,
    colorbar=false,
)
GMT.scatter!(
    ΔDVVdf[:,:LON],
    ΔDVVdf[:,:LAT],
    zcolor=-ΔDVVdf[:,:SLOPE],
    markeredgecolor=:black,
    marker=:c,
    transparency= 15,
    # colorbar=true,
    cmap=CDVV,
    markersize="6p",
)
GMT.colorbar!(
    C=CDVV,
    frame=(annot=:auto, ticks=:auto, xlabel="dv/v [%]"),
    pos=(anchor=:BL,horizontal=true,offset=(-5.5,-2),length=5),
)
GMT.colorbar!(
    C=CGRACE,
    frame=(annot=:auto, ticks=:auto, xlabel="LWE [cm]"), 
    pos=(anchor=:BL, horizontal=true,offset=(-5.5,-2),move_annot=true,length=5), 
    W=-1,
)
rect = [LAminlon LAminlat; LAminlon LAmaxlat; LAmaxlon LAmaxlat; LAmaxlon LAminlat; LAminlon LAminlat];
GMT.plot!(
    rect, 
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed
)
GMT.plot!(
    [LAminlon LAmaxlat; -119.65 38.9],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed,
    alpha=15,
)
GMT.plot!(
    [LAmaxlon LAmaxlat; -114.55 38.9],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed, 
    alpha=15,
)
# inset map 
tt = mapproject(region=(LAminlon,LAmaxlon,LAminlat,LAmaxlat), proj=:merc, figsize=6, map_size=true);
mapW = tt[1].data[1];   mapH = tt[1].data[2]
GMT.basemap!(inset=(size=(mapW, mapH), anchor=:TR, width=1, offset=(-6.75, -10.25), save="xx000"))
t = readdlm("xx000")
GMT.grdimage!(
    "@srtm_relief_03s",
    cmap=gray, 
    J=:merc,
    shade=true,
    region=(LAminlon,LAmaxlon,LAminlat,LAmaxlat),
    coast=false,
    colorbar=false,
    xshift=t[1], yshift=t[2],
)
GMT.coast!(
    shore=true, 
    ocean=:white,
    N="a",
)
GMT.scatter!(
    ΔDVVdf[:,:LON],
    ΔDVVdf[:,:LAT],
    zcolor=-ΔDVVdf[:,:SLOPE],
    markeredgecolor=:black,
    marker=:c,
    transparency= 15,
    # colorbar=true,
    cmap=CDVV,
    markersize="4p",
    show=1,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/DVV-2004-2005-$freqmin-$freqmax.png"),
)

# 2012 - 2016 
maxgap = 180 # maximum data gap [days]
mindays = 3 * 365 # minimum number of days needed for analysis 
arrowfiles = glob("*",joinpath(@__DIR__,"../data/DVV-90-DAY-COMP/$freqmin-$freqmax/"))
ΔDVVdf = DataFrame()

for jj in 1:length(arrowfiles)
    # read dv/v for each station 
    DVV = Arrow.Table(arrowfiles[jj]) |> Arrow.columntable |> DataFrame
    ind = findall(Date(2011,9,1) .<= DVV[:,:DATE] .<= Date(2016,9,1))
    # only keep stations with > 2 years of data after 2012 
    if length(ind) < mindays
        continue 
    end

    # make sure observations start and end in October 

    # check for gaps in the data 
    tDVV = (DVV[ind,:DATE] .- DVV[ind[1],:DATE]) ./ Day(1) .+ 1
    if maximum(diff(tDVV)) > maxgap
        continue 
    end

    # get station name and location 
    netsta = replace(basename(arrowfiles[jj]),".arrow"=>"")
    net, sta = split(netsta,".")
    staind = findfirst(CAdf[!,:Station] .== sta)
    if isnothing(staind)
        println("No lat/lon for $netsta")
        continue 
    end
    stalat = CAdf[staind,:Latitude]
    stalon = CAdf[staind,:Longitude]

    println("Calculating Δdv/v $netsta $jj")

    # remove temperature effects
    fitind = findfirst(fitdf[:,:NETSTA] .== netsta)
    lonind = argmin(abs.(pptlon .- fitdf[fitind,:LON]))
    latind = argmin(abs.(pptlat .- fitdf[fitind,:LAT]))
    temp = tmean[latind,lonind,:]
    temp .-= mean(temp)
    tempfilt = smooth_withfiltfilt(temp,window_len=45)
    tempind = findall(DVV[1,:DATE] .<= ttmean .<= DVV[end,:DATE])
    # ttemp = ttmean[tempind]

    # filter temperature data 
    dvvT = tempfilt[tempind .- round(Int,fitdf[fitind,:E5])] 
    dvvT .-= mean(dvvT)
    dvvT ./= std(dvvT)
    dvvT .*= fitdf[fitind,:E4]

    # remove temperature from dv/v 
    dvvtempind = findall(in(ttmean[tempind]),DVV[ind,:DATE])
    DVV[ind,:DVV] .-= dvvT[dvvtempind] 
    
    # create linear model 
    data = DataFrame(X=tDVV ./ tDVV[end],Y=DVV[ind,:DVV])
    ols = lm(@formula(Y ~ X), data, wts=DVV[ind,:CC] .^ 3)

    # get relevant data 
    inter, slope = coef(ols)
    rtwo = r2(ols)
    interSTD, slopeSTD = stderror(ols)
    mindate = DVV[ind[1],:DATE]
    maxdate = DVV[ind[end],:DATE]

    # create dataframe
    df = DataFrame(
        MINDATE=mindate,
        MAXDATE=maxdate,
        NOBS=size(data,1),
        INTER=inter,
        SLOPE=slope,
        R2=rtwo,
        INTERSTD=interSTD,
        SLOPESTD=slopeSTD,
        NETSTA=netsta,
        NET=net,
        STA=sta,
        LON=stalon,
        LAT=stalat,
    )
    append!(ΔDVVdf,df)
end

# plot with GMT 
minlon = minimum(ΔDVVdf[:,:LON]) - 0.5
maxlon = maximum(ΔDVVdf[:,:LON]) + 0.5
minlat = minimum(ΔDVVdf[:,:LAT]) - 0.5
maxlat = maximum(ΔDVVdf[:,:LAT]) + 0.5

# plot Δdv/v 
colorlimit=0.3
CDVV = GMT.makecpt(T=(-colorlimit,colorlimit,colorlimit / 101), cmap=:vik, reverse=true);
CGRACE = GMT.makecpt(T=(-15,15,15 / 101), cmap=:vik);
GMT.grdimage(
    "@earth_relief_15s",
    cmap=gray, 
    J=:guess,
    shade=true,
    region=(minlon,maxlon,minlat,maxlat),
    coast=true,
    colorbar=false,
    )
GMT.coast!(
    N="a",
    ocean=:white,
)
GMT.grdimage!(
    -(G2016 .- G2011),
    cmap=CGRACE,
    alpha=50,
    colorbar=false,
)
GMT.scatter!(
    ΔDVVdf[:,:LON],
    ΔDVVdf[:,:LAT],
    zcolor=-ΔDVVdf[:,:SLOPE],
    markeredgecolor=:black,
    marker=:c,
    transparency= 15,
    # colorbar=true,
    cmap=CDVV,
    markersize="6p",
)
GMT.colorbar!(
    C=CDVV,
    frame=(annot=:auto, ticks=:auto, xlabel="dv/v [%]"),
    pos=(anchor=:BL,horizontal=true,offset=(-5.5,-2),length=5),
)
GMT.colorbar!(
    C=CGRACE,
    frame=(annot=:auto, ticks=:auto, xlabel="LWE [cm]"), 
    pos=(anchor=:BL, horizontal=true,offset=(-5.5,-2),move_annot=true,length=5), 
    W=-1,
)
rect = [LAminlon LAminlat; LAminlon LAmaxlat; LAmaxlon LAmaxlat; LAmaxlon LAminlat; LAminlon LAminlat];
GMT.plot!(
    rect, 
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed
)
GMT.plot!(
    [LAminlon LAmaxlat; -119.65 38.9],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed,
    alpha=15,
)
GMT.plot!(
    [LAmaxlon LAmaxlat; -114.55 38.9],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed, 
    alpha=15,
)
# inset map 
tt = mapproject(region=(LAminlon,LAmaxlon,LAminlat,LAmaxlat), proj=:merc, figsize=6, map_size=true);
mapW = tt[1].data[1];   mapH = tt[1].data[2]
GMT.basemap!(inset=(size=(mapW, mapH), anchor=:TR, width=1, offset=(-6.75, -10.25), save="xx000"))
t = readdlm("xx000")
GMT.grdimage!(
    "@srtm_relief_03s",
    cmap=gray, 
    J=:merc,
    shade=true,
    region=(LAminlon,LAmaxlon,LAminlat,LAmaxlat),
    coast=false,
    colorbar=false,
    xshift=t[1], yshift=t[2],
)
GMT.coast!(
    shore=true, 
    ocean=:white,
    N="a",
)
GMT.scatter!(
    ΔDVVdf[:,:LON],
    ΔDVVdf[:,:LAT],
    zcolor=-ΔDVVdf[:,:SLOPE],
    markeredgecolor=:black,
    marker=:c,
    transparency= 15,
    # colorbar=true,
    cmap=CDVV,
    markersize="4p",
    show=1,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/DVV-2011-2016-$freqmin-$freqmax.png"),
)


# 2005 - 2020
maxgap = 2*365 # maximum data gap [days]
mindays = 13 * 365 # minimum number of days needed for analysis 

arrowfiles = glob("*",joinpath(@__DIR__,"../data/DVV-90-DAY-COMP/$freqmin-$freqmax/"))
ΔDVVdf = DataFrame()

for jj in 1:length(arrowfiles)
    # read dv/v for each station 
    DVV = Arrow.Table(arrowfiles[jj]) |> Arrow.columntable |> DataFrame
    ind = findall(Date(2005,10,1) .< DVV[:,:DATE] .< Date(2020,10,1)) 
    # only keep stations with > 2 years of data after 2012 
    if length(ind) < mindays
        continue 
    end

    # check for gaps in the data 
    tDVV = (DVV[ind,:DATE] .- DVV[ind[1],:DATE]) ./ Day(1) .+ 1
    if maximum(diff(tDVV)) > maxgap
        continue 
    end

    # get station name and location 
    netsta = replace(basename(arrowfiles[jj]),".arrow"=>"")
    net, sta = split(netsta,".")
    staind = findfirst(CAdf[!,:Station] .== sta)
    if isnothing(staind)
        println("No lat/lon for $netsta")
        continue 
    end
    stalat = CAdf[staind,:Latitude]
    stalon = CAdf[staind,:Longitude]

    println("Calculating Δdv/v $netsta $(now())")

    # remove temperature effects
    fitind = findfirst(fitdf[:,:NETSTA] .== netsta)
    lonind = argmin(abs.(pptlon .- fitdf[fitind,:LON]))
    latind = argmin(abs.(pptlat .- fitdf[fitind,:LAT]))
    temp = tmean[latind,lonind,:]
    temp .-= mean(temp)
    tempfilt = smooth_withfiltfilt(temp,window_len=45)
    tempind = findall(DVV[1,:DATE] .<= ttmean .<= DVV[end,:DATE])
    # ttemp = ttmean[tempind]

    # filter temperature data 
    dvvT = tempfilt[tempind .- round(Int,fitdf[fitind,:E5])] 
    dvvT .-= mean(dvvT)
    dvvT ./= std(dvvT)
    dvvT .*= fitdf[fitind,:E4]

    # remove temperature from dv/v 
    dvvtempind = findall(in(ttmean[tempind]),DVV[ind,:DATE])
    DVV[ind,:DVV] .-= dvvT[dvvtempind] 
    
    # create linear model 
    data = DataFrame(X=tDVV ./ tDVV[end],Y=DVV[ind,:DVV])
    ols = lm(@formula(Y ~ X), data, wts=DVV[ind,:CC] .^ 3)

    # get relevant data 
    inter, slope = coef(ols)
    rtwo = r2(ols)
    interSTD, slopeSTD = stderror(ols)
    mindate = DVV[ind[1],:DATE]
    maxdate = DVV[ind[end],:DATE]

    # create dataframe
    df = DataFrame(
        MINDATE=mindate,
        MAXDATE=maxdate,
        NOBS=size(data,1),
        INTER=inter,
        SLOPE=slope,
        R2=rtwo,
        INTERSTD=interSTD,
        SLOPESTD=slopeSTD,
        NETSTA=netsta,
        NET=net,
        STA=sta,
        LON=stalon,
        LAT=stalat,
    )
    append!(ΔDVVdf,df)
end

# plot with GMT 
minlon = minimum(ΔDVVdf[:,:LON]) - 0.5
maxlon = maximum(ΔDVVdf[:,:LON]) + 0.5
minlat = minimum(ΔDVVdf[:,:LAT]) - 0.5
maxlat = maximum(ΔDVVdf[:,:LAT]) + 0.5

# plot Δdv/v 
colorlimit=0.3
CDVV = GMT.makecpt(T=(-colorlimit,colorlimit,colorlimit / 101), cmap=:vik, reverse=true);
CGRACE = GMT.makecpt(T=(-10,10,10 / 101), cmap=:vik);
GMT.grdimage(
    "@earth_relief_15s",
    cmap=:gray, 
    J=:guess,
    shade=true,
    region=(minlon,maxlon,minlat,maxlat),
    coast=true,
    colorbar=false,
)
GMT.grdimage!(
    -(G2020 .- G2005b),
    cmap=CGRACE,
    alpha=50,
    colorbar=false,
)
GMT.coast!(
    N="a",
    ocean=:white,
)
GMT.scatter!(
    ΔDVVdf[:,:LON],
    ΔDVVdf[:,:LAT],
    zcolor=-ΔDVVdf[:,:SLOPE],
    markeredgecolor=:black,
    marker=:c,
    transparency= 15,
    colorbar=false,
    cmap=CDVV,
)
GMT.colorbar!(
    C=CDVV,
    frame=(annot=:auto, ticks=:auto, xlabel="dv/v [%]"),
    pos=(anchor=:BL,horizontal=true,offset=(-5.5,-2),length=5),
)
GMT.colorbar!(
    C=CGRACE,
    frame=(annot=:auto, ticks=:auto, xlabel="LWE [cm]"), 
    pos=(anchor=:BL, horizontal=true,offset=(-5.5,-2),move_annot=true,length=5), 
    W=-1,
    show=1,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/DVV-2005-2020-$freqmin-$freqmax.png"),
)

# Ridgecrest 
maxgap = 10 # maximum data gap [days]
mindays = 75 # minimum number of days needed for analysis 
arrowfiles = glob("*",joinpath(@__DIR__,"../data/FIT-DVV-SSE/90-DAY/"))
ΔDVVdf = DataFrame()

for jj in 1:length(arrowfiles)
    # read dv/v for each station 
    DVV = Arrow.Table(arrowfiles[jj]) |> Arrow.columntable |> DataFrame
    dropmissing!(DVV)
    indafter = findall(Date(2019,7,6) .<= DVV[:,:DATE] .< Date(2019,7,6) + Day(90))
    # only keep stations with > 90 days of data after Ridgecrest 
    if length(indafter) < mindays
        continue 
    end

    tafter = (DVV[indafter,:DATE] .- DVV[indafter[1],:DATE]) ./ Day(1) .+ 1
    if maximum(diff(tafter)) > maxgap
        continue 
    end

    # get station name and location 
    netsta = replace(basename(arrowfiles[jj]),".arrow"=>"")
    net, sta = split(netsta,".")
    staind = findfirst(CAdf[!,:Station] .== sta)
    if isnothing(staind)
        println("No lat/lon for $netsta")
        continue 
    end
    stalat = CAdf[staind,:Latitude]
    stalon = CAdf[staind,:Longitude]

    println("Calculating Δdv/v $netsta $(now())")

    
    # create linear model before 
    DVV[:,:DVV] .-= DVV[:,:ELASTIC]
    interA, slopeA, rtwoA = DVVslope(DVV,indafter)

    # create EQ drop 
    DVV[:,:EQ] = DVV[:,:DVV] .- DVV[:,:ELASTIC]
    diffdvv = DVV[indafter[end],:EQ] - DVV[indafter[1],:EQ]

    # create dataframe
    df = DataFrame(
        INTERA=interA,
        SLOPEA=slopeA,
        R2A=rtwoA,
        DIFFDVV=diffdvv,
        NETSTA=netsta,
        NET=net,
        STA=sta,
        LON=stalon,
        LAT=stalat,
    )
    append!(ΔDVVdf,df)
end

# plot with GMT 
minlon = -122
maxlon = -113
minlat = 32
maxlat = 40
EQlon = -117.599
EQlat = 35.778

# plot Δdv/v 
colorlimit = 1.0
C = GMT.makecpt(T=(-colorlimit,colorlimit,colorlimit / 101), cmap=:vik);
GMT.grdimage(
    "@earth_relief_15s",
    cmap=:gray, 
    J=:guess,
    shade=true,
    region=(minlon,maxlon,minlat,maxlat),
    coast=true,
    colorbar=false,
)
GMT.coast!(
    N="a",
    ocean=:white,
)
GMT.plot!(
    GMT.circgeo(
        [EQlon,EQlat], 
        radius=50, 
        unit=:k,
        np=500,
    ),
    lc=:dodgerblue,
    lw=1,
    linestyle="--",
)
GMT.plot!(
    GMT.circgeo(
        [EQlon,EQlat], 
        radius=100, 
        unit=:k,
        np=500,
    ),
    lc=:deepskyblue,
    lw=1.25,
    linestyle="--",
)
GMT.plot!(
    GMT.circgeo(
        [EQlon,EQlat], 
        radius=200, 
        unit=:k,
        np=500,
    ),
    lc=:skyblue,
    lw=1.5,
    linestyle="--",
)
GMT.plot!(
    GMT.circgeo(
        [EQlon,EQlat], 
        radius=400, 
        unit=:k,
        np=500,
    ),
    lc=:lightskyblue,
    lw=2,
    linestyle="--",
)
GMT.scatter!(
    ΔDVVdf[:,:LON],
    ΔDVVdf[:,:LAT],
    zcolor=ΔDVVdf[:,:SLOPEA],
    markeredgecolor=:black,
    marker=:c,
    transparency= 15,
    colorbar=false,
    cmap=C,
    markersize="6p",
)
GMT.colorbar!(
    C=C,
    frame=(annot=:auto, ticks=:auto, xlabel="dv/v [%]"),
    pos=(anchor=:BL,horizontal=true,offset=(-5.5,-2),length=5),
)
GMT.psmeca!(
    [EQlon EQlat 8. 322 81 -173 1 0 0], 
    aki=true, 
    fill=:red, 
    show=true,
    fmt=:png,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/DVV-Ridgecrest-$freqmin-$freqmax.png"),
)