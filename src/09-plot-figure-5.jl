using Arrow
using CSV
using Dates
using DataFrames
using DelimitedFiles
using DSP
using GLM
using Glob
using Interpolations
using LaTeXStrings
using NetCDF
using Statistics
using StatsModels
import GMT  
import Plots


####### Some functions
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
function smooth_withfiltfilt(A::AbstractArray; window_len::Int=11, window::Symbol=:rect)
    w = getfield(DSP.Windows, window)(window_len)
    A = DSP.filtfilt(w ./ sum(w), A)
    return A
end

######################### station data ################## 
#load station locations 
NCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/NCstations.csv")))
SCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/CIstations.csv")))
CAdf = vcat(NCdf, SCdf)

# lat/lon for LA map 
LAminlon = -119.5
LAmaxlon = -116.5 
LAminlat = 32.75 
LAmaxlat = 35 

# lat/lon for CA
minlon = -125
maxlon = -114
minlat = 32
maxlat = 42

# # plot maps 2-4 Hz 
freqmin = 2.0
freqmax = 4.0


######################### play with temperature data ##################
# load temperature values 
filename2 = joinpath(@__DIR__,"../data/tmean.nc")
pptlon = ncread(filename2,"lon")
pptlat = ncread(filename2,"lat")
tmean = ncread(filename2,"tmean")
ttmean = ncread(filename2,"t")
ttmean = Date.(Dates.unix2datetime.(ttmean))
ttmeanday = (ttmean .- ttmean[1]) ./ Day(1)

######################### play with precipitation data ##################
# load precip data 
filename = joinpath(@__DIR__,"../data/ppt.nc")
plon = ncread(filename,"lon")
plat = ncread(filename,"lat")
tppt = ncread(filename,"t")
tppt = Date.(Dates.unix2datetime.(tppt))
ppt = ncread(filename,"ppt")
ppt[ppt .== -9999.0] .= 0 
tpptday = (tppt .- tppt[1]) ./ Day(1)


# wateryear precipitation 
years = 1986:2020
cummprecip=zeros(size(ppt[:,:,eachindex(years)]))
for year in years
    ind = findall(Date(year - 1,10,1) .<= tppt.< Date(year,10,1))
    cummprecip[:,:,year-first(years) + 1] .= sum(ppt[:,:,ind], dims=3)
end

# average water year precipitation
precipmean = dropdims(mean(cummprecip, dims=3), dims=3)

# mean 2012-2016 data 
precipmean_drought = dropdims(mean(cummprecip[:,:, 2012 .<= years .<= 2016], dims=3), dims=3)
########################## seismic dv/v stuff #############

# load fitted values dvv 
fitdf = Arrow.Table(joinpath(@__DIR__,"../data/hydro-model-90-day.arrow")) |> Arrow.columntable |> DataFrame
arrowfiles = glob("*",joinpath(@__DIR__,"../data/DVV-90-DAY-COMP/$freqmin-$freqmax/"))

########################## Oct 2004 - May 2005  #################
maxgap = 50 # maximum data gap [days]
mindays = 300 # minimum number of days needed for analysis 
ΔDVVdf = DataFrame()
for jj in 1:length(arrowfiles)
    # read dv/v for each station 
    DVV = Arrow.Table(arrowfiles[jj]) |> Arrow.columntable |> DataFrame
    ind = findall(Date(2004,10,1) .<= DVV[:,:DATE] .< Date(2005,10,1))
    # remove station without data
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

    println("\nCalculating Δdv/v $netsta $jj")

    # remove temperature effects
    fitind = findfirst(fitdf[:,:NETSTA] .== netsta)
    lonind = argmin(abs.(pptlon .- fitdf[fitind,:LON]))
    latind = argmin(abs.(pptlat .- fitdf[fitind,:LAT]))
    temp = tmean[latind,lonind,:]
    temp .-= mean(temp)
    if sum(abs.(temp))<1
        continue
    end
    tempfilt = smooth_withfiltfilt(temp,window_len=45)
    tempind = findall(DVV[1,:DATE] .<= ttmean .<= DVV[end,:DATE])

    # filter temperature data 
    dvvT = tempfilt[tempind .- round(Int,fitdf[jj,:D5])] 
    dvvT .-= mean(dvvT)
    dvvT ./= std(dvvT)
    dvvT .*= fitdf[fitind,:D4]

    # get precip data at that station
        # get precip nearest station lat, lon 
    staind = findfirst(CAdf[!,:Station] .== sta)
    if isnothing(staind)
        continue
    end
    lonind = argmin(abs.(plon .- stalon))
    latind = argmin(abs.(plat .- stalat))
    precip =  cummprecip[latind,lonind,findall(years .== 2005)]./precipmean[latind,lonind]

    # remove temperature from dv/v 
    dvvtempind = findall(in(ttmean[tempind]),DVV[ind,:DATE])
    dvvptind = findall(in(tppt),DVV[ind,:DATE])
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
    correction = (maxdate - mindate) / Day(1) / 364
    slope /= correction 

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
        cumprecip=precip
    )
    append!(ΔDVVdf,df)
end


# plot sensitivity of dv/v to preciptation
colorlimit=1.0
CDVV = GMT.makecpt(T=(-colorlimit,colorlimit,colorlimit / 101), cmap=:vik);
preciplimit = 250 
CPRECIP = GMT.makecpt(T=(-preciplimit,preciplimit, preciplimit / 101), cmap=:vik);
land = GMT.makecpt(color=:terra, range=(-5000,5000), no_bg=:true);
precip_2005 = GMT.mat2grid(
    (cummprecip[:,:, findall(years .== 2005)] .- precipmean) ./ precipmean .* 100,
    x = plon, 
    y = plat, 
)


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
    -precip_2005,
    cmap=CPRECIP,
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
    zcolor=ΔDVVdf[:,:SLOPE],
    markeredgecolor=:black,
    marker=:c,
    transparency= 15,
    colorbar=false,
    cmap=CDVV,
    markersize="5p",
)
GMT.colorbar!(
    C=CDVV,
    frame=(annot=:auto, ticks=:auto, xlabel="dv/v [%]"),
    pos=(anchor=:BL,horizontal=true,offset=(-5.5,-2),length=-5),
)
GMT.colorbar!(
    C=CPRECIP,
    frame=(annot=:auto, ticks=:auto, xlabel="@~D@~ Precipitation [%]"), 
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
    [LAminlon LAmaxlat; -119.04 38.98],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed,
    alpha=15,
)
GMT.plot!(
    [LAmaxlon LAmaxlat; -114.8 38.8],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed, 
    alpha=15,
)
# inset map 
tt = GMT.mapproject(
    region=(LAminlon,LAmaxlon,LAminlat,LAmaxlat), 
    proj=:merc, 
    figsize=5, 
    map_size=true,
);
mapW = tt.data[1]
mapH = tt.data[2]
GMT.basemap!(
    inset=(size=(mapW, mapH), 
    anchor=:TR, 
    width=1, 
    offset=(-7.5, -10.25), 
    save="xx000"),
)
t = readdlm("xx000")
GMT.grdimage!(
    "@earth_relief_15s",
    cmap=land, 
    J=:merc,
    shade=true,
    region=(LAminlon,LAmaxlon,LAminlat,LAmaxlat),
    coast=:white,
    colorbar=false,
    shore=:true,
    xshift=t[1], 
    yshift=t[2],
)
GMT.coast!(
    shore=true, 
    ocean=:white,
    N="a",
)
GMT.scatter!(
    ΔDVVdf[:,:LON],
    ΔDVVdf[:,:LAT],
    zcolor=ΔDVVdf[:,:SLOPE],
    markeredgecolor=:black,
    marker=:c,
    transparency= 15,
    # colorbar=true,
    cmap=CDVV,
    markersize="4p",
    show=1,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/DVV-2004-2005-$freqmin-$freqmax.png"),
)

######################################### 2011 - 2016 ###########################
maxgap = 180 # maximum data gap [days]
mindays = 3 * 365 # minimum number of days needed for analysis 
ΔDVVdf = DataFrame()
for jj in 1:length(arrowfiles)
    # read dv/v for each station 
    DVV = Arrow.Table(arrowfiles[jj]) |> Arrow.columntable |> DataFrame
    ind = findall(Date(2011,10,1) .<= DVV[:,:DATE] .< Date(2016,10,1))
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

    println("Calculating Δdv/v $netsta $jj")

    # remove temperature effects
    fitind = findfirst(fitdf[:,:NETSTA] .== netsta)
    lonind = argmin(abs.(pptlon .- fitdf[fitind,:LON]))
    latind = argmin(abs.(pptlat .- fitdf[fitind,:LAT]))
    temp = tmean[latind,lonind,:]
    if maximum(temp)<1
        continue
    end
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
    inter, slope = coef(ols)
    rtwo = r2(ols)
    interSTD, slopeSTD = stderror(ols)
    mindate = DVV[ind[1],:DATE]
    maxdate = DVV[ind[end],:DATE]

    # correct for data coverage / annualized data 
    correction = (maxdate - mindate) / Day(1) / 1826
    slope /= correction 
    slope /= 5 

    # create dataframe
    df = DataFrame(
        MINDATE=mindate,
        MAXDATE=maxdate,
        NOBS=size(data,1),
        INTER=inter,
        SLOPE= slope,
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

# plot topography + precipitation 
colorlimit=0.15
CDVV = GMT.makecpt(T=(-colorlimit,colorlimit,colorlimit / 101), cmap=:vik);
preciplimit = 50 
CPRECIP = GMT.makecpt(T=(-preciplimit,preciplimit,preciplimit / 101), cmap=:vik);
precip_2012_2016 = GMT.mat2grid(
    (precipmean_drought .- precipmean) ./ precipmean .* 100,
    x = plon, 
    y = plat, 
)


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
    -precip_2012_2016,
    cmap=CPRECIP,
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
    zcolor=ΔDVVdf[:,:SLOPE],
    markeredgecolor=:black,
    marker=:c,
    transparency= 15,
    colorbar=false,
    cmap=CDVV,
    markersize="5p",
)
GMT.colorbar!(
    C=CDVV,
    frame=(annot=:auto, ticks=:auto, xlabel="dv/v [%]"),
    pos=(anchor=:BL,horizontal=true,offset=(-5.5,-2),length=-5),
)
GMT.colorbar!(
    C=CPRECIP,
    frame=(annot=:auto, ticks=:auto, xlabel="@~D@~ Precipitation [%]"), 
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
    [LAminlon LAmaxlat; -119.04 38.98],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed,
    alpha=15,
)
GMT.plot!(
    [LAmaxlon LAmaxlat; -114.8 38.8],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed, 
    alpha=15,
)
# inset map 
tt = GMT.mapproject(
    region=(LAminlon,LAmaxlon,LAminlat,LAmaxlat), 
    proj=:merc, 
    figsize=5, 
    map_size=true,
);
mapW = tt.data[1]
mapH = tt.data[2]
GMT.basemap!(
    inset=(size=(mapW, mapH), 
    anchor=:TR, 
    width=1, 
    offset=(-7.5, -10.25), 
    save="xx000"),
)
t = readdlm("xx000")
GMT.grdimage!(
    "@earth_relief_15s",
    cmap=land, 
    J=:merc,
    shade=true,
    region=(LAminlon,LAmaxlon,LAminlat,LAmaxlat),
    coast=:white,
    colorbar=false,
    shore=:true,
    xshift=t[1], 
    yshift=t[2],
)
GMT.coast!(
    shore=true, 
    ocean=:white,
    N="a",
)
GMT.scatter!(
    ΔDVVdf[:,:LON],
    ΔDVVdf[:,:LAT],
    zcolor=ΔDVVdf[:,:SLOPE],
    markeredgecolor=:black,
    marker=:c,
    transparency= 15,
    # colorbar=true,
    cmap=CDVV,
    markersize="4p",
    show=1,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/DVV-2011-2016-$freqmin-$freqmax.png"),
)

############################ LA Area  # 2005 - 2020 ###################################

maxgap = 2*365 # maximum data gap [days]
mindays = 13 * 365 # minimum number of days needed for analysis 
# arrowfiles = glob("*",joinpath(@__DIR__,"../data/DVV-90-DAY-COMP/$freqmin-$freqmax/"))
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