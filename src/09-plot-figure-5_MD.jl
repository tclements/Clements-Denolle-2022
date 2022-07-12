using Arrow ,CSV,  Dates, NetCDF
using DataFrames ,DelimitedFiles, Glob
using DSP , GLM, Interpolations, Statistics, StatsModels
using GMT  
using Plots

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
ind2005e = findfirst(tlwe .== Date(2005,5,16))
ind2005f = findfirst(tlwe .== Date(2005,10,16))
ind2011 = findfirst(tlwe .== Date(2011,9,16))
ind2016 = findfirst(tlwe .== Date(2016,8,21))
ind2020 = findfirst(tlwe .== Date(2020,10,16))
GG=ncread(filename,"lwe_thickness")

# calculate the standard deviation
Gstd=zeros(size(GG[:,:,1]))
for i in 1:length(GG[:,1,1])
    for j in 1:length(GG[1,:,1])
        Gstd[i,j]=std(GG[i,j,:])
    end
end
# standard deviation 
Plots.heatmap(lon,lat,Gstd',clim=(0,20),
xlim=(360-125,360-115),ylim=(32,40))

# standard deviation 
Plots.heatmap(lon,lat,(abs.(Gmax)./Gstd)',clim=(0,2),
xlim=(360-125,360-115),ylim=(32,40))

Gmax=GG[:,:,ind2004]-GG[:,:,ind2005e]

G2004a=sum(GG[:,:,ind2004:ind2005e],dims=3)
G2012a=sum(GG[:,:,ind2011:ind2016],dims=3)
# p=heatmap(lat,lon,G2004a[:,:,1],clim=(-100,100),xlim=(31,40),ylim=(-125+360,-114+360),c =:redsblues)
G2004 = gmtread(filename,layer=ind2004,varname="lwe_thickness")
G2005e = gmtread(filename,layer=ind2005e,varname="lwe_thickness")
G2005f = gmtread(filename,layer=ind2005f,varname="lwe_thickness")
G2011 = gmtread(filename,layer=ind2011,varname="lwe_thickness")
G2016 = gmtread(filename,layer=ind2016,varname="lwe_thickness")
G2020 = gmtread(filename,layer=ind2020,varname="lwe_thickness")


# load fitted values 
fitdf = Arrow.Table(joinpath(@__DIR__,"../data/hydro-model-90-day.arrow")) |> Arrow.columntable |> DataFrame

# load temperature values 
filename2 = joinpath(@__DIR__,"../data/tmean.nc")
pptlon = ncread(filename2,"lon")
pptlat = ncread(filename2,"lat")
tmean = ncread(filename2,"tmean")
ttmean = ncread(filename2,"t")
ttmean = Date.(Dates.unix2datetime.(ttmean))
ttmeanday = (ttmean .- ttmean[1]) ./ Day(1)


# load precip data 
filename = joinpath(@__DIR__,"../data/ppt.nc")
plon = ncread(filename,"lon")
plat = ncread(filename,"lat")
tppt = ncread(filename,"t")
tppt = Date.(Dates.unix2datetime.(tppt))
ppt = ncread(filename,"ppt")
tpptday = (tppt .- tppt[1]) ./ Day(1)
# create time series of yearly precipitation centered on each winter
cummprecip=zeros(size(ppt[:,:,1:21]))
for year in 2001:2022
    println(year)
    ind = findall(Date(year,10,1) .< tppt.< Date(year+1,6,1))
    for i in 1:length(ppt[:,1,1])
        for j in 1:length(ppt[1,:,1])
            crap=ppt[i,j,ind]
            # print(i,j,sum(crap[crap.>0]))
            cummprecip[i,j,year-2000]=sum(crap[crap.>0])
        end
    end
end
# cool plot of rain in 2004-2005 winter
Plots.heatmap(plon,plat,cummprecip[:,:,4]/1000,xlim=(-125,-115),clim=(0,2))

# get variability in precipitation
precipmean=zeros(size(cummprecip[:,:,1]))
crap=0
for i in 1:length(cummprecip[:,1,1])
    for j in 1:length(cummprecip[1,:,1])
        crap=cummprecip[i,j,:]
        precipmean[i,j]=mean(crap)
    end
end
Plots.heatmap(plon,plat,precipmean/1000,xlim=(-125,-115),clim=(0,2))


Plots.heatmap(plon,plat,-log10.(cummprecip[:,:,4]./precipmean),
xlim=(-125,-115),clim=(-0.7,0.7),color=:bluesreds)
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
    ind = findall(Date(2004,10,1) .< DVV[:,:DATE] .< Date(2005,5,1))
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

    println("")
    println("Calculating Δdv/v $netsta $jj")

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
    # ttemp = ttmean[tempind]

    println("got some temp")
    # filter temperature data 
    dvvT = tempfilt[tempind .- round(Int,fitdf[jj,:D5])] 
    dvvT .-= mean(dvvT)
    dvvT ./= std(dvvT)
    dvvT .*= fitdf[fitind,:D4]
    println("fit the temp")

    # get precip data at that station
        # get precip nearest station lat, lon 
    staind = findfirst(CAdf[!,:Station] .== sta)
    if isnothing(staind)
        continue
    end
    lonind = argmin(abs.(plon .- stalon))
    latind = argmin(abs.(plat .- stalat))
    precip = ppt[latind,lonind,:]
    tpminind = findfirst(tppt .== DVV[1,:DATE])
    tpmaxind = findlast(tppt .== DVV[end, :DATE])

    println("got some precip")
    octs = Date(2004,10,1):Year(1):Date(2005,10,1)
    juns = Date(2005,6,1):Year(1):Date(2006,6,1)
    cumprecip = zeros(eltype(precip),length(octs))
    for dd = 1:length(octs)-1
        ind3 = findall((tppt .> octs[dd]) .& (tppt .< juns[dd]))
        if length(ind3)==0
            continue
        else
         cumprecip[dd] = sum(precip[ind])/1000
        end
    end
    println(cumprecip[1])


    # println("")
    # println("")
    # display(plot(DVV[ind,:DATE],dvvT[dvvtempind] ))
    # display(plot!(DVV[ind,:DATE],DVV[ind,:DVV]))
    # remove temperature from dv/v 
    dvvtempind = findall(in(ttmean[tempind]),DVV[ind,:DATE])
    dvvptind = findall(in(tppt),DVV[ind,:DATE])
    
    # display(Plots.plot(DVV[ind,:DATE],dvvT[dvvtempind] ))
    # display(Plots.plot!(DVV[ind,:DATE],DVV[ind,:DVV]))
    DVV[ind,:DVV] .-= dvvT[dvvtempind]  

    # display(Plots.plot!(DVV[ind,:DATE],DVV[ind,:DVV]))
    # display(Plots.plot!(DVV[ind,:DATE],tempfilt[dvvtempind]/10))
    # display(Plots.plot!(DVV[ind,:DATE],precip[dvvptind]))
    # display(Plots.plot!(DVV[ind,:DATE],precip[dvvptind]/10))
    # create linear model 
    data = DataFrame(X=tDVV ./ tDVV[end],Y=DVV[ind,:DVV])
    ols = lm(@formula(Y ~ X), data, wts=DVV[ind,:CC] .^ 3)

    # get relevant data 
    inter, slope = coef(ols)
    rtwo = r2(ols)
    interSTD, slopeSTD = stderror(ols)
    mindate = DVV[ind[1],:DATE]
    maxdate = DVV[ind[end],:DATE]

    ind2=Int(floor(length(ind)/2))
    p2p=-maximum(DVV[ind[1:ind2],:DVV])+minimum(DVV[ind[ind2+1:end],:DVV])
    println("p2p $p2p")
    println("compare this to std in dv/v")
    dvvstd = std(DVV[:,:DVV])
    println(p2p/dvvstd)
    # display(Plots.plot!(DVV[ind,:DATE],(DVV[ind,:DATE]-mindate)./Day(1)./365.0 .*slope.+inter))
    println("std $dvvstd")
    println("Fit between $mindate and $maxdate is $slope and peak-to-peak is $p2p")
    # create dataframe
    df = DataFrame(
        MINDATE=mindate,
        MAXDATE=maxdate,
        NOBS=size(data,1),
        INTER=inter,
        SLOPE=slope,
        P2P=p2p,
        dvvstd=dvvstd,
        R2=rtwo,
        INTERSTD=interSTD,
        SLOPESTD=slopeSTD,
        NETSTA=netsta,
        NET=net,
        STA=sta,
        LON=stalon,
        LAT=stalat,
        cumprecip=cumprecip[1]
    )
    append!(ΔDVVdf,df)
end

# plot with GMT 
minlon = minimum(ΔDVVdf[:,:LON]) - 0.5
maxlon = maximum(ΔDVVdf[:,:LON]) + 0.5
minlat = minimum(ΔDVVdf[:,:LAT]) - 0.5
maxlat = maximum(ΔDVVdf[:,:LAT]) + 0.5


# plot variability in dv/v
crap =Array{Float64}(undef, length(arrowfiles))
slat =Array{Float64}(undef, length(arrowfiles))
slon =Array{Float64}(undef, length(arrowfiles))
for jj in 1:length(arrowfiles)
    # read dv/v for each station 
    DVV = Arrow.Table(arrowfiles[jj]) |> Arrow.columntable |> DataFrame
    crap[jj]=std(DVV[:,:DVV])

    # get station name and location 
    netsta = replace(basename(arrowfiles[jj]),".arrow"=>"")
    net, sta = split(netsta,".")
    staind = findfirst(CAdf[!,:Station] .== sta)
    if isnothing(staind)
        println("No lat/lon for $netsta")
        continue 
    end
    slat[jj] = CAdf[staind,:Latitude]
    slon[jj] = CAdf[staind,:Longitude]
end
Plots.scatter(slon,slat,zcolor=crap,
title="STD(dv/v)",color=:bilbao,markeralpha=1,
colorbar_title="",legend=false,colorbar=true,clim=(0,0.5),
xlim=(-125,-115),ylim=(32,42))
savefig("../data/FINAL-FIGURES/scatter_dvvstd.png")


# plot sensitivity of dv/v to preciptation
Plots.scatter(ΔDVVdf[:,:LON],ΔDVVdf[:,:LAT],zcolor=-ΔDVVdf[:,:P2P]./ΔDVVdf[:,:dvvstd],
title="Peak2Peak / STD ",color=:bilbao,markeralpha=1,
colorbar_title="",legend=false,colorbar=true,clim=(1,5))
savefig("../data/FINAL-FIGURES/scatter_sensitivity_of_dv.png")

# plot peak2peak
Plots.scatter(ΔDVVdf[:,:LON],ΔDVVdf[:,:LAT],zcolor=-2*ΔDVVdf[:,:P2P],
    title="Peak2Peak drop vs total precip",color=:bilbao,markeralpha=2*ΔDVVdf[:,:cumprecip],
    colorbar_title="",legend=false,colorbar=true,clim=(-1,1))
savefig("../data/FINAL-FIGURES/scatter_plot_p2p_precip_winter2005.png")






########################################

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

#####################
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