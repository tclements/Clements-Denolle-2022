using ArchGDAL
using Arrow
using DataFrames
using Dates  
using DelimitedFiles
using Glob
using Statistics
import GMT 
import Plots 

# lat/lon for LA map 
LAminlon = -119.0
LAmaxlon = -115.5 
LAminlat = 32.75 
LAmaxlat = 34.75 

# load dv/v 
fitdf = Arrow.Table(joinpath(@__DIR__,"../data/hydro-model-90-day_L12.arrow")) |> Arrow.columntable |> DataFrame
minlon = floor(minimum(fitdf[:,:LON]))
maxlon = ceil(maximum(fitdf[:,:LON]))
minlat = floor(minimum(fitdf[:,:LAT]))
maxlat = ceil(maximum(fitdf[:,:LAT])) - 0.75 

# add VS30 values to PGV 
# load VS30 data from USGS https://earthquake.usgs.gov/data/vs30/
ds = ArchGDAL.readraster(joinpath(@__DIR__,"../data/VS30/global_vs30.tif"))
geotransform = ArchGDAL.getgeotransform(ds)
vs30lon = range(geotransform[1],step=geotransform[2],length=ds.size[1])
vs30lat = range(geotransform[4],step=geotransform[6],length=ds.size[2])

# subset to area with PGV measurements 
lonind = findall(minlon - 0.5 .< vs30lon .< maxlon + 0.5)
latind = findall(minlat - 0.5 .< vs30lat .< maxlat + 0.5)
vs30lon = vs30lon[lonind]
vs30lat = vs30lat[latind]
dataset = ds[lonind,latind,1]

# set oceans values to zero 
dataset[dataset .== 600.0] .= 0.0

# find nearest VS30 values 
VS30 = zeros(size(fitdf,1))
for ii in 1:length(VS30)
    jj = argmin(abs.(fitdf[ii,:LON] .- vs30lon))
    kk = argmin(abs.(fitdf[ii,:LAT] .- vs30lat))
    VS30[ii] = dataset[jj,kk,1]
end

# find which model is beats
label = ["drained","elastic","fully coupled","cdm","ssw"]
AAmax = zeros(length(VS30))
for i in 1:length(fitdf[:,:LAT])
    AA = [fitdf[i,:r2DL1], fitdf[i,:r2EL1],fitdf[i,:r2FCL1],fitdf[i,:r2CDML1],fitdf[i,:r2SSWL1]]
    AAmax[i]=argmax(AA)
    println(AAmax[i])
end

# number of stations that have the drained model as the best model:
sum(AAmax.==1.0)

# first ignore the stations where the diffusivity was found to reach the bounds.
Xlat = convert(Vector,fitdf.LAT)
Xlon = convert(Vector,fitdf.LON)
Xr2 = convert(Vector,fitdf.r2DL1)
XC = convert(Vector,fitdf.D3)
XT = convert(Vector,fitdf.D5)
Sind = findall(fitdf[:,:D3] .< 10.)
Sind2 = findall(fitdf[:,:D3] .> 10.)

# plot diffusivity
Plots.scatter(Xlon[Sind],Xlat[Sind],zcolor=log10.(XC[Sind]),
    title="Diffusivity in Drained model",markeralpha=fitdf[Sind,:r2DL1],
    colorbar_title="(log_{10} (m^2/s)",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_plot_diffusivity.png")



# plot diffusivity
Plots.scatter(Xlon[Sind2],Xlat[Sind2],zcolor=log10.(fitdf[Sind2,:E3]),
    title="Diffusivity in Drained model",markeralpha=fitdf[Sind2,:r2EL1],
    colorbar_title="(log_{10} (m^2/s)",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_plot_diffusivity_E.png")


# Are the diffusivity values reasonable?
#seismic wavelength ^2 / day 
# Plots.scatter(Xlon[Sind],Xlat[Sind],zcolor=((1.5.*VS30[Sind]./3).^2/(3600*24)),
#     title="\lambda^2/day",color=:jet,
#     colorbar_title="(days)",legend=false,colorbar=true)
Plots.histogram(sqrt.(fitdf[Sind,:D3].*3600*24))
Plots.scatter(Xlon[Sind],Xlat[Sind],zcolor=log10.(sqrt.(fitdf[Sind,:D3].*3600*24)),
    title="Diffusion distance per day",markeralpha=fitdf[Sind,:r2DL1],
    colorbar_title="(log_{10} (m)",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_plot_diffusion_distance_per_day.png")



# plot temp delay
Plots.scatter(Xlon[Sind2],Xlat[Sind2],zcolor=log10.(fitdf[Sind2,:E3]),
    title="Diffusivity for elastic model",color=:bilbao,markeralpha=fitdf[Sind2,:r2EL1],
    colorbar_title="(days)",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_diff_E.png")

# plot temp delay
Plots.scatter(Xlon[Sind],Xlat[Sind],zcolor=(XT[Sind]),
    title="Phase delay in temp model",color=:jet,
    colorbar_title="(days)",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_plot_temp_delay.png")


# plot R2
Plots.scatter(Xlon[Sind],Xlat[Sind],zcolor=(fitdf[Sind,:r2DL1]),
title="R2 for drained model",color=:bilbao,
colorbar_title="(R2)",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_r2_drained.png")


Plots.scatter(Xlon[Sind],Xlat[Sind],zcolor=(fitdf[Sind,:r2EL1]),
title="R2 for elastic model",color=:bilbao,
colorbar_title="(R2)",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_r2_elastic.png")


Plots.scatter(Xlon[Sind],Xlat[Sind],zcolor=(fitdf[Sind,:r2FCL1]),
title="R2 for FC model",color=:bilbao,
colorbar_title="(R2)",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_r2_fc.png")


Plots.scatter(Xlon[Sind],Xlat[Sind],zcolor=(fitdf[Sind,:r2CDML1]),
title="R2 for CDM model",color=:bilbao,
colorbar_title="(R2)",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_r2_cdm.png")


Plots.scatter(Xlon[Sind],Xlat[Sind],zcolor=(fitdf[Sind,:r2SSWL1]),
title="R2 for SSW model",color=:bilbao,
colorbar_title="(R2)",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_r2_ssw.png")


# 


# plot temp vs hydro
tempcomp = abs.(fitdf[:,:D4]) ./ (abs.(fitdf[:,:D4]) .+  abs.(fitdf[:,:D2]))
tt=(tempcomp[Sind] .- 0.5) .*2
Plots.scatter(Xlon[Sind],Xlat[Sind],zcolor=tt,
    title="Temp vs hydro",color=:bluesreds,markeralpha=fitdf[Sind,:r2EL1],
    colorbar_title="",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_plot_temp_vs_hydro.png")



##################### USING GMT ################
# plot station map first 
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
    shore=true, 
    ocean=:white,
    N="a",
)
GMT.scatter!(
    fitdf[!,:LON],
    fitdf[!,:LAT],
    mc=:yellow,
    markeredgecolor=:black,
    show=true,
    fmt=:png,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/CI-station-map.png"),
)

# show fitted dv/v and diffusion coefficient 
Sind = findall((fitdf[:,:r2D3L1] .> 0.5) .& (fitdf[:,:D3] .< 1.))
hotmap = GMT.makecpt(range=(-4.0,0.0,20,:number),log=true,cmap=:hot)
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
    shore=true, 
    ocean=:white,
    N="a",
)
GMT.scatter!(
    fitdf[Sind,:LON],
    fitdf[Sind,:LAT],
    zcolor=fitdf[Sind,:D3],
    markeredgecolor=:black,
    title="SSW Model Decay Parameter",
    colorbar=false,
    cmap=hotmap,
    show=false,
)
GMT.colorbar!(cmap=hotmap,log=true,show=true,savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/SSW-Decay.png"),
)

# show fitted dv/v and diffusion coefficient 
Eind = findall((fitdf[:,:corE] .> 0.5) .& (fitdf[:,:E3] .< 10.))
hotmap = GMT.makecpt(range=(-4.0,-1.0,20,:number),log=true,cmap=:hot)
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
    shore=true, 
    ocean=:white,
    N="a",
)
GMT.scatter!(
    fitdf[Eind,:LON],
    fitdf[Eind,:LAT],
    zcolor=fitdf[Eind,:E3],
    markeredgecolor=:black,
    title="Elastic Model Hydraulic Diffusivity",
    colorbar=false,
    cmap=hotmap,
    show=false,
)
GMT.colorbar!(
    cmap=hotmap,
    log=true,
    show=true,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/Elastic-Diffusivity.png"),
)

# plot temperature delay 
tempmap = GMT.makecpt(range=(0,180,20),cmap=:hot)
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
    shore=true, 
    ocean=:white,
    N="a",
)
GMT.scatter!(
    fitdf[Eind,:LON],
    fitdf[Eind,:LAT],
    zcolor=fitdf[Eind,:E5],
    markeredgecolor=:black,
    title="Elastic Model Temperature Delay [Days]",
    colorbar=true,
    cmap=:hot,
    show=true,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/Elastic-Temp.png"),
)

# plot SSW temp delay 
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
    shore=true, 
    ocean=:white,
    N="a",
)
GMT.scatter!(
    fitdf[Sind,:LON],
    fitdf[Sind,:LAT],
    zcolor=fitdf[Sind,:SSW5],
    markeredgecolor=:black,
    title="SSW Model Temperature Delay [Days]",
    colorbar=true,
    cmap=:hot,
    show=true,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/SSW-Temp.png"),
)

# scatter ratio of hydro to temperature 
tempcomp = abs.(fitdf[:,:E4]) ./ (abs.(fitdf[:,:E4]) .+  abs.(fitdf[:,:E2]))
Eind = findall((fitdf[:,:E2] .< 0.0) .& (fitdf[:,:E4] .> 0.))
C = GMT.makecpt(T=(-100.0,100.0,1.0), cmap=:polar)
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
    shore=true, 
    region=(minlon,maxlon,minlat,maxlat),
    ocean=:white,
    N="a",
)
GMT.scatter!(
    fitdf[Eind,:LON],
    fitdf[Eind,:LAT],
    markeredgecolor=:black,
    zcolor=(tempcomp[Eind] .- 0.5) .* 200,
    colormap=C,
    alpha=5,
    markersize="6p",
)
rect = [LAminlon LAminlat; LAminlon LAmaxlat; LAmaxlon LAmaxlat; LAmaxlon LAminlat; LAminlon LAminlat];
GMT.plot!(
    rect, 
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed
)
GMT.colorbar!(
    C=C,
    frame=(annot=:auto, ticks=:auto, xlabel="dv/v Hydro    <--->   dv/v Thermal", suffix=" %"),
    pos=(anchor=:BL,horizontal=true,offset=(-9,-2),length=8),
)
GMT.plot!(
    [LAminlon LAmaxlat; -119.75 38.9],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed,
    alpha=15,
)
GMT.plot!(
    [LAmaxlon LAmaxlat; -114.57 38.88],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed, 
    alpha=15,
)
# inset map 
GMT.basemap!(inset=(anchor=:TR, width=2, offset=(5.33, 3.1), save="xx000"))
t = readdlm("xx000")
# GMT.coast!(region=[LAminlon LAmaxlon LAminlat LAmaxlat], proj=:merc,
#            land=:lightgray, area=5000, shore=:faint, ocean=:dodgerblue,
#            x_offset=t[1], y_offset=t[2],
#            N="a",
#            frame=:bare,
#            DCW=(state="CA", fill=:bisque),
# )
GMT.grdimage!(
    "@srtm_relief_03s",
    cmap=:elevation, 
    J=:merc,
    shade=true,
    region=(LAminlon,LAmaxlon,LAminlat,LAmaxlat),
    coast=false,
    colorbar=false,
    x_offset=t[1], y_offset=t[2],
)
GMT.coast!(
    shore=true, 
    ocean=:white,
    N="a",
)
GMT.scatter!(
    fitdf[Eind,:LON],
    fitdf[Eind,:LAT],
    markeredgecolor=:black,
    zcolor=(tempcomp[Eind] .- 0.5) .* 200,
    colormap=C,
    alpha=0,
    show=true,
    markersize="4p",
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/Hydro-Thermal-Mixing-map.png"),
)
