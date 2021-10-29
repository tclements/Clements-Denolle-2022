using ArchGDAL
using Arrow
using DataFrames
using Dates  
using DelimitedFiles
using Glob

import GMT 
import Plots 

# lat/lon for LA map 
LAminlon = -118.75 
LAmaxlon = -117.5 
LAminlat = 33.6 
LAmaxlat = 34.52 

# load dv/v 
fitdf = Arrow.Table("/media/FOUR/data/hydro-model-90-day.arrow") |> Arrow.columntable |> DataFrame
minlon = floor(minimum(fitdf[:,:LON]))
maxlon = ceil(maximum(fitdf[:,:LON]))
minlat = floor(minimum(fitdf[:,:LAT]))
maxlat = ceil(maximum(fitdf[:,:LAT])) - 0.75 

# add VS30 values to PGV 
# load VS30 data from USGS https://earthquake.usgs.gov/data/vs30/
ds = ArchGDAL.readraster("/media/FOUR/data/VS30/global_vs30.tif")
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

# find nearest PGV values 
VS30 = zeros(size(fitdf,1))
for ii in 1:length(VS30)
    jj = argmin(abs.(fitdf[ii,:LON] .- vs30lon))
    kk = argmin(abs.(fitdf[ii,:LAT] .- vs30lat))
    VS30[ii] = dataset[jj,kk,1]
end


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
    savefig="/media/FOUR/data/FINAL-FIGURES/CI-station-map.png",
)

# show fitted dv/v and diffusion coefficient 
Sind = findall((fitdf[:,:corSSW] .> 0.5) .& (fitdf[:,:SSW3] .< 1.))
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
    zcolor=fitdf[Sind,:SSW3],
    markeredgecolor=:black,
    title="SSW Model Decay Parameter",
    colorbar=false,
    cmap=hotmap,
    show=false,
)
GMT.colorbar!(cmap=hotmap,log=true,show=true,savefig="/media/FOUR/data/FINAL-FIGURES/SSW-Decay.png",)

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
GMT.colorbar!(cmap=hotmap,log=true,show=true,savefig="/media/FOUR/data/FINAL-FIGURES/Elastic-Diffusivity.png",)

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
    savefig="/media/FOUR/data/FINAL-FIGURES/Elastic-Temp.png",
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
    savefig="/media/FOUR/data/FINAL-FIGURES/SSW-Temp.png",
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
    markersize="7p",
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
    [LAminlon LAmaxlat; -119.05 39.0],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed,
    alpha=15,
)
GMT.plot!(
    [LAmaxlon LAmaxlat; -114.75 39.0],
    region=(minlon,maxlon,minlat,maxlat), 
    lw=1,
    ls=:dashed, 
    alpha=15,
)
# inset map 
GMT.basemap!(inset=(anchor=:TR, width=2, offset=(4.5, 3), save="xx000"))
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
    cmap=:terra, 
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
    markersize="7p",
    savefig="/media/FOUR/data/FINAL-FIGURES/Hydro-Thermal-Mixing-map.png",
)
