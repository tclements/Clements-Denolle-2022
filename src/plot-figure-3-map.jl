using CSV
using DataFrames
using Dates  
using DelimitedFiles
using NetCDF
using Shapefile
using Statistics
using GMT

# get CIstation location 
cistations = DataFrame(CSV.File("/home/timclements/CALI/CIstations.csv"))
LJRind = findfirst(cistations[!,:Station] .== "LJR")
LJRlat = cistations[LJRind,:Latitude]
LJRlon = cistations[LJRind,:Longitude]
minlon, maxlon = LJRlon - 0.2, LJRlon + 0.2
minlat, maxlat = LJRlat - 0.2, LJRlat + 0.2 
CAminlon, CAmaxlon = [-125,-113.5]
CAminlat, CAmaxlat = [32,43]

# GWL location 
GWLlon = [-118.866718, -118.8161, -118.820786, -118.8624] 
GWLlat = [34.824446, 34.850942, 34.848739, 34.83199] 

# get precip location 
filename = "/media/FOUR/data/ppt.nc"
pptlon = ncread(filename,"lon")
pptlat = ncread(filename,"lat")
pptlonind = argmin(abs.(pptlon .- LJRlon))
pptlatind = argmin(abs.(pptlat .- LJRlat))
pptlonminus = mean(pptlon[pptlonind-1:pptlonind])
pptlonplus =  mean(pptlon[pptlonind:pptlonind+1])
pptlatminus = mean(pptlat[pptlatind-1:pptlatind])
pptlatplus =  mean(pptlat[pptlatind:pptlatind+1])
pptrect = [pptlonminus pptlatminus;
        pptlonminus pptlatplus;
        pptlonplus pptlatplus;
        pptlonplus pptlatminus;
        pptlonminus pptlatminus]

# load GRACE data 
# from http://www2.csr.utexas.edu/grace/RL06_mascons.html
filename = "/media/FOUR/data/CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02.nc"
GRACElon = ncread(filename,"lon")
GRACElat = ncread(filename,"lat")
GRACElon[GRACElon .> 180] .-= 360
GRACElonind = argmin(abs.(GRACElon .- LJRlon))
GRACElatind = argmin(abs.(GRACElat .- LJRlat))
GRACElonminus = mean(GRACElon[GRACElonind-1:GRACElonind])
GRACElonplus =  mean(GRACElon[GRACElonind:GRACElonind+1])
GRACElatminus = mean(GRACElat[GRACElatind-1:GRACElatind])
GRACElatplus =  mean(GRACElat[GRACElatind:GRACElatind+1])
GRACErect = [GRACElonminus GRACElatminus;
        GRACElonminus GRACElatplus;
        GRACElonplus GRACElatplus;
        GRACElonplus GRACElatminus;
        GRACElonminus GRACElatminus]

# get location of I5 
table = Shapefile.Table("/media/FOUR/data/interstate/intrstat.shp")
I5 = Shapefile.shapes(table)[findall(table.ROUTE_NUM .== "I5")[1]]
I5lon = [point.x for point in I5.points]
I5lat = [point.y for point in I5.points]
I5lonind = findall((I5lon .> minlon - 0.01) .& (I5lon .< maxlon + 0.01))
I5latind = findall((I5lat .> minlat - 0.01) .& (I5lat .< maxlat + 0.01))
I5ind = intersect(I5lonind,I5latind)
I5lon = I5lon[I5ind]
I5lat = I5lat[I5ind]

# GPS station location 
WNAMloc = DataFrame(
    CSV.File(
        "/media/FOUR/data/WNAM-LOC.tsv",
        delim=' ',
        ignorerepeated=true,
    )
)
LJRNind = findfirst(WNAMloc[!,:STATION] .== "LJRN")
LJRNlon = WNAMloc[LJRNind,:LON]
LJRNlat = WNAMloc[LJRNind,:LAT]

# plot station map first 
topo = makecpt(color=:geo, range=(0,5000,500), continuous=true);
grdimage(
    "@srtm_relief_01s", 
    J=:guess,
    shade=true,
    region=(minlon,maxlon,minlat,maxlat),
    coast=false,
    colorbar=false,
    color=topo,
    conf=(MAP_FRAME_TYPE="fancy+", MAP_GRID_PEN_PRIMARY="thinnest,black,.",FONT_ANNOT_PRIMARY="+16"),
    dpi=250,
)
plot!(
    I5lon,
    I5lat,
    linecolor=:black,
    lw=2,
)
plot!(
    pptrect,
    lw=2,
    linecolor=:chartreuse,
)
plot!(
    GRACErect,
    lw=2,
    linecolor=:gold,
    linestyle="--",
)
# scatter!(
#     [LJRNlon],
#     [LJRNlat],
#     marker=:pentagon, 
#     markerfacecolor=:purple,
#     markeredgecolor=:black,
#     markersize=0.5,
#     ml=1,
# )
scatter!(
    [LJRlon],
    [LJRlat],
    marker=:circle, 
    markersize="3.5k",
    ml=(1.25,:dodgerblue), 
    linestyle="--",
)
scatter!(
    [LJRlon],
    [LJRlat],
    marker=:circle, 
    markersize="7.5k",
    ml=(1.25,:dodgerblue), 
    linestyle="--",
)

scatter!(
    [LJRlon],
    [LJRlat],
    mc=:gold,
    markeredgecolor=:black,
    marker="t",
    markersize=0.5,
    # show=true,
    # fmt=:png,
    # savefig="/media/FOUR/data/FINAL-FIGURES/CI-station-map.png",
)
scatter!(
    GWLlon,
    GWLlat,
    mc=:lightskyblue,
    markeredgecolor=:black,
    marker="circle",
    markersize=0.3,
    # show=true,
    # fmt=:png,
    # savefig="/media/FOUR/data/FINAL-FIGURES/CI-station-map.png",
)
basemap!(inset=(anchor=:BL, width=0.1, offset=(0.5, 0.5), save="xx000"))
t = readdlm("xx000")
coast!(region=[CAminlon CAmaxlon CAminlat CAmaxlat], proj=:merc,
           land=:lightgray, area=5000, shore=:faint, ocean=:dodgerblue,
           x_off=t[1], y_off=t[2],
           N="a",
           figsize=4,
           frame=:bare,
           DCW=(state="CA", fill=:bisque),
)
scatter!(
     [LJRlon],
     [LJRlat],
     marker=:triangle,
     markeredgecolor=:black,
     fill=:gold,
     markersize=0.5,
     show=true,
     fmt=:png,
     savefig=expanduser("/media/FOUR/data/FINAL-FIGURES/LJR-map.png"),
)
