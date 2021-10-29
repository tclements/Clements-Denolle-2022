using Arrow 
using CSV
using DataFrames
using Dates 
using Glob 
try  
    using GMT
catch
end
using GMT 

NCdf = DataFrame(CSV.File("/home/timclements/CALI/NCstations.csv"))
SCdf = DataFrame(CSV.File("/home/timclements/CALI/CIstations.csv"))
CAdf = vcat(NCdf, SCdf)
CAdf[:,:NETSTA] = CAdf[:,:Network] .* "." .* CAdf[:,:Station]
CAdf = CAdf[findall(in(["NJQ"]),CAdf[:,:Station]),:]


# map boundaries 
minlon, maxlon = [-120.26,-120.0]
minlat, maxlat = [34.45,34.6]


grdimage(
    "@srtm_relief_01s",
    cmap=:dem3, 
    J=:guess,
    shade=true,
    region=(minlon,maxlon,minlat,maxlat),
    coast=true,
    colorbar=false,
)
GMT.coast!(
    shore=true, 
    ocean=:dodgerblue,
    N="a",
)
GMT.text!(
    GMT.mat2ds(
        hcat(
            CAdf[:,:Longitude] .- 0.01,
            CAdf[:,:Latitude] .+ 0.01 ,
        ),
        CAdf[:,:NETSTA],
    ),
    font=(10,)
)
GMT.text!(
    GMT.text_record(
        [-120.22 34.46],
        "Pacific Ocean",
    ),
)
GMT.text!(
    GMT.text_record(
        [-120.15 34.5],
        "Santa Ynez Mountains",
    ),
    angle=15,
)
GMT.basemap!(
    map_scale=(inside=true, anchor=:BC, scale_at_lat="34N", fancy=true,
    length="4k", units=true, offset=0.5),
)
scatter!(
    CAdf[!,:Longitude],
    CAdf[!,:Latitude],
    markeredgecolor=:black,
    marker=:triangle,
    color=:gold,
    transparency= 15,
    markersize=[5.0,5.0],
    show=true,
    fmt=:png,
    savefig="/media/FOUR/data/FINAL-FIGURES/NJQ-map.png",
)



