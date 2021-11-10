using Arrow 
using CSV
using DataFrames
using Dates 
using Glob 
using Images
try  
    using GMT
catch
end
using GMT 

NCdf = DataFrame(CSV.File("/home/timclements/CALI/NCstations.csv"))
SCdf = DataFrame(CSV.File("/home/timclements/CALI/CIstations.csv"))
CAdf = vcat(NCdf, SCdf)
CAdf[:,:NETSTA] = CAdf[:,:Network] .* "." .* CAdf[:,:Station]
CAdf = CAdf[findall(in(["RXH","SAL"]),CAdf[:,:Station]),:]

# map boundaries 
minlon, maxlon = [-116.33,-115.2]
minlat, maxlat = [33.0,33.66]

# Obsidian Butte Swarm 
obEQlat = 33.153
obEQlon = -115.646

# 2014 M4.2 EQ
M42EQlat = 33.185
M42EQlon = -115.609

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
            CAdf[:,:Longitude] .+ 0.05,
            CAdf[:,:Latitude] .- 0.025 ,
        ),
        CAdf[:,:NETSTA],
    ),
    font=(10,)
)
GMT.text!(
    GMT.text_record(
        [-115.85 33.35],
        "Salton Sea",
    ),
    angle=-45,
)
# plot focal mechanisms with lines between them 
offset = 0.05
GMT.plot!(
    [obEQlon obEQlat; obEQlon obEQlat - offset ],
    lw=0.5,
    lc=:red,
)
GMT.psmeca!(
        [obEQlon obEQlat - offset  5.0 332 80 178 1.0 0 0], 
        aki=true, 
        fill=:red, 
)
GMT.plot!(
    [M42EQlon M42EQlat; M42EQlon + offset M42EQlat + offset ],
    lw=0.5,
    lc=:red,
)
GMT.psmeca!(
        [M42EQlon + offset M42EQlat + offset 5.0 204 30 -25 0.75 0 0], 
        aki=true, 
        fill=:red, 
)
GMT.basemap!(
    map_scale=(inside=true, anchor=:RB, scale_at_lat="33N", fancy=true,
    length="10k", units=true, offset=1.0),
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
    savefig="/media/FOUR/data/FINAL-FIGURES/SALTON-map.png",
)


# plot salton mapind
img1 = load("/media/FOUR/data/FINAL-FIGURES/SALTON-map.png")
p1 = plot(img1,border=:none)
annotate!((-160,100,text("A",20)))

# plot dv/v at salton city 
img2 = load("/media/FOUR/data/FINAL-FIGURES/CISAL-DVV-EQ.png")
p2 = plot(img2,border=:none)
annotate!((157,100,text("B",20)))

# plot dv/v at CI.RXH 
img3 = load("/media/FOUR/data/FINAL-FIGURES/CIRXH-DVV.png")
p3 = plot(img3,border=:none)
annotate!((-100,100,text("C",20)))
# l1 = @layout  [a [b{0.95w,0.5h}; c{0.9w, 0.5h}]]
l1 = @layout  [a{0.36h}; b{0.34h}; c{0.3h}]
plot(p1,p2,p3,layout=l1,size=(1600,1600),dpi=250)
savefig("/media/FOUR/data/FINAL-FIGURES/SALTON-SAL-RXH.png")
