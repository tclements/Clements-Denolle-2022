using Arrow 
using CSV
using DataFrames
using Dates 
using Glob 
using GMT 

NCdf = DataFrame(CSV.File("/home/timclements/CALI/NCstations.csv"))
SCdf = DataFrame(CSV.File("/home/timclements/CALI/CIstations.csv"))
CAdf = vcat(NCdf, SCdf)
CAdf[:,:NETSTA] = CAdf[:,:Network] .* "." .* CAdf[:,:Station]

# grab number of days per netsta 
dvvfiles = glob("*","/media/FOUR/data/DVV-90-DAY-COMP/2.0-4.0/")
N = length(dvvfiles)
DVVdf = DataFrame()
netstas = Array{String}(undef,N)
days = zeros(Int,N)
for ii in 1:N
    dvvarrow = Arrow.Table(dvvfiles[ii]) |> DataFrame
    netstas[ii] = replace(basename(dvvfiles[ii]),".arrow"=>"")
    days[ii] = round(Int,size(dvvarrow,1) / 365.25)
end
DVVdf = DataFrame(NETSTA = netstas, DAYS = days)
DAYdf = innerjoin(CAdf, DVVdf, on=:NETSTA)

minlon, maxlon = [-125,-113.5]
minlat, maxlat = [31.5,42.5]

# plot station map first
C = GMT.makecpt(T=(0,maximum(DAYdf[:,:DAYS]) + 4,0.25), cmap=:lajolla);
grdimage(
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
scatter!(
    DAYdf[!,:Longitude],
    DAYdf[!,:Latitude],
    cmap=C,
    markeredgecolor=:black,
    marker=:triangle,
    zcolor=DAYdf[:,:DAYS],
    colorbar=false,
)
GMT.colorbar!(
    C=C,
    frame=(annot=:auto, ticks=:auto, xlabel="Years of data"),
    pos=(anchor=:BL,horizontal=true,offset=(-7.5,-1.75),length=6),
    show=true,
    fmt=:png,
    savefig="/media/FOUR/data/FINAL-FIGURES/CA-station-map.png",
)
    