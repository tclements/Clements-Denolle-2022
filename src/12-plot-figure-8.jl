using ArchGDAL
using Arrow 
using CSV
using Dates
using DataFrames 
using GeoInterface
using Glob 
using GLM 
using GMT
using Interpolations 
using Statistics
using StatsModels
const AG = ArchGDAL

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

#load station locations 
NCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/NCstations.csv")))
SCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/CIstations.csv")))
CAdf = vcat(NCdf, SCdf)

# load PGV values 
PGV = Arrow.Table(joinpath(@__DIR__,"../data/ridgecrest-PGV.arrow")) |> Arrow.columntable |> DataFrame 

# plot maps for 3 frequency octaves 
freqmin = 2.0 .^ (0:2)
freqmax = 2.0 .^ (1:3)

# plot PGV 
for ii in 1:length(freqmin)
    fmin = freqmin[ii]
    fmax = freqmax[ii]

    # subset frequencies 
    PGVFREQ = PGV[PGV[:,:FREQMIN] .== fmin,:]

    # plot with GMT 
    minlon = -122
    maxlon = -113
    minlat = 32
    maxlat = 40
    EQlon = -117.599
    EQlat = 35.778

    # plot PGV
    colorlimit = 20
    C = GMT.makecpt(T=(0,colorlimit,colorlimit / 101), cmap=:oslo,reverse=true);
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
        PGVFREQ[:,:LON],
        PGVFREQ[:,:LAT],
        zcolor=PGVFREQ[:,:PGV],
        markeredgecolor=:black,
        marker=:c,
        transparency= 15,
        colorbar=false,
        cmap=C,
    )
    GMT.colorbar!(
        C=C,
        frame=(annot=:auto, ticks=:auto, xlabel="PGV [cm/s]"),
        pos=(anchor=:BL,horizontal=true,offset=(-5.5,-2),length=5),
    )
    GMT.psmeca!(
        [EQlon EQlat 8. 322 81 -173 1 0 0], 
        aki=true, 
        fill=:red, 
        show=true,
        fmt=:png,
        savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/PGV-Ridgecrest-$fmin-$fmax.png"),
    )
end

# plot PGA
for ii in 1:length(freqmin)
    fmin = freqmin[ii]
    fmax = freqmax[ii]

    # subset frequencies 
    PGVFREQ = PGV[PGV[:,:FREQMIN] .== fmin,:]

    # plot with GMT 
    minlon = -122
    maxlon = -113
    minlat = 32
    maxlat = 40
    EQlon = -117.599
    EQlat = 35.778

    # plot PGA
    colorlimit = 250.
    C = GMT.makecpt(T=(0,colorlimit,colorlimit / 101), cmap=:oslo,reverse=true);
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
        PGVFREQ[:,:LON],
        PGVFREQ[:,:LAT],
        zcolor=PGVFREQ[:,:PGA],
        markeredgecolor=:black,
        marker=:c,
        transparency= 15,
        colorbar=false,
        cmap=C,
    )
    GMT.colorbar!(
        C=C,
        frame=(annot=:auto, ticks=:auto, xlabel="PGA [cm/s^2]"),
        pos=(anchor=:BL,horizontal=true,offset=(-5.5,-2),length=5),
    )
    GMT.psmeca!(
        [EQlon EQlat 8. 322 81 -173 1 0 0], 
        aki=true, 
        fill=:red, 
        show=true,
        fmt=:png,
        savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/PGA-Ridgecrest-$fmin-$fmax.png"),
    )
end


# Baja
maxgap = 10 # maximum data gap [days]
mindays = 75 # minimum number of days needed for analysis 
arrowfiles = glob("*",joinpath(@__DIR__,"../data/FIT-DVV-SSE/90-DAY/"))
ΔDVVdf = DataFrame()

for jj in 1:length(arrowfiles)
    # read dv/v for each station 
    DVV = Arrow.Table(arrowfiles[jj]) |> Arrow.columntable |> DataFrame
    dropmissing!(DVV)
    indafter = findall(Date(2010,4,4) .<= DVV[:,:DATE] .<= Date(2010,4,4) + Day(90))
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


# subset good observations 
# ΔDVVdf[:,:SLOPE] = ΔDVVdf[:,:SLOPEA] .* 10 
# ΔDVVdf = ΔDVVdf[ΔDVVdf[:,:R2A] .> 0.5,:]

# plot with GMT 
# plot with GMT 
minlon = -121
maxlon = -113
minlat = 31
maxlat = 38
EQlon = -115.295
EQlat = 32.286

# plot Δdv/v 
colorlimit = 1.5
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
    [EQlon EQlat 10. 222 59 15 1 0 0], 
    aki=true, 
    fill=:red, 
    show=true,
    fmt=:png,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/DVV-Baja-$freqmin-$freqmax.png"),
)