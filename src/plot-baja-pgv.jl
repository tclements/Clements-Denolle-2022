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
const AG = ArchGDAL

#load station locations 
NCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/NCstations.csv")))
SCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/CIstations.csv")))
CAdf = vcat(NCdf, SCdf)

# load PGV values 
PGV = Arrow.Table(joinpath(@__DIR__,"../data/baja-PGV.arrow")) |> Arrow.columntable |> DataFrame 

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
    minlon = -121
    maxlon = -113
    minlat = 31
    maxlat = 38
    EQlon = -115.295
    EQlat = 32.286

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
        [EQlon EQlat 10. 222 59 15 1 0 0], 
        aki=true, 
        fill=:red, 
        show=true,
        fmt=:png,
        savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/PGV-Baja-$fmin-$fmax.png"),
    )
end

# plot PGA
for ii in 1:length(freqmin)
    fmin = freqmin[ii]
    fmax = freqmax[ii]

    # subset frequencies 
    PGVFREQ = PGV[PGV[:,:FREQMIN] .== fmin,:]

    # plot with GMT 
    minlon = -121
    maxlon = -113
    minlat = 31
    maxlat = 38
    EQlon = -115.295
    EQlat = 32.286

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
        [EQlon EQlat 10. 222 59 15 1 0 0], 
        aki=true, 
        fill=:red, 
        show=true,
        fmt=:png,
        savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/PGA-Baja-$fmin-$fmax.png"),
    )
end
