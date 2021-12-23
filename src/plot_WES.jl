using Arrow 
using CSV
using DataFrames
using Dates
using DSP 
using HDF5
using JSON
using LaTeXStrings
using LinearAlgebra
using Optim 
using Plots 
using Plots.Measures
using Statistics 
try
    import GMT
catch
end
import GMT

function smooth_withfiltfilt(A::AbstractArray; window_len::Int=11, window::Symbol=:rect)
    w = getfield(DSP.Windows, window)(window_len)
    A = DSP.filtfilt(w ./ sum(w), A)
    return A
end

# read data from CI.WES 
WES = Arrow.Table(joinpath(@__DIR__,"../data/FIT-DVV-SSE/90-DAY/CI.WES.arrow")) |> DataFrame

# get WES lat, lon from CI station locations 
SCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/CIstations.csv")))
NCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/NCstations.csv")))
CAdf = vcat(NCdf, SCdf)
CAdf = CAdf[CAdf[:,:StartTime] .< Date(2010,4,4),:]
CAdf = CAdf[CAdf[:,:EndTime] .> Date(2010,4,4),:]
CAdf = CAdf[CAdf[:,:Longitude] .>= -116.108,:]
# remove SGL station - only has EHZ
CAdf = CAdf[CAdf[:,:Station] .!= "SGL",:]
WESind = findfirst(CAdf[!,:Station] .== "WES")
WESlat = CAdf[WESind,:Latitude]
WESlon = CAdf[WESind,:Longitude]

# GPS location 
p494lon = -115.732065
p494lat = 32.759655

# M 7.2 Baja California Sierra El Mayor earthquake
EQlon = -115.295
EQlat = 32.286

# load EQ locations 
EQ = DataFrame(
    CSV.File(
    joinpath(@__DIR__,"../data/baja.txt"),
    delim=" ",
    ignorerepeated=true,
    )
)
EQ = EQ[EQ[:,:MAG] .> 3.0,:]

# find groundwater wells nearby 
GWL = CSV.File(joinpath(@__DIR__,"../data/324603115480501.tsv"),comment="#",skipto=3) |> DataFrame
dropmissing!(GWL,:lev_va) # lev_va is ft below surface
GWL[!,:lev_va] ./= 3.28 # convert to meters 
GWL = GWL[GWL[!,:lev_dt] .>= Date(2001),:]
GWL = GWL[GWL[:,:lev_va ] .> 14.9,:]
GWL[!,:lev_va] .*= -1

# read GPS data for station P494 
GPS = CSV.File(
    joinpath(@__DIR__,"../data/WNAMdetrend/p494FilterDetrend.neu"), 
    comment="#",
    header=[
        "DATE",
        "YEAR",
        "DAY",
        "E",
        "N",
        "V",
        "Esig",
        "Nsig",
        "Vsig",
        "ENcor",
        "NVcor",
        "EVcor"
    ],
    ignorerepeated=true,
     delim=" "
) |> DataFrame
GPS[!,:DATE] = Date.(GPS[!,:YEAR]) .+ Day.(GPS[!,:DAY] .- 1)
# select after Baja earthquake
GPS = GPS[(GPS[!,:DATE] .>= Date(2010,4,5)) .& (GPS[!,:DATE] .<= WES[end,:DATE]),:]
smoothEN = smooth_withfiltfilt(sqrt.( (GPS[!,:E].^2 .+ GPS[!,:N] .^2) ./ 2),window_len=45)
smoothEN .-= minimum(smoothEN)

# read PGA from USGS https://earthquake.usgs.gov/earthquakes/eventpage/ci9108652/shakemap/intensity
baja = h5open(joinpath(@__DIR__,"../data/Baja-EQ-shake.hdf"),"r") 
arrays = read(baja,"arrays")
PGA = arrays["imts"]["GREATER_OF_TWO_HORIZONTAL"]["PGA"]["mean"]
PGV = arrays["imts"]["GREATER_OF_TWO_HORIZONTAL"]["PGV"]["mean"]

# rotate array for plotting 
PGV = Array(transpose(reverse(PGV,dims=2)))
PGA = Array(transpose(reverse(PGA,dims=2)))

# transform from log(PGV) to PGV 
PGV = exp.(PGV)
PGA = exp.(PGA)

# get PGA grid spacing 
info = read(baja,"dictionaries/info.json") |> JSON.parse
Nlon = parse(Int,info["output"]["map_information"]["grid_points"]["longitude"])
Nlat = parse(Int,info["output"]["map_information"]["grid_points"]["latitude"])
latmax = parse(Float64,info["output"]["map_information"]["max"]["latitude"])
latmin = parse(Float64,info["output"]["map_information"]["min"]["latitude"])
lonmax = parse(Float64,info["output"]["map_information"]["max"]["longitude"])
lonmin = parse(Float64,info["output"]["map_information"]["min"]["longitude"])
PGAlon = range(lonmin,lonmax,length=Nlon)
PGAlat = range(latmin,latmax,length=Nlat)

# convert to GMT grid 
GPGV = GMT.mat2grid(PGV,x=PGAlon,y=PGAlat)
GPGA = GMT.mat2grid(PGA,x=PGAlon,y=PGAlat)
C = GMT.makecpt(T=(0,ceil(maximum(PGV)),0.25), cmap=:lajolla)
region = (-116.33,-114.5,31.66,33.33)

# plot HEC and LDES station map 
GMT.basemap(
    region=region,
    coast=true,
    colorbar=false,
    J=:guess,
)
GMT.coast!(
    N="a",
    ocean=:dodgerblue
    )
GMT.grdimage!(
    "@srtm_relief_03s",
    region=region,
    cmap=:gray,
    shade=true,
    colorbar=false,
    alpha=50,
)
GMT.grdimage!(
    GPGV,
    cmap=C,
    region=region,
    colorbar=true,
    alpha=50,
)
GMT.colorbar!(
    C=C,
    frame=(annot=:auto, ticks=:auto, xlabel="Peak Ground Velocity [cm/s]"),
)
GMT.basemap!(
    map_scale=(inside=true, anchor=:MB, scale_at_lat="32N", fancy=true,
    length="50k", units=true, offset=1.0),
)
GMT.psxy!(
    joinpath(@__DIR__,"../data/CFM/obj/traces/gmt/Baja_traces.lonLat"),
    pen=1,
)
GMT.scatter!(
    EQ[:,:LON],
    EQ[:,:LAT],
    markersize=EQ[:,:MAG] ./ 30.0,
    markeredgecolor=:blue,
)
GMT.scatter!(
    [p494lon],
    [p494lat],
    mc=:red,
    markeredgecolor=:black,
    marker=:circle,
    markersize=0.5,
)
GMT.text!(
    GMT.text_record(
        [p494lon .- 0.03 p494lat .+ 0.07],
        "p494",
    ),
)
GMT.scatter!(
    CAdf[:,:Longitude],
    CAdf[:,:Latitude],
    mc=:dodgerblue,
    markeredgecolor=:black,
    marker="t",
    markersize=0.5,
)
GMT.text!(
    GMT.text_record(
        [CAdf[:,:Longitude] .+ 0.12 CAdf[:,:Latitude] .+ 0.02],
        CAdf[:,:Network] .* "." .* CAdf[:,:Station],
    ),
    font=(10,)
)
# GMT.text!(
#     GMT.text_record(
#         [WESlon WESlat .+ 0.07],
#         "CI.WES",
#     ),
# )
GMT.text!(
    GMT.text_record(
        [-115.15 32.73],
        "California",
    ),
    angle=5,
)
GMT.text!(
    GMT.text_record(
        [-115.15 32.65],
        "Mexico",
    ),
    angle=5,
)      
GMT.text!(
    GMT.text_record(
        [-115.7 32.4],
        "Borrego Fault",
    ),
    angle=-45,
)     
GMT.psmeca!(
        [EQlon EQlat 10.0 222 59 -15 2 0 0], 
        aki=true, 
        fill=:red, 
        show=true,
        fmt=:png,
        savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/Baja-map.png"),
)

# fit model to dv/v decay 
eqind = argmin(WES[:,:DVV])

# get residual 
resid = WES[eqind:end,:DVV] .- WES[eqind:end,:ELASTIC]
resid1 = resid[1]
resid0 = resid .- resid1 
dtfit = (WES[eqind:end,:DATE] .- WES[eqind,:DATE]) ./ Day(1) ./ 365.25
preddate = WES[eqind,:DATE] :Day(1) : Date(2026,1,1)
dtpred = (preddate .- preddate[1]) ./ Day(1) ./ 365.25

# fitting functions 
# log model 
function logmodel(p)
    model = p[1] .* log.(1.0 .+ dtfit ./ p[2])
    fit = resid0 .- model 
    return sum(fit .^ 2)
end

# exp model 
function expmodel(p)
    model =  p[1] .* (1.0 .- exp.(-dtfit ./ p[2]))
    fit = resid0 .- model 
    return sum(fit .^ 2)
end

# fit log-log model 
# solve for log model 
reslog = optimize(
    logmodel, 
    [-Inf,5e-5],
    [Inf,Inf],
    [0.0,3.0],
    Fminbox(LBFGS()),
    Optim.Options(time_limit=60.0),
)
plog = Optim.minimizer(reslog)
l = plog[1] .* log.(1.0 .+ dtpred ./ plog[2])

# solve for exp model 
resexp = optimize(
    expmodel, 
    [-Inf,5e-5],
    [Inf,Inf],
    [0.0,3.0],
    Fminbox(LBFGS()),
    Optim.Options(time_limit=60.0),
)
pexp = Optim.minimizer(resexp)
e = pexp[1] .* (1.0 .- exp.(-dtpred ./ pexp[2]))


# plot WES dv/V
Ylims = (-2.0, 1.0)
Xlims = (WES[1,:DATE],preddate[end])
scatter(
    WES[!,:DATE],
    WES[!,:DVV] .- WES[:,:ELASTIC],
    seriescolor=:Blues_9,
    marker_z=WES[!,:CC],
    alpha=WES[!,:CC] ./ 10,
    label="",
    colorbar=false,
    ylabel="          dv/v [%]",
    # yaxis=:flip,
    ylims=Ylims,
    xlims=Xlims,
    # xlims=(WES[1,:DATE],WES[end,:DATE]),
    right_margin = 1.5cm,
    left_margin = 2cm,
    yguidefontcolor=:blue3,
    ytickfont=font(10,:blue3),
    yforeground_color_axis=:blue3,
    ytick_direction=:out,
    dpi=500,
    size=(800,400),
)
# # plot line for earthquake 
# vline!(
#     [Date(2010,4,4)],
#     label="Mw 7.2 Baja EQ",
#     linewidth=2,
#     c=:black,
#     linestyle=:dash,
#     legend=:topleft,
#     legendfontsize=9,
# )
# plot arrow for earthquake 
xdate = Date(2005) : Day(1) : Date(2010,1,1)
x = (xdate .- xdate[1]) ./ Day(1) ./ length(xdate) .* (π /3) .+ π / 6
Plots.plot!(
    xdate,
    sin.(x) .* 0.62 .- 1.3,
    arrow=true,
    lw=2,
    color=:black,
    label="",
)
Plots.annotate!(
    (
        Date(2005,6,1), 
        -1.3 ,
        Plots.text(
            "Baja CA Earthquake\n M7.2",
            12,
            :bold,
            :black,
        ),
    )
)
Plots.plot!(
    preddate,
    l .+ resid1,
    lw=3,
    c=:purple,
    alpha=0.75,
    ls=:solid,
    # label="$(round(plog[1],digits=2)) * log(t[years] / $(round(plog[2],digits=2)))",
    label="",
    legend=:bottomright,
    ylims=Ylims,
)
Plots.plot!(
    preddate,
    e .+ resid1,
    lw=3,
    c=:orange,
    alpha=0.75,
    ls=:solid,
    # label="$(round(pexp[1],digits=2)) * (1 - exp(-t[years] / $(round(pexp[2],digits=2))))",
    label="",
    legend=:bottomright,
    ylims=Ylims,
)
Plots.annotate!(
    (
        Date(2024,4,1), 
        0.57 ,
        Plots.text(
            "log",
            12,
            rotation = 10,
            :bold,
            :purple,
        ),
    ),
    subplot=1,
)
Plots.annotate!(
    (
        Date(2019,1,1), 
        -1.3 ,
        Plots.text(
            L"dv/v = %$(round(plog[1],digits=2)) * log(t[years] / %$(round(plog[2],digits=2)))",
            16,
            :bold,
            :purple,
        ),
    ),
    subplot=1,
)
Plots.annotate!(
    (
        Date(2019,1,1), 
        -1.7,
        Plots.text(
            L"dv/v = %$(round(pexp[1],digits=2)) * (1 - e^{-t[years] / %$(round(pexp[2],digits=2))})",
            16,
            :bold,
            :orange,
        ),
    ),
    subplot=1,
)

Plots.annotate!(
    (
        Date(2024,4,1), 
        0.1 ,
        Plots.text(
            "exp",
            12,
            :bold,
            :orange,
        ),
    ),
    subplot=1,
)

# plot GPS after earthquake 
plot!(
    twinx(), 
    GPS[!,:DATE], 
    smoothEN, 
    sharex=true, 
    xticks=false,
    grid=:off, 
    ylims=(0,90),
    linewidth=2.5,
    c=:red,
    alpha=0.95,
    label="",
    ylabel="GPS H. Displacement [mm]",
    # yaxis=:flip,
    # xlims=(WES[1,:DATE],WES[end,:DATE]),
    yguidefontcolor=:red,
    ytickfont=font(10,:red),
    yforeground_color_axis=:red,
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CIWES-DVV.png"))
