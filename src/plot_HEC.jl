using Arrow 
using CSV
using DataFrames
using Dates
using DSP 
using HDF5
using Images
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

# get CIstation location 
SCdf = DataFrame(CSV.File("/home/timclements/CALI/CIstations.csv"))
NCdf = DataFrame(CSV.File("/home/timclements/CALI/NCstations.csv"))
CAdf = vcat(NCdf, SCdf)
CAdf = CAdf[CAdf[:,:StartTime] .< Date(1999,10,16),:]
CAdf = CAdf[CAdf[:,:EndTime] .> Date(1999,10,16),:]
# remove CI.RAG - no data until 2015 
CAdf = CAdf[CAdf[:,:Station] .!= "RAG",:]
CAdf = CAdf[CAdf[:,:Longitude] .> -117.16,:]
CAdf = CAdf[CAdf[:,:Latitude] .> 33.85,:]
HECind = findfirst(CAdf[!,:Station] .== "HEC")
HEClat = CAdf[HECind,:Latitude]
HEClon = CAdf[HECind,:Longitude]

# GPS location 
LDESlon = -116.432803
LDESlat = 34.267338

# Hector Mine Location 
EQlon = -116.268
EQlat = 34.60

# load EQ locations 
EQ = DataFrame(
    CSV.File(
    "/media/FOUR/data/hectormine.txt",
    delim=" ",
    ignorerepeated=true,
    )
)

# read PGA from USGS https://earthquake.usgs.gov/earthquakes/eventpage/ci9108652/shakemap/intensity
hmine = h5open("/media/FOUR/data/Hector-Mine-EQ.hdf","r") 
arrays = read(hmine,"arrays")
PGA = arrays["imts"]["GREATER_OF_TWO_HORIZONTAL"]["PGA"]["mean"]
PGV = arrays["imts"]["GREATER_OF_TWO_HORIZONTAL"]["PGV"]["mean"]

# rotate array for plotting 
PGV = Array(transpose(reverse(PGV,dims=2)))
PGA = Array(transpose(reverse(PGA,dims=2)))

# transform from log(PGV) to PGV 
PGV = exp.(PGV)
PGA = exp.(PGA)

# get PGA grid spacing 
info = read(hmine,"dictionaries/info.json") |> JSON.parse
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
region = (-117.0,-115.5,33.85,35.15)

# plot HEC and LDES station map 
GMT.basemap(
    region=region,
    coast=true,
    colorbar=false,
    J=:guess,
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
    map_scale=(inside=true, anchor=:MB, scale_at_lat="34N", fancy=true,
    length="40k", units=true, offset=1.0),
)
GMT.psxy!("/media/FOUR/data/CFM/obj/traces/gmt/HectorMine_traces.lonLat",pen=1)
GMT.scatter!(
    EQ[:,:LON],
    EQ[:,:LAT],
    markersize=EQ[:,:MAG] ./ 30.0,
    markeredgecolor=:blue,
)
GMT.scatter!(
    [LDESlon],
    [LDESlat],
    mc=:red,
    markeredgecolor=:black,
    marker=:circle,
    markersize=0.5,
)
GMT.text!(
    GMT.text_record(
        [LDESlon .+ 0.1 LDESlat],
        "LDES",
    ),
    font=(10,),
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
        [CAdf[:,:Longitude] .+ 0.02 CAdf[:,:Latitude] .+ 0.07],
        CAdf[:,:Network] .* "." .* CAdf[:,:Station],
    ),
    font=(10,)
)

# GMT.scatter!(
#     [HEClon],
#     [HEClat],
#     mc=:dodgerblue,
#     markeredgecolor=:black,
#     marker="t",
#     markersize=0.5,
#     # fmt=:png,
#     # savefig="/media/FOUR/data/FINAL-FIGURES/CI-station-map.png",
# )
# # GMT.text!(
#     GMT.text_record(
#         [HEClon .+ 0.1 HEClat + 0.02],
#         "CI.HEC",
#     ),
# )
GMT.text!(
    GMT.text_record(
        [-116.47 34.67],
        "Lavic Lake Fault",
    ),
    angle=-55,
)
GMT.text!(
    GMT.text_record(
        [-116.05 34.36],
        "Bullion Fault",
    ),
    angle=-55,
)      
GMT.psmeca!(
        [-116.268 34.60 13.7 157 83 -162 1.5 0 0], 
        aki=true, 
        fill=:red, 
        show=true,
        fmt=:png,
        savefig="/media/FOUR/data/FINAL-FIGURES/Hector-Mine-map.png",
)


# load dv/v 
HEC = Arrow.Table("/media/FOUR/data/FIT-DVV-SSE/90-DAY/CI.HEC.arrow") |> DataFrame

# load GPS from station LDES 
GPS = CSV.File(
    "/media/FOUR/data/WNAMtrend/ldesFilterTrend.neu", 
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
# select after Hector Mine earthquake
GPS = GPS[(GPS[!,:DATE] .>= Date(1999,10,15)),:]
smoothV = smooth_withfiltfilt(GPS[:,:V],window_len=45)

# get day of Hector Mine EQ 
eqind = argmin(HEC[:,:DVV])

# get residual 
resid = HEC[eqind:end,:DVV] .- HEC[eqind:end,:ELASTIC]
resid1 = resid[1]
resid0 = resid .- resid1 
dtfit = (HEC[eqind:end,:DATE] .- HEC[eqind,:DATE]) ./ Day(1) ./ 365.25
preddate = HEC[eqind,:DATE] :Day(1) : Date(2025,1,1)
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
    [0.0,20.0],
    Fminbox(LBFGS()),
    Optim.Options(time_limit=60.0),
)
pexp = Optim.minimizer(resexp)
e = pexp[1] .* (1.0 .- exp.(-dtpred ./ pexp[2]))

# plot 
Ylims = (-1.1,0.5)
Xlims = (HEC[1,:DATE],preddate[end])
scatter(
    HEC[!,:DATE],
    HEC[!,:DVV] .- HEC[:,:ELASTIC],
    seriescolor=:Blues_9,
    marker_z=HEC[!,:CC],
    alpha=HEC[!,:CC] ./ 10,
    label="",
    colorbar=false,
    ylabel="          dv/v [%]",
    # yaxis=:flip,
    ylims=Ylims,
    xlims=Xlims,
    # xlims=(HEC[1,:DATE],HEC[end,:DATE]),
    right_margin = 1.5cm,
    left_margin = 2cm,
    yguidefontcolor=:blue3,
    ytickfont=font(10,:blue3),
    yforeground_color_axis=:blue3,
    ytick_direction=:out,
    dpi=500,
    size=(800,400),
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
    xlims=Xlims,
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
    xlims=Xlims,
)
Plots.annotate!(
    (
        Date(2023,4,1), 
        0.33 ,
        Plots.text(
            "log",
            12,
            rotation = 10,
            :bold,
            :purple,
        ),
    )
)
Plots.annotate!(
    (
        Date(2016,6,1), 
        -0.6 ,
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
        Date(2016,6,1), 
        -0.8 ,
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
        Date(2023,4,1), 
        0.17 ,
        Plots.text(
            "exp",
            12,
            :bold,
            :orange,
        ),
    ),
    subplot=1,
)
# plot arrow for earthquake 
xdate = Date(2005) : Day(-1) : Date(2000,2,1)
x = (xdate .- xdate[end]) ./ Day(1) ./ length(xdate) .* (π /3) .+ π / 6
Plots.plot!(
    xdate,
    sin.(x) .* 0.62 .- 0.3,
    arrow=true,
    lw=2,
    color=:black,
    label="",
    xlims=Xlims,
)
Plots.annotate!(
    (
        Date(2006,4,1), 
        0.35 ,
        Plots.text(
            "Hector Mine Earthquake\n M7.1",
            12,
            :bold,
            :black,
        ),
    )
)

# plot LDES vertical displacment 
Plots.plot!(
    twinx(),
    GPS[:,:DATE],
    smoothV,
    xticks=false,
    # sharex=true,
    xlims=Xlims,
    grid=:off, 
    ylims=(-5,85),
    lw=2.5,
    c=:red,
    alpha=0.95,
    label="",
    ylabel="GPS V. Displacement [mm]",
    yguidefontcolor=:red,
    ytickfont=font(10,:red),
    yforeground_color_axis=:red,
)
Plots.savefig("/media/FOUR/data/FINAL-FIGURES/CIHEC-DVV.png")
