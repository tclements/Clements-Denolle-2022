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

# get WES lat, lon from CI station locations 
SCdf = DataFrame(CSV.File("/home/timclements/CALI/CIstations.csv"))
NCdf = DataFrame(CSV.File("/home/timclements/CALI/NCstations.csv"))
CAdf = vcat(NCdf, SCdf)
CAdf = CAdf[findall(in(["WES","HEC","JRC2"]),CAdf[:,:Station]),:]

PGVthresh = 10.0 # cm/s 

# read PGV from USGS https://earthquake.usgs.gov/earthquakes/eventpage/ci9108652/shakemap/intensity
baja = h5open("/media/FOUR/data/Baja-EQ-shake.hdf","r") 
arrays = read(baja,"arrays")
PGVbaja = arrays["imts"]["GREATER_OF_TWO_HORIZONTAL"]["PGV"]["mean"]

# rotate array for plotting 
PGVbaja = Array(transpose(reverse(PGVbaja,dims=2)))

# transform from log(PGV) to PGV 
PGVbaja = exp.(PGVbaja)
PGVbaja[PGVbaja .< PGVthresh] .= NaN

# get PGA grid spacing 
info = read(baja,"dictionaries/info.json") |> JSON.parse
Nlon = parse(Int,info["output"]["map_information"]["grid_points"]["longitude"])
Nlat = parse(Int,info["output"]["map_information"]["grid_points"]["latitude"])
latmax = parse(Float64,info["output"]["map_information"]["max"]["latitude"])
latmin = parse(Float64,info["output"]["map_information"]["min"]["latitude"])
lonmax = parse(Float64,info["output"]["map_information"]["max"]["longitude"])
lonmin = parse(Float64,info["output"]["map_information"]["min"]["longitude"])
bajalon = range(lonmin,lonmax,length=Nlon)
bajalat = range(latmin,latmax,length=Nlat)

# convert to GMT grid 
GPGVbaja = GMT.mat2grid(PGVbaja,x=bajalon,y=bajalat)

# read PGV from USGS https://earthquake.usgs.gov/earthquakes/eventpage/ci9108652/shakemap/intensity
hec = h5open("/media/FOUR/data/Hector-Mine-EQ.hdf","r") 
arrays = read(hec,"arrays")
PGVhec = arrays["imts"]["GREATER_OF_TWO_HORIZONTAL"]["PGV"]["mean"]

# rotate array for plotting 
PGVhec = Array(transpose(reverse(PGVhec,dims=2)))

# transform from log(PGV) to PGV 
PGVhec = exp.(PGVhec)
PGVhec[PGVhec .< PGVthresh] .= NaN


# get PGA grid spacing 
info = read(hec,"dictionaries/info.json") |> JSON.parse
Nlon = parse(Int,info["output"]["map_information"]["grid_points"]["longitude"])
Nlat = parse(Int,info["output"]["map_information"]["grid_points"]["latitude"])
latmax = parse(Float64,info["output"]["map_information"]["max"]["latitude"])
latmin = parse(Float64,info["output"]["map_information"]["min"]["latitude"])
lonmax = parse(Float64,info["output"]["map_information"]["max"]["longitude"])
lonmin = parse(Float64,info["output"]["map_information"]["min"]["longitude"])
heclon = range(lonmin,lonmax,length=Nlon)
heclat = range(latmin,latmax,length=Nlat)
# convert to GMT grid 
GPGVhec = GMT.mat2grid(PGVhec,x=heclon,y=heclat)

# read PGV from USGS https://earthquake.usgs.gov/earthquakes/eventpage/ci9108652/shakemap/intensity
rc = h5open("/media/FOUR/data/Ridgecrest-EQ.hdf","r") 
arrays = read(rc,"arrays")
PGVrc = arrays["imts"]["GREATER_OF_TWO_HORIZONTAL"]["PGV"]["mean"]

# rotate array for plotting 
PGVrc = Array(transpose(reverse(PGVrc,dims=2)))

# transform from log(PGV) to PGV 
PGVrc = exp.(PGVrc)
PGVrc[PGVrc .< PGVthresh] .= NaN

# get PGA grid spacing 
info = read(rc,"dictionaries/info.json") |> JSON.parse
Nlon = parse(Int,info["output"]["map_information"]["grid_points"]["longitude"])
Nlat = parse(Int,info["output"]["map_information"]["grid_points"]["latitude"])
latmax = parse(Float64,info["output"]["map_information"]["max"]["latitude"])
latmin = parse(Float64,info["output"]["map_information"]["min"]["latitude"])
lonmax = parse(Float64,info["output"]["map_information"]["max"]["longitude"])
lonmin = parse(Float64,info["output"]["map_information"]["min"]["longitude"])
rclon = range(lonmin,lonmax,length=Nlon)
rclat = range(latmin,latmax,length=Nlat)

# convert to GMT grid 
GPGVrc = GMT.mat2grid(PGVrc,x=rclon,y=rclat)
 
# GPS locations 
PBO = CSV.File(
    "/media/FOUR/data/PBO.csv",
) |> DataFrame
PBO = PBO[findall(in(["P494","BEPK","LDSW"]),PBO[:,:pnum]),:]

# M 7.2 Baja California Sierra El Mayor earthquake
bajaEQlon = -115.295
bajaEQlat = 32.286

# M7.1 Hector Mine Location 
hecEQlon = -116.268
hecEQlat = 34.60

# M7.1 Ridgecrest Location 
rcEQlon = -117.599
rcEQlat = 35.77

# load Baja Earthquake locations 
bajaEQ = DataFrame(
    CSV.File(
    "/media/FOUR/data/baja.txt",
    delim=" ",
    ignorerepeated=true,
    )
)
bajaEQ = bajaEQ[bajaEQ[:,:MAG] .> 3.0,:]

# load Hector Mine Earthquake locations 
hecEQ = DataFrame(
    CSV.File(
    "/media/FOUR/data/hectormine.txt",
    delim=" ",
    ignorerepeated=true,
    )
)
hecEQ = hecEQ[hecEQ[:,:MAG] .> 3.0, :]

# load Ridgecrest Earthquake locations 
rcEQ = DataFrame(
    CSV.File(
    "/media/FOUR/data/ridgecrest.txt",
    delim=" ",
    ignorerepeated=true,
    )
)
rcEQ = rcEQ[rcEQ[:,:MAG] .> 3.0, :]

# setup colormap 
C = GMT.makecpt(T=(PGVthresh,ceil(maximum(x->isnan(x) ? -Inf : x,PGVrc)),0.25), cmap=:lajolla)
region = (-118.5,-114.0,31.0,36.66)

# plot HEC and LDES station map 
GMT.basemap(
    region=region,
    coast=true,
    colorbar=false,
    J=:guess,
)
GMT.grdimage!(
    "@earth_relief_15s",
    region=region,
    cmap=:gray,
    shade=true,
    colorbar=false,
    alpha=50,
)
GMT.grdimage!(
    GPGVbaja,
    cmap=C,
    region=region,
    colorbar=false,
    alpha=50,
    Q=true,
)
GMT.grdimage!(
    GPGVhec,
    cmap=C,
    region=region,
    colorbar=false,
    alpha=50,
    Q=true,
)
GMT.grdimage!(
    GPGVrc,
    cmap=C,
    region=region,
    colorbar=false,
    alpha=50,
    Q=true,
)
GMT.coast!(
    N="a",
    ocean=:white,
)
GMT.basemap!(
    map_scale=(inside=true, anchor=:RB, scale_at_lat="32N", fancy=true,
    length="100k", units=true, offset=1.0),
)
GMT.psxy!(
    "/media/FOUR/data/CFM/obj/traces/gmt/Baja_traces.lonLat",
    pen=1,
)
GMT.psxy!(
    "/media/FOUR/data/CFM/obj/traces/gmt/HectorMine_traces.lonLat",
    pen=1,
)
GMT.psxy!(
    "/media/FOUR/data/CFM/obj/traces/gmt/Ridgecrest_traces.lonLat",
    pen=1,
)
GMT.scatter!(
    PBO[:,:lon],
    PBO[:,:lat],
    mc=:red,
    markeredgecolor=:black,
    marker=:circle,
    markersize=0.5,
)
GMT.scatter!(
    CAdf[:,:Longitude],
    CAdf[:,:Latitude],
    mc=:dodgerblue,
    markeredgecolor=:black,
    marker="t",
    markersize=0.5,
)
GMT.colorbar!(
    C=C,
    frame=(annot=:auto, ticks=:auto, xlabel="PGV [cm/s]"),
    pos=(anchor=:BL,horizontal=true,offset=(-5.5,-2),length=5),
)
# plot focal mechanisms with lines between them 
offset = 0.25
GMT.plot!(
    [hecEQlon hecEQlat; hecEQlon + offset hecEQlat - offset / 2],
    lw=0.5,
    lc=:red,
)
GMT.psmeca!(
        [hecEQlon + offset hecEQlat - offset / 2 13.7 157 83 -162 1.5 0 0], 
        aki=true, 
        fill=:red, 
)
GMT.plot!(
    [rcEQlon rcEQlat; rcEQlon + offset rcEQlat + offset],
    lw=0.5,
    lc=:red,
)
GMT.psmeca!(
        [rcEQlon + offset rcEQlat + offset 8.0 230 81 -6 1.5 0 0], 
        aki=true, 
        fill=:red, 
)
GMT.plot!(
    [bajaEQlon bajaEQlat; bajaEQlon + offset bajaEQlat + offset],
    lw=0.5,
    lc=:red,
)
GMT.psmeca!(
        [bajaEQlon + offset bajaEQlat + offset 10.0 222 59 -15 1.5 0 0], 
        aki=true, 
        fill=:red, 
        show=true,
        fmt=:png,
        savefig="/media/FOUR/data/FINAL-FIGURES/EQ-map.png",
)

##### plot dv/v and GPS ####
# load dv/v 
WES = Arrow.Table("/media/FOUR/data/FIT-DVV-SSE/90-DAY/CI.WES.arrow") |> DataFrame
HEC = Arrow.Table("/media/FOUR/data/FIT-DVV-SSE/90-DAY/CI.HEC.arrow") |> DataFrame
JRC2 = Arrow.Table("/media/FOUR/data/FIT-DVV-SSE/90-DAY/CI.JRC2.arrow") |> DataFrame

# load GPS data 
BEPK = CSV.File(
    "/media/FOUR/data/WNAMdetrend/bepkFilterDetrend.neu", 
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
BEPK[!,:DATE] = Date.(BEPK[!,:YEAR]) .+ Day.(BEPK[!,:DAY] .- 1)
LDSW = CSV.File(
    "/media/FOUR/data/WNAMdetrend/ldswFilterDetrend.neu", 
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
LDSW[!,:DATE] = Date.(LDSW[!,:YEAR]) .+ Day.(LDSW[!,:DAY] .- 1)
P494 = CSV.File(
    "/media/FOUR/data/WNAMdetrend/p494FilterDetrend.neu", 
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
P494[!,:DATE] = Date.(P494[!,:YEAR]) .+ Day.(P494[!,:DAY] .- 1)

# select GPS after day of each earthquake 
# Baja earthquake
P494 = P494[(P494[!,:DATE] .>= Date(2010,4,5)) .& (P494[!,:DATE] .<= WES[end,:DATE]),:]

# Hector Mine earthquake 
LDSW = LDSW[(LDSW[!,:DATE] .>= Date(1999,10,15)),:]

# Ridgecrest earthquake 
BEPK = BEPK[BEPK[:,:DATE] .>= Date(2019,7,7),:]

# smooth GPS stations 
smoothBEPK = smooth_withfiltfilt(
    sqrt.( 
        (BEPK[!,:E].^2 .+ BEPK[!,:N] .^2) ./ 2,
    ),
    window_len=45,
)
smoothP494 = smooth_withfiltfilt(
    sqrt.( 
        (P494[!,:E].^2 .+ P494[!,:N] .^2) ./ 2,
    ),
    window_len=45,
)
smoothLDSW = smooth_withfiltfilt(
    sqrt.( 
        (LDSW[!,:E].^2 .+ LDSW[!,:N] .^2) ./ 2,
    ),
    window_len=45,
)

# remove minimum value 
smoothBEPK .-= minimum(smoothBEPK)
smoothP494 .-= minimum(smoothP494)
smoothLDSW .-= minimum(smoothLDSW)

# last date for fitting 
lastdate = Date(2026,1,1)

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

#### fit model to Ridgecrest earthquake dv/v ####
RCdate = Date(2019,7,4)
beforeRC = findlast(JRC2[:,:DATE] .< RCdate)
eqind = argmin(JRC2[:,:DVV] .- JRC2[:,:ELASTIC])

# get residual 
resid = JRC2[eqind:end,:DVV] .- JRC2[eqind:end,:ELASTIC]
resid1 = resid[1]
resid0 = resid .- resid1 
dtfit = (JRC2[eqind:end,:DATE] .- JRC2[eqind,:DATE]) ./ Day(1) ./ 365.25
preddate = JRC2[eqind,:DATE] :Day(1) : lastdate
dtJRC2 = (preddate .- preddate[1]) ./ Day(1) ./ 365.25

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
plogJRC2 = Optim.minimizer(reslog)
lJRC2 = plogJRC2[1] .* log.(1.0 .+ dtJRC2 ./ plogJRC2[2])
println("Ridgecrest Log Relaxation: $(round(plogJRC2[2],digits=2)) years")

# solve for exp model 
resexp = optimize(
    expmodel, 
    [-Inf,5e-5],
    [Inf,Inf],
    [0.0,3.0],
    Fminbox(LBFGS()),
    Optim.Options(time_limit=60.0),
)
pexpJRC2 = Optim.minimizer(resexp)
eJRC2 = pexpJRC2[1] .* (1.0 .- exp.(-dtJRC2 ./ pexpJRC2[2]))
println("Ridgecrest Exp Relaxation: $(round(pexpJRC2[2],digits=2)) years")

Xlims = (HEC[1,:DATE],preddate[end])

# plot JRC2 dv/v
JRC2lims = (-1.62, 0.9)
p1 = scatter(
    JRC2[!,:DATE],
    JRC2[!,:DVV] .- JRC2[:,:ELASTIC] .- JRC2[beforeRC,:DVV],
    seriescolor=:Blues_9,
    marker_z=JRC2[!,:CC],
    alpha=JRC2[!,:CC] ./ 10,
    label="",
    colorbar=false,
    ylabel="dv/v [%]",
    # yaxis=:flip,
    ylims=JRC2lims,
    xlims=Xlims,
    # xlims=(WES[1,:DATE],WES[end,:DATE]),
    right_margin = 2.0cm,
    left_margin = 2cm,
    yguidefontcolor=:blue3,
    ytickfont=font(10,:blue3),
    yforeground_color_axis=:blue3,
    ytick_direction=:out,
    dpi=500,
    size=(800,400),
)
annotate!((Date(1994,1,1),0.82,text("B",20)))
# # plot arrow for earthquake 
# xdate = Date(2005) : Day(1) : Date(2010,1,1)
# x = (xdate .- xdate[1]) ./ Day(1) ./ length(xdate) .* (π /3) .+ π / 6
# Plots.plot!(
#     xdate,
#     sin.(x) .* 0.62 .- 1.3,
#     arrow=true,
#     lw=2,
#     color=:black,
#     label="",
# )
# Plots.annotate!(
#     (
#         Date(2005,6,1), 
#         -1.3 ,
#         Plots.text(
#             "Baja CA Earthquake\n M7.2",
#             12,
#             :bold,
#             :black,
#         ),
#     )
# )
Plots.plot!(
    preddate,
    lJRC2 .+ resid1 .- JRC2[beforeRC,:DVV],
    lw=3,
    c=:purple,
    alpha=0.75,
    ls=:solid,
    # label="$(round(plog[1],digits=2)) * log(t[years] / $(round(plog[2],digits=2)))",
    label="",
    legend=:bottomright,
    ylims=JRC2lims,
)
Plots.plot!(
    preddate,
    eJRC2 .+ resid1 .- JRC2[beforeRC,:DVV],
    lw=3,
    c=:orange,
    alpha=0.75,
    ls=:solid,
    # label="$(round(pexp[1],digits=2)) * (1 - exp(-t[years] / $(round(pexp[2],digits=2))))",
    label="",
    legend=:bottomright,
    ylims=JRC2lims,
)
Plots.annotate!(
    (
        Date(2023,4,1), 
        -0.5,
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
# Plots.annotate!(
#     (
#         Date(2019,1,1), 
#         -1.3 ,
#         Plots.text(
#             L"dv/v = %$(round(plog[1],digits=2)) * log(t[years] / %$(round(plog[2],digits=2)))",
#             16,
#             :bold,
#             :purple,
#         ),
#     ),
#     subplot=1,
# )
# Plots.annotate!(
#     (
#         Date(2019,1,1), 
#         -1.7,
#         Plots.text(
#             L"dv/v = %$(round(pexp[1],digits=2)) * (1 - e^{-t[years] / %$(round(pexp[2],digits=2))})",
#             16,
#             :bold,
#             :orange,
#         ),
#     ),
#     subplot=1,
# )
Plots.annotate!(
    (
        Date(2023,4,1), 
        -1.02,
        Plots.text(
            "exp",
            12,
            :bold,
            :orange,
        ),
    ),
    subplot=1,
)
Plots.annotate!(
    (
        Date(2010,1,1), 
        -1.25,
        Plots.text(
            "M7.1 Ridgecrest Earthquake",
            12,
            :bold,
            :black,
        ),
    ),
    subplot=1,
)
# plot GPS after earthquake 
plot!(
    twinx(), 
    BEPK[:,:DATE],
    smoothBEPK,
    xticks=false,
    grid=:off, 
    ylims=(0,45),
    linewidth=2.5,
    xlims=Xlims,
    c=:red,
    alpha=0.95,
    label="",
    ylabel="GPS Displ. [mm]",
    # yaxis=:flip,
    # xlims=(WES[1,:DATE],WES[end,:DATE]),
    yguidefontcolor=:red,
    ytickfont=font(10,:red),
    yforeground_color_axis=:red,
)

#### fit model to Hector Mine earthquake dv/v #### 
# get day of Hector Mine EQ 
HMdate = Date(1999,10,16)
beforeHM = findlast(HEC[:,:DATE] .< HMdate)
eqind = argmin(HEC[:,:DVV])

# get residual 
resid = HEC[eqind:end,:DVV] .- HEC[eqind:end,:ELASTIC]
resid1 = resid[1]
resid0 = resid .- resid1 
dtfit = (HEC[eqind:end,:DATE] .- HEC[eqind,:DATE]) ./ Day(1) ./ 365.25
preddate = HEC[eqind,:DATE] :Day(1) : lastdate
dtHEC = (preddate .- preddate[1]) ./ Day(1) ./ 365.25

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
plogHEC = Optim.minimizer(reslog)
lHEC = plogHEC[1] .* log.(1.0 .+ dtHEC ./ plogHEC[2])
println("Hector Mine Log Relaxation: $(round(plogHEC[2],digits=2)) years")

# solve for exp model 
resexp = optimize(
    expmodel, 
    [-Inf,5e-5],
    [Inf,Inf],
    [0.0,20.0],
    Fminbox(LBFGS()),
    Optim.Options(time_limit=60.0),
)
pexpHEC = Optim.minimizer(resexp)
eHEC = pexpHEC[1] .* (1.0 .- exp.(-dtHEC ./ pexpHEC[2]))
println("Hector Mine Exp Relaxation: $(round(pexpHEC[2],digits=2)) years")

HEClims = (-1.75,0.3)
p2 = scatter(
    HEC[!,:DATE],
    HEC[!,:DVV] .- HEC[:,:ELASTIC] .- HEC[beforeHM,:DVV],
    seriescolor=:Blues_9,
    marker_z=HEC[!,:CC],
    alpha=HEC[!,:CC] ./ 10,
    label="",
    colorbar=false,
    ylabel="dv/v [%]",
    # yaxis=:flip,
    ylims=HEClims,
    xlims=Xlims,
    # xlims=(HEC[1,:DATE],HEC[end,:DATE]),
    right_margin = 2.0cm,
    left_margin = 2cm,
    yguidefontcolor=:blue3,
    ytickfont=font(10,:blue3),
    yforeground_color_axis=:blue3,
    ytick_direction=:out,
    dpi=500,
    size=(800,400),
)
annotate!((Date(1994,1,1),0.24,text("C",20)))
Plots.plot!(
    preddate,
    lHEC .+ resid1 .- HEC[beforeHM,:DVV],
    lw=3,
    c=:purple,
    alpha=0.75,
    ls=:solid,
    # label="$(round(plog[1],digits=2)) * log(t[years] / $(round(plog[2],digits=2)))",
    label="",
    legend=:bottomright,
    ylims=HEClims,
    xlims=Xlims,
)
Plots.plot!(
    preddate,
    eHEC .+ resid1 .- HEC[beforeHM,:DVV],
    lw=3,
    c=:orange,
    alpha=0.75,
    ls=:solid,
    # label="$(round(pexp[1],digits=2)) * (1 - exp(-t[years] / $(round(pexp[2],digits=2))))",
    label="",
    legend=:bottomright,
    ylims=HEClims,
    xlims=Xlims,
)
Plots.annotate!(
    (
        Date(2023,4,1), 
        -0.05,
        Plots.text(
            "log",
            12,
            rotation = 10,
            :bold,
            :purple,
        ),
    )
)
# Plots.annotate!(
#     (
#         Date(2016,6,1), 
#         -0.6 ,
#         Plots.text(
#             L"dv/v = %$(round(plog[1],digits=2)) * log(t[years] / %$(round(plog[2],digits=2)))",
#             16,
#             :bold,
#             :purple,
#         ),
#     ),
#     subplot=1,
# )
# Plots.annotate!(
#     (
#         Date(2016,6,1), 
#         -0.8 ,
#         Plots.text(
#             L"dv/v = %$(round(pexp[1],digits=2)) * (1 - e^{-t[years] / %$(round(pexp[2],digits=2))})",
#             16,
#             :bold,
#             :orange,
#         ),
#     ),
#     subplot=1,
# )
Plots.annotate!(
    (
        Date(2023,4,1), 
        -0.45,
        Plots.text(
            "exp",
            12,
            :bold,
            :orange,
        ),
    ),
    subplot=1,
)
# # plot arrow for earthquake 
# xdate = Date(2005) : Day(-1) : Date(2000,2,1)
# x = (xdate .- xdate[end]) ./ Day(1) ./ length(xdate) .* (π /3) .+ π / 6
# Plots.plot!(
#     xdate,
#     sin.(x) .* 0.62 .- 0.3,
#     arrow=true,
#     lw=2,
#     color=:black,
#     label="",
#     xlims=Xlims,
# )
Plots.annotate!(
    (
        Date(2010,1,1), 
        -1.25 ,
        Plots.text(
            "M7.1 Hector Mine Earthquake",
            12,
            :bold,
            :black,
        ),
    ),
    subplot=1,
)
# plot LDSW vertical displacment 
Plots.plot!(
    twinx(),
    LDSW[:,:DATE],
    smoothLDSW,
    xticks=false,
    # sharex=true,
    xlims=Xlims,
    grid=:off, 
    ylims=(0,22),
    lw=2.5,
    c=:red,
    alpha=0.95,
    label="",
    ylabel="GPS Displ. [mm]",
    yguidefontcolor=:red,
    ytickfont=font(10,:red),
    yforeground_color_axis=:red,
)


#### fit model to Baja earthquake dv/v ####
BCdate = Date(2010,4,4)
beforeBC = findlast(WES[:,:DATE] .< BCdate)
eqind = argmin(WES[:,:DVV])

# get residual 
resid = WES[eqind:end,:DVV] .- WES[eqind:end,:ELASTIC]
resid1 = resid[1]
resid0 = resid .- resid1 
dtfit = (WES[eqind:end,:DATE] .- WES[eqind,:DATE]) ./ Day(1) ./ 365.25
preddate = WES[eqind,:DATE] :Day(1) : lastdate
dtWES = (preddate .- preddate[1]) ./ Day(1) ./ 365.25

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
plogWES = Optim.minimizer(reslog)
lWES = plogWES[1] .* log.(1.0 .+ dtWES ./ plogWES[2])
println("Baja CA Log Relaxation: $(round(plogWES[2],digits=2)) years")

# solve for exp model 
resexp = optimize(
    expmodel, 
    [-Inf,5e-5],
    [Inf,Inf],
    [0.0,3.0],
    Fminbox(LBFGS()),
    Optim.Options(time_limit=60.0),
)
pexpWES = Optim.minimizer(resexp)
eWES = pexpWES[1] .* (1.0 .- exp.(-dtWES ./ pexpWES[2]))
println("Baja CA Exp Relaxation: $(round(pexpWES[2],digits=2)) years")

# plot WES dv/v
WESlims = (-2.61, 0.39)
p3 = scatter(
    WES[!,:DATE],
    WES[!,:DVV] .- WES[:,:ELASTIC] .- WES[beforeBC,:DVV],
    seriescolor=:Blues_9,
    marker_z=WES[!,:CC],
    alpha=WES[!,:CC] ./ 10,
    label="",
    colorbar=false,
    ylabel="dv/v [%]",
    # yaxis=:flip,
    ylims=WESlims,
    xlims=Xlims,
    # xlims=(WES[1,:DATE],WES[end,:DATE]),
    right_margin = 2.0cm,
    left_margin = 2cm,
    yguidefontcolor=:blue3,
    ytickfont=font(10,:blue3),
    yforeground_color_axis=:blue3,
    ytick_direction=:out,
    dpi=500,
    size=(800,400),
)
annotate!((Date(1994,1,1),0.3,text("D",20)))
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
# # plot arrow for earthquake 
# xdate = Date(2005) : Day(1) : Date(2010,1,1)
# x = (xdate .- xdate[1]) ./ Day(1) ./ length(xdate) .* (π /3) .+ π / 6
# Plots.plot!(
#     xdate,
#     sin.(x) .* 0.62 .- 1.3,
#     arrow=true,
#     lw=2,
#     color=:black,
#     label="",
# )
# Plots.annotate!(
#     (
#         Date(2005,6,1), 
#         -1.3 ,
#         Plots.text(
#             "Baja CA Earthquake\n M7.2",
#             12,
#             :bold,
#             :black,
#         ),
#     )
# )
Plots.plot!(
    preddate,
    lWES .+ resid1 .- WES[beforeBC,:DVV],
    lw=3,
    c=:purple,
    alpha=0.75,
    ls=:solid,
    # label="$(round(plog[1],digits=2)) * log(t[years] / $(round(plog[2],digits=2)))",
    label="",
    legend=:bottomright,
    ylims=WESlims,
)
Plots.plot!(
    preddate,
    eWES .+ resid1 .- WES[beforeBC,:DVV],
    lw=3,
    c=:orange,
    alpha=0.75,
    ls=:solid,
    # label="$(round(pexp[1],digits=2)) * (1 - exp(-t[years] / $(round(pexp[2],digits=2))))",
    label="",
    legend=:bottomright,
    ylims=WESlims,
)
Plots.annotate!(
    (
        Date(2023,4,1), 
        0.05,
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
# Plots.annotate!(
#     (
#         Date(2019,1,1), 
#         -1.3 ,
#         Plots.text(
#             L"dv/v = %$(round(plog[1],digits=2)) * log(t[years] / %$(round(plog[2],digits=2)))",
#             16,
#             :bold,
#             :purple,
#         ),
#     ),
#     subplot=1,
# )
# Plots.annotate!(
#     (
#         Date(2019,1,1), 
#         -1.7,
#         Plots.text(
#             L"dv/v = %$(round(pexp[1],digits=2)) * (1 - e^{-t[years] / %$(round(pexp[2],digits=2))})",
#             16,
#             :bold,
#             :orange,
#         ),
#     ),
#     subplot=1,
# )
Plots.annotate!(
    (
        Date(2023,4,1), 
        -0.53,
        Plots.text(
            "exp",
            12,
            :bold,
            :orange,
        ),
    ),
    subplot=1,
)
Plots.annotate!(
    (
        Date(2005,1,1), 
        -1.75 ,
        Plots.text(
            "M7.2 Baja California\n Earthquake ",
            12,
            :bold,
            :black,
        ),
    ),
    subplot=1,
)
# plot GPS after earthquake 
plot!(
    twinx(), 
    P494[!,:DATE], 
    smoothP494, 
    xlims=Xlims, 
    xticks=false,
    grid=:off, 
    ylims=(0,90),
    linewidth=2.5,
    c=:red,
    alpha=0.95,
    label="",
    ylabel="GPS Displ. [mm]",
    # yaxis=:flip,
    # xlims=(WES[1,:DATE],WES[end,:DATE]),
    yguidefontcolor=:red,
    ytickfont=font(10,:red),
    yforeground_color_axis=:red,
)

# plot final figure with map 
img = load("/media/FOUR/data/FINAL-FIGURES/EQ-map.png")
p4 = plot(img,border=:none,dpi=250,size=(400,600))
annotate!((-100,100,text("A",20)))
l = Plots.@layout [a{0.5w} [b; c; d;]]
Plots.plot(p4,p1,p2,p3,layout=l,size=(1200,700))
Plots.savefig("/media/FOUR/data/FINAL-FIGURES/relaxation-map.png")
