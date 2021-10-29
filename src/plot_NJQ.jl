using Arrow 
using CSV
using DataFrames
using Dates 
using LaTeXStrings
using LsqFit
using Plots 
using Plots.Measures
using NetCDF
using Statistics

function CDM(A::AbstractArray,k::Int)
    Amean = [cumsum(A[1:k-1]) ./ (1:k-1) ;rolling_mean(A,k)]
    return cumsum(A .- Amean)
end

function rolling_mean(A::AbstractArray,k::Int)
    B = cumsum([zero(eltype(A)); A])
    return (B[k+1:end] - B[1:end-k]) ./ k
end

# read data from CI.NJQ 
NJQ = Arrow.Table("/media/FOUR/data/FIT-DVV-SSE/90-DAY/CI.NJQ.arrow") |> DataFrame

# Alisal reservoir elevation 
ALI = CSV.File("/media/FOUR/data/alisalreservoir.tsv",comment="#",skipto=3) |> DataFrame
ALI[!,:WSE] ./= 3.28 # convert water surface elevation to m 
dropmissing!(ALI,:WSE)
ALI[:,:WSE] .-= mean(ALI[:,:WSE])

# daily discharge in nearby river
SYD = CSV.File("/media/FOUR/data/santaynezdischarge.tsv",comment="#",skipto=3) |> DataFrame
dropmissing!(SYD,:DISCHARGE) # DISCHARGE is in ft^3 / s (every 15 minute)
SYD[!,:DISCHARGE] ./= 35.314 # convert to meters^3 / s
SYD[!,:DISCHARGE] .*= 900 # convert m^3 / s to m^3
SYD[!,:TIME] = DateTime.(SYD[!,:DATE],"yyyy-mm-dd HH:MM")
SYD[!,:DATE] = Date.(SYD[!,:DATE],"yyyy-mm-dd HH:MM")

SYDg = groupby(SYD,:DATE)
SYDdaily = combine(SYDg,:DISCHARGE => sum)
tdis = SYDdaily[!,:DATE]
discharge = Array(SYDdaily[!,:DISCHARGE_sum])
cumdischarge = zeros(eltype(discharge),size(discharge,1))
octs = Date(2001,10,1):Year(1):Date(2020,10,1)
juns = Date(2002,6,1):Year(1):Date(2021,6,1)
for dd = 1:length(octs)
    ind = findall((tdis .> octs[dd]) .& (tdis .< juns[dd]))
    cumdischarge[ind] .= cumsum(discharge[ind])
end
cumdischarge[iszero.(cumdischarge)] .= NaN

# subset t with dv/v 
tind = findall((tdis .>= NJQ[1,:DATE]) .& (tdis .<= NJQ[end,:DATE]))
tdis = tdis[tind]
cumdischarge = cumdischarge[tind]

# get NJQ lat, lon from CI station locations 
cistations = DataFrame(CSV.File("/home/timclements/CALI/CIstations.csv"))
NJQind = findfirst(cistations[!,:Station] .== "NJQ")
NJQlat = cistations[NJQind,:Latitude]
NJQlon = cistations[NJQind,:Longitude]

# load precip data
filename = "/media/FOUR/data/ppt.nc"
lon = ncread(filename,"lon")
lat = ncread(filename,"lat")
lonind = argmin(abs.(lon .- NJQlon))
latind = argmin(abs.(lat .- NJQlat))
precip = ncread(filename,"ppt",start=[latind,lonind,1],count=[1,1,-1])[1,1,:]
t = ncread(filename,"t")
t = Date.(Dates.unix2datetime.(t))

cumprecip = zeros(eltype(precip),length(precip))
octs = Date(2001,10,1):Year(1):Date(2020,10,1)
aprs = Date(2002,4,1):Year(1):Date(2021,4,1)
for dd = 1:length(octs)
    ind = findall((t .> octs[dd]) .& (t .< aprs[dd]))
    cumprecip[ind] .= cumsum(precip[ind])
end
cumprecip ./= 1000

# subset t with dv/v 
tind = findall((t .>= NJQ[1,:DATE]) .& (t .<= NJQ[end,:DATE]))
tppt = t[tind]
cumprecip = cumprecip[tind]

# find best fitting precip model against dv/v 
days = 365 : 15 * 365
N = length(days)
cors = zeros(N)
tminind = findfirst(t .== NJQ[1,:DATE])
tmaxind = findfirst(t .== NJQ[end,:DATE])
# find best CDM date range 
for ii in 1:N
    CDMii = CDM(precip[tminind - days[ii]:end],days[ii])
    tind = findall(in(NJQ[:,:DATE]),t[tminind - days[ii]:end])
    cors[ii] = cor(NJQ[:,:DVV],CDMii[tind])
end
bestday = days[argmin(cors)]
CDMbest = CDM(precip[tminind - bestday:end],bestday) ./ 1000
CDMbest = CDMbest[bestday:bestday+(tmaxind-tminind)]
CDMbest .-= CDMbest[1]
dvvprecipind = tminind:tmaxind
modelind = findall(in(NJQ[:,:DATE]),t[dvvprecipind])

# fit precip model to dv/v
@. model(x,p) = p[1] + p[2] * x
dvvprecipfit = curve_fit(
    model,
    CDMbest[modelind],
    NJQ[!,:DVV],
    [0.,-1.],
)

# fit ALI to dv/v 
dvvaliind = findall(in(ALI[:,:DATE]),NJQ[:,:DATE])
alidvvind = findall(in(NJQ[:,:DATE]),ALI[:,:DATE])
alidvvfit = curve_fit(
    model,
    NJQ[dvvaliind,:DVV],
    ALI[alidvvind,:WSE],
    [0,-1.],
)

# plot dv/v vs precip and water level 
NJQylims = (-1.5,1.5)
scatter(
    NJQ[!,:DATE],
    NJQ[!,:DVV],
    alpha=NJQ[!,:CC] ./ 10,
    label="",
    seriescolor=:Reds_9,
    marker_z=NJQ[!,:CC],
    ylims=NJQylims,
    colorbar=false,
    ylabel="                        dv/v [%]",
    right_margin = 2cm,
    left_margin = 3cm,
    yaxis=:flip,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    dpi=500,
    size=(800,400),
)
plot!(
    twinx(),
    ALI[!,:DATE],
    ALI[!,:WSE],
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    ylims=reverse(NJQylims .* coef(alidvvfit)[2] .+ coef(alidvvfit)[1]),
    seriescolor=:royalblue3,
    linewidth=2.5,
    alpha=0.85,
    label="",
    ylabel="\\Delta Alisal Reservoir Surface Elevation [m]",
    yguidefontcolor=:royalblue3,
    ytickfont=font(10,:royalblue3),
    yforeground_color_axis=:royalblue3,
    grid=false,
)
plot!(
    tppt, 
    cumprecip, 
    fillrange = 0.0, 
    fillcolor = :dodgerblue,
    fillalpha=0.3,
    inset = (1, bbox(0.0, 0.0, 1.0, 0.5, :bottom, :left)), 
    subplot = 3,
    label = "",
    grid=:off,
    bg_inside = nothing,
    ylabel="Cumulative Annual\n Rainfall [m]\n\n",
    # ylims = (0,2.5e8),
    left_margin = 5mm,
    yticks = (collect(0. : 0.5 : 1.5),["0.0        ","0.5        ","1.0        ","1.5        "]),
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    xtick=:off,
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    ylims=(0.,1.5)
)
savefig("/media/FOUR/data/FINAL-FIGURES/CINJQ-DVV-PRECIP.svg")
savefig("/media/FOUR/data/FINAL-FIGURES/CINJQ-DVV-PRECIP.png")

# plot dv/v vs precip and water level 
scatter(
    NJQ[!,:DATE],
    NJQ[!,:DVV],
    alpha=NJQ[!,:CC] ./ 10,
    label="",
    seriescolor=:Reds_9,
    marker_z=NJQ[!,:CC],
    ylims=NJQylims,
    colorbar=false,
    ylabel="                        dv/v [%]",
    right_margin = 2cm,
    left_margin = 3cm,
    yaxis=:flip,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    dpi=500,
    size=(800,400),
)
plot!(
    twinx(),
    ALI[!,:DATE],
    ALI[!,:WSE],
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    ylims=reverse(NJQylims .* coef(alidvvfit)[2] .+ coef(alidvvfit)[1]),
    seriescolor=:royalblue3,
    linewidth=2.5,
    alpha=0.85,
    label="",
    ylabel="\\Delta Alisal Reservoir Surface Elevation [m]",
    yguidefontcolor=:royalblue3,
    ytickfont=font(10,:royalblue3),
    yforeground_color_axis=:royalblue3,
    grid=false,
)
plot!(
    tdis, 
    cumdischarge, 
    fillrange = 100.0, 
    fillcolor = :dodgerblue,
    fillalpha=0.3,
    inset = (1, bbox(0.0, 0.0, 1.0, 0.5, :bottom, :left)), 
    subplot = 3,
    label = "",
    grid=:off,
    bg_inside = nothing,
    ylabel="Alisal Creek\nCumulative Annual\n Streamflow [m^3]\n\n",
    # ylims = (0,2.5e8),
    left_margin = 5mm,
    yticks = ([100,1e8,2e8,3e8],["0.0        ","1.0x10^{8}        ","2.0x10^{8}        ","3.0x10^{8}        "]),
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    xtick=:off,
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    ylims=(100.,3e8)
)
savefig("/media/FOUR/data/FINAL-FIGURES/CINJQ-DVV-STREAMFLOW.svg")
savefig("/media/FOUR/data/FINAL-FIGURES/CINJQ-DVV-STREAMFLOW.png")

# plot dv/v vs CDM
scatter(
    NJQ[!,:DATE],
    NJQ[!,:DVV],
    alpha=NJQ[!,:CC] ./ 10,
    label="",
    seriescolor=:Reds_9,
    marker_z=NJQ[!,:CC],
    ylims=NJQylims,
    colorbar=false,
    ylabel="                        dv/v [%]",
    right_margin = 2cm,
    left_margin = 3cm,
    yaxis=:flip,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    dpi=500,
    size=(800,400),
)
plot!(
    twinx(),
    t[dvvprecipind][modelind],
    CDMbest[modelind],
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    ylims=reverse((NJQylims .- coef(dvvprecipfit)[1]) ./ coef(dvvprecipfit)[2]),
    seriescolor=:royalblue3,
    linewidth=2.5,
    alpha=0.85,
    label="",
    ylabel="CDM [m]",
    yguidefontcolor=:royalblue3,
    ytickfont=font(10,:royalblue3),
    yforeground_color_axis=:royalblue3,
    grid=false,
)
plot!(
    tppt, 
    cumprecip, 
    fillrange = 0.0, 
    fillcolor = :dodgerblue,
    fillalpha=0.3,
    inset = (1, bbox(0.0, 0.0, 1.0, 0.5, :bottom, :left)), 
    subplot = 3,
    label = "",
    grid=:off,
    bg_inside = nothing,
    ylabel="Cumulative Annual\n Rainfall [m]\n\n",
    # ylims = (0,2.5e8),
    left_margin = 5mm,
    yticks = (collect(0. : 0.5 : 1.5),["0.0        ","0.5        ","1.0        ","1.5        "]),
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    xtick=:off,
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    ylims=(0.,1.5)
)
savefig("/media/FOUR/data/FINAL-FIGURES/CINJQ-DVV-CDM.svg")
savefig("/media/FOUR/data/FINAL-FIGURES/CINJQ-DVV-CDM.png")

# plot dv/v vs Baseflow model 
scatter(
    NJQ[!,:DATE],
    NJQ[!,:DVV],
    alpha=NJQ[!,:CC] ./ 10,
    label="",
    seriescolor=:Reds_9,
    marker_z=NJQ[!,:CC],
    ylims=NJQylims,
    colorbar=false,
    ylabel="                        dv/v [%]",
    right_margin = 2cm,
    left_margin = 3cm,
    yaxis=:flip,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    dpi=500,
    size=(800,400),
)
plot!(
    NJQ[!,:DATE],
    NJQ[!,:SSW],
    lw=3,
    c=:gold,
    alpha=0.85,
    label="dv/v Baseflow Model"
)
plot!(
    twinx(),
    ALI[!,:DATE],
    ALI[!,:WSE],
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    ylims=reverse(NJQylims .* coef(alidvvfit)[2] .+ coef(alidvvfit)[1]),
    seriescolor=:royalblue3,
    linewidth=2.5,
    alpha=0.85,
    label="",
    ylabel="\\Delta Alisal Reservoir Surface Elevation [m]",
    yguidefontcolor=:royalblue3,
    ytickfont=font(10,:royalblue3),
    yforeground_color_axis=:royalblue3,
    grid=false,
)
plot!(
    tppt, 
    cumprecip, 
    fillrange = 0.0, 
    fillcolor = :dodgerblue,
    fillalpha=0.3,
    inset = (1, bbox(0.0, 0.0, 1.0, 0.5, :bottom, :left)), 
    subplot = 3,
    label = "",
    grid=:off,
    bg_inside = nothing,
    ylabel="Cumulative Annual\n Rainfall [m]\n\n",
    # ylims = (0,2.5e8),
    left_margin = 5mm,
    yticks = (collect(0. : 0.5 : 1.5),["0.0        ","0.5        ","1.0        ","1.5        "]),
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    xtick=:off,
    xlims=(NJQ[1,:DATE],NJQ[end,:DATE]),
    ylims=(0.,1.5)
)
savefig("/media/FOUR/data/FINAL-FIGURES/CINJQ-DVV-SSW.svg")
savefig("/media/FOUR/data/FINAL-FIGURES/CINJQ-DVV-SSW.png")