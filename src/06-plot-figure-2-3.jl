using Arrow 
using CSV
using DataFrames 
using Dates 
using DSP
using Images
using Interpolations 
using JLD2
using LaTeXStrings
using LsqFit
using Plots
using Plots.PlotMeasures
using NetCDF
using SeisIO
using SeisNoise
using Serialization
using SpecialFunctions
using Statistics
import SeisNoise: smooth!, smooth

function CDM(A::AbstractArray,k::Int)
    Amean = [cumsum(A[1:k-1]) ./ (1:k-1) ;rolling_mean(A,k)]
    return cumsum(A .- Amean)
end

function rolling_mean(A::AbstractArray,k::Int)
    B = cumsum([zero(eltype(A)); A])
    return (B[k+1:end] - B[1:end-k]) ./ k
end


function elastic(precip,c,α,r,δt)
    Δprecip = precip .- mean(precip)
    P = zeros(size(precip,1))

    # compute P_i for each day 
    for ii in 1:length(P)
        for jj in 1:ii 
            P[ii] += Δprecip[jj] * α * erf(r / sqrt(4. * c * Float64(ii - jj) * δt ))
        end
    end
    return P
end

function talwani(precip,c,α,r,δt)
    Δprecip = precip .- mean(precip)           
    P = zeros(size(precip,1))

    # compute P_i for each day            
    for ii in 1:length(P)
        for jj in 1:ii 
            P[ii] += Δprecip[jj] * erfc(r / sqrt(4. * c * Float64(ii - jj) * δt )) + Δprecip[jj] * α * erf(r / sqrt(4. * c * Float64(ii - jj) * δt ))
        end
    end
    return P
end

function SSW06(precip,ϕ,a)
    GWL = zeros(length(precip))

    for ii in 1:length(GWL)
        for jj in 1:ii 
            GWL[ii] += precip[jj] / ϕ * exp(-a * (ii-jj))
        end
    end
    return GWL
end

function smooth!(df::DataFrame,col::Symbol;interval::Period=Day(1))
    Nrows = size(df,1)
    smooth_out = deepcopy(df[:,col])

    for ii in 1:Nrows
        firstind = findfirst(df[1:ii,:DATE] .+ interval .> df[ii,:DATE])
        smooth_out[ii] = mean(df[firstind:ii,col])
    end
    df[:,col] = smooth_out
    return nothing
end
smooth(df::DataFrame,col::Symbol;interval::Period=Day(1)) = (
    U = deepcopy(df); 
    smooth!(U,col,interval=interval);
    return U
)

"""
    smooth_withfiltfilt(A::AbstractArray; window_len::Int=7, window::Symbol=:bartlett)
Apply DSP.filtfilt() to smooth the waveform.
# Arguments
- "A::AbstractArray": periodogram of spectrum
- "window_len::Int=7": window length
- "window::Symbol=:bartlett": window type: (see https://juliadsp.github.io/DSP.jl/stable/windows/)
"""
function smooth_withfiltfilt(A::AbstractArray; window_len::Int=11, window::Symbol=:rect)
    w = getfield(DSP.Windows, window)(window_len)
    A = DSP.filtfilt(w ./ sum(w), A)
    return A
end


# read CI.LJR dv/v
LJR = Arrow.Table(joinpath(@__DIR__,"../data/DVV-90-DAY-COMP/2.0-4.0/CI.LJR.arrow")) |> DataFrame
mindate = Date(2000)

# remove days with instrument problems 
ind = findall(LJR[:,:CC] .>= 0.875)
LJR = LJR[ind,:]

# read nearest groundwater well and convert to m 
GWL = CSV.File(joinpath(@__DIR__,"../data/344614118454101.tsv"),comment="#",skipto=3) |> DataFrame
dropmissing!(GWL,:lev_va) # lev_va is ft below surface
GWL[!,:lev_va] ./= 3.28 # convert to meters 
GWL = GWL[GWL[!,:lev_dt] .>= mindate,:]
GWL[!,:lev_va] .*= -1
GWL[:,:lev_va] .-= mean(GWL[:,:lev_va])
MW1 = CSV.File(joinpath(@__DIR__,"../data/TRC-MW1.csv"),dateformat="Y/m/d") |> DataFrame
MW1[:,:GWL] .-= mean(MW1[:,:GWL])
MW1[:,:GWL] ./= 3.28

# get LJR lat, lon from CI station locations 
cistations = DataFrame(CSV.File(joinpath(@__DIR__,"../data/CIstations.csv")))
LJRind = findfirst(cistations[!,:Station] .== "LJR")
LJRlat = cistations[LJRind,:Latitude]
LJRlon = cistations[LJRind,:Longitude]

# get LJRN GPS data  
WNAMloc = DataFrame(
    CSV.File(
        joinpath(@__DIR__,"../data/WNAM-LOC.tsv"),
        delim=' ',
        ignorerepeated=true,
    )
)
LJRNind = findfirst(WNAMloc[!,:STATION] .== "LJRN")
LJRNlon = WNAMloc[LJRNind,:LON]
LJRNlat = WNAMloc[LJRNind,:LAT]
GPS = DataFrame(
    CSV.File(
        joinpath(@__DIR__,"../data/WNAMdetrend/ljrnFilterDetrend.neu"), 
        comment="#",
        header=[
            "DECDATE",
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
    )
)
GPS[:,:DATE] = Date.(GPS[:,:YEAR]) .+ Day.(GPS[:,:DAY] .- 1)
GPS = GPS[GPS[:,:DATE] .> Date(2013,11,6),:]
GPSLJRind = findall(in(LJR[:,:DATE]),GPS[:,:DATE])
LJRGPSind = findall(in(GPS[:,:DATE]),LJR[:,:DATE])
smooth!(GPS,:V,interval=Day(30))
smooth!(GPS,:E,interval=Day(30))
smooth!(GPS,:N,interval=Day(30))

# load GRACE data 
# from http://www2.csr.utexas.edu/grace/RL06_mascons.html
filename = joinpath(@__DIR__,"../data/GRACE-LJR.nc")
lwe = ncread(filename,"lwe_thickness")[1,1,:]
tlwe = Date(2002,1,1) .+ Day.(floor.(ncread(filename,"time")))
lwe .-= mean(lwe) # change baseline to entire signal 
lwe ./= 100 # convert to m 

# load precip data
filename = joinpath(@__DIR__,"../data/PRECIP-LJR.nc")
t = ncread(filename,"t")
tppt = Date.(Dates.unix2datetime.(t))
precip = ncread(filename,"ppt")[1,1,:]
precip ./= 1000

cumprecip = zeros(eltype(precip),length(precip))
octs = Date(1999,10,1):Year(1):Date(2020,10,1)
juns = Date(2000,6,1):Year(1):Date(2021,6,1)
for dd = 1:length(octs)
    ind = findall((tppt .> octs[dd]) .& (tppt .< juns[dd]))
    cumprecip[ind] .= cumsum(precip[ind])
end

# find best fitting precip model against dv/v 
days = 365 : 15 * 365
N = length(days)
cors = zeros(N)
tminind = findfirst(tppt .== LJR[1,:DATE])
tmaxind = findfirst(tppt .== LJR[end,:DATE])
# find best CDM date range 
for ii in 1:N
    CDMii = CDM(precip[tminind - days[ii]:end],days[ii])
    tind = findall(in(LJR[:,:DATE]),tppt[tminind - days[ii]:end])
    cors[ii] = cor(LJR[:,:DVV],CDMii[tind])
end
bestday = days[argmin(cors)]
CDMbest = CDM(precip[tminind - bestday:end],bestday) 
CDMbest = CDMbest[bestday:bestday+(tmaxind-tminind)]
# CDMbest .-= CDMbest[1]
dvvprecipind = tminind:tmaxind
modelind = findall(in(LJR[:,:DATE]),tppt[dvvprecipind])


## Elastic model of Talwani 
r = 500.
δt = 86400.
α = (1 + 0.27) / 3 / (1 - 0.27)
c = 0.0038 
PE = elastic(precip,c,α,r,δt)

## Fully-coupled model of Talwani 
cfc = 4.
PFC = talwani(precip,cfc,α,r,δt)

## model of Sens-Schonfelder and Wegler, 2006 
ϕ = 0.25   # porosity 
a = 0.0008  # exponential constant 
PSSW06 = SSW06(precip,ϕ,a)


# model of dv/v <-> groundwater 
g = 9.81 # m/s^2 
ρ = 1000 # kg / m^3 
ν = 0.27 # undrained poisson's ratio 
B = 1.0 # Skempton's coefficient 
G = 2e10 # Shear Modulus GPa
C = 1.5 * ρ * g / G / B * (1 - 2 * ν) / ( 1 + ν)
βs = [5e3 7.5e3 1e4]
S_y = 0.15 

# fit precip model to dv/v
@. model(x,p) = p[1] + p[2] * x
dvvprecipfit = curve_fit(
    model,
    CDMbest[modelind],
    LJR[!,:DVV],
    [0.,-1.],
)

# fit fully-coupled model to dv/v 
dvvPFCfit = curve_fit(
    model,
    PFC[dvvprecipind][modelind],
    LJR[!,:DVV],
    [0.,-1.],
)


# fit elastic model to dv/v 
dvvPEfit = curve_fit(
    model,
    PE[dvvprecipind][modelind],
    LJR[!,:DVV],
    [0.,-1.],
)

# fit exponential model to dv/v 
dvvSSWfit = curve_fit(
    model,
    PSSW06[dvvprecipind][modelind],
    LJR[!,:DVV],
    [0.,-1.],
)

# fit GRACE to dv/v 
dvvgracefit = curve_fit(
    model,
    lwe[findall(in(LJR[:,:DATE]),tlwe)],
    LJR[findall(in(tlwe),LJR[:,:DATE]),:DVV],
    [0,-.1],
)

# fit GRACE to precip 
precipgracefit = curve_fit(
    model,
    lwe[findall(in(tppt[dvvprecipind]),tlwe)] ./ 100,
    CDMbest[findall(in(tlwe),tppt[dvvprecipind])],
    [0,-.1],
)

# first create interpolated groundwater vector
tGWLdays = MW1[1,:DATE]:Day(1):MW1[end,:DATE]
tGWL = (MW1[:,:DATE] .- MW1[1,:DATE]) ./ Day(1) .+ 1 
tGWLinterp = tGWL[1]:tGWL[end]
interp_linear = LinearInterpolation(tGWL,MW1[:,:GWL])
GWLinterp = interp_linear.(tGWLinterp)
gracegwlind = findall(in(Date.(tlwe)),tGWLdays)
gwlgraceind = findall(in(tGWLdays),Date.(tlwe))

# now fit groundwater model to dv/v 
gwlprecipind = findall(in(tppt[dvvprecipind]),tGWLdays)
precipgwlind = findall(in(tGWLdays),tppt[dvvprecipind])
gwlprecipfit = curve_fit(
    model,
    CDMbest[precipgwlind],
    GWLinterp[gwlprecipind],
    [0,1.],
)

# fit dvv to gwl
gwldvvind = findall(in(LJR[:,:DATE]),tGWLdays)
dvvgwlind = findall(in(tGWLdays),LJR[:,:DATE])
dvvgwlfit = curve_fit(
    model,
    LJR[dvvgwlind,:DVV],
    GWLinterp[gwldvvind],
    [0,-1.],
)

# fit gwl to dvv
gwldvvfit = curve_fit(
    model,
    GWLinterp[gwldvvind],
    LJR[dvvgwlind,:DVV],
    [0,-1.],
)

# fit GPS to dvv
GPSLJRind = findall(in(LJR[:,:DATE]),GPS[:,:DATE])
LJRGPSind = findall(in(GPS[:,:DATE]),LJR[:,:DATE])
dvvgpsfit = curve_fit(
    model,
    GPS[GPSLJRind,:V],
    LJR[LJRGPSind,:DVV],
    [0,1.],
)

# fit GRACE to GWL
gracegwlfit = curve_fit(
    model,
    GWLinterp[gracegwlind],
    lwe[gwlgraceind],
    [0,1.],
)

###### CI-LJR DVV FIGURE #####
# plot dv/v first  
LJRylims = (-2.0,1.75)
LJRticks = ceil(LJR[1,:DATE],Year):Year(1):floor(LJR[end,:DATE],Year)
# scatter(
#     LJR[!,:DATE],
#     LJR[!,:DVV],
#     alpha=LJR[!,:CC],
#     label="",
#     seriescolor=:Reds_9,
#     marker_z=LJR[!,:CC],
#     markersize=2,
#     markerstrokewidth=0.1,
#     ylims=LJRylims,
#     xlims=(LJR[1,:DATE],LJR[end,:DATE]),
#     colorbar=false,
#     ylabel="dv/v [%]",
#     right_margin = 2cm,
#     left_margin = 1.5cm,
#     yaxis=:flip,
#     yguidefontcolor=:darkred,
#     ytickfont=font(10,:darkred),
#     ytick_direction=:out,
#     yforeground_color_axis=:darkred,
#     dpi=500,
#     size=(800,400)
# )
# plot scaled GRACE / GRACE-FO 
indbefore = findall(tlwe .< Date(2018))
indafter = findall(tlwe .> Date(2018))
# GRACE 2003-2017
plot(
    tlwe[indbefore],
    lwe[indbefore] .* coef(dvvgracefit)[2] .+ coef(dvvgracefit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    label="Scaled GRACE",
    c=:gold,
    linewidth=2,
    ylims=LJRylims,
    alpha=0.7,
    ylabel="dv/v [%]",
    right_margin = 2cm,
    left_margin = 1.5cm,
    yaxis=:flip,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
    dpi=500,
    size=(800,400),
)
# GRACE 2018- now 
plot!(
    tlwe[indafter],
    lwe[indafter] .* coef(dvvgracefit)[2] .+ coef(dvvgracefit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    label="",
    c=:gold,
    linewidth=2,
    ylims=LJRylims,
    alpha=0.7,
)
# scatter dv/v
scatter!(
    LJR[!,:DATE],
    LJR[!,:DVV],
    alpha=LJR[!,:CC],
    label="",
    seriescolor=:Reds_9,
    marker_z=LJR[!,:CC],
    markersize=3,
    markerstrokewidth=0.1,
    ylims=LJRylims,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    colorbar=false,
    yaxis=:flip,
)
# plot scaled CDM
plot!(
    tppt[dvvprecipind],
    PE[dvvprecipind] .* coef(dvvPEfit)[2] .+ coef(dvvPEfit)[1],
    c=:chartreuse,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    linewidth=2,
    alpha=0.7,
    # yguidefontcolor=:dodgerblue,
    # ytickfont=font(10,:dodgerblue),
    # yforeground_color_axis=:dodgerblue,
    # ylabel="Cumulative Deviation of\n Precipitation from Mean [m]",
    label="Scaled Elastic Model",
)

# plot groundwater level 
plot!(
    twinx(),
    MW1[!,:DATE],
    MW1[!,:GWL],
    sharex=true,
    ylim=reverse(LJRylims .* coef(dvvgwlfit)[2] .+ coef(dvvgwlfit)[1]),
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    linewidth=2,
    alpha=0.7,
    seriescolor=:royalblue3,
    label="",
    ylabel="\\Delta Groundwater Level [m]",
    grid=:off,
    yguidefontcolor=:royalblue3,
    ytickfont=font(10,:royalblue3),
    yforeground_color_axis=:royalblue3,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
# plot cumulative precipiation 
# subset t with dv/v 
cumind = findall((tppt .>= LJR[1,:DATE]) .& (tppt .<= LJR[end,:DATE]))
plot!(
    Date.(tppt[cumind]), 
    cumprecip[cumind], 
    fillrange = 0.0, 
    fillcolor = :dodgerblue,
    fillalpha=0.3,
    inset = (1, bbox(0.0, 0.0, 1.0, 0.33, :bottom, :left)), 
    subplot = 3,
    label = "",
    grid=:off,
    bg_inside = nothing,
    ylabel="Cumulative Annual\n Precipitation [m]\n\n",
    ylims = (0,1),
    left_margin = 5mm,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    yticks = (collect(0. : 0.2 : 1.),["0.0    ","0.2    ","0.4    ","0.6    ","0.8    ","1.0    "]),
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-DVV.png"))

##### DV/V - GRACE #####
# plot dv/v and GRACE
scatter(
    Date.(LJR[!,:DATE]),
    LJR[!,:DVV],
    alpha=LJR[!,:CC],
    label="",
    seriescolor=:Reds_9,
    marker_z=LJR[!,:CC],
    markersize=2,
    markerstrokewidth=0.1,
    ylims=LJRylims,
    colorbar=false,
    ylabel="dv/v [%]",
    right_margin = 2cm,
    left_margin = 1.5cm,
    yaxis=:flip,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    dpi=500,
    size=(800,400),
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
plot!(
    twinx(),
    tlwe[indbefore],
    lwe[indbefore],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    label="",
    c=:dodgerblue2,
    linewidth=2.5,
    ylabel="Liquid Water Equivalent [m]",
    ylims= reverse((LJRylims .- coef(dvvgracefit)[1]) ./ coef(dvvgracefit)[2]),
    yguidefontcolor=:dodgerblue2,
    ytickfont=font(10,:dodgerblue2),
    ytick_direction=:out,
    yforeground_color_axis=:dodgerblue2,
    grid=false,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
plot!(
    twinx(),
    tlwe[indafter],
    lwe[indafter],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    label="",
    c=:dodgerblue2,
    linewidth=2.5,
    ylabel="Liquid Water Equivalent [m]",
    ylims= reverse((LJRylims .- coef(dvvgracefit)[1]) ./ coef(dvvgracefit)[2]),
    yguidefontcolor=:dodgerblue2,
    ytickfont=font(10,:dodgerblue2),
    ytick_direction=:out,
    yforeground_color_axis=:dodgerblue2,
    grid=false,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-GRACE-DVV.png"))

##### DV/V - Fully Coupled #####
scatter(
    Date.(LJR[!,:DATE]),
    LJR[!,:DVV],
    alpha=LJR[!,:CC],
    label="",
    seriescolor=:Reds_9,
    marker_z=LJR[!,:CC],
    markersize=2,
    markerstrokewidth=0.1,
    ylims=LJRylims,
    colorbar=false,
    ylabel="dv/v [%]",
    right_margin = 2cm,
    left_margin = 1.5cm,
    yaxis=:flip,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    dpi=500,
    size=(800,400),
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
plot!(
    tppt[dvvprecipind],
    CDMbest .* coef(dvvprecipfit)[2] .+ coef(dvvprecipfit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    c=:chartreuse,
    linewidth=1.5,
    alpha=0.75,
    label="CDM",
    ylims= LJRylims,
    # yguidefontcolor=:dodgerblue2,
    # ytickfont=font(10,:dodgerblue2),
    # ytick_direction=:out,
    # yforeground_color_axis=:dodgerblue2,
    # grid=false,
)
plot!(
    tppt[dvvprecipind],
    PFC[dvvprecipind] .* coef(dvvPFCfit)[2] .+ coef(dvvPFCfit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    c=:orange,
    linewidth=1.5,
    alpha=0.75,
    label="Fully-Coupled Model",
    ylims= LJRylims,
    # yguidefontcolor=:dodgerblue2,
    # ytickfont=font(10,:dodgerblue2),
    # ytick_direction=:out,
    # yforeground_color_axis=:dodgerblue2,
    # grid=false,
)
plot!(
    tppt[dvvprecipind],
    PSSW06[dvvprecipind] .* coef(dvvSSWfit)[2] .+ coef(dvvSSWfit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    c=:dodgerblue2,
    linewidth=1.5,
    alpha=0.75,
    label="SSW06 Model",
    ylims= LJRylims,
    # yguidefontcolor=:dodgerblue2,
    # ytickfont=font(10,:dodgerblue2),
    # ytick_direction=:out,
    # yforeground_color_axis=:dodgerblue2,
    # grid=false,
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-FC-MODELS.png"))

##### DV/V - CDM #####
scatter(
    Date.(LJR[!,:DATE]),
    LJR[!,:DVV],
    alpha=LJR[!,:CC],
    label="",
    seriescolor=:Reds_9,
    marker_z=LJR[!,:CC],
    markersize=2,
    markerstrokewidth=0.1,
    ylims=LJRylims,
    colorbar=false,
    ylabel="dv/v [%]",
    right_margin = 2cm,
    left_margin = 1.5cm,
    yaxis=:flip,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    dpi=500,
    size=(800,400),
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
plot!(
    tppt[dvvprecipind],
    CDMbest .* coef(dvvprecipfit)[2] .+ coef(dvvprecipfit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    c=:chartreuse,
    linewidth=1.5,
    alpha=0.75,
    label="CDM",
    ylims= LJRylims,
    # yguidefontcolor=:dodgerblue2,
    # ytickfont=font(10,:dodgerblue2),
    # ytick_direction=:out,
    # yforeground_color_axis=:dodgerblue2,
    # grid=false,
)
plot!(
    tppt[dvvprecipind],
    PE[dvvprecipind] .* coef(dvvPEfit)[2] .+ coef(dvvPEfit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    c=:orange,
    linewidth=1.5,
    alpha=0.75,
    label="Elastic Model",
    ylims= LJRylims,
    # yguidefontcolor=:dodgerblue2,
    # ytickfont=font(10,:dodgerblue2),
    # ytick_direction=:out,
    # yforeground_color_axis=:dodgerblue2,
    # grid=false,
)
plot!(
    tppt[dvvprecipind],
    PSSW06[dvvprecipind] .* coef(dvvSSWfit)[2] .+ coef(dvvSSWfit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    c=:dodgerblue2,
    linewidth=1.5,
    alpha=0.75,
    label="SSW06 Model",
    ylims= LJRylims,
    # yguidefontcolor=:dodgerblue2,
    # ytickfont=font(10,:dodgerblue2),
    # ytick_direction=:out,
    # yforeground_color_axis=:dodgerblue2,
    # grid=false,
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-MODELS.png"))

# hysteresis for two time periods 
scatter(
    LJR[:,:DVV][1997:5000],
    CDMbest[modelind][1997:5000],
    alpha=0.75,
    markerz=(ind[1997:5000] .- ind[1997]) ./365.25 .+ 2008.5,
    ylabel="CDMk [m]",
    xlabel="dv/v [%]",
    right_margin=1.5cm,
    left_margin=1.5cm,
    label="",
    c=:vik,
    dpi=250,
    ylims=(-0.4,0.3),
    xlims=(-0.9,1.2),
)
annotate!((-0.49,-0.22,text("Drought",14,:left,:top,:indianred)))
annotate!((0.07,0.18,text("High groundwater\n         level",14,:left,:top,:dodgerblue)))
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-DVV-CDM-HYST-1"))

scatter(
    LJR[:,:DVV][3993:end-100],
    CDMbest[modelind][3993:end-100],
    alpha=0.75,
    markerz=(ind[3993:end-100] .- ind[3993]) ./365.25 .+ 2014.,
    ylabel="CDMk [m]",
    xlabel="dv/v [%]",
    right_margin=1.5cm,
    left_margin=1.5cm,
    label="",
    c=cgrad(:vik,rev=true,scale=:log),
    dpi=250,
    ylims=(-0.4,0.3),
    xlims=(-0.9,1.2),
)
annotate!((0.24,0.12,text("   Increasing\n groundwater\n        level",14,:left,:top,:dodgerblue)))
annotate!((-0.49,-0.22,text("Drought",14,:left,:top,:indianred)))
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-DVV-CDM-HYST-2"))

# plot predicted dv/v 
scatter(
    Date.(LJR[!,:DATE]),
    LJR[!,:DVV],
    alpha=LJR[!,:CC],
    label="",
    seriescolor=:Reds_9,
    marker_z=LJR[!,:CC],
    markersize=2,
    markerstrokewidth=0.1,
    ylims=LJRylims,
    colorbar=false,
    ylabel="dv/v [%]",
    right_margin = 3cm,
    left_margin = 1.5cm,
    yaxis=:flip,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    dpi=500,
    size=(800,400),
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
plot!(
    tppt[dvvprecipind],
    CDMbest .* coef(dvvprecipfit)[2] .+ coef(dvvprecipfit)[1],
    c=:dodgerblue,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    linewidth=2,
    alpha=0.85,
    # yguidefontcolor=:dodgerblue,
    # ytickfont=font(10,:dodgerblue),
    # yforeground_color_axis=:dodgerblue,
    # ylabel="Cumulative Deviation of\n Precipitation from Mean [m]",
    label="Predicted dv/v",
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-PREDICT-DVV.png"))

# plot dvv vs precip 
startdate = (LJR[1,:DATE] - Date(year(LJR[1,:DATE]))).value / 365.25 + year(LJR[1,:DATE])
cticks = 2003:2:2021
scatter(
    CDMbest[findall(in(LJR[:,:DATE]),
    tppt[tminind:tmaxind])],
    LJR[:,:DVV],
    alpha=LJR[:,:CC] ./ 3,
    ylims=LJRylims,
    label="",
    zcolor=(datetime2unix.(DateTime.(LJR[:,:DATE])) .- datetime2unix(DateTime(LJR[1,:DATE]))) ./ 86400 ./ 365.25 .+ startdate,
    xlabel="CDM [m]",
    ylabel="dv/v [%]",
    right_margin = 1cm,
    dpi=250,
    colorbar_ticks=string.(cticks),
)
pt = -0.4:0.1:0.65
plot!(
    pt,
    pt .* coef(dvvprecipfit)[2] .+ coef(dvvprecipfit)[1],
    color=:red,
    linewidth=2,
    linestyle=:dash,
    label="Linear Model",
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-HYST.png"))

# plot predicted dv/v and hysteresis 
p1 = scatter(
    Date.(LJR[!,:DATE]),
    LJR[!,:DVV],
    alpha=LJR[!,:CC] ./ 10,
    label="",
    # seriescolor=:Reds_9,
    zcolor=(LJR[:,:DATE] .- LJR[1,:DATE])./Day(1),
    ylims=LJRylims,
    colorbar=false,
    ylabel="dv/v [%]",
    legend=:topleft,
    dpi=500,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
plot!(
    p1,
    tppt[dvvprecipind],
    CDMbest .* coef(dvvprecipfit)[2] .+ coef(dvvprecipfit)[1],
    c=:dodgerblue,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    linewidth=2,
    alpha=0.85,
    # yguidefontcolor=:dodgerblue,
    # ytickfont=font(10,:dodgerblue),
    # yforeground_color_axis=:dodgerblue,
    # ylabel="Cumulative Deviation of\n Precipitation from Mean [m]",
    label="Scaled CDM",
)

# plot dvv vs precip 
Ndvv = size(LJR,1)
Nstep = 2
p2 = scatter(
    CDMbest[findall(
        in(LJR[:,:DATE]),
        tppt[tminind:tmaxind],
    )][1:Nstep:Ndvv],
    LJR[:,:DVV][1:Nstep:Ndvv],
    alpha=LJR[:,:CC][1:Nstep:Ndvv] ./ 3,
    ylims=LJRylims,
    label="",
    zcolor=((LJR[:,:DATE] .- LJR[1,:DATE])./Day(1))[1:Nstep:Ndvv],
    xlabel="CDM [m]",
    ylabel="dv/v [%]",
    colorbar=false,
    dpi=500,
)
pt = -0.4:0.1:0.65
plot!(
    p2,
    pt,
    pt .* coef(dvvprecipfit)[2] .+ coef(dvvprecipfit)[1],
    color=:dodgerblue,
    linewidth=2.5,
    linestyle=:dash,
    label="Linear Model",
)
l = @layout [a b]
plot(p1,p2,layout=l,size=(800,300),margin=5mm)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-DVV-HYST.png"))

# plot predicted groundwater 
cs = palette(:Reds_3,6)
p = plot(    
    dpi=500,
    size=(800,400),
    right_margin = 2cm,
    left_margin = 1.5cm,
)
for ii = 1:length(βs)
    exponent = floor(Int,log10(βs[ii]))
    value = βs[ii] / 10^exponent
    if value % 1 == 0 
        value = Int(value)
    end
    plot!(
        p,
        LJR[:,:DATE],
        -LJR[:,:DVV]  ./ C ./ βs[ii] ./ 100,
        # c=:dodgerblue,
        alpha=0.85,
        xlims=(LJR[1,:DATE],LJR[end,:DATE]),
        linewidth=1,
        label="\\beta = -$value x 10^$exponent",
        c=cs[ii * 2],
        xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
        xrot=60,
    )
end
plot!(
    MW1[!,:DATE],
    MW1[!,:GWL],
    c=:dodgerblue,
    label="GWL measured",
    ylabel="Measured \\Delta Groundwater Level [m]",
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    linewidth=2,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
display(p)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-PREDICT-GWL.png"))

#### plot groundwater and GRACE ####
plot(
    tlwe[indbefore],
    lwe[indbefore],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    label="",
    c=:gold,
    linewidth=2.5,
    ylabel="Liquid Water Equivalent [m]",
    # ylims= reverse((LJRylims .- coef(gracegwlfit)[1]) ./ coef(gracegwlfit)[2]),
    yguidefontcolor=:gold,
    ytickfont=font(10,:gold),
    ytick_direction=:out,
    yforeground_color_axis=:gold,
    dpi=250,
    right_margin=2cm,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
plot!(
    tlwe[indafter],
    lwe[indafter],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    label="",
    c=:gold,
    linewidth=2.5,
    ylabel="Liquid Water Equivalent [m]",
    # ylims= reverse((LJRylims .- coef(gracegwlfit)[1]) ./ coef(gracegwlfit)[2]),
    yguidefontcolor=:gold,
    ytickfont=font(10,:gold),
    ytick_direction=:out,
    yforeground_color_axis=:gold,
)
plot!(
    twinx(),
    MW1[!,:DATE],
    MW1[!,:GWL],
    c=:dodgerblue,
    label="",
    ylabel="\\Delta Groundwater Level [m]",
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    # ylims=reverse((LJRylims .- coef(gracegwlfit)[1]) ./ coef(gracegwlfit)[2]),
    linewidth=3,
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    grid=false,
    alpha=0.75,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-GRACE-GWL.png"))

## plot dv/v + GWL ## 
scatter(
    LJR[!,:DATE],
    LJR[!,:DVV],
    alpha=LJR[!,:CC],
    label="",
    seriescolor=:Reds_9,
    marker_z=LJR[!,:CC],
    markersize=3,
    markerstrokewidth=0.1,
    ylims=LJRylims,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    colorbar=false,
    yaxis=:flip,
    ylabel="                    dv/v [%]",
    right_margin = 2.5cm,
    left_margin = 2.5cm,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    dpi=500,
    size=(800,400),
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
# plot groundwater level 
plot!(
    twinx(),
    MW1[!,:DATE],
    MW1[!,:GWL],
    sharex=true,
    ylims=reverse(LJRylims .* coef(dvvgwlfit)[2] .+ coef(dvvgwlfit)[1]) ,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    linewidth=4,
    alpha=0.8,
    seriescolor=:royalblue3,
    label="",
    ylabel="\\Delta Groundwater Level [m]",
    grid=:off,
    yguidefontcolor=:royalblue3,
    ytickfont=font(10,:royalblue3),
    yforeground_color_axis=:royalblue3,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-MW1.png"))

#### plot autocorrletion for LJR ####
NZ = deserialize(joinpath(@__DIR__,"../data/ONECORR/CI.LJR..BHN.CI.LJR..BHZ"))
remove_nan!(NZ)
ind = setdiff(1:length(NZ.t),5297:5359)
NZ = NZ[ind]
clean_up!(NZ,2.,4.)
SeisNoise.smooth!(NZ,Day(10))
# abs_max!(C)

# load PSD data 
@load joinpath(@__DIR__,"../data/LJR-psd.jld2") psd tpsd fpsd
psd_range = tpsd[1] : Day(1) : tpsd[end]
PSD = zeros(eltype(psd),size(psd,1),size(psd_range,1))
psdind = findall(in(tpsd),psd_range)
for ii in 1:size(psdind,1)
    PSD[:,psdind[ii]] = psd[:,ii]
end
# remove sensor problems 
psdind = findall(median(PSD,dims=1)[:] .> 2e-14)
PSD[:,psdind] .= 0
PSD[iszero.(PSD)] .= NaN

# custom corelation plotting based SeisNoise recipe 
minlag = 1.5 
maxlag = 10.
lags = -NZ.maxlag:1/NZ.fs:NZ.maxlag
posind = findall( minlag .<= lags .<= maxlag)
negind = findall( -maxlag .<= lags .<= -minlag)
Cdates = Date.(round.(u2d.(NZ.t),Dates.Day))
date_range = Cdates[1] : Day(1) : Cdates[end]
times = Dates.format.(date_range,"yyyy/m/d")
Apos = zeros(eltype(NZ.corr),size(posind,1),size(date_range,1))
Aneg = zeros(eltype(NZ.corr),size(negind,1),size(date_range,1))
mapind = findall(in(Cdates),date_range)
# for ii in 1:size(mapind,1)
#     A[:,mapind[ii]] = (reverse(C.corr[negind,ii]) .+ C.corr[posind,ii]) ./ 2
# end
for ii in 1:size(mapind,1)
    Apos[:,mapind[ii]] =  NZ.corr[posind,ii]
end
for ii in 1:size(mapind,1)
    Aneg[:,mapind[ii]] =  reverse(NZ.corr[negind,ii])
end

# scale by t^2
Apos .*= lags[posind] 
Aneg .*= lags[posind] 

# normalize 
Apos ./= maximum(abs.(Apos))

# plot spetrogram, correlations and dv/v 
p1 = heatmap(
    psd_range,
    fpsd, 
    10 .* log10.(PSD),
    ylabel = "Frequency [Hz]",
    # legend=:none,
    colorbar_title=" \n\n\n  10 * log10 m^2/s^4/Hz dB",
    clims=(-170,-130),
    xlims=(Cdates[1],tpsd[end]),
    grid=false,
    left_margin=5mm,
    right_margin=5mm,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
hline!(
    [2.],
    color=:black,
    linestyle=:dash,
    xlims=(Cdates[1],tpsd[end]),
    linewidth=1.5,
    label="",
)
hline!(
    [4.],
    color=:black,
    linestyle=:dash,
    xlims=(Cdates[1],tpsd[end]),
    linewidth=1.5,
    label="",
)
annotate!((Date(2001,3,1),9,text("A",20)))
p2 = heatmap(
    date_range,
    lags[posind],
    Apos,
    c=cgrad(:diverging_bwr_40_95_c42_n256,rev=true),
    # legend=:none,
    xlims=(Cdates[1],tpsd[end]),
    ylabel = "Lag [s]",
    left_margin=5mm,
    colorbar_title=" \n\n\n  Amplitude",
    # clims=(-0.025,0.025),
    clims=(-0.5,0.5),
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
hline!(
    [2.],
    color=:black,
    linestyle=:dash,
    xlims=(Cdates[1],tpsd[end]),
    linewidth=1.5,
    label="",)
hline!(
    [8.],
    color=:black,
    linestyle=:dash,
    xlims=(Cdates[1],tpsd[end]),
    linewidth=1.5,
    label="",
)
annotate!((Date(2001,3,1),9,text("B",20)))
plot(p1,p2,layout = grid(2,1,heights=[0.5,0.5]), 
    dpi=500,
    size=(800,600),
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-CORR.png"))

## MAP + DVV  
# load the map from GMT 
img = load(joinpath(@__DIR__,"../data/FINAL-FIGURES/LJR-map.png"))
p1 = plot(img,border=:none,dpi=250,margins = -2Plots.cm)
annotate!((-100,100,text("A",20)))

indbefore = findall(tlwe .< Date(2018))
indafter = findall(tlwe .> Date(2018))
# scatter dv/v
p2 = scatter(
    LJR[!,:DATE],
    LJR[!,:DVV],
    alpha=LJR[!,:CC],
    label="",
    seriescolor=:Reds_9,
    marker_z=LJR[!,:CC],
    markersize=3,
    markerstrokewidth=0.1,
    ylims=LJRylims,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    colorbar=false,
    yaxis=:flip,
    ylabel="                    dv/v [%]",
    right_margin = 2.5cm,
    left_margin = 2.5cm,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    legendfontsize=10,
    dpi=500,
    size=(800,400),
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
    margins=3Plots.mm
)
# GRACE 2003-2017
plot!(
    tlwe[indbefore],
    lwe[indbefore] .* coef(dvvgracefit)[2] .+ coef(dvvgracefit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    label="Scaled GRACE",
    c=:gold,
    linewidth=3,
    ylims=LJRylims,
    alpha=0.9,
)
# GRACE 2018- now 
plot!(
    tlwe[indafter],
    lwe[indafter] .* coef(dvvgracefit)[2] .+ coef(dvvgracefit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    label="",
    c=:gold,
    linewidth=3,
    ylims=LJRylims,
    alpha=0.9,
)
# plot GPS 
# plot!(
#     GPS[GPSLJRind,:DATE],
#     GPS[GPSLJRind,:V] .* coef(dvvgpsfit)[2] .+ coef(dvvgpsfit)[1], 
#     label="Scaled GPS",
#     c=:purple,
#     linewidth=3,
#     ylims=LJRylims,
#     alpha=0.9,
# )
annotate!((Date(2000,6,1),-1.85,text("B",20)))
# plot scaled Elastic model 
plot!(
    tppt[dvvprecipind],
    PE[dvvprecipind] .* coef(dvvPEfit)[2] .+ coef(dvvPEfit)[1],
    c=:chartreuse,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    linewidth=3,
    alpha=0.7,
    # yguidefontcolor=:dodgerblue,
    # ytickfont=font(10,:dodgerblue),
    # yforeground_color_axis=:dodgerblue,
    # ylabel="Cumulative Deviation of\n Precipitation from Mean [m]",
    label="Scaled Elastic Model",
)
# plot groundwater level 
plot!(
    twinx(),
    MW1[!,:DATE],
    MW1[!,:GWL],
    sharex=true,
    ylim=reverse(LJRylims .* coef(dvvgwlfit)[2] .+ coef(dvvgwlfit)[1]),
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    linewidth=3,
    alpha=0.7,
    seriescolor=:royalblue3,
    label="",
    ylabel="\\Delta Groundwater Level [m]",
    grid=:off,
    yguidefontcolor=:royalblue3,
    ytickfont=font(10,:royalblue3),
    yforeground_color_axis=:royalblue3,
    # xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xticks=false,
    xrot=60,
)
# plot cumulative precipiation 
cumind = findall((tppt .>= LJR[1,:DATE]) .& (tppt .<= LJR[end,:DATE]))
plot!(
    Date.(tppt[cumind]), 
    cumprecip[cumind], 
    fillrange = 0.0, 
    fillcolor = :dodgerblue,
    fillalpha=0.3,
    inset = (1, bbox(0.0, 0.0, 1.0, 0.45, :bottom, :left)), 
    subplot = 3,
    label = "",
    grid=:off,
    bg_inside = nothing,
    ylabel="Cumulative Annual\n Precipitation [m]\n\n",
    ylims = (0,1),
    left_margin = 5mm,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    yticks = (collect(0. : 0.2 : 1.),["0.0     ","0.2     ","0.4     ","0.6     ","0.8     ","1.0     "]),
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    # xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xticks = false,
    xrot=60,
)
l1 = @layout [a b{0.6w}]
plot(p1,p2,layout=l1,size=(1200,400),dpi=500)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-MAP-DVV.png"))
 

#### MULTI-FIGURE #### 
p1 = heatmap(
    psd_range,
    fpsd, 
    10 .* log10.(PSD),
    ylabel = "Frequency [Hz]",
    # legend=:none,
    colorbar_title=" \n\n\n     10 * log10 m^2/s^4/Hz dB",
    clims=(-170,-130),
    xlims=(Cdates[1],tpsd[end]),
    grid=false,
    left_margin=5mm,
    right_margin=5mm,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
hline!(
    [2.],
    color=:black,
    linestyle=:dash,
    xlims=(Cdates[1],tpsd[end]),
    linewidth=1.5,
    label="",
)
hline!(
    [4.],
    color=:black,
    linestyle=:dash,
    xlims=(Cdates[1],tpsd[end]),
    linewidth=1.5,
    label="",
)
annotate!((Date(2001,1,1),9,text("B",20)))
p2 = heatmap(
    date_range,
    lags[posind],
    Apos,
    c=cgrad(:diverging_bwr_40_95_c42_n256,rev=true),
    # legend=:none,
    xlims=(Cdates[1],tpsd[end]),
    colorbar_title=" \n\n\n     Amplitude",
    ylabel = "Lag [s]",
    left_margin=5mm,
    clims=(-0.025,0.025),
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
hline!(
    [2.],
    color=:black,
    linestyle=:dash,
    xlims=(Cdates[1],tpsd[end]),
    linewidth=1.5,
    label="",)
hline!(
    [8.],
    color=:black,
    linestyle=:dash,
    xlims=(Cdates[1],tpsd[end]),
    linewidth=1.5,
    label="",
)
annotate!((Date(2001,1,1),9,text("C",20)))
indbefore = findall(tlwe .< Date(2018))
indafter = findall(tlwe .> Date(2018))
# scatter dv/v
p3 = scatter(
    LJR[!,:DATE],
    LJR[!,:DVV],
    alpha=LJR[!,:CC],
    label="",
    seriescolor=:Reds_9,
    marker_z=LJR[!,:CC],
    markersize=3,
    markerstrokewidth=0.1,
    ylims=LJRylims,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    colorbar=false,
    yaxis=:flip,
    ylabel="                    dv/v [%]",
    right_margin = 2.5cm,
    left_margin = 2.5cm,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    ytick_direction=:out,
    yforeground_color_axis=:darkred,
    legendfontsize=16,
    dpi=500,
    size=(800,400),
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
# GRACE 2003-2017
plot!(
    tlwe[indbefore],
    lwe[indbefore] .* coef(dvvgracefit)[2] .+ coef(dvvgracefit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    label="Scaled GRACE",
    c=:gold,
    linewidth=4,
    ylims=LJRylims,
    alpha=0.9,
)
# GRACE 2018- now 
plot!(
    tlwe[indafter],
    lwe[indafter] .* coef(dvvgracefit)[2] .+ coef(dvvgracefit)[1],
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    label="",
    c=:gold,
    linewidth=4,
    ylims=LJRylims,
    alpha=0.9,
)
# plot GPS 
# plot!(
#     GPS[GPSLJRind,:DATE],
#     GPS[GPSLJRind,:V] .* coef(dvvgpsfit)[2] .+ coef(dvvgpsfit)[1], 
#     label="Scaled GPS",
#     c=:purple,
#     linewidth=3,
#     ylims=LJRylims,
#     alpha=0.9,
# )
annotate!((Date(2001,1,1),-2.25,text("D",20)))
# plot scaled Elastic model 
plot!(
    tppt[dvvprecipind],
    PE[dvvprecipind] .* coef(dvvPEfit)[2] .+ coef(dvvPEfit)[1],
    c=:chartreuse,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    linewidth=3,
    alpha=0.7,
    # yguidefontcolor=:dodgerblue,
    # ytickfont=font(10,:dodgerblue),
    # yforeground_color_axis=:dodgerblue,
    # ylabel="Cumulative Deviation of\n Precipitation from Mean [m]",
    label="Scaled Elastic Model",
)
# plot groundwater level 
plot!(
    twinx(),
    MW1[!,:DATE],
    MW1[!,:GWL],
    sharex=true,
    ylim=reverse(LJRylims .* coef(dvvgwlfit)[2] .+ coef(dvvgwlfit)[1]),
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    linewidth=3,
    alpha=0.7,
    seriescolor=:royalblue3,
    label="",
    ylabel="\\Delta Groundwater Level [m]",
    grid=:off,
    yguidefontcolor=:royalblue3,
    ytickfont=font(10,:royalblue3),
    yforeground_color_axis=:royalblue3,
    # xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xticks=false,
    xrot=60,
)
# plot cumulative precipiation 
cumind = findall((tppt .>= LJR[1,:DATE]) .& (tppt .<= LJR[end,:DATE]))
plot!(
    Date.(tppt[cumind]), 
    cumprecip[cumind], 
    fillrange = 0.0, 
    fillcolor = :dodgerblue,
    fillalpha=0.3,
    inset = (1, bbox(0.0, 0.0, 1.0, 0.45, :bottom, :left)), 
    subplot = 3,
    label = "",
    grid=:off,
    bg_inside = nothing,
    ylabel="Cumulative Annual\n Precipitation [m]\n\n",
    ylims = (0,1),
    left_margin = 5mm,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    yticks = (collect(0. : 0.2 : 1.),["0.0     ","0.2     ","0.4     ","0.6     ","0.8     ","1.0     "]),
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    # xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xticks = false,
    xrot=60,
)

# plot(p1,p2,p3,layout = grid(3,1,heights=[0.3,0.3,0.4]), 
#     dpi=500,
#     size=(800,800),
# )
# load the map from GMT 
img = load(joinpath(@__DIR__,"../data/FINAL-FIGURES/LJR-map.png"))
p4 = plot(img,border=:none)
annotate!((-100,100,text("A",20)))

l1 = @layout [ [a [b{0.95w,0.5h}; c{0.9w, 0.5h}]]; d{1.0w, 0.5h}]
plot(p4,p1,p2,p3,layout=l1,size=(1600,1600),dpi=250)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-ALL.png"))
