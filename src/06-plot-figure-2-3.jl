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
LJR = Arrow.Table(joinpath(@__DIR__,"../data/FIT-DVV-SSE/90-DAY/CI.LJR.arrow")) |> DataFrame
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
dvvprecipind = tminind:tmaxind
modelind = findall(in(LJR[:,:DATE]),tppt[dvvprecipind])

@. model(x,p) = p[1] + p[2] * x
dvvCDMfit = curve_fit(
    model,
    LJR[!,:CDM],
    LJR[!,:DVV],
    [0.,-1.],
)

# fit fully-coupled model to dv/v 
dvvFCfit = curve_fit(
    model,
    LJR[!,:FC],
    LJR[!,:DVV],
    [0.,-1.],
)

# fit elastic model to dv/v 
dvvPEfit = curve_fit(
    model,
    LJR[:,:ELASTIC],
    LJR[!,:DVV],
    [0.,-1.],
)

# fit exponential model to dv/v 
dvvSSWfit = curve_fit(
    model,
    LJR[!,:SSW],
    LJR[!,:DVV],
    [0.,-1.],
)

# fit drained model to dv/v 
dvvDRAINfit = curve_fit(
    model,
    LJR[!,:DRAINDED], # typo 
    LJR[!,:DVV],
    [0.,-1.],
)


# model of dv/v <-> groundwater 
g = 9.81 # m/s^2 
ρ = 1000 # kg / m^3 
ν = 0.27 # undrained poisson's ratio 
B = 1.0 # Skempton's coefficient 
G = 2e10 # Shear Modulus GPa
C = 1.5 * ρ * g / G / B * (1 - 2 * ν) / ( 1 + ν)
βs = [5e3 7.5e3 1e4]
S_y = 0.15 

# first create interpolated groundwater vector
tGWLdays = MW1[1,:DATE]:Day(1):MW1[end,:DATE]
tGWL = (MW1[:,:DATE] .- MW1[1,:DATE]) ./ Day(1) .+ 1 
tGWLinterp = tGWL[1]:tGWL[end]
interp_linear = LinearInterpolation(tGWL,MW1[:,:GWL])
GWLinterp = interp_linear.(tGWLinterp)
gracegwlind = findall(in(Date.(tlwe)),tGWLdays)
gwlgraceind = findall(in(tGWLdays),Date.(tlwe))

###### CI-LJR DVV FIGURE #####
# plot dv/v first  
LJRylims = (-2.0,1.75)
LJRticks = ceil(LJR[1,:DATE],Year):Year(1):floor(LJR[end,:DATE],Year)

#### plot autocorrletion for LJR ####
NZ = deserialize(joinpath(@__DIR__,"../data/ONECORR/CI.LJR..BHN.CI.LJR..BHZ"))
remove_nan!(NZ)
ind = setdiff(1:length(NZ.t),5297:5359)
NZ = NZ[ind]
clean_up!(NZ,2.,4.)
SeisNoise.smooth!(NZ,Day(10))
# abs_max!(C)

# load PSD data 
psdfile = joinpath(@__DIR__,"../data/LJR-psd.jld2")
@load psdfile psd tpsd fpsd
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

# plot map, spetrogram, correlations and dv/v 
# load the map from GMT 
img = load(joinpath(@__DIR__,"../data/FINAL-FIGURES/LJR-map.png"))
p1 = plot(img,border=:none,dpi=500)
annotate!((-100,100,text("A",20)))

p2 = heatmap(
    psd_range,
    fpsd, 
    10 .* log10.(PSD),
    ylabel = "Frequency [Hz]",
    # legend=:none,
    colorbar_title=" \n  10 * log10 m^2/s^4/Hz dB",
    clims=(-170,-130),
    xlims=(Cdates[1],tpsd[end]),
    grid=false,
    left_margin=5mm,
    right_margin=15mm,
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
annotate!((Date(2001,3,1),9,text("B",20)))
p3 = heatmap(
    date_range,
    lags[posind],
    Apos,
    c=cgrad(:diverging_bwr_40_95_c42_n256,rev=true),
    # legend=:none,
    xlims=(Cdates[1],tpsd[end]),
    ylabel = "Lag [s]",
    left_margin=15mm,
    colorbar_title=" \n  Amplitude",
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
annotate!((Date(2001,3,1),9,text("C",20)))
l1 = @layout [a{0.4w} [b; c]]
plot(p1,p2,p3,layout=l1,size=(1260,640),dpi=500,margins=5Plots.mm)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-MAP-SPEC-CORR.png")) 

# A. dv/v + groundwater B. modeled groundwater level 
indbefore = findall(tlwe .< Date(2018))
indafter = findall(tlwe .> Date(2018))
# scatter dv/v
p1 = scatter(
    LJR[!,:DATE],
    LJR[!,:DVV],
    alpha=LJR[!,:CC] ./ 5,
    label="dv/v",
    seriescolor=:Reds_9,
    marker_z=LJR[!,:CC],
    markersize=4,
    markerstrokewidth=0.1,
    ylims=LJRylims,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    colorbar=false,
    yaxis=:flip,
    ylabel="        dv/v [%]",
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
# plot!(
#     tlwe[indbefore],
#     lwe[indbefore] .* coef(dvvgracefit)[2] .+ coef(dvvgracefit)[1],
#     xlims=(LJR[1,:DATE],LJR[end,:DATE]),
#     label="Scaled GRACE",
#     c=:gold,
#     linewidth=3,
#     ylims=LJRylims,
#     alpha=0.9,
# )
# # GRACE 2018- now 
# plot!(
#     tlwe[indafter],
#     lwe[indafter] .* coef(dvvgracefit)[2] .+ coef(dvvgracefit)[1],
#     xlims=(LJR[1,:DATE],LJR[end,:DATE]),
#     label="",
#     c=:gold,
#     linewidth=3,
#     ylims=LJRylims,
#     alpha=0.9,
# )
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
annotate!((Date(2000,6,1),-1.85,text("A",20)))
# plot scaled Elastic model 
plot!(
    LJR[:,:DATE],
    LJR[:,:ELASTIC] .* coef(dvvPEfit)[2] .+ coef(dvvPEfit)[1],
    c=:darkcyan,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    linewidth=3,
    alpha=0.75,
    # yguidefontcolor=:dodgerblue,
    # ytickfont=font(10,:dodgerblue),
    # yforeground_color_axis=:dodgerblue,
    # ylabel="Cumulative Deviation of\n Precipitation from Mean [m]",
    label="Elastic Model",
)
plot!(
    LJR[:,:DATE],
    LJR[:,:DRAINDED] .* coef(dvvDRAINfit)[2] .+ coef(dvvDRAINfit)[1],
    c=:red,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    linewidth=3,
    alpha=0.75,
    # yguidefontcolor=:dodgerblue,
    # ytickfont=font(10,:dodgerblue),
    # yforeground_color_axis=:dodgerblue,
    # ylabel="Cumulative Deviation of\n Precipitation from Mean [m]",
    label="Drained Model",
)
# plot cumulative precipiation 
cumind = findall((tppt .>= LJR[1,:DATE]) .& (tppt .<= LJR[end,:DATE]))
plot!(
    twinx(),
    Date.(tppt[cumind]), 
    cumprecip[cumind], 
    fillrange = 0.0, 
    fillcolor = :dodgerblue,
    fillalpha=0.3,
    # inset = (1, bbox(0.0, 0.0, 1.0, 0.45, :bottom, :left)), 
    # subplot = 3,
    label = "",
    grid=:off,
    bg_inside = nothing,
    ylabel="Cumulative Annual                         \n Precipitation [m]                         ",
    ylims = (0,2),
    left_margin = 10Plots.mm,
    right_margin=2.5cm,
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    yticks = (collect(0. : 0.2 : 1.),["0.0","0.2","0.4","0.6","0.8","1.0"]),
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    # xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xticks = false,
    xrot=60,
)
cs = palette(:Reds_9,6)
p2 = plot(    
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
        p2,
        LJR[:,:DATE],
        -LJR[:,:DVV]  ./ C ./ βs[ii] ./ 100,
        # c=:dodgerblue,
        alpha=0.9,
        xlims=(LJR[1,:DATE],LJR[end,:DATE]),
        linewidth=3,
        label="\\beta = -$value x 10^$exponent",
        c=cs[ii * 2],
        xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
        xrot=60,
    )
end
plot!(
    p2,
    MW1[!,:DATE],
    MW1[!,:GWL],
    label="",
    ylabel="\\Delta Groundwater Level [m]   ",
    xlims=(LJR[1,:DATE],LJR[end,:DATE]),
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    linewidth=3,
    alpha=0.9,
    c=:dodgerblue,
    xticks = (LJRticks, Dates.format.(LJRticks,"yyyy")),
    xrot=60,
)
xdate = Date(2007,7,1) : Day(1) : Date(2009,10,1)
revdate = Date(2013,5,1) : -Day(1) : Date(2012,6,1)
x = (xdate .- xdate[1]) ./ Day(1) ./ length(xdate) .* (π /3) .+ π / 6
y = (revdate .- revdate[end]) ./ Day(1) ./ length(revdate) .* (π /3) .+ π / 6
Plots.plot!(
    xdate,
    sin.(x) .* 4 .- 5.4,
    arrow=true,
    lw=2,
    color=:black,
    label="",
)
Plots.plot!(
    revdate,
    -cos.(y) .* 6 .+ 7.4,
    arrow=true,
    lw=2,
    color=:black,
    label="",
)
annotate!((Date(2013,6,1),10,text("      Modeled    \n Groundwater Level",12, :red)))
annotate!((Date(2006,6,1),-4,text("     Measured    \n Groundwater Level",12, :dodgerblue)))
annotate!((Date(2000,6,1),15,text("B",20)))
l1 = @layout [a; b]
plot(p1,p2,layout=l1,size=(1200,800),dpi=500,margins=5Plots.mm)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-DVV-GWL.png"))
