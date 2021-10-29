using Arrow 
using CSV
using DataFrames
using Dates
using DSP 
using Impute
using Interpolations
using LaTeXStrings
using LsqFit
using Plots 
using Plots.Measures
using Statistics

function smooth_withfiltfilt(A::AbstractArray; window_len::Int=11, window::Symbol=:rect)
    w = getfield(DSP.Windows, window)(window_len)
    A = DSP.filtfilt(w ./ sum(w), A)
    return A
end

# read data from CI.RXH 
RXH = Arrow.Table("/media/FOUR/data/COMP-ONE-DVV/1.0-2.0/CI.RXH.arrow") |> DataFrame
EN = Arrow.Table("/media/FOUR/data/ONE-DVV/1.0-2.0/CI.RXH..BHE.CI.RXH..BHN.arrow") |> DataFrame
EZ = Arrow.Table("/media/FOUR/data/ONE-DVV/1.0-2.0/CI.RXH..BHE.CI.RXH..BHZ.arrow") |> DataFrame
NZ = Arrow.Table("/media/FOUR/data/ONE-DVV/1.0-2.0/CI.RXH..BHN.CI.RXH..BHZ.arrow") |> DataFrame

# get RXH lat, lon from CI station locations 
cistations = DataFrame(CSV.File("/home/timclements/CALI/CIstations.csv"))
RXHind = findfirst(cistations[!,:Station] .== "RXH")
RXHlat = cistations[RXHind,:Latitude]
RXHlon = cistations[RXHind,:Longitude]

# find groundwater wells nearby 
LAKE = CSV.File("/media/FOUR/data/saltonelevation.tsv",comment="#",skipto=3) |> DataFrame
dropmissing!(LAKE,:WSE) # lev_va is ft below surface
LAKE[!,:WSE] ./= 3.28 # convert to meters 
LAKE[:,:WSE] .-= mean(LAKE[:,:WSE])

# fit surface elevation model to dv/v 
@. model(x,p) = p[1] + p[2] * x
DVV = (EN[!,:DVVNEG] .+ EZ[!,:DVVNEG] .+ NZ[!,:DVVNEG]) ./ 3
CC = (EN[!,:CCNEG] .+ EZ[!,:CCNEG] .+ NZ[!,:CCNEG]) ./ 3

# subset CC over 0.5
ind = findall(CC .> 0.62)
DVV = DVV[ind]
CC = CC[ind]
RXH = RXH[ind,:]
dvvgwlind = findall(in(LAKE[!,:DATE]),RXH[:,:DATE])
gwldvvind = findall(in(RXH[:,:DATE]),LAKE[:,:DATE])
dvvgwlfit = curve_fit(
    model,
    LAKE[gwldvvind,:WSE],
    DVV[dvvgwlind],
    [0,-1.],
)

resid = DVV[dvvgwlind] .- (coef(dvvgwlfit)[1] .+ coef(dvvgwlfit)[2] .*  LAKE[gwldvvind,:WSE])

# plot groundwater level 
RXHylims = (-1.0,1.25)
p1 = scatter(
    RXH[!,:DATE],
    DVV,
    alpha=CC ./ 5,
    marker_z=(EZ[!,:CCNEG] .+ NZ[!,:CCNEG]) ./ 2,
    seriescolor=:Reds_9,
    xlims=(RXH[1,:DATE],LAKE[end,:DATE]),
    ylims=RXHylims,
    label="",
    yaxis=:flip,
    ylabel="dv/v [%]",
    colorbar=false,
    right_margin = 2cm,
    left_margin = 1cm,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    yforeground_color_axis=:darkred,
    dpi=500,
    size=(800,400),
)
plot!(
    twinx(),
    LAKE[!,:DATE],
    LAKE[!,:WSE],
    xlims=(RXH[1,:DATE],LAKE[end,:DATE]),
    ylims=reverse((RXHylims .- coef(dvvgwlfit)[1]) ./ coef(dvvgwlfit)[2]),
    seriescolor=:dodgerblue,
    linewidth=4,
    alpha=0.85,
    grid=false,
    label="",
    xticks=false,
    ylabel="\\Delta Salton Sea [m]",
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
)
# plot residual 
p2 = vline(
    [Date(2005,8,31)],
    label="2005 Obsidian Butte Swarm",
    ylabel="dv/v residual [%]",
    xlims=(RXH[1,:DATE],LAKE[end,:DATE]),
    legend=:bottom,
    alpha=0.85,
    linewidth=1.5,
    c=:black,
    linestyle=:dash,
    right_margin = 2cm,
    left_margin = 1cm,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    yforeground_color_axis=:darkred,
    dpi=500,
    size=(800,400),
)
vline!(
    [Date(2014,12,24)],
    label="2014 M4.2",
    alpha=0.85,
    linewidth=1.5,
    c=:red,
    linestyle=:dash,
)
scatter!(
    RXH[dvvgwlind,:DATE],
    resid,
    alpha=CC ./ 5,
    label="",
    seriescolor=:Reds_9,
    marker_z=CC,
    colorbar=false,
)
l = Plots.@layout [a; b]
plot(p1,p2,layout=l)
savefig("/media/FOUR/data/FINAL-FIGURES/CIRXH-DVV.svg")
savefig("/media/FOUR/data/FINAL-FIGURES/CIRXH-DVV.png")

# plot groundwater level 
RXHylims = (-1.0,1.25)
p1 = scatter(
    RXH[1:1,:DATE],
    DVV[1:1],
    alpha=CC ./ 5,
    marker_z=(EZ[1:1,:CCNEG] .+ NZ[1:1,:CCNEG]) ./ 1000,
    seriescolor=:Reds_9,
    xlims=(RXH[1,:DATE],LAKE[end,:DATE]),
    ylims=RXHylims,
    label="",
    yaxis=:flip,
    ylabel="dv/v [%]",
    colorbar=false,
    right_margin = 2cm,
    left_margin = 1cm,
    yguidefontcolor=:darkred,
    ytickfont=font(10,:darkred),
    yforeground_color_axis=:darkred,
    dpi=500,
    size=(800,400),
)
plot!(
    twinx(),
    LAKE[!,:DATE],
    LAKE[!,:WSE],
    xlims=(RXH[1,:DATE],LAKE[end,:DATE]),
    ylims=reverse((RXHylims .- coef(dvvgwlfit)[1]) ./ coef(dvvgwlfit)[2]),
    seriescolor=:dodgerblue,
    linewidth=2.5,
    alpha=0.85,
    grid=false,
    label="",
    xticks=false,
    ylabel="\\Delta Salton Sea [m]",
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
)
l = Plots.@layout [a; b]
plot(p1,p2,layout=l)
savefig("/media/FOUR/data/FINAL-FIGURES/SALTON-ELEV.svg")
savefig("/media/FOUR/data/FINAL-FIGURES/SALTON-ELEV.png")

