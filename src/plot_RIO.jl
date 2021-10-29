using Arrow 
using CSV
using DataFrames 
using Dates 
using Interpolations
using LaTeXStrings
using LsqFit
using Plots
using Plots.PlotMeasures
using NetCDF
using Statistics


# read CI.RIO dv/v
RIO = Arrow.Table("/media/FOUR/data/DVV-90-DAY-COMP/2.0-4.0/CI.RIO.arrow") |> DataFrame
EN = Arrow.Table("/media/FOUR/data/DVV-90-DAY/2.0-4.0/CI.RIO..BHE.CI.RIO..BHN.arrow") |> DataFrame
EZ = Arrow.Table("/media/FOUR/data/DVV-90-DAY/2.0-4.0/CI.RIO..BHE.CI.RIO..BHZ.arrow") |> DataFrame
NZ = Arrow.Table("/media/FOUR/data/DVV-90-DAY/2.0-4.0/CI.RIO..BHN.CI.RIO..BHZ.arrow") |> DataFrame

mindate = Date(1999,10,1)

# read nearest groundwater well and convert to m 
GWL = Arrow.Table("/media/FOUR/data/baldwinpark.arrow") |> Arrow.columntable |> DataFrame
dropmissing!(GWL,:WSE) # WSE is GWL elevation 
GWL[!,:WSE] ./= 3.28 # convert to meters 
GWL = GWL[GWL[!,:DATE] .>= mindate,:]
GWL = unique(GWL,:DATE) # take one measurement per day 
GWL[:,:WSE] .-= mean(GWL[:,:WSE])
# GWL[!,:lev_va] .*= -1

# get RIO lat, lon from CI station locations 
cistations = DataFrame(CSV.File("/home/timclements/CALI/CIstations.csv"))
RIOind = findfirst(cistations[!,:Station] .== "RIO")
RIOlat = cistations[RIOind,:Latitude]
RIOlon = cistations[RIOind,:Longitude]

# load precip data
filename = "/media/FOUR/data/ppt.nc"
lon = ncread(filename,"lon")
lat = ncread(filename,"lat")
t = ncread(filename,"t")
t = Dates.unix2datetime.(t)
lonind = argmin(abs.(lon .- RIOlon))
latind = argmin(abs.(lat .- RIOlat))
precip = ncread(filename,"ppt",start=[latind,lonind,1],count=[1,1,-1])[1,1,:]

cumprecip = zeros(eltype(precip),length(precip))
octs = Date(1999,10,1):Year(1):Date(2020,10,1)
mays = Date(2000,5,1):Year(1):Date(2021,5,1)
for dd = 1:length(octs)
    ind = findall((t .> octs[dd]) .& (t .< mays[dd]))
    cumprecip[ind] .= cumsum(precip[ind])
end
cumprecip ./= 1000

# subset t with dv/v 
tind = findall((t .>= RIO[1,:DATE]) .& (t .<= RIO[end,:DATE]))
t = t[tind]
cumprecip = cumprecip[tind]

tGWL = (GWL[!,:DATE] .- GWL[1,:DATE]) ./ Day(1) .+ 1 
tGWLinterp = tGWL[1]:tGWL[end]
interp_linear = LinearInterpolation(tGWL,GWL[:,:WSE])
GWLinterp = interp_linear.(tGWLinterp)

# now fit groundwater model to dv/v 
@. model(x,p) = p[1] + p[2] * x
dvvgwlind = findall(in(RIO[!,:DATE]),GWL[1,:DATE]:Day(1):GWL[end,:DATE])
gwldvvind = findall(in(GWL[1,:DATE]:Day(1):GWL[end,:DATE]),RIO[:,:DATE])
dvvgwlfit = curve_fit(
    model,
    RIO[dvvgwlind,:DVV],
    GWLinterp[gwldvvind],
    [0,-1.],
)

# plot dv/v first 
RIOylims = (-4,4) 
scatter(
    RIO[!,:DATE],
    RIO[!,:DVV],
    alpha=RIO[!,:CC] ./ 10,
    label="",
    seriescolor=:Reds_9,
    marker_z=RIO[!,:CC],
    ylims=RIOylims,
    xlims=(RIO[1,:DATE],RIO[end,:DATE]),
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
    size=(800,400)
)
# plot groundwater level 
plot!(
    twinx(),
    GWL[!,:DATE],
    GWL[!,:WSE] .- mean(GWL[!,:WSE]),
    sharex=true,
    ylim=reverse(RIOylims .* coef(dvvgwlfit)[2] .+ coef(dvvgwlfit)[1]),
    xlims=(RIO[1,:DATE],RIO[end,:DATE]),
    linewidth=2.5,
    alpha=0.85,
    seriescolor=:royalblue3,
    label="",
    ylabel="\\Delta Groundwater Level Change [m]",
    grid=:off,
    yguidefontcolor=:royalblue3,
    ytickfont=font(10,:royalblue3),
    yforeground_color_axis=:royalblue3,
)

# plot cumulative precipiation 
plot!(
    Date.(t), 
    cumprecip, 
    fillrange = 0.0, 
    fillcolor = :dodgerblue,
    fillalpha=0.3,
    inset = (1, bbox(0.0, 0.0, 1.0, 0.33, :bottom, :left)), 
    subplot = 3,
    label = "",
    grid=:off,
    bg_inside = nothing,
    ylabel="Cumulative Annual\n Precipitation [m]\n\n",
    ylims = (0,1.5),
    xlims=(RIO[1,:DATE],RIO[end,:DATE]),
    left_margin = 5mm,
    # yticks = (collect(0. : 0.2 : 1.),["0.0    ","0.2    ","0.4    ","0.6    ","0.8    ","1.0    "]),
    yguidefontcolor=:dodgerblue,
    ytickfont=font(10,:dodgerblue),
    yforeground_color_axis=:dodgerblue,
    xtick=:off,
)
savefig("/media/FOUR/data/FINAL-FIGURES/CIRIO-DVV.svg")
savefig("/media/FOUR/data/FINAL-FIGURES/CIRIO-DVV.png")
