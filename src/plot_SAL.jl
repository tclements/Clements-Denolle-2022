using Arrow 
using CSV
using DataFrames
using Dates 
using DSP
using LinearAlgebra
using LsqFit
using Plots 
using Plots.Measures
using NetCDF
using Statistics

function smooth_withfiltfilt(A::AbstractArray; window_len::Int=11, window::Symbol=:rect)
    w = getfield(DSP.Windows, window)(window_len)
    A = DSP.filtfilt(w ./ sum(w), A)
    return A
end

# define a function that returns a Plots.Shape
rectangle(w, h, x, y) = Shape(x + [0,w,w,0], y + [0,0,h,h])

# read dv/v for CI.SAL 
SAL = Arrow.Table("/media/FOUR/data/DVV-10-DAY-COMP/2.0-4.0/CI.SAL.arrow") |> DataFrame

# get SAL lat, lon from CI station locations 
SCdf = DataFrame(CSV.File("/home/timclements/CALI/CIstations.csv"))
SALind = findfirst(SCdf[!,:Station] .== "SAL")
SALlat = SCdf[SALind,:Latitude]
SALlon = SCdf[SALind,:Longitude]

# read temperature
filename = "/media/FOUR/data/tmean.nc"
ttmean = ncread(filename,"t")
lon = ncread(filename,"lon")
lat = ncread(filename,"lat")
ttmean = Date.(Dates.unix2datetime.(ttmean))
lonind = argmin(abs.(lon .- SALlon))
latind = argmin(abs.(lat .- SALlat))
temp = ncread(filename,"tmean",start=[latind,lonind,1],count=[1,1,-1])[1,1,:]
temp .-= mean(temp)
smoothtemp = smooth_withfiltfilt(temp,window_len=45)

# load precip data
filename = "/media/FOUR/data/ppt.nc"
lon = ncread(filename,"lon")
lat = ncread(filename,"lat")
t = ncread(filename,"t")
t = Dates.unix2datetime.(t)
lonind = argmin(abs.(lon .- SALlon))
latind = argmin(abs.(lat .- SALlat))
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
tind = findall((t .>= SAL[1,:DATE]) .& (t .<= SAL[end,:DATE]))
t = t[tind]
cumprecip = cumprecip[tind]

# before EQ indices 
before2010 = findall(SAL[:,:DATE] .< Date(2010,4,1))

# fit model to temp + dv/v 
@. model(x,p) = p[1] + p[2] * x
dvvtempind = findall(in(ttmean),SAL[before2010,:DATE])
tempdvvind = findall(in(SAL[before2010,:DATE]),ttmean)
dvvtempfit = curve_fit(
    model,
    SAL[before2010,:DVV],
    smoothtemp[tempdvvind],
    [0,1.],
)

tempdvvfit = curve_fit(
    model,
    smoothtemp[tempdvvind],
    SAL[before2010,:DVV],
    [0,1.],
)

# plot SAL dv/V
SALylims = (-1.0,1.5)
scatter(
    SAL[!,:DATE],
    SAL[!,:DVV],
    seriescolor=:Blues_9,
    marker_z=SAL[!,:CC],
    alpha=SAL[!,:CC] ./ 2,
    markersize=2,
    label="",
    colorbar=false,
    ylabel="dv/v [%]",
    ylims=SALylims,
    xlims=(SAL[1,:DATE],SAL[end,:DATE]),
    right_margin = 2cm,
    left_margin = 2cm,
    yguidefontcolor=:darkblue,
    ytickfont=font(10,:darkblue),
    yforeground_color_axis=:darkblue,
    ytick_direction=:out,
    dpi=500,
    size=(800,400),
)
# plot temperature
plot!(
    twinx(), 
    ttmean, 
    smoothtemp, 
    sharex=true, 
    xticks=false,
    grid=:off, 
    ylims=SALylims .* coef(dvvtempfit)[2] .+ coef(dvvtempfit)[1],
    linewidth=1.5,
    c=:red,
    alpha=0.75,
    label="",
    ylabel="\\Delta Surface Temperature [C]",
    seriescolor=:orange,
    xlims=(SAL[1,:DATE],SAL[end,:DATE]),
    yguidefontcolor=:red,
    ytickfont=font(10,:red),
    yforeground_color_axis=:red,
    colorbar=false,
)

savefig("/media/FOUR/data/FINAL-FIGURES/CISAL-DVV.svg")
savefig("/media/FOUR/data/FINAL-FIGURES/CISAL-DVV.png")

SALylims = (-1.0,1.5)
scatter(
    SAL[!,:DATE],
    SAL[!,:DVV],
    seriescolor=:Blues_9,
    marker_z=SAL[!,:CC],
    alpha=SAL[!,:CC] ./ 2,
    markersize=2,
    label="",
    colorbar=false,
    ylabel="dv/v [%]",
    ylims=SALylims,
    xlims=(SAL[1,:DATE],SAL[end,:DATE]),
    right_margin = 2cm,
    left_margin = 2cm,
    yguidefontcolor=:darkblue,
    ytickfont=font(10,:darkblue),
    yforeground_color_axis=:darkblue,
    ytick_direction=:out,
    dpi=500,
    size=(800,400),
)
# plot inset for 2010 EQ
ind = findall((SAL[:,:DATE] .>= Date(2009,11,1)) .& (SAL[:,:DATE] .< Date(2010,9,1)))
# plot!(rectangle(SAL[ind[end],:DATE] - SAL[ind[1],:DATE],-0.9,SAL[ind[1],:DATE],1.5), opacity=.5)
vline!(
    [Date(2010,4,4)],
    label="2010 M7.2 Baja, CA",
    legend=:bottomleft,
    alpha=0.85,
    linewidth=1.5,
    c=:black,
    linestyle=:dash,
)
# plot temperature
plot!(
    twinx(), 
    ttmean, 
    smoothtemp, 
    sharex=true, 
    xticks=false,
    grid=:off, 
    ylims=SALylims .* coef(dvvtempfit)[2] .+ coef(dvvtempfit)[1],
    linewidth=1.5,
    c=:red,
    alpha=0.75,
    label="",
    ylabel="\\Delta Surface Temperature [C]",
    seriescolor=:orange,
    xlims=(SAL[1,:DATE],SAL[end,:DATE]),
    yguidefontcolor=:red,
    ytickfont=font(10,:red),
    yforeground_color_axis=:red,
    colorbar=false,
)

# inset for velocity drop 
scatter!(
    SAL[ind,:DATE], 
    SAL[ind,:DVV], 
    inset = (1, bbox(0.165, 0.0, 0.66, 0.23, :top)), 
    subplot = 3,
    markersize=2,
    seriescolor=:Blues_9,
    marker_z=SAL[!,:CC],
    alpha=SAL[!,:CC] ./ 2,
    colorbar=false,
    label = "",
    grid=:off,
    bg_inside = nothing,
    ylabel="dv/v [%]",
    ylims = (-0.9,0.6),
    left_margin = 5mm,
    xlims=(SAL[ind[1],:DATE],SAL[ind[end],:DATE]),
    # yticks = (collect(0. : 0.2 : 1.),["0.0    ","0.2    ","0.4    ","0.6    ","0.8    ","1.0    "]),
    # yguidefontcolor=:dodgerblue,
    # ytickfont=font(10,:dodgerblue),
    # yforeground_color_axis=:dodgerblue,
    xtick=:off,
)

savefig("/media/FOUR/data/FINAL-FIGURES/CISAL-DVV-EQ.svg")
savefig("/media/FOUR/data/FINAL-FIGURES/CISAL-DVV-EQ.png")
