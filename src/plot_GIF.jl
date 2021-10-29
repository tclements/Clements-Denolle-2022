using Arrow 
using Dates
using DataFrames
using DSP
using GLM 
using Glob 
using LsqFit
using NetCDF
using Plots
using Plots.Measures
using Statistics
using StatsModels

function CDM(A::AbstractArray,k::Int)
    Amean = [cumsum(A[1:k-1]) ./ (1:k-1) ;rolling_mean(A,k)]
    return cumsum(A .- Amean)
end

function rolling_mean(A::AbstractArray,k::Int)
    B = cumsum([zero(eltype(A)); A])
    return (B[k+1:end] - B[1:end-k]) ./ k
end

@userplot CDMPlot
@recipe function f(cp::CDMPlot)
    precip, t, days = cp.args
    linewidth --> 2.5
    seriesalpha --> 0.85
    seriescolor --> :dodgerblue
    label --> false
    xlims --> (tppt[pptminind],tppt[pptmaxind])
    ylims --> (-1000, 1000)
    yguide --> "Cumulative Deviation from Mean [m]"
    guide --> "$days day lagging mean"
    dpi --> 250
    t, precip .- mean(precip)
end

@userplot hystPlot
@recipe function f(hp::hystPlot)
    precip, dvv, cc, t, minind, maxind, dvvprecipfit = hp.args 
    N = length(precip)

    layout := @layout [a{0.75h}; b]

    # hysteresis plot 
    @series begin 
        subplot := 1
        seriestype --> :path
        linestyle --> :dash
        alpha --> 0.5
        color --> :red
        label --> ""
        xguide --> "Cumulative Deviation from Mean [cm]"
        yguide --> "dv/v [%]"
        dpi --> 250 
        xlims --> (-0.43,0.7)
        ylims --> (-2.5,2)
        -0.43:0.01:0.7, coef(dvvprecipfit)[2] .* collect(-0.43:0.01:0.7) .+ coef(dvvprecipfit)[1]
    end 
    
    @series begin 
        subplot := 1
        seriestype --> :scatter
        seriesalpha --> LinRange(0,1,N)
        seriescolor --> :magma
        marker_z --> dvv[minind:maxind]
        label --> last(t[minind:maxind])
        xlims --> (-0.43,0.7)
        ylims --> (-2.5,2)
        clims --> (-2.5,2)
        precip[minind:maxind], dvv[minind:maxind]
    end 

    # dv/v plot 
    @series begin
        subplot := 2
        seriestype --> :scatter
        seriesalpha --> cc ./ 10
        yguide --> "dv/v [%]"
        label --> ""
        ylims --> (-2.5,2)
        t, dvv 
    end
    @series begin
        subplot := 2
        seriestype --> :scatter
        seriesalpha --> 1.
        label --> ""
        yguide --> "dv/v [%]"
        ylims --> (-2.5,2)
        seriescolor --> :orange
        markersize --> 6
        t[maxind:maxind], dvv[maxind:maxind] 
    end
end

# read CI.LJR dv/v
LJR = Arrow.Table("/media/FOUR/data/DVV-10-DAY-COMP/2.0-4.0/CI.LJR.arrow") |> DataFrame
LJRlon = -118.86775
LJRlat = 34.80762

# load precip data 
filename = "/media/FOUR/data/ppt.nc"
pptlon = ncread(filename,"lon")
pptlat = ncread(filename,"lat")
t = ncread(filename,"t")
tppt = Date.(Dates.unix2datetime.(t))
lonind = argmin(abs.(pptlon .- LJRlon))
latind = argmin(abs.(pptlat .- LJRlat))
precip = ncread(filename,"ppt",start=[latind,lonind,1],count=[1,1,-1])[1,1,:]


days = [1:29;30: 30 : 8 * 365]
N = length(days)
pptminind = findfirst(tppt .== Date(2002,11,12))
pptmaxind = findfirst(tppt .== Date(2021,2,18))

# gif of CDM 
anim = @animate for jj in 1:N
    CDMjj = CDM(precip[pptminind - days[jj]:end],days[jj])
    cdmplot(CDMjj[days[jj]+1:end - (length(tppt) - pptmaxind)],tppt[pptminind:pptmaxind],days[jj])
end
gif(anim, "/media/FOUR/data/FINAL-FIGURES/CDM-LJR_fps20.gif", fps = 20)

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
CDMbest = CDM(precip[tminind - bestday:end],bestday) ./ 1000
CDMbest = CDMbest[bestday:bestday+(tmaxind-tminind)]
# CDMbest .-= CDMbest[1]
dvvprecipind = tminind:tmaxind
modelind = findall(in(LJR[:,:DATE]),tppt[dvvprecipind])

# fit precip model to dv/v
@. model(x,p) = p[1] + p[2] * x
dvvprecipfit = curve_fit(
    model,
    CDMbest[modelind],
    LJR[!,:DVV],
    [0.,-1.],
)


# gif of dv/v vs CDM 
days = 30:30:size(LJR,1)
anim = @animate for maxind in days
    minind = max(1,maxind-365*3)
    println(maxind / 30)
    hystplot(
        CDMbest[modelind],
        LJR[:,:DVV],
        LJR[:,:CC],
        LJR[:,:DATE],
        minind,
        maxind,
        dvvprecipfit,
    )
end
gif(anim, "/media/FOUR/data/FINAL-FIGURES/HYST-DVV-PRECIP-fps5.gif", fps = 5)

# plot daily precip 
plot(
    tppt,
    precip,
    seriestype=:sticks,
    xlims=(tppt[pptminind],tppt[pptmaxind]),
    ylims=(0,80),
    ylabel="Daily Precipitation [mm]",
    label="",
    dpi=250,
)
savefig("/media/FOUR/data/FINAL-FIGURES/LJR-PRECIP.png")
savefig("/media/FOUR/data/FINAL-FIGURES/LJR-PRECIP.svg")

# plot 4 different lag times
bestday = 2789 
days = [30, 30*6, 30 * 24, bestday]
N = length(days)
pptminind = findfirst(tppt .== Date(2002,11,12))
pptmaxind = findfirst(tppt .== Date(2021,2,18))
ps = Array{Plots.Plot}(undef,N)
for jj in 1:N
    CDMjj = CDM(precip[pptminind - days[jj]:end],days[jj]) ./ 1000
    CDMjj = CDMjj[days[jj]+1:end - (length(tppt) - pptmaxind)]
    CDMjj .-=CDMjj[1]
    months = round(Int,days[jj] / 30.4375)
    ps[jj] = plot(
        tppt[pptminind:pptmaxind],
        CDMjj,
        label="$months month" * "s" ^ Int(months > 1),
        linewidth=1.,
        alpha=0.75,
        linecolor=:dodgerblue,
        ylabel="CDM [m]",
        dpi=250,
    )
end
l = @layout [a b; c d]
plot(
    ps...,
    layout=l,
    sharex=true,
    sharey=true,
    xformatter = x -> Dates.format(Date(Dates.UTD(x)), "yyyy"),
)
savefig("/media/FOUR/data/FINAL-FIGURES/CDM-MONTHS.png")
savefig("/media/FOUR/data/FINAL-FIGURES/CDM-MONTHS.svg")

# load precipitation data 
filename = "/media/FOUR/data/PRECIP-LJR.nc"
t = ncread(filename,"t")
tppt = Date.(Dates.unix2datetime.(t))
precip = ncread(filename,"ppt")[1,1,:]

# subset data from Oct 1985 - Oct 2020 
startdate = Date(1985,10,1)
enddate = Date(2020,9,30)
ind = findall(startdate .<= tppt .<= enddate)
tppt = tppt[ind]
precip = precip[ind]

# find days with 90% rainfall 
wetind = findall(precip .> 0.0)
p90 = quantile(precip[wetind],0.90)
wetind10 = findall(precip .> p90)
wetind90 = setdiff(wetind,wetind10)

# cumulative precip 
octs = Date(1985,10,1):Year(1):Date(2019,10,1)
seps = Date(1986,9,30):Year(1):Date(2020,9,30)
yearprecip = zeros(length(octs))
yearprecip10 = zeros(length(octs))
yearprecip90 = zeros(length(octs))
cumprecip = zeros(length(tppt))
cumprecip10 = zeros(length(tppt))
cumprecip90 = zeros(length(tppt))
for ii = 1:length(seps)
    ind = findall((tppt .>= octs[ii]) .& (tppt .<= seps[ii]))
    year10 = intersect(ind,wetind10)
    year90 = intersect(ind,wetind90)
    yr = year(tppt[ind][1])
    println("$yr-$(yr+1) $(length(year10))")

    # yearly precip 
    yearprecip[ii] = sum(precip[ind])
    yearprecip10[ii] = sum(precip[year10])
    yearprecip90[ii] = sum(precip[year90])

    # daily cumulative 
    cumprecip[ind] .= cumsum(precip[ind])
    cumprecip10[year10] = precip[year10]
    cumprecip90[year90] = precip[year90]
    cumprecip10[ind] .= cumsum(cumprecip10[ind])
    cumprecip90[ind] .= cumsum(cumprecip90[ind])
end


ind92 = findall(octs .> Date(1992))
tDVV = (octs[ind92] .- octs[ind92][1]) ./ Day(1) .+ 1
data = DataFrame(X=tDVV ./ tDVV[end],Y=yearprecip90[ind92])
ols = lm(@formula(Y ~ X), data)

# fit exponential model to precip 
h = fit(Histogram,precip[wetind],0.0:1.0:80.0)
fitexp = fit(Exponential,h.weights)
dataexp = DataFrame(
    X=log.(h.edges[1][1:end-1] ./2 .+ h.edges[1][2:end] ./2),
    Y=-log(1/scale(fitexp)) .- 1/scale(fitexp) .* h.weights,
    )
olsexp = lm(@formula(Y ~ X), dataexp)


# plot cumulative annual precip 
plot(
    tppt, 
    cumprecip, 
    fillrange = cumprecip90, 
    fillalpha = 0.35, 
    c = :red, 
    label = "Wettest 10%",
    ylabel="Cumulative Annual Precipitation [mm]",
    rigth_margin=1.5cm,
    left_margin=1cm,
    size=(800,400),
    dpi=250,
)
plot!(
    tppt, 
    cumprecip90, 
    fillrange = 0., 
    fillalpha = 0.35, 
    c = 1, 
    label = "Remaining 90%",
)
savefig("/media/FOUR/data/FINAL-FIGURES/LJR-CUMPRECIP.png")