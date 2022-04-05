using Arrow
using CSV
using DataFrames
using Dates
using DSP 
using Glob
using Images
using Impute
using NetCDF
using SpecialFunctions
using Statistics
try 
    import GMT
catch
end
import GMT 
import Plots
import Plots.Measures
import SeisNoise

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

function CDM(A::AbstractArray,k::Int)
    Amean = [cumsum(A[1:k-1]) ./ (1:k-1) ;rolling_mean(A,k)]
    return cumsum(A .- Amean)
end

function rolling_mean(A::AbstractArray,k::Int)
    B = cumsum([zero(eltype(A)); A])
    return (B[k+1:end] - B[1:end-k]) ./ k
end

function SSW06(precip,ϕ,a)
    expij = @fastmath(exp.(-a .* (0:length(precip))))
    GWL = conv(expij,precip .- mean(precip))[1:end÷2] ./ ϕ
    return GWL
end

function elastic(precip,c,α,r,δt)
    erfij = erf.(r ./ sqrt.(4. .* c .* Float64.(0:length(precip)) .* δt ))
    P = conv(erfij,precip .- mean(precip))[1:end÷2] .* α
    return P
end

function fullycoupled(precip,c,α,r,δt)
    erfij = erf.(r ./ sqrt.(4. .* c .* Float64.(0:length(precip)) .* δt ))
    erfcij = erfc.(r ./ sqrt.(4. .* c .* Float64.(0:length(precip)) .* δt ))
    P = conv(erfij,precip .- mean(precip))[1:end÷2] .* α .+ conv(erfcij,precip .- mean(precip))[1:end÷2]
    return P
end

function drained(precip,c,r,δt)
    erfcij = erfc.(r ./ sqrt.(4. .* c .* Float64.(0:length(precip)) .* δt ))
    P = conv(erfcij,precip .- mean(precip))[1:end÷2]
    return P
end

# create lowpass butterworth filter 
responsetype = Lowpass(0.01,fs=1.)
designmethod = Butterworth(4)

# load fitting DataFrame 
fitdf = Arrow.Table(joinpath(@__DIR__,"../data/hydro-model-90-day.arrow")) |> Arrow.columntable |> DataFrame

# constants for models 
r = 500.
δt = 86400.
α = (1 + 0.27) / 3 / (1 - 0.27)
ϕ = 0.15   # porosity 


# load precip data 
filename = joinpath(@__DIR__,"../data/ppt.nc")
pptlon = ncread(filename,"lon")
pptlat = ncread(filename,"lat")
tppt = ncread(filename,"t")
tppt = Date.(Dates.unix2datetime.(tppt))
ppt = ncread(filename,"ppt")
tpptday = (tppt .- tppt[1]) ./ Day(1)

# load tmean data 
filename = joinpath(@__DIR__,"../data/tmean.nc")
tmean = ncread(filename,"tmean")
ttmean = ncread(filename,"t")
ttmean = Date.(Dates.unix2datetime.(ttmean))
ttmeanday = (ttmean .- ttmean[1]) ./ Day(1)

# load GRACE data 
filename = joinpath(@__DIR__,"../data/CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02.nc")
tlwe = Date(2002,1,1) .+ Day.(floor.(ncread(filename,"time")))
lwe = ncread(filename, "lwe_thickness")
lwelon = ncread(filename,"lon")
lwelat = ncread(filename,"lat")
lwelon[end ÷ 2 + 1: end]  .-= 360

# load files 
LA1 = CSV.File(joinpath(@__DIR__,"../data/LA-stations-1.txt"),delim="|") |> DataFrame
LA2 = CSV.File(joinpath(@__DIR__,"../data/LA-stations-2.txt"),delim="|") |> DataFrame
LA = vcat(LA1,LA2)
LA[:,:NETSTA] = LA[:,:Network] .* "." .* LA[:,:Station]
sort!(LA,:NETSTA)
LA = innerjoin(LA,fitdf,on=:NETSTA)
LAminlat = minimum(LA[:,:Latitude])
LAmaxlat = maximum(LA[:,:Latitude])
LAminlon = minimum(LA[:,:Longitude])
LAmaxlon = maximum(LA[:,:Longitude])
LAminlat = floor(LAminlat)
LAmaxlat = ceil(LAmaxlat)
LAminlon = round(LAminlon + (-0.5 - LAminlon % 0.5),digits=2)
LAmaxlon = ceil(LAmaxlon)

# load dv/v data 
DVVDIR = joinpath(@__DIR__,"../data/DVV-90-DAY-COMP/2.0-4.0")
files = glob("*",DVVDIR)
LAfiles = [f for f in files if in(replace(basename(f),".arrow"=>""),LA[:,:NETSTA])]
LA[:,:FILES] = LAfiles

# join dataframes 
LA1 = innerjoin(LA1,LA,on=[:Station,:Network,:Latitude,:Longitude,:Elevation,:Sitename,:StartTime,:EndTime])
LA2 = innerjoin(LA2,LA,on=[:Station,:Network,:Latitude,:Longitude,:Elevation,:Sitename,:StartTime,:EndTime])

# sort by distance from LJR
LJRind = findall(LA[:,:Station] .== "LJR")[1]
LJRlat = LA[LJRind,:Latitude]
LJRlon = LA[LJRind,:Longitude]
LJRdist1 = SeisNoise.surface_distance.(LJRlon,LJRlat,LA1[:,:Longitude],LA1[:,:Latitude],true)
ind1 = sortperm(LJRdist1,rev=true)
LJRdist2 = SeisNoise.surface_distance.(LJRlon,LJRlat,LA2[:,:Longitude],LA2[:,:Latitude],true)
ind2 = sortperm(LJRdist2,rev=true)
LA1 = LA1[ind1,:]
LA2 = LA2[ind2,:]

# read LWE for LA 
ind2006 = findfirst(tlwe .>= Date(2006,1,1))
ind2021 = length(tlwe)
G2006 = GMT.gmtread(filename,layer=ind2006,varname="lwe_thickness")
G2021 = GMT.gmtread(filename,layer=ind2021,varname="lwe_thickness")

# plot with GMT 
gray = GMT.makecpt(color=150, range=(-10000,10000), no_bg=:true);
CGRACE = GMT.makecpt(T=(-30,30,30 / 101), cmap=:vik, reverse=true);
GMT.grdimage(
    "@srtm_relief_01s",
    cmap=gray, 
    J=:guess,
    shade=true,
    region=(LAminlon,LAmaxlon,LAminlat,LAmaxlat),
    coast=true,
    colorbar=false,
)
GMT.grdimage!(
    (G2021 .- G2006),
    cmap=CGRACE,
    alpha=50,
    colorbar=false,
)
GMT.coast!(
    N=2,
    ocean=:white,
)
GMT.colorbar!(
    C=CGRACE,
    frame=(annot=:auto, ticks=:auto, xlabel="LWE Change [cm]"), 
)
GMT.scatter!(
    LA[:,:Longitude],
    LA[:,:Latitude],
    # zcolor=-ΔDVVdf[:,:SLOPE],
    markeredgecolor=:black,
    marker=:triangle,
    color=:gold,
    transparency= 15,
    markersize=15 .* ones(size(LA,1)),
    # colorbar=true,
)
GMT.text!(
    GMT.mat2ds(
        hcat(
            LA[:,:Longitude] .+ 0.05,
            LA[:,:Latitude] .+ 0.05 ,
        ),
        LA[:,:NETSTA],
    ),
    show=true,
    savefig=joinpath(@__DIR__,"../data/FINAL-FIGURES/LA-map.png"),
)

# plot dv/v in LA 
Xticks = Date(2006):Year(2):Date(2021)
p2 = Plots.plot(
    dpi=500,
    size=(200,400),
    xlims=(Date(2006),Date(2021)),
    left_margin=5Plots.mm,
    xticks = (Xticks, Dates.format.(Xticks,"yyyy")),
    xrot=60,
)
# plot first dashed line 
Plots.hline!(
        p2,
        [size(LA1,1) * 10 + 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
for ii in 1:size(LA1,1)
    Plots.hline!(
        p2,
        [10 * ii],
        xlims=(Date(2006),Date(2021)),
        color=:black,
        label="",
    )
    Plots.hline!(
        p2,
        [10 * ii - 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
end
Plots.yticks!(p2,(10:10:10*size(LA1,1),LA1[:,:NETSTA]))
Plots.plot!(
    p2,
    [Date(2021),Date(2021)],
    [size(LA1,1) * 10 + 5, size(LA1,1) * 10 - 5],
    color=:black,
    lw=5,
    label="",
)
Plots.plot!(
    p2,
    [Date(2009),Date(2012,6,1)],
    [2.5, 2.5],
    color=:dodgerblue,
    lw=2,
    label="",
)
Plots.annotate!((Date(2021,6,1),size(LA1,1) * 10 + 5.5,Plots.text("1",10)))
Plots.annotate!((Date(2021,6,1),size(LA1,1) * 10 - 5.5,Plots.text("-1",10)))
Plots.annotate!((Date(2022,2,1),size(LA1,1) * 10 ,Plots.text("dv/v [%]",10,rotation = 90,:bold)))
Plots.annotate!((Date(2016,4,1),2.5,Plots.text("Scaled LWE",10)))
for ii in 1:size(LA1,1)
    df = Arrow.Table(LA1[ii,:FILES]) |> DataFrame 
    df = df[df[:,:DATE] .> Date(2006),:]
    Plots.scatter!(
        p2,
        df[:,:DATE],
        df[:,:DVV] .* 5 .+ 10 * ii,
        # alpha=df[:,:CC] .^ 3 ./ 20,
        alpha = 0.04,
        # markerz=df[:,:CC] .^ 3,
        c=Plots.theme_palette(:auto)[ii],
        label="",
        colorbar=false,
        markersize=3,
    )

    # get lwe time series 
    lonind = argmin(abs.(lwelon .- LA1[ii,:LON]))
    latind = argmin(abs.(lwelat .- LA1[ii,:LAT]))
    stalwe = lwe[lonind,latind,:] 
    # remove mean from 2002-2017 and 2018 - 2021
    stalwe .-= mean(stalwe)
    scaling = quantile(abs.(stalwe),0.95) ./ quantile(abs.(df[:,:DVV]),0.95)
    stalwe ./= scaling

    # plot 2002-2017
    Plots.plot!(
        p2,
        tlwe[1:163], 
        -stalwe[1:163] .* 5 .+ 10 * ii,
        lw=2, 
        alpha=0.85,
        label="",
        color=:dodgerblue,
    )

    # plot 2018-2021
    Plots.plot!(
        p2,
        tlwe[164:end], 
        -stalwe[164:end] .* 5 .+ 10 * ii,
        lw=2, 
        alpha=0.85,
        label="",
        color=:dodgerblue,
    )
end

Xticks = Date(2006):Year(2):Date(2021)
p3 = Plots.plot(
    dpi=500,
    size=(200,400),
    xlims=(Date(2006),Date(2021)),
    left_margin=5Plots.mm,
    xticks = (Xticks, Dates.format.(Xticks,"yyyy")),
    xrot=60,
)
# plot first dashed line 
Plots.hline!(
        p3,
        [size(LA2,1) * 10 + 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
for ii in 1:size(LA2,1)
    Plots.hline!(
        p3,
        [10 * ii],
        xlims=(Date(2006),Date(2021)),
        color=:black,
        label="",
    )
    Plots.hline!(
        p3,
        [10 * ii - 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
end
for ii in 1:size(LA2,1)
    df = Arrow.Table(LA2[ii,:FILES]) |> DataFrame 
    df = df[df[:,:DATE] .> Date(2006),:]

    Plots.scatter!(
        p3,
        df[:,:DATE],
        df[:,:DVV] .* 5 .+ 10 * ii,
        # alpha=df[:,:CC] .^ 3 ./ 20,
        alpha = 0.04,
        # markerz=df[:,:CC] .^ 3,
        c=Plots.theme_palette(:auto)[ii + 8],
        label="",
        colorbar=false,
        markersize=3,
    )

    # get lwe time series 
    lonind = argmin(abs.(lwelon .- LA2[ii,:LON]))
    latind = argmin(abs.(lwelat .- LA2[ii,:LAT]))
    stalwe = lwe[lonind,latind,:] 
    # remove mean from 2002-2017 and 2018 - 2021
    stalwe .-= mean(stalwe)
    scaling = quantile(abs.(stalwe),0.95) ./ quantile(abs.(df[:,:DVV]),0.95)
    stalwe ./= scaling

    # plot 2002-2017
    Plots.plot!(
        p3,
        tlwe[1:163], 
        -stalwe[1:163] .* 5 .+ 10 * ii,
        lw=2, 
        alpha=0.85,
        label="",
        color=:dodgerblue,
    )

    # plot 2018-2021
    Plots.plot!(
        p3,
        tlwe[164:end], 
        -stalwe[164:end] .* 5 .+ 10 * ii,
        lw=2, 
        alpha=0.85,
        label="",
        color=:dodgerblue,
    )
end
Plots.yticks!(p3,(10:10:10*size(LA2,1),LA2[:,:NETSTA]))

l1 = Plots.@layout [a b]
Plots.plot(p2,p3,layout=l1,size=(600,800),dpi=500)
Plots.savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/LA-dvv-lwe.png"))

# plot dv/v & Elastic Component in LA 
Xticks = Date(2006):Year(2):Date(2021)
p2 = Plots.plot(
    dpi=500,
    size=(200,400),
    xlims=(Date(2006),Date(2021)),
    left_margin=5Plots.mm,
    xticks = (Xticks, Dates.format.(Xticks,"yyyy")),
    xrot=60,
)
# plot first dashed line 
Plots.hline!(
        p2,
        [size(LA1,1) * 10 + 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
for ii in 1:size(LA1,1)
    Plots.hline!(
        p2,
        [10 * ii],
        xlims=(Date(2006),Date(2021)),
        color=:black,
        label="",
    )
    Plots.hline!(
        p2,
        [10 * ii - 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
end
Plots.yticks!(p2,(10:10:10*size(LA1,1),LA1[:,:NETSTA]))
Plots.plot!(
    p2,
    [Date(2021),Date(2021)],
    [size(LA1,1) * 10 + 5, size(LA1,1) * 10 - 5],
    color=:black,
    lw=5,
    label="",
)
Plots.annotate!((Date(2021,6,1),size(LA1,1) * 10 + 5.5,Plots.text("1",10)))
Plots.annotate!((Date(2021,6,1),size(LA1,1) * 10 - 5.5,Plots.text("-1",10)))
Plots.annotate!((Date(2022,2,1),size(LA1,1) * 10 ,Plots.text("dv/v [%]",10,rotation = 90,:bold)))
for ii in 1:size(LA1,1)
    df = Arrow.Table(LA1[ii,:FILES]) |> Arrow.columntable |> DataFrame 
    # df = df[df[:,:DATE] .> Date(2006),:]
    sort!(df,:DATE)
    dt = df[1,:DATE]:Day(1):df[end,:DATE]
    dtdf = DataFrame(DATE = dt)
    DVVmiss = outerjoin(dtdf,df,on=:DATE)
    sort!(DVVmiss,:DATE)
    DVVimpute = Impute.interp(DVVmiss)

    # find nearest precip grid cell 
    lonind = argmin(abs.(pptlon .- LA1[ii,:Longitude]))
    latind = argmin(abs.(pptlat .- LA1[ii,:Latitude]))
    precip = ppt[latind,lonind,:]
    precipind = findall(df[1,:DATE] .<= tppt .<= df[end,:DATE])
    tprecip = tppt[precipind]

    dvvE = elastic(precip[precipind],LA1[ii,:E3],α,r,δt)
    dvvE .-= mean(dvvE)
    dvvE./= std(dvvE)
    dvvE .*= LA1[ii,:E2]

    Plots.scatter!(
        p2,
        tprecip,
        dvvE .* 5 .+ 10 * ii,
        alpha=DVVimpute[:,:CC] .^ 3 ./ 20,
        # markerz=df[:,:CC] .^ 3,
        c=Plots.theme_palette(:auto)[ii],
        markersize=3,
        label="",
        colorbar=false,
        xlims=(Date(2006),Date(2021)),
    )
    # Plots.plot!(
    #     p2,
    #     pptfit[precipind][dvvshiftind],
    #     tempfilt[shiftind] .* coef(dvvfit)[3] .+ CDMbest[precipind][dvvshiftind] .* coef(dvvfit)[2] .+ coef(dvvfit)[1] ,
    #     alpha=0.75,
    #     label="",
    #     c=:dodgerblue,
    #     lw=3,
    # )
end

Xticks = Date(2006):Year(2):Date(2021)
p3 = Plots.plot(
    dpi=500,
    size=(200,400),
    xlims=(Date(2006),Date(2021)),
    left_margin=5Plots.mm,
    xticks = (Xticks, Dates.format.(Xticks,"yyyy")),
    xrot=60,
)
# plot first dashed line 
Plots.hline!(
        p3,
        [size(LA2,1) * 10 + 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
for ii in 1:size(LA2,1)
    Plots.hline!(
        p3,
        [10 * ii],
        xlims=(Date(2006),Date(2021)),
        color=:black,
        label="",
    )
    Plots.hline!(
        p3,
        [10 * ii - 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
end
for ii in 1:size(LA2,1)
    df = Arrow.Table(LA2[ii,:FILES]) |> Arrow.columntable |> DataFrame 
    df = df[df[:,:DATE] .> Date(2006),:]
    sort!(df,:DATE)
    dt = df[1,:DATE]:Day(1):df[end,:DATE]
    dtdf = DataFrame(DATE = dt)
    DVVmiss = outerjoin(dtdf,df,on=:DATE)
    sort!(DVVmiss,:DATE)
    DVVimpute = Impute.interp(DVVmiss)

    # find nearest precip grid cell 
    lonind = argmin(abs.(pptlon .- LA2[ii,:Longitude]))
    latind = argmin(abs.(pptlat .- LA2[ii,:Latitude]))
    precip = ppt[latind,lonind,:]
    precipind = findall(df[1,:DATE] .<= tppt .<= df[end,:DATE])
    tprecip = tppt[precipind]

    dvvE = elastic(precip[precipind],LA2[ii,:E3],α,r,δt)
    dvvE .-= mean(dvvE)
    dvvE./= std(dvvE)
    dvvE .*= LA2[ii,:E2]

    Plots.scatter!(
        p3,
        tprecip,
        dvvE .* 5 .+ 10 * ii,
        alpha=DVVimpute[:,:CC] .^ 3 ./ 20,
        # markerz=df[:,:CC] .^ 3,
        c=Plots.theme_palette(:auto)[ii + 8],
        markersize=3,
        label="",
        colorbar=false,
        xlims=(Date(2006),Date(2021)),
    )
    # Plots.plot!(
    #     p2,
    #     pptfit[precipind][dvvshiftind],
    #     tempfilt[shiftind] .* coef(dvvfit)[3] .+ CDMbest[precipind][dvvshiftind] .* coef(dvvfit)[2] .+ coef(dvvfit)[1] ,
    #     alpha=0.75,
    #     label="",
    #     c=:dodgerblue,
    #     lw=3,
    # )
end
Plots.yticks!(p3,(10:10:10*size(LA2,1),LA2[:,:NETSTA]))

l1 = Plots.@layout [a b]
Plots.plot(p2,p3,layout=l1,size=(600,800),dpi=500)
Plots.savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/LA-dvv-elastic.png"))

# plot dv/v & temp Component in LA 
Xticks = Date(2006):Year(2):Date(2021)
p2 = Plots.plot(
    dpi=500,
    size=(200,400),
    xlims=(Date(2006),Date(2021)),
    left_margin=5Plots.mm,
    xticks = (Xticks, Dates.format.(Xticks,"yyyy")),
    xrot=60,
)
# plot first dashed line 
Plots.hline!(
        p2,
        [size(LA1,1) * 10 + 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
for ii in 1:size(LA1,1)
    Plots.hline!(
        p2,
        [10 * ii],
        xlims=(Date(2006),Date(2021)),
        color=:black,
        label="",
    )
    Plots.hline!(
        p2,
        [10 * ii - 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
end
Plots.yticks!(p2,(10:10:10*size(LA1,1),LA1[:,:NETSTA]))
Plots.plot!(
    p2,
    [Date(2021),Date(2021)],
    [size(LA1,1) * 10 + 5, size(LA1,1) * 10 - 5],
    color=:black,
    lw=5,
    label="",
)
Plots.annotate!((Date(2021,6,1),size(LA1,1) * 10 + 5.5,Plots.text("1",10)))
Plots.annotate!((Date(2021,6,1),size(LA1,1) * 10 - 5.5,Plots.text("-1",10)))
Plots.annotate!((Date(2022,2,1),size(LA1,1) * 10 ,Plots.text("dv/v [%]",10,rotation = 90,:bold)))
for ii in 1:size(LA1,1)
    df = Arrow.Table(LA1[ii,:FILES]) |> Arrow.columntable |> DataFrame 
    # df = df[df[:,:DATE] .> Date(2006),:]
    sort!(df,:DATE)
    dt = df[1,:DATE]:Day(1):df[end,:DATE]
    dtdf = DataFrame(DATE = dt)
    DVVmiss = outerjoin(dtdf,df,on=:DATE)
    sort!(DVVmiss,:DATE)
    DVVimpute = Impute.interp(DVVmiss)

    # find nearest precip grid cell 
    lonind = argmin(abs.(pptlon .- LA1[ii,:Longitude]))
    latind = argmin(abs.(pptlat .- LA1[ii,:Latitude]))
    temp = tmean[latind,lonind,:]
    temp .-= mean(temp)
    tempfilt = smooth_withfiltfilt(temp,window_len=45)
    tempind = findall(df[1,:DATE] .<= ttmean .<= df[end,:DATE])
    ttemp = ttmean[tempind]

    dvvT = tempfilt[tempind .- round(Int,LA1[ii,:E5])] 
    dvvT .-= mean(dvvT)
    dvvT ./= std(dvvT)
    dvvT .*= LA1[ii,:E4]

    Plots.scatter!(
        p2,
        ttemp,
        dvvT .* 5 .+ 10 * ii,
        alpha=DVVimpute[:,:CC] .^ 3 ./ 20,
        # markerz=df[:,:CC] .^ 3,
        c=Plots.theme_palette(:auto)[ii],
        markersize=3,
        label="",
        colorbar=false,
    )
    # Plots.plot!(
    #     p2,
    #     pptfit[precipind][dvvshiftind],
    #     tempfilt[shiftind] .* coef(dvvfit)[3] .+ CDMbest[precipind][dvvshiftind] .* coef(dvvfit)[2] .+ coef(dvvfit)[1] ,
    #     alpha=0.75,
    #     label="",
    #     c=:dodgerblue,
    #     lw=3,
    # )
end
Plots.yticks!(p2,(10:10:10*size(LA1,1),LA1[:,:NETSTA]))

Xticks = Date(2006):Year(2):Date(2021)
p3 = Plots.plot(
    dpi=500,
    size=(200,400),
    xlims=(Date(2006),Date(2021)),
    left_margin=5Plots.mm,
    xticks = (Xticks, Dates.format.(Xticks,"yyyy")),
    xrot=60,
)
# plot first dashed line 
Plots.hline!(
        p3,
        [size(LA2,1) * 10 + 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
for ii in 1:size(LA2,1)
    Plots.hline!(
        p3,
        [10 * ii],
        xlims=(Date(2006),Date(2021)),
        color=:black,
        label="",
    )
    Plots.hline!(
        p3,
        [10 * ii - 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
end
for ii in 1:size(LA2,1)
    df = Arrow.Table(LA2[ii,:FILES]) |> DataFrame 
    df = df[df[:,:DATE] .> Date(2006),:]
    sort!(df,:DATE)
    dt = df[1,:DATE]:Day(1):df[end,:DATE]
    dtdf = DataFrame(DATE = dt)
    DVVmiss = outerjoin(dtdf,df,on=:DATE)
    sort!(DVVmiss,:DATE)
    DVVimpute = Impute.interp(DVVmiss)

    # find nearest precip grid cell 
    lonind = argmin(abs.(pptlon .- LA2[ii,:Longitude]))
    latind = argmin(abs.(pptlat .- LA2[ii,:Latitude]))
    temp = tmean[latind,lonind,:]
    temp .-= mean(temp)
    tempfilt = smooth_withfiltfilt(temp,window_len=45)
    tempind = findall(df[1,:DATE] .<= ttmean .<= df[end,:DATE])
    ttemp = ttmean[tempind]

    dvvT = tempfilt[tempind .- round(Int,LA2[ii,:E5])] 
    dvvT .-= mean(dvvT)
    dvvT ./= std(dvvT)
    dvvT .*= LA2[ii,:E4]

    Plots.scatter!(
        p3,
        ttemp,
        dvvT .* 5 .+ 10 * ii,
        alpha=DVVimpute[:,:CC] .^ 3 ./ 20,
        # markerz=df[:,:CC] .^ 3,
        c=Plots.theme_palette(:auto)[ii + 8],
        markersize=3,
        label="",
        colorbar=false,
    )
    # Plots.plot!(
    #     p2,
    #     pptfit[precipind][dvvshiftind],
    #     tempfilt[shiftind] .* coef(dvvfit)[3] .+ CDMbest[precipind][dvvshiftind] .* coef(dvvfit)[2] .+ coef(dvvfit)[1] ,
    #     alpha=0.75,
    #     label="",
    #     c=:dodgerblue,
    #     lw=3,
    # )
end
Plots.yticks!(p3,(10:10:10*size(LA2,1),LA2[:,:NETSTA]))

l1 = Plots.@layout [a b]
Plots.plot(p2,p3,layout=l1,size=(600,800),dpi=500)
Plots.savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/LA-dvv-temp.png"))


# plot dv/v & residual Component in LA 
Xticks = Date(2006):Year(2):Date(2021)
p2 = Plots.plot(
    dpi=500,
    size=(200,400),
    xlims=(Date(2006),Date(2021)),
    left_margin=5Plots.mm,
    xticks = (Xticks, Dates.format.(Xticks,"yyyy")),
    xrot=60,
)
# plot first dashed line 
Plots.hline!(
        p2,
        [size(LA1,1) * 10 + 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
for ii in 1:size(LA1,1)
    Plots.hline!(
        p2,
        [10 * ii],
        xlims=(Date(2006),Date(2021)),
        color=:black,
        label="",
    )
    Plots.hline!(
        p2,
        [10 * ii - 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
end
Plots.yticks!(p2,(10:10:10*size(LA1,1),LA1[:,:NETSTA]))
Plots.plot!(
    p2,
    [Date(2021),Date(2021)],
    [size(LA1,1) * 10 + 5, size(LA1,1) * 10 - 5],
    color=:black,
    lw=5,
    label="",
)
Plots.annotate!((Date(2021,6,1),size(LA1,1) * 10 + 5.5,Plots.text("1",10)))
Plots.annotate!((Date(2021,6,1),size(LA1,1) * 10 - 5.5,Plots.text("-1",10)))
Plots.annotate!((Date(2022,2,1),size(LA1,1) * 10 ,Plots.text("dv/v [%]",10,rotation = 90,:bold)))
for ii in 1:size(LA1,1)
    df = Arrow.Table(LA1[ii,:FILES]) |>  Arrow.columntable |> DataFrame 
    # df = df[df[:,:DATE] .> Date(2006),:]
    sort!(df,:DATE)
    dt = df[1,:DATE]:Day(1):df[end,:DATE]
    dtdf = DataFrame(DATE = dt)
    DVVmiss = outerjoin(dtdf,df,on=:DATE)
    sort!(DVVmiss,:DATE)
    DVVimpute = Impute.interp(DVVmiss)

    # find nearest precip grid cell 
    lonind = argmin(abs.(pptlon .- LA1[ii,:Longitude]))
    latind = argmin(abs.(pptlat .- LA1[ii,:Latitude]))
    precip = ppt[latind,lonind,:]
    precipind = findall(df[1,:DATE] .<= tppt .<= df[end,:DATE])
    tprecip = tppt[precipind]
    temp = tmean[latind,lonind,:]
    temp .-= mean(temp)
    temp ./= std(temp)
    tempfilt = smooth_withfiltfilt(temp,window_len=45)
    tempind = findall(df[1,:DATE] .<= ttmean .<= df[end,:DATE])
    ttemp = ttmean[tempind]

    # fitted dvv 
    tminind = findfirst(tppt .== DVVimpute[1,:DATE])
    tmaxind = findlast(tppt .== DVVimpute[end, :DATE])
    cdmk = CDM(precip[tminind-ceil(Int,LA1[ii,:CDM3]):tmaxind],ceil(Int,LA1[ii,:CDM3]))[ceil(Int,LA1[ii,:CDM3])+1:end]
    cdmk .-= mean(cdmk)
    cdmk ./= std(cdmk)
    
    dvvE = elastic(precip[precipind],LA1[ii,:E3],α,r,δt)
    dvvE .-= mean(dvvE)
    dvvE ./= std(dvvE)
    dvvT = tempfilt[tempind .- round(Int,LA1[ii,:E5])] 
    dvvET = LA1[ii,:E1] .+  LA1[ii,:E2] .* dvvE .+ LA1[ii,:E4] .* dvvT
    dvvCT = LA1[ii,:CDM1] .+  LA1[ii,:CDM2] .* cdmk .+ LA1[ii,:CDM4] .* dvvT
    dvvind = findall(in(df[:,:DATE]),ttemp)

    Plots.scatter!(
        p2,
        ttemp[dvvind],
        (dvvCT[dvvind] .- df[:,:DVV]) .* 5 .+ 10 * ii,
        alpha=df[:,:CC] .^ 3 ./ 20,
        # markerz=df[:,:CC] .^ 3,
        c=Plots.theme_palette(:auto)[ii],
        markersize=3,
        label="",
        colorbar=false,
        xlims=(Date(2006),Date(2021)),
    )
    # Plots.plot!(
    #     p2,
    #     pptfit[precipind][dvvshiftind],
    #     tempfilt[shiftind] .* coef(dvvfit)[3] .+ CDMbest[precipind][dvvshiftind] .* coef(dvvfit)[2] .+ coef(dvvfit)[1] ,
    #     alpha=0.75,
    #     label="",
    #     c=:dodgerblue,
    #     lw=3,
    # )
end
Plots.yticks!(p2,(10:10:10*size(LA1,1),LA1[:,:NETSTA]))

Xticks = Date(2006):Year(2):Date(2021)
p3 = Plots.plot(
    dpi=500,
    size=(200,400),
    xlims=(Date(2006),Date(2021)),
    left_margin=5Plots.mm,
    xticks = (Xticks, Dates.format.(Xticks,"yyyy")),
    xrot=60,
)
# plot first dashed line 
Plots.hline!(
        p3,
        [size(LA2,1) * 10 + 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
for ii in 1:size(LA2,1)
    Plots.hline!(
        p3,
        [10 * ii],
        xlims=(Date(2006),Date(2021)),
        color=:black,
        label="",
    )
    Plots.hline!(
        p3,
        [10 * ii - 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
end
for ii in 1:size(LA2,1)
    df = Arrow.Table(LA2[ii,:FILES]) |> Arrow.columntable |> DataFrame 
    # df = df[df[:,:DATE] .> Date(2006),:]
    sort!(df,:DATE)
    dt = df[1,:DATE]:Day(1):df[end,:DATE]
    dtdf = DataFrame(DATE = dt)
    DVVmiss = outerjoin(dtdf,df,on=:DATE)
    sort!(DVVmiss,:DATE)
    DVVimpute = Impute.interp(DVVmiss)

     # find nearest precip grid cell 
     lonind = argmin(abs.(pptlon .- LA2[ii,:Longitude]))
     latind = argmin(abs.(pptlat .- LA2[ii,:Latitude]))
     precip = ppt[latind,lonind,:]
     precipind = findall(df[1,:DATE] .<= tppt .<= df[end,:DATE])
     tprecip = tppt[precipind]
     temp = tmean[latind,lonind,:]
     temp .-= mean(temp)
     temp ./= std(temp)
     tempfilt = smooth_withfiltfilt(temp,window_len=45)
     tempind = findall(df[1,:DATE] .<= ttmean .<= df[end,:DATE])
     ttemp = ttmean[tempind]
 
     # fitted dvv 
     tminind = findfirst(tppt .== DVVimpute[1,:DATE])
    tmaxind = findlast(tppt .== DVVimpute[end, :DATE])
    cdmk = CDM(precip[tminind-ceil(Int,LA2[ii,:CDM3]):tmaxind],ceil(Int,LA2[ii,:CDM3]))[ceil(Int,LA2[ii,:CDM3])+1:end]
    cdmk .-= mean(cdmk)
    cdmk ./= std(cdmk)

     dvvE = elastic(precip[precipind],LA2[ii,:E3],α,r,δt)
     dvvE .-= mean(dvvE)
     dvvE ./= std(dvvE)
     dvvTE = tempfilt[tempind .- round(Int,LA2[ii,:E5])] 
     dvvTC = tempfilt[tempind .- round(Int,LA2[ii,:CDM5])] 
     dvvET = LA2[ii,:E1] .+  LA2[ii,:E2] .* dvvE .+ LA2[ii,:E4] .* dvvTE
     dvvCT = LA2[ii,:CDM1] .+  LA2[ii,:CDM2] .* dvvE .+ LA2[ii,:CDM4] .* dvvTC
     dvvind = findall(in(df[:,:DATE]),ttemp)
 
     Plots.scatter!(
         p3,
         ttemp[dvvind],
         (dvvCT[dvvind] .- df[:,:DVV]) .* 5 .+ 10 * ii,
         alpha=df[:,:CC] .^ 3 ./ 20,
         # markerz=df[:,:CC] .^ 3,
         c=Plots.theme_palette(:auto)[ii],
         markersize=3,
         label="",
         colorbar=false,
         xlims=(Date(2006),Date(2021)),
     )
    # Plots.plot!(
    #     p2,
    #     pptfit[precipind][dvvshiftind],
    #     tempfilt[shiftind] .* coef(dvvfit)[3] .+ CDMbest[precipind][dvvshiftind] .* coef(dvvfit)[2] .+ coef(dvvfit)[1] ,
    #     alpha=0.75,
    #     label="",
    #     c=:dodgerblue,
    #     lw=3,
    # )
end
Plots.yticks!(p3,(10:10:10*size(LA2,1),LA2[:,:NETSTA]))

l1 = Plots.@layout [a b]
Plots.plot(p2,p3,layout=l1,size=(600,800),dpi=500)
Plots.savefig(joinpath(@__DIR,"../data/FINAL-FIGURES/LA-dvv-residual.png"))


# plot dv/v with CDMk 
Xticks = Date(2006):Year(2):Date(2021)
p2 = Plots.plot(
    dpi=500,
    size=(200,400),
    xlims=(Date(2006),Date(2021)),
    left_margin=5Plots.mm,
    xticks = (Xticks, Dates.format.(Xticks,"yyyy")),
    xrot=60,
)
# plot first dashed line 
Plots.hline!(
        p2,
        [size(LA1,1) * 10 + 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
for ii in 1:size(LA1,1)
    Plots.hline!(
        p2,
        [10 * ii],
        xlims=(Date(2006),Date(2021)),
        color=:black,
        label="",
    )
    Plots.hline!(
        p2,
        [10 * ii - 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
end
for ii in 1:size(LA1,1)
    df = Arrow.Table(LA1[ii,:FILES]) |> DataFrame 
    df = df[df[:,:DATE] .> Date(2006),:]

    # find nearest precip grid cell 
    lonind = argmin(abs.(pptlon .- LA1[ii,:Longitude]))
    latind = argmin(abs.(pptlat .- LA1[ii,:Latitude]))
    precip = ppt[latind,lonind,:]
    temp = tmean[latind,lonind,:]
    temp .-= mean(temp)
    tempfilt = filtfilt(digitalfilter(responsetype, designmethod), temp)

    # find best fitting temperature 
    maxlag = isodd(length(df[!,:DVV])) ? size(df[!,:DVV],1) : size(df[!,:DVV],1) - 1
    lags = -maxlag:maxlag
    lagind = findall((lags .<= 360) .& (lags .>= 0))
    tempind = findall(in(df[!,:DATE]),ttmean)
    tempdvvcorr = xcorr(tempfilt[tempind],df[!,:DVV])[lagind]
    tempcor = cor(df[!,:DVV],tempfilt[tempind])
    tempshift = argmax(tempdvvcorr) 

    # find common times between DV/V and shifted temperature 
    shiftind = findall(in(df[!,:DATE]),ttmean .- Day(tempshift))
    dvvshiftind = findall(in(ttmean .- Day(tempshift)),df[:,:DATE])
    tempcor = cor(df[dvvshiftind,:DVV],tempfilt[shiftind])
    println("$(LA2[ii,:NETSTA]) delay=$(lags[lagind][tempshift]), CC=$tempcor")

    # find best fitting precip model against dv/v 
    days = 180 : 10 : 20 * 365
    N = length(days)
    cors = zeros(N)
    pptminind = findfirst(tppt .== df[1,:DATE])
    pptmaxind = findfirst(tppt .== df[end,:DATE])
    # find best CDM date range 
    for jj in 1:N
        CDMii = CDM(precip[pptminind - days[jj]:end],days[jj])
        pptind = findall(in(df[!,:DATE]),tppt[pptminind - days[jj]:end])
        cors[jj] = cor(df[!,:DVV],CDMii[pptind])
    end
    CDMbest = CDM(precip[pptminind - days[argmin(cors)]:end],days[argmin(cors)]) ./ 1000
    CDMbest .-= CDMbest[1]
    pptfit = tppt[pptminind - days[argmin(cors)]:end]
    precipind = findall(in(df[!,:DATE]),pptfit)

    # fit multimodel to dv/v 
    @. multimodel(x, p) = p[1] + p[2] * x[:,1] + p[3] * x[:,2]
    dvvfit = curve_fit(
        multimodel,
        [CDMbest[precipind][dvvshiftind] tempfilt[shiftind]],
        df[dvvshiftind,:DVV] .* 5 .+ 10 * ii,
        [0.,-1.,1.],
    )
    println("$(LA1[ii,:NETSTA]) CDMk=$(coef(dvvfit)[2]), TEMP=$(coef(dvvfit)[3])")

    Plots.scatter!(
        p2,
        df[:,:DATE],
        df[:,:DVV] .* 5 .+ 10 * ii,
        alpha=df[:,:CC] .^ 3 ./ 20,
        # markerz=df[:,:CC] .^ 3,
        c=Plots.theme_palette(:auto)[ii],
        markersize=3,
        label="",
        colorbar=false,
    )
    Plots.plot!(
        p2,
        pptfit[precipind][dvvshiftind],
        tempfilt[shiftind] .* coef(dvvfit)[3] .+ CDMbest[precipind][dvvshiftind] .* coef(dvvfit)[2] .+ coef(dvvfit)[1] ,
        alpha=0.75,
        label="",
        c=:dodgerblue,
        lw=3,
    )
end
Plots.yticks!(p2,(10:10:10*size(LA1,1),LA1[:,:NETSTA]))

Xticks = Date(2006):Year(2):Date(2021)
p3 = Plots.plot(
    dpi=500,
    size=(200,400),
    xlims=(Date(2006),Date(2021)),
    left_margin=5Plots.mm,
    xticks = (Xticks, Dates.format.(Xticks,"yyyy")),
    xrot=60,
)
# plot first dashed line 
Plots.hline!(
        p3,
        [size(LA2,1) * 10 + 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
for ii in 1:size(LA2,1)
    Plots.hline!(
        p3,
        [10 * ii],
        xlims=(Date(2006),Date(2021)),
        color=:black,
        label="",
    )
    Plots.hline!(
        p3,
        [10 * ii - 5],
        xlims=(Date(2006),Date(2021)),
        color=:gray,
        linestyle=:dash,
        label="",
    )
end
for ii in 1:size(LA2,1)
    df = Arrow.Table(LA2[ii,:FILES]) |> DataFrame 
    df = df[df[:,:DATE] .> Date(2006),:]

    # find nearest precip grid cell 
    lonind = argmin(abs.(pptlon .- LA2[ii,:Longitude]))
    latind = argmin(abs.(pptlat .- LA2[ii,:Latitude]))
    precip = ppt[latind,lonind,:]
    temp = tmean[latind,lonind,:]
    temp .-= mean(temp)
    tempfilt = filtfilt(digitalfilter(responsetype, designmethod), temp)

    # find best fitting temperature 
    maxlag = isodd(length(df[!,:DVV])) ? size(df[!,:DVV],1) : size(df[!,:DVV],1) - 1
    lags = -maxlag:maxlag
    lagind = findall((lags .<= 360) .& (lags .>= 0))
    tempind = findall(in(df[!,:DATE]),ttmean)
    tempdvvcorr = xcorr(tempfilt[tempind],df[!,:DVV])[lagind]
    tempcor = cor(df[!,:DVV],tempfilt[tempind])
    tempshift = argmax(tempdvvcorr) 

    # find common times between DV/V and shifted temperature 
    shiftind = findall(in(df[!,:DATE]),ttmean .- Day(tempshift))
    dvvshiftind = findall(in(ttmean .- Day(tempshift)),df[:,:DATE])
    tempcor = cor(df[dvvshiftind,:DVV],tempfilt[shiftind])
    println("$(LA2[ii,:NETSTA]) delay=$(lags[lagind][tempshift]), CC=$tempcor")

    # find best fitting precip model against dv/v 
    days = 180 : 10 : 20 * 365
    N = length(days)
    cors = zeros(N)
    pptminind = findfirst(tppt .== df[1,:DATE])
    pptmaxind = findfirst(tppt .== df[end,:DATE])
    # find best CDM date range 
    for jj in 1:N
        CDMii = CDM(precip[pptminind - days[jj]:end],days[jj])
        pptind = findall(in(df[!,:DATE]),tppt[pptminind - days[jj]:end])
        cors[jj] = cor(df[!,:DVV],CDMii[pptind])
    end
    CDMbest = CDM(precip[pptminind - days[argmin(cors)]:end],days[argmin(cors)]) ./ 1000
    CDMbest .-= CDMbest[1]
    pptfit = tppt[pptminind - days[argmin(cors)]:end]
    precipind = findall(in(df[!,:DATE]),pptfit)

    # fit multimodel to dv/v 
    @. multimodel(x, p) = p[1] + p[2] * x[:,1] + p[3] * x[:,2]
    dvvfit = curve_fit(
        multimodel,
        [CDMbest[precipind][dvvshiftind] tempfilt[shiftind]],
        df[dvvshiftind,:DVV] .* 5 .+ 10 * ii,
        [0.,-1.,1.],
    )
    println("$(LA2[ii,:NETSTA]) CDMk=$(coef(dvvfit)[2]), TEMP=$(coef(dvvfit)[3])")

    Plots.scatter!(
        p3,
        df[:,:DATE],
        df[:,:DVV] .* 5 .+ 10 * ii,
        alpha=df[:,:CC] .^ 3 ./ 20,
        # markerz=df[:,:CC] .^ 3,
        c=Plots.theme_palette(:auto)[ii + 8],
        label="",
        markersize=3,
        colorbar=false,
    )
    Plots.plot!(
        p3,
        pptfit[precipind][dvvshiftind],
        tempfilt[shiftind] .* coef(dvvfit)[3] .+ CDMbest[precipind][dvvshiftind] .* coef(dvvfit)[2] .+ coef(dvvfit)[1] ,
        alpha=0.75,
        label="",
        c=:dodgerblue,
        lw=3,
    )
end
Plots.yticks!(p3,(10:10:10*size(LA2,1),LA2[:,:NETSTA]))

l1 = Plots.@layout [a b]
Plots.plot(p2,p3,layout=l1,size=(600,800),dpi=500)
Plots.savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/LA-dvv-fit.png"))