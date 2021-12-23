# this script: 
# 1. Defines the five hydrology models 
# 2. Defines the thermo-elastic model 
# 3. fits a combined thermo+hydro-elastic/poroelastic model to each dv/v time series 

using Arrow 
using CSV
using DataFrames 
using Dates 
using DSP 
using Glob 
using Impute 
using NetCDF
using Optim 
using Plots
using SpecialFunctions
using Statistics

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

# number of days 
days = 90 

# load dv/v data 
files = glob("*",joinpath(@__DIR__,"../data/DVV-$days-DAY-COMP/2.0-4.0"))
FITDIR = joinpath(@__DIR__,"../data/FIT-DVV-SSE/$days-DAY")
FIGDIR = joinpath(@__DIR__,"../data/FIGURES-DVV-SSE/$days-DAY")
if !isdir(FITDIR)
    mkpath(FITDIR)
end
if !isdir(FIGDIR)
    mkpath(FIGDIR)
end

# load precip data 
filename = joinpath(@__DIR__,"../data/ppt.nc")
lon = ncread(filename,"lon")
lat = ncread(filename,"lat")
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

# get CI station locations 
SCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/CIstations.csv")))
NCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../data/NCstations.csv")))
CAdf = vcat(NCdf, SCdf)

# df to hold parameters 
fitdf = DataFrame()

# constants for models 
r = 500.
δt = 86400.
α = (1 + 0.27) / 3 / (1 - 0.27)
ϕ = 0.15   # porosity 

for kk in 1:length(files)
    file = files[kk]
    DVVdf = Arrow.Table(files[kk]) |> Arrow.columntable |> DataFrame
    sort!(DVVdf,:DATE)
    dt = DVVdf[1,:DATE]:Day(1):DVVdf[end,:DATE]
    DVVdt = DataFrame(DATE = dt)
    DVVmiss = outerjoin(DVVdt,DVVdf,on=:DATE)
    sort!(DVVmiss,:DATE)
    DVVimpute = Impute.interp(DVVmiss)

    netsta = replace(basename(file),".arrow"=>"")
    net, sta = split(netsta,'.')
    
    if size(DVVdf,1) < 365
        continue
    end

    # get precip nearest station lat, lon 
    staind = findfirst(CAdf[!,:Station] .== sta)
    if isnothing(staind)
        continue
    end
    stalat = CAdf[staind,:Latitude]
    stalon = CAdf[staind,:Longitude]
    staele = CAdf[staind,:Elevation]
    lonind = argmin(abs.(lon .- stalon))
    latind = argmin(abs.(lat .- stalat))
    precip = ppt[latind,lonind,:]
    temp = tmean[latind,lonind,:]
    temp .-= mean(temp)
    temp ./= std(temp)
    smoothtemp = smooth_withfiltfilt(temp,window_len=days÷2)
    tempind = findall(in(DVVimpute[:,:DATE]),ttmean)
    tminind = findfirst(tppt .== DVVimpute[1,:DATE])
    tmaxind = findlast(tppt .== DVVimpute[end, :DATE])

    # skip stations that are outside PRISM for now 
    if iszero(temp)
        continue
    end

    # define functions for fitting 
    function modelCDMk(p)
        cdmk = CDM(precip[tminind-ceil(Int,p[3]):tmaxind],ceil(Int,p[3]))[ceil(Int,p[3])+1:end]
        cdmk .-= mean(cdmk)
        cdmk ./= std(cdmk)
        ypred = p[1] .+ p[2] .* cdmk .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
        return  sum((ypred .- DVVimpute[:,:DVV]) .^ 2)
    end

    function modelSSW(p)
        ssw = SSW06(precip[tminind:tmaxind],ϕ,p[3])
        ssw .-= mean(ssw)
        ssw ./= std(ssw)
        ypred = p[1] .+ p[2] .* ssw .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
        return  sum((ypred .- DVVimpute[:,:DVV]) .^ 2)
    end
    
    function modelElastic(p)
        E = elastic(precip[tminind:tmaxind],p[3],α,r,δt)
        E .-= mean(E)
        E ./= std(E)
        ypred = p[1] .+ p[2] .* E .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
        return  sum((ypred .- DVVimpute[:,:DVV]) .^ 2)
    end
    
    function modelFC(p)
        FC = fullycoupled(precip[tminind:tmaxind],p[3],α,r,δt)
        FC .-= mean(FC)
        FC ./= std(FC)
        ypred = p[1] .+ p[2] .* FC .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
        return  sum((ypred .- DVVimpute[:,:DVV]) .^ 2)
    end
    
    function modelDrained(p)
        D = drained(precip[tminind:tmaxind],p[3],r,δt)
        D .-= mean(D)
        D ./= std(D)
        ypred = p[1] .+ p[2] .* D .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
        return  sum((ypred .- DVVimpute[:,:DVV]) .^ 2)
    end    

    # solve for elastic model 
    resElastic = optimize(
        modelElastic, 
        [-Inf,-Inf,5e-5,-Inf,0],
        [Inf,Inf,Inf,Inf,200],
        [0.0, -0.001,0.01,1.0,55],
        Fminbox(LBFGS()),
        Optim.Options(time_limit=60.0),
    )
    pE = Optim.minimizer(resElastic)
    E = elastic(precip[tminind:tmaxind],pE[3],α,r,δt)
    E .-= mean(E)
    E ./= std(E)
    ypredElastic = pE[1] .+ pE[2] .* E .+ pE[4] .* smoothtemp[tempind .- round(Int,pE[5])]
    corE = cor(ypredElastic,DVVimpute[:,:DVV])
    r2E = corE ^ 2

    # solve for fully-coupled model 
    resFC = optimize(
        modelFC, 
        [-Inf,-Inf,5e-5,-Inf,0],
        [Inf,Inf,Inf,Inf,200],
        [0.0, -0.001,0.01,1.0,55],
        Fminbox(LBFGS()),
        Optim.Options(time_limit=60.0),
    )
    pFC = Optim.minimizer(resFC)
    FC = fullycoupled(precip[tminind:tmaxind],pFC[3],α,r,δt)
    FC .-= mean(FC)
    FC ./= std(FC)
    ypredFC = pFC[1] .+ pFC[2] .* FC .+ pFC[4] .* smoothtemp[tempind .- round(Int,pFC[5])]
    corFC = cor(ypredFC,DVVimpute[:,:DVV])
    r2FC = corFC ^ 2

    # solve for drained model 
    resDrained = optimize(
        modelDrained, 
        [-Inf,-Inf,5e-5,-Inf,0],
        [Inf,Inf,100.0,Inf,200],
        [0.0, -0.001,0.01,1.0,55],
        Fminbox(LBFGS()),
        Optim.Options(time_limit=60.0),
    )
    pD = Optim.minimizer(resDrained)
    D = drained(precip[tminind:tmaxind],pD[3],r,δt)
    D .-= mean(D)
    D ./= std(D)
    ypredDrained = pD[1] .+ pD[2] .* D .+ pD[4] .* smoothtemp[tempind .- round(Int,pD[5])]
    corD = cor(ypredDrained,DVVimpute[:,:DVV])
    r2D = corD ^ 2

    # solve for baseflow model 
    resSSW = optimize(
        modelSSW, 
        [-Inf,-Inf,5e-5,-Inf,0],
        [Inf,Inf,1.,Inf,200],
        [0.0, -0.001,0.01,1.0,55],
        Fminbox(LBFGS()),
        Optim.Options(time_limit=60.0),
    )
    pS = Optim.minimizer(resSSW)
    ssw = SSW06(precip[tminind:tmaxind],ϕ,pS[3])
    ssw .-= mean(ssw)
    ssw ./= std(ssw)
    ypredSSW = pS[1] .+ pS[2] .* ssw .+ pS[4] .* smoothtemp[tempind .- round(Int,pS[5])]
    corSSW = cor(ypredSSW,DVVimpute[:,:DVV])
    r2SSW = corSSW ^ 2

    # solve for CDMk model 
    resCDM = optimize(
        modelCDMk, 
        [-Inf,-Inf,365. ,-Inf,0],
        [Inf,Inf,365.0 * 14,Inf,200],
        [0.0, -0.001,365.0 * 8.,1.0,55],
        Fminbox(LBFGS()),
        Optim.Options(time_limit=60.0),
    )
    pC = Optim.minimizer(resCDM)
    cdmk = CDM(precip[tminind-ceil(Int,pC[3]):tmaxind],ceil(Int,pC[3]))[ceil(Int,pC[3])+1:end]
    cdmk .-= mean(cdmk)
    cdmk ./= std(cdmk)
    ypredCDM = pC[1] .+ pC[2] .* cdmk .+ pC[4] .* smoothtemp[tempind .- round(Int,pC[5])]
    corCDM = cor(ypredCDM,DVVimpute[:,:DVV])
    r2CDM = corCDM ^ 2

    # fitted dvv 
    stationdf = DataFrame(
        Dict(
            :DATE=>DVVimpute[:,:DATE],
            :DVV=>DVVimpute[:,:DVV],
            :CC=>DVVimpute[:,:CC],
            :CDM=>ypredCDM,
            :ELASTIC=>ypredElastic,
            :DRAINED=>ypredDrained,
            :FC=>ypredFC,
            :SSW=>ypredSSW,
        )
    )
    fitpath = joinpath(FITDIR,netsta*".arrow")
    Arrow.write(fitpath,stationdf)
    
    # write to disk 
    append!(
        fitdf,
        DataFrame(
            Dict(
                :NETSTA=>netsta,
                :CDM1=>pC[1],
                :CDM2=>pC[2],
                :CDM3=>pC[3],
                :CDM4=>pC[4],
                :CDM5=>pC[5],
                :D1=>pD[1],
                :D2=>pD[2],
                :D3=>pD[3],
                :D4=>pD[4],
                :D5=>pD[5],
                :E1=>pE[1],
                :E2=>pE[2],
                :E3=>pE[3],
                :E4=>pE[4],
                :E5=>pE[5],
                :FC1=>pFC[1],
                :FC2=>pFC[2],
                :FC3=>pFC[3],
                :FC4=>pFC[4],
                :FC5=>pFC[5],
                :SSW1=>pS[1],
                :SSW2=>pS[2],
                :SSW3=>pS[3],
                :SSW4=>pS[4],
                :SSW5=>pS[5],
                :corE=>corE,
                :corD=>corD,
                :corFC=>corFC,
                :corCDM=>corCDM,
                :corSSW=>corSSW,
                :r2E=>r2E,
                :r2D=>r2D,
                :r2FC=>r2FC,
                :r2CDM=>r2CDM,
                :r2SSW=>r2SSW,
                :LAT=>stalat,
                :LON=>stalon,
                :ELEV=>staele,
            )
        )
    )
    println("$(netsta)")
    println("--------------------")
    println("Elastic = $(round(corE,digits=2))")
    println("Drained = $(round(corD,digits=2))")
    println("Fully-Coupled = $(round(corFC,digits=2))")
    println("CDM = $(round(corCDM,digits=2))")
    println("SSW = $(round(corSSW,digits=2))\n\n")

    # plot each model 
    plot(DVVimpute[:,:DATE],
        DVVimpute[:,:DVV],
        legend=:outerbottom,
        label="dv/v",
        ylabel="dv/v [%]",
        legendfontsize=4,
        figsize=(800,400),
        dpi=250,
    )
    plot!(DVVimpute[:,:DATE],ypredElastic,label="Elastic = $(round(corE,digits=2))")
    plot!(DVVimpute[:,:DATE],ypredDrained,label="Drained = $(round(corD,digits=2))")
    plot!(DVVimpute[:,:DATE],ypredFC,label="Fully-coupled = $(round(corFC,digits=2))")
    plot!(DVVimpute[:,:DATE],ypredSSW,label="SSW = $(round(corSSW,digits=2))")
    plot!(DVVimpute[:,:DATE],ypredCDM,label="CDMk = $(round(corCDM,digits=2))")
    figpath = joinpath(FIGDIR,netsta*".png")
    savefig(figpath)

end
# write to arrow 
Arrow.write(joinpath(@__DIR__,"../data/hydro-model-$days-day.arrow"),fitdf)
