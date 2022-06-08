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
using StatsPlots
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


# constants for models 
r = 500.
δt = 86400.
α = (1 + 0.27) / 3 / (1 - 0.27)
ϕ = 0.15   # porosity 


# define functions for fitting 
function modelCDMk(p)
    cdmk = CDM(precip[tminind-ceil(Int,p[3]):tmaxind],ceil(Int,p[3]))[ceil(Int,p[3])+1:end]
    cdmk .-= mean(cdmk)
    cdmk ./= std(cdmk)
    ypred = p[1] .+ p[2] .* cdmk .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum((ypred .- DVVimpute[:,:DVV]) .^ 2)
    # p1 is the mean, p2 is scalar in front of CDMk, p4 is scalar in front of T, p5 is phase shift in temp, p3 is the lag time in CDMk
end

# define functions for fitting 
function modelCDMkL1(p)
    cdmk = CDM(precip[tminind-ceil(Int,p[3]):tmaxind],ceil(Int,p[3]))[ceil(Int,p[3])+1:end]
    cdmk .-= mean(cdmk)
    cdmk ./= std(cdmk)
    ypred = p[1] .+ p[2] .* cdmk .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum(abs.(ypred .- DVVimpute[:,:DVV]) )
end

function modelSSW(p)
    ssw = SSW06(precip[tminind:tmaxind],ϕ,p[3])
    ssw .-= mean(ssw)
    ssw ./= std(ssw)
    ypred = p[1] .+ p[2] .* ssw .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum((ypred .- DVVimpute[:,:DVV]) .^ 2)
    # p1 is mean, p2 is scalar in front of hydro, p3 is 1/decay time, p4 - p5 are temp.q
end

function modelSSWL1(p)
    ssw = SSW06(precip[tminind:tmaxind],ϕ,p[3])
    ssw .-= mean(ssw)
    ssw ./= std(ssw)
    ypred = p[1] .+ p[2] .* ssw .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum(abs.(ypred .- DVVimpute[:,:DVV]) )
end

function modelElastic(p)
    E = elastic(precip[tminind:tmaxind],p[3],α,r,δt)
    E .-= mean(E)
    E ./= std(E)
    ypred = p[1] .+ p[2] .* E .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum((ypred .- DVVimpute[:,:DVV]) .^ 2)
end

function modelElasticL1(p)
    E = elastic(precip[tminind:tmaxind],p[3],α,r,δt)
    E .-= mean(E)
    E ./= std(E)
    ypred = p[1] .+ p[2] .* E .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum(abs.((ypred .- DVVimpute[:,:DVV]) ))
end

function modelFC(p)
    FC = fullycoupled(precip[tminind:tmaxind],p[3],α,r,δt)
    FC .-= mean(FC)
    FC ./= std(FC)
    ypred = p[1] .+ p[2] .* FC .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum((ypred .- DVVimpute[:,:DVV]) .^ 2)
end

function modelFCL1(p)
    FC = fullycoupled(precip[tminind:tmaxind],p[3],α,r,δt)
    FC .-= mean(FC)
    FC ./= std(FC)
    ypred = p[1] .+ p[2] .* FC .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum(abs.(ypred .- DVVimpute[:,:DVV]) )
end


function modelDrained(p)
    D = drained(precip[tminind:tmaxind],p[3],r,δt)
    D .-= mean(D)
    D ./= std(D)
    ypred = p[1] .+ p[2] .* D .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum((ypred .- DVVimpute[:,:DVV]) .^ 2)
end    

function modelDrainedL1(p)
    D = drained(precip[tminind:tmaxind],p[3],r,δt)
    D .-= mean(D)
    D ./= std(D)
    ypred = p[1] .+ p[2] .* D .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum(abs.(ypred .- DVVimpute[:,:DVV]) )
end    

# df to hold parameters 
fitdf = DataFrame()

for kk in 328:length(files)
    print(kk)
    try
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

        # skip  if nan or constant vallues (out of bounds for PRISM)
        if any(isnan.(temp)) || std(temp)==0 || any(isnan.(precip)) || std(precip)==0
            continue
        end

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


        # Now we will solve L2 and L2 norms for all models and save the parameters into a dataframe


        # solve for elastic model L2 norm
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
        corE = cor(ypredElastic,DVVimpute[:,:DVV]) # pearson correlation
        r2E = corE ^ 2 # explained variance

        # solve for elastic model L1 norm
        resElastic = optimize(
            modelElasticL1, 
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
        ypredElastic1 = pE[1] .+ pE[2] .* E .+ pE[4] .* smoothtemp[tempind .- round(Int,pE[5])]
        corEL1 = cor(ypredElastic1,DVVimpute[:,:DVV]) # pearson correlation
        r2EL1 = corEL1 ^ 2# explained variance

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
        corFC = cor(ypredFC,DVVimpute[:,:DVV]) # pearson correlation
        r2FC = corFC ^ 2


        # solve for fully-coupled model 
        resFC = optimize(
            modelFCL1, 
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
        ypredFC1 = pFC[1] .+ pFC[2] .* FC .+ pFC[4] .* smoothtemp[tempind .- round(Int,pFC[5])]
        corFCL1 = cor(ypredFC1,DVVimpute[:,:DVV]) # pearson correlation
        r2FCL1 = corFCL1 ^ 2

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
        corD = cor(ypredDrained,DVVimpute[:,:DVV]) # pearson correlation
        r2D = corD ^ 2

        # solve for drained model 
        resDrained = optimize(
            modelDrainedL1, 
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
        ypredDrained1 = pD[1] .+ pD[2] .* D .+ pD[4] .* smoothtemp[tempind .- round(Int,pD[5])]
        corDL1 = cor(ypredDrained1,DVVimpute[:,:DVV]) # pearson correlation
        r2DL1 = corDL1 ^ 2

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
        corSSW = cor(ypredSSW,DVVimpute[:,:DVV]) # pearson correlation
        r2SSW = corSSW ^ 2


        # solve for baseflow model 
        resSSW = optimize(
            modelSSWL1, 
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
        ypredSSW1 = pS[1] .+ pS[2] .* ssw .+ pS[4] .* smoothtemp[tempind .- round(Int,pS[5])]
        corSSWL1 = cor(ypredSSW1,DVVimpute[:,:DVV]) # pearson correlation
        r2SSWL1 = corSSWL1 ^ 2


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
        corCDM = cor(ypredCDM,DVVimpute[:,:DVV]) # pearson correlation
        r2CDM = corCDM ^ 2

        # solve for CDMk model 
        resCDM = optimize(
            modelCDMkL1, 
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
        ypredCDM1 = pC[1] .+ pC[2] .* cdmk .+ pC[4] .* smoothtemp[tempind .- round(Int,pC[5])]
        corCDML1 = cor(ypredCDM1,DVVimpute[:,:DVV]) # pearson correlation
        r2CDML1 = corCDML1 ^ 2

        # fitted dvv , save the time series
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
                :CDM1=>ypredCDM1,
                :ELASTIC1=>ypredElastic1,
                :DRAINED1=>ypredDrained1,
                :FC1=>ypredFC1,
                :SSW1=>ypredSSW1,
            )
        )
        fitpath = joinpath(FITDIR,netsta*"_L2_L1.arrow")
        Arrow.write(fitpath,stationdf)
        
        # write to disk 
        append!(
            fitdf,
            DataFrame(
                Dict(
                    :NETSTA=>netsta,
                    :CDM1=>pC[1], 
                    # p1 is the mean, p2 is scalar in front of CDMk, p3 is scalar in front of T, p5 is phase shift in temp.
                    :CDM2=>pC[2], # scalar in front of CDMk
                    :CDM3=>pC[3], # time lag in CDMk
                    :CDM4=>pC[4], # amp in front of temp
                    :CDM5=>pC[5], # phase shift in Temp
                    :D1=>pD[1],# scalar mean
                    :D2=>pD[2], # amplitude for hydro term
                    :D3=>pD[3], # diffusivity
                    :D4=>pD[4],# amp in front of temp
                    :D5=>pD[5], # phase shift in T
                    :E1=>pE[1],# scalar mean
                    :E2=>pE[2], # amplitude for hydro term
                    :E3=>pE[3], # diffusivity
                    :E4=>pE[4],# amp in front of temp
                    :E5=>pE[5], # phase shift in T
                    :FC1=>pFC[1], # scalar mean
                    :FC2=>pFC[2], # amplitude for hydro term
                    :FC3=>pFC[3], # diffusivity
                    :FC4=>pFC[4], # amplitude for temp
                    :FC5=>pFC[5], # time shift for temp
                    :SSW1=>pS[1], # mean amp
                    :SSW2=>pS[2], # amp of hydro
                    :SSW3=>pS[3], # time scale of decay
                    :SSW4=>pS[4], # amp of temp
                    :SSW5=>pS[5], # phase shift in T
                    :r2EL1=>r2EL1, # explained variance for elastic model, L1 norm
                    :r2DL1=>r2DL1, # explained variance for drained model, L1 norm
                    :r2FCL1=>r2FCL1, # explained variance for fully coupled model, L1 norm
                    :r2CDML1=>r2CDML1, # explained variance for CDM model, L1 norm
                    :r2SSWL1=>r2SSWL1, # explained variance for baseflow model, L1 norm
                    :r2E=>r2E,# explained variance for elastic model, L2 norm
                    :r2D=>r2D,# explained variance for drained model, L2 norm
                    :r2FC=>r2FC,# explained variance for fully coupled model, L2 norm
                    :r2CDM=>r2CDM,# explained variance for CDM model, L2 norm
                    :r2SSW=>r2SSW,# explained variance for baseflow model, L2 norm
                    :LAT=>stalat,
                    :LON=>stalon,
                    :ELEV=>staele,
                )
            )
        )

        if isnan(r2E)
            break
        end
        println("$(netsta)")
        println("--------------------")
        println("Elastic = $(round(r2E,digits=2)),$(round(r2EL1,digits=2)) ")
        println("Drained = $(round(r2D,digits=2)), $(round(r2DL1,digits=2)) ")
        println("Fully-Coupled = $(round(r2FC,digits=2)),$(round(r2FCL1,digits=2)) ")
        println("CDM = $(round(r2CDM,digits=2)),$(round(r2CDML1,digits=2)) ")
        println("SSW = $(round(r2SSW,digits=2)), $(round(r2SSWL1,digits=2)) \n\n")

        # plot each model 
        plot(DVVimpute[:,:DATE],
            DVVimpute[:,:DVV],
            legend=:outerbottom,
            label="dv/v",
            ylabel="dv/v [%]",
            linewidth=4,
            legendfontsize=4,
            figsize=(800,400),
            dpi=250,
        )
        plot!(DVVimpute[:,:DATE],ypredElastic,label="Elastic = $(round(r2E,digits=2))",color=theme_palette(:auto).colors.colors[2],linewidth=2)
        plot!(DVVimpute[:,:DATE],ypredDrained,label="Drained = $(round(r2D,digits=2))",color=theme_palette(:auto).colors.colors[3],linewidth=2)
        plot!(DVVimpute[:,:DATE],ypredFC,label="Fully-coupled = $(round(r2FC,digits=2))",color=theme_palette(:auto).colors.colors[4],linewidth=2)
        plot!(DVVimpute[:,:DATE],ypredSSW,label="SSW = $(round(r2SSW,digits=2))",color=theme_palette(:auto).colors.colors[5],linewidth=2)
        plot!(DVVimpute[:,:DATE],ypredCDM,label="CDMk = $(round(r2CDM,digits=2))",color=theme_palette(:auto).colors.colors[6],linewidth=2)
        # figpath = joinpath(FIGDIR,netsta*"_L2.png")
        # savefig(figpath)
        
        # # plot each model 
        # plot(DVVimpute[:,:DATE],
        #     DVVimpute[:,:DVV],
        #     legend=:outerbottom,
        #     label="dv/v",
        #     ylabel="dv/v [%]",
        #     legendfontsize=4,
        #     figsize=(800,400),
        #     dpi=250,
        # )
        plot!(DVVimpute[:,:DATE],ypredElastic1,label="L1 Elastic = $(round(r2EL1,digits=2))",color=theme_palette(:auto).colors.colors[2])
        plot!(DVVimpute[:,:DATE],ypredDrained1,label="L1 Drained = $(round(r2DL1,digits=2))",color=theme_palette(:auto).colors.colors[3])
        plot!(DVVimpute[:,:DATE],ypredFC1,label="L1 Fully-coupled = $(round(r2FCL1,digits=2))",color=theme_palette(:auto).colors.colors[4])
        plot!(DVVimpute[:,:DATE],ypredSSW1,label="L1 SSW = $(round(r2SSWL1,digits=2))",color=theme_palette(:auto).colors.colors[5])
        plot!(DVVimpute[:,:DATE],ypredCDM1,label="L1 CDMk = $(round(r2CDML1,digits=2))",color=theme_palette(:auto).colors.colors[6])
        figpath = joinpath(FIGDIR,netsta*"_L1_L2.png")
        savefig(figpath)
    catch e
        println("something was wrong")
    end


end
# write to arrow 
Arrow.write(joinpath(@__DIR__,"../data/hydro-model-$days-day_L12.arrow"),fitdf)

# write out overall goodness of fit 
println("L2 norm fit")
println("median R2  Drained model "*"$(mean(fitdf[:,:r2D]))")
println("median R2  Elastic model "*"$(mean(fitdf[:,:r2E]))")
println("median R2  Fully-Coupled model "*"$(mean(fitdf[:,:r2FC]))")
println("median R2  Empirical model "*"$(mean(fitdf[:,:r2CDM]))")
println("median R2  baseflow model "*"$(mean(fitdf[:,:r2SSW]))")

# write out overall goodness of fit 
println("  ")
println("L1 norm fit")
println("median R2  Drained model "*"$(mean(fitdf[:,:r2DL1]))")
println("median R2  Elastic model "*"$(mean(fitdf[:,:r2EL1]))")
println("median R2  Fully-Coupled model "*"$(mean(fitdf[:,:r2FCL1]))")
println("median R2  Empirical model "*"$(mean(fitdf[:,:r2CDML1]))")
println("median R2  baseflow model "*"$(mean(fitdf[:,:r2SSWL1]))")


# is L1 better than L2?
println("is L1 better than L2?")
println("median R2  Drained model "*"$(mean(fitdf[:,:r2DL1].-fitdf[:,:r2D]))")
println("median R2  Elastic model "*"$(mean(fitdf[:,:r2EL1]-fitdf[:,:r2E]))")
println("median R2  Fully-Coupled model "*"$(mean(fitdf[:,:r2FCL1]-fitdf[:,:r2FC]))")
println("median R2  Empirical model "*"$(mean(fitdf[:,:r2CDML1]-fitdf[:,:r2CDM]))")
println("median R2  baseflow model "*"$(mean(fitdf[:,:r2SSWL1]-fitdf[:,:r2SSW]))")

#We conclude that the drained model is favored.
# is it better at most stations?
println("number of stations where Drain beats elastic  "*"$(sum((fitdf[:,:r2D] -fitdf[:,:r2E] ).>0 )/sum((fitdf[:,:r2D]  ).>0 ))")
println("number of stations where Drain beats Fully coupled  "*"$(sum((fitdf[:,:r2D] -fitdf[:,:r2FC] ).>0 )/sum((fitdf[:,:r2D]  ).>0 ))")
println("number of stations where Drain beats CDM  "*"$(sum((fitdf[:,:r2D] -fitdf[:,:r2CDM] ).>0 )/sum((fitdf[:,:r2D]  ).>0 ))")
println("number of stations where Drain beats SSW  "*"$(sum((fitdf[:,:r2D] -fitdf[:,:r2SSW] ).>0 )/sum((fitdf[:,:r2D]  ).>0 ))")

#some Plots


Plots.scatter(fitdf[:,:LON],fitdf[:,:LAT],zcolor=(fitdf[:,:r2DL1]),
title="R2 for drained model",color=:bilbao,
colorbar_title="(R2)",legend=false,colorbar=true)
savefig("../data/FINAL-FIGURES/scatter_r2_drained.png")


# Plots.histogram(fitdf[:,:r2D]-fitdf[:,:r2E])
# @df fitdf plot(:r2E)
