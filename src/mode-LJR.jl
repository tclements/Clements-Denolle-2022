using Arrow 
using DataFrames 
using Dates 
using DSP 
using Impute 
using NetCDF
using Optim 
using Plots
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

# load dv/v and impute values 
LJR = Arrow.Table("/media/FOUR/data/DVV-10-DAY-COMP/2.0-4.0/CI.LJR.arrow") |> DataFrame
ind = findall(LJR[:,:CC] .>= 0.875)
LJR = LJR[ind,:]
dt = LJR[1,:DATE]:Day(1):LJR[end,:DATE]
LJRdt = DataFrame(DATE = dt)
LJRmiss = outerjoin(LJRdt,LJR,on=:DATE)
sort!(LJRmiss,:DATE)
LJRimpute = Impute.interp(LJRmiss)

# get temperature data 
filename = "/media/FOUR/data/TEMP-LJR.nc"
temp = ncread(filename,"tmean")[1,1,:]
t = ncread(filename,"t")
ttemp = Date.(Dates.unix2datetime.(t))
tempind = findall(in(LJRimpute[:,:DATE]),ttemp)
smoothtemp = smooth_withfiltfilt(temp,window_len=5)
# temp = temp[tempind]
smoothtemp .-= mean(smoothtemp)
smoothtemp ./= std(smoothtemp)

# get precip data 
filename = "/media/FOUR/data/PRECIP-LJR.nc"
precip = ncread(filename,"ppt")[1,1,:]
t = ncread(filename,"t")
tppt = Date.(Dates.unix2datetime.(t))

# constants for models 
r = 500.
δt = 86400.
α = (1 + 0.27) / 3 / (1 - 0.27)
ϕ = 0.15   # porosity 

# solve for elastic model 
function modelElastic(p)
    ypred = p[1] .+ p[2] .* elastic(precip[tminind:tmaxind],p[3],α,r,δt) .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum((ypred .- LJRimpute[:,:DVV]) .^ 2)
end
resElastic = optimize(modelElastic, [-Inf,-Inf,0.0,-Inf,0],[Inf,Inf,Inf,Inf,150],[0.0, -0.001,0.001,1.0,1])
pE = Optim.minimizer(resElastic)
ypredElastic = pE[1] .+ pE[2] .* elastic(precip[tminind:tmaxind],pE[3],α,r,δt) .+ pE[4] .* smoothtemp[tempind .- round(Int,pE[5])]

# solve for fully-coupled model 
function modelFC(p)
    ypred = p[1] .+ p[2] .* fullycoupled(precip[tminind:tmaxind],p[3],α,r,δt) .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum((ypred .- LJRimpute[:,:DVV]) .^ 2)
end
resFC = optimize(modelFC, [-Inf,-Inf,0.0,-Inf,0],[Inf,Inf,Inf,Inf,150],[0.0, -0.001,0.001,1.0,1])
pFC = Optim.minimizer(resFC)
ypredFC = pFC[1] .+ pFC[2] .* fullycoupled(precip[tminind:tmaxind],pFC[3],α,r,δt) .+ pFC[4] .* smoothtemp[tempind .- round(Int,pFC[5])]

# solve for drained model 
function modelDrained(p)
    ypred = p[1] .+ p[2] .* drained(precip[tminind:tmaxind],p[3],r,δt) .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum((ypred .- LJRimpute[:,:DVV]) .^ 2)
end
resDrained = optimize(modelDrained, [-Inf,-Inf,0.0,-Inf,0],[Inf,Inf,100.0,Inf,150],[0.0, -0.001,1.0,1.0,1])
pD = Optim.minimizer(resDrained)
ypredDrained = pD[1] .+ pD[2] .* drained(precip[tminind:tmaxind],pD[3],r,δt) .+ pD[4] .* smoothtemp[tempind .- round(Int,pD[5])]

# solve for baseflow model 
function modelSSW(p)
    ypred = p[1] .+ p[2] .* SSW06(precip[tminind:tmaxind],ϕ,p[3]) .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum((ypred .- LJRimpute[:,:DVV]) .^ 2)
end
resSSW = optimize(modelSSW, [-Inf,-Inf,0.,-Inf,0],[Inf,Inf,1.,Inf,150],[0.0, -0.001,0.001,1.0,1])
pS= Optim.minimizer(resSSW)
ypredSSW = pS[1] .+ pS[2] .* SSW06(precip[tminind:tmaxind],ϕ,pS[3]) .+ pS[4] .* smoothtemp[tempind .- round(Int,pS[5])]

# solve for CDMk model 
function modelCDMk(p)
    ypred = p[1] .+ p[2] .* CDM(precip[tminind-ceil(Int,p[3]):tmaxind],ceil(Int,p[3]))[ceil(Int,p[3])+1:end] .+ p[4] .* smoothtemp[tempind .- round(Int,p[5])]
    return  sum((ypred .- LJRimpute[:,:DVV]) .^ 2)
end
resCDM = optimize(modelCDMk, [-Inf,-Inf,365. ,-Inf,0],[Inf,Inf,365.0 * 15,Inf,150],[0.0, -0.001,365.0 * 8.,1.0,1])
pC = Optim.minimizer(resCDM)
ypredCDM = pC[1] .+ pC[2] .* CDM(precip[tminind-ceil(Int,pC[3]):tmaxind],ceil(Int,pC[3]))[ceil(Int,pC[3])+1:end] .+ pC[4] .* smoothtemp[tempind .- round(Int,pC[5])]

# plot each model 
plot(LJRimpute[:,:DATE],LJRimpute[:,:DVV],legend=:topleft,label="dv/v")
plot!(LJRimpute[:,:DATE],ypredSSW,label="SSW")
plot!(LJRimpute[:,:DATE],ypredFC,label="Fully-coupled")
plot!(LJRimpute[:,:DATE],ypredCDM,label="CDMk")
plot!(LJRimpute[:,:DATE],ypredDrained,label="Drained")