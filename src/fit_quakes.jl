using Arrow , CSV, DataFrames, Dates,  DSP , HDF5, JSON
using LaTeXStrings, LinearAlgebra, Optim ,Plots , Plots.Measures
using Statistics , QuadGK



cd(joinpath(@__DIR__,"Dropbox/RESEARCH_GROUP/TIM_MARINE_PROJEcTS/Clements-Denolle-2022/src/"))
# read data from the three stations
WES = Arrow.Table(joinpath(@__DIR__,"../../data/FIT-DVV-SSE/90-DAY/CI.WES.arrow")) |> DataFrame
HEC = Arrow.Table(joinpath(@__DIR__,"../../data/FIT-DVV-SSE/90-DAY/CI.HEC.arrow")) |> DataFrame
JRC2 = Arrow.Table(joinpath(@__DIR__,"../../data/FIT-DVV-SSE/90-DAY/CI.JRC2.arrow")) |> DataFrame

# dates: El Cucapha
QuakeDate0=Date(2010,4,4) # add 90 days due to the smoothing of the dv/v time series
QuakeDate=Date(2010,4,4)+Day(90) # add 90 days due to the smoothing of the dv/v time series
ind_postquake = findall(WES[:,:DATE].>QuakeDate)
tyear_wes=(WES[ind_postquake,:DATE]-QuakeDate)./Day(1)./365.0
indq = findall(WES[:,:DATE].>QuakeDate0)[1]
yb_wes=mean(WES[indq-15:indq+15,:DVV] - WES[indq-15:indq+15,:DRAINDED])
# Hectore Mine
QuakeDate0_Hectore=Date(1999,10,16) # add 90 days due to the smoothing of the dv/v time series
QuakeDate_Hectore=Date(1999,10,16)+Day(90) # add 90 days due to the smoothing of the dv/v time series
ind_postquake_hectore = findall(HEC[:,:DATE].>QuakeDate_Hectore)
tyear_hec=(HEC[ind_postquake_hectore,:DATE]-QuakeDate_Hectore)./Day(1)./365.0
indq = findall(HEC[:,:DATE].>QuakeDate0_Hectore)[1]
yb_hec=mean(HEC[indq-15:indq+15,:DVV] - HEC[indq-15:indq+15,:DRAINDED])
# Ridgecrest
QuakeDate0_ridge=Date(2019,7,5) # add 90 days due to the smoothing of the dv/v time series
QuakeDate_ridge=Date(2019,7,5)+Day(90) # add 90 days due to the smoothing of the dv/v time series
ind_postquake_ridge = findall(JRC2[:,:DATE].>QuakeDate_ridge)
tyear_jrc=(JRC2[ind_postquake_ridge,:DATE]-QuakeDate_ridge)./Day(1)./365.0
indq = findall(JRC2[:,:DATE].>QuakeDate0_ridge)[1]
yb_jrc=mean(JRC2[indq-15:indq+15,:DVV] - JRC2[indq-15:indq+15,:DRAINDED])


# get WES lat, lon from CI station locations 
SCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../../data/CIstations.csv")))
NCdf = DataFrame(CSV.File(joinpath(@__DIR__,"../../data/NCstations.csv")))
CAdf = vcat(NCdf, SCdf)
CAdf = CAdf[CAdf[:,:StartTime] .< QuakeDate,:]
CAdf = CAdf[CAdf[:,:EndTime] .> QuakeDate,:]
CAdf = CAdf[CAdf[:,:Longitude] .>= -116.108,:]


# fitting functions 


function smooth_withfiltfilt(A::AbstractArray; window_len::Int=11, window::Symbol=:rect)
    w = getfield(DSP.Windows, window)(window_len)
    A = DSP.filtfilt(w ./ sum(w), A)
    return A
end

# log model 
function logmodel(p)
    model = p[1] .* log.(1.0 .+ dtfit ./ p[2])
    fit = resid0 .- model 
    return sum(fit .^ 2)
end

# exp model 
function expmodel(p)
    model =  p[1] .* (1.0 .- exp.(-dtfit ./ p[2]))
    fit = resid0 .- model 
    return sum(fit .^ 2)
end

# healing function
function healing_s(t::AbstractArray,tmin1::Float64,tmax1::Float64)
    c=zeros( length(t))
    for x = 1:length(t)
        fexp(tau)= exp(-t[x]/tau)./tau
        c[x],e=QuadGK.quadgk(fexp, tmin1, tmax1;maxevals=10_000) # here tau varies
    end
    return c
end

function modelhealingWES(p::AbstractArray)
    ypred = p[1].*healing_s(tyear_wes,p[2],p[3])#.+p[3]
    crap=sum(abs.(ypred .- ( WES[ind_postquake,:DVV]
    -WES[ind_postquake,:DRAINDED].-yb_wes)).^2 )
    return  crap
end   


function modelhealingHEC(p::AbstractArray)
    ypred = p[1].*healing_s(tyear_hec,p[2],p[3])#.+p[3]
    crap=sum(abs.(ypred .- (HEC[ind_postquake_hectore,:DVV]
    -HEC[ind_postquake_hectore,:DRAINDED].-yb_hec)).^2 )
    return  crap
end    


function modelhealingJRC(p::AbstractArray)
    ypred = p[1].*healing_s(tyear_jrc,p[2],p[3])#.+p[3]
    crap=   sum(abs.(ypred .- ( JRC2[ind_postquake_ridge,:DVV]
    -JRC2[ind_postquake_ridge,:DRAINDED].-yb_jrc)).^2 )
        return crap
end    
# fit log-log model 
# solve for log model 
reslog = optimize(
    logmodel, 
    [-Inf,5e-5],
    [Inf,Inf],
    [0.0,7.0],
    Fminbox(LBFGS()),
    Optim.Options(time_limit=60.0),
)
plog = Optim.minimizer(reslog)
l = plog[1] .* log.(1.0 .+ dtpred ./ plog[2])

# solve for exp model 
resexp = optimize(
    expmodel, 
    [-Inf,5e-5],
    [Inf,Inf],
    [0.0,7.0],
    Fminbox(LBFGS()),
    Optim.Options(time_limit=60.0),
)
pexp = Optim.minimizer(resexp)
e = pexp[1] .* (1.0 .- exp.(-dtpred ./ pexp[2]))

# optimize Snieder's model
# El Cucapah Mayor
reshealing = optimize(
        modelhealingWES, 
        [-3.,0.,0.1], # lower
        [0.,10.,30.], # upper
        [-1.5,0.1,7.], # initial values
        Fminbox(LBFGS()),
        Optim.Options(time_limit=60.0))
p = Optim.minimizer(reshealing)
ypredWES =  p[1].*healing_s(tyear_wes,p[2],p[3])#.+p[3]  
println("Fit for El Cucapah Mayor returns "*string(p[2])*" "*string(p[3])*" years")
# 13.8 years     
plot(WES[ind_postquake,:DATE],WES[ind_postquake,:DVV]
-WES[ind_postquake,:DRAINDED].-yb_wes)
plot!(WES[ind_postquake,:DATE],ypredWES)
tmin_wes=p[2]
tmax_wes=p[3]
# Hectore Mine
reshealing = optimize(
        modelhealingHEC, 
        [-3.,0.,0.1],
        [0.,10.,30.],
        [-1.5,3.,7. ],
        Fminbox(LBFGS()),
        Optim.Options(time_limit=60.0))
p = Optim.minimizer(reshealing)
ypredHEC =  p[1].*healing_s(tyear_hec,p[2],p[3]) #.+p[3]
println("Fit for Hector Mine returns "*string(p[2])*" years")
tmin_hec=p[2]
tmax_hec=p[3]   
plot(HEC[ind_postquake_hectore,:DATE],HEC[ind_postquake_hectore,:DVV]
-HEC[ind_postquake_hectore,:DRAINDED].-yb_hec)
plot!(HEC[ind_postquake_hectore,:DATE],ypredHEC)
plot!(HEC[:,:DATE],HEC[:,:DVV]-HEC[:,:DRAINDED].-yb_hec)

# Ridgecrest
reshealing = optimize(
        modelhealingJRC, 
        [-3.,0.,0.1],
        [0.,10.,30.],
        [-1.5,3.,7.],
        Fminbox(LBFGS()),
        Optim.Options(time_limit=60.0))
p = Optim.minimizer(reshealing)
ypredJRC =  p[1].*healing_s(tyear_jrc,p[2],p[3])#.+p[3]        
println("Fit for Ridgecrest returns "*string(p[2])*" years")
plot(JRC2[ind_postquake_ridge,:DATE],JRC2[ind_postquake_ridge,:DVV]
-JRC2[ind_postquake_ridge,:DRAINDED].-yb_jrc)
plot!(JRC2[ind_postquake_ridge,:DATE],ypredJRC)
plot!(JRC2[:,:DATE],JRC2[:,:DVV]
-JRC2[:,:DRAINDED].-yb_jrc)
tmin_jrc=p[2]
tmax_jrc=p[3]

# plot WES dv/V
Ylims = (-2.0, 1.0)
Xlims = (WES[1,:DATE],preddate[end])
scatter(
    WES[!,:DATE],
    WES[!,:DVV] .- WES[:,:DRAINDED],
    seriescolor=:Blues_9,
    marker_z=WES[!,:CC],
    alpha=WES[!,:CC] ./ 10,
    label="",
    colorbar=false,
    ylabel="          dv/v [%]",
    # yaxis=:flip,
    ylims=Ylims,
    xlims=Xlims,
    # xlims=(WES[1,:DATE],WES[end,:DATE]),
    right_margin = 1.5cm,
    left_margin = 2cm,
    yguidefontcolor=:blue3,
    # ytickfont=font(10,:blue3),
    yforeground_color_axis=:blue3,
    ytick_direction=:out,
    dpi=500,
    size=(800,400),
)
# # plot line for earthquake 
# vline!(
#     [Date(2010,4,4)],
#     label="Mw 7.2 Baja EQ",
#     linewidth=2,
#     c=:black,
#     linestyle=:dash,
#     legend=:topleft,
#     legendfontsize=9,
# )
# plot arrow for earthquake 
xdate = Date(2005) : Day(1) : Date(2010,1,1)
x = (xdate .- xdate[1]) ./ Day(1) ./ length(xdate) .* (π /3) .+ π / 6
Plots.plot!(
    xdate,
    sin.(x) .* 0.62 .- 1.3,
    arrow=true,
    lw=2,
    color=:black,
    label="",
)
Plots.annotate!(
    (
        Date(2005,6,1), 
        -1.3 ,
        Plots.text(
            "Baja CA Earthquake\n M7.2",
            12,
            :bold,
            :black,
        ),
    )
)
Plots.plot!(
    preddate,
    l .+ resid1,
    lw=3,
    c=:purple,
    alpha=0.75,
    ls=:solid,
    # label="$(round(plog[1],digits=2)) * log(t[years] / $(round(plog[2],digits=2)))",
    label="",
    legend=:bottomright,
    ylims=Ylims,
)
Plots.plot!(
    preddate,
    e .+ resid1,
    lw=3,
    c=:orange,
    alpha=0.75,
    ls=:solid,
    # label="$(round(pexp[1],digits=2)) * (1 - exp(-t[years] / $(round(pexp[2],digits=2))))",
    label="",
    legend=:bottomright,
    ylims=Ylims,
)
Plots.annotate!(
    (
        Date(2024,4,1), 
        0.57 ,
        Plots.text(
            "log",
            12,
            rotation = 10,
            :bold,
            :purple,
        ),
    ),
    subplot=1,
)
Plots.annotate!(
    (
        Date(2019,1,1), 
        -1.3 ,
        Plots.text(
            L"dv/v = %$(round(plog[1],digits=2)) * log(t[years] / %$(round(plog[2],digits=2)))",
            16,
            :bold,
            :purple,
        ),
    ),
    subplot=1,
)
Plots.annotate!(
    (
        Date(2019,1,1), 
        -1.7,
        Plots.text(
            L"dv/v = %$(round(pexp[1],digits=2)) * (1 - e^{-t[years] / %$(round(pexp[2],digits=2))})",
            16,
            :bold,
            :orange,
        ),
    ),
    subplot=1,
)

Plots.annotate!(
    (
        Date(2024,4,1), 
        0.1 ,
        Plots.text(
            "exp",
            12,
            :bold,
            :orange,
        ),
    ),
    subplot=1,
)

# plot GPS after earthquake 
plot!(
    twinx(), 
    GPS[!,:DATE], 
    smoothEN, 
    sharex=true, 
    xticks=false,
    grid=:off, 
    ylims=(0,90),
    linewidth=2.5,
    c=:red,
    alpha=0.95,
    label="",
    ylabel="GPS H. Displacement [mm]",
    # yaxis=:flip,
    # xlims=(WES[1,:DATE],WES[end,:DATE]),
    yguidefontcolor=:red,
    ytickfont=font(10,:red),
    yforeground_color_axis=:red,
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CIWES-DVV.png"))




# find groundwater wells nearby 
# 344614118454101.tsv
GWL = CSV.File(joinpath(@__DIR__,"../../data/344614118454101.tsv"),comment="#",skipto=3) |> DataFrame
# GWL = CSV.File(joinpath(@__DIR__,"../../data/324603115480501.tsv"),comment="#",skipto=3) |> DataFrame
dropmissing!(GWL,:lev_va) # lev_va is ft below surface
GWL[!,:lev_va] ./= 3.28 # convert to meters 
GWL = GWL[GWL[!,:lev_dt] .>= Date(2001),:]
GWL = GWL[GWL[:,:lev_va ] .> 14.9,:]
GWL[!,:lev_va] .*= -1


# 2010 quake
# read PGA from USGS https://earthquake.usgs.gov/earthquakes/eventpage/ci9108652/shakemap/intensity
baja = h5open(joinpath(@__DIR__,"../data/Baja-EQ-shake.hdf"),"r") 
arrays = read(baja,"arrays")
PGA = arrays["imts"]["GREATER_OF_TWO_HORIZONTAL"]["PGA"]["mean"]
PGV = arrays["imts"]["GREATER_OF_TWO_HORIZONTAL"]["PGV"]["mean"]
# rotate array for plotting 
PGV = Array(transpose(reverse(PGV,dims=2)))
PGA = Array(transpose(reverse(PGA,dims=2)))
# transform from log(PGV) to PGV 
PGV = exp.(PGV)
PGA = exp.(PGA)

# get PGA grid spacing 
info = read(baja,"dictionaries/info.json") |> JSON.parse
Nlon = parse(Int,info["output"]["map_information"]["grid_points"]["longitude"])
Nlat = parse(Int,info["output"]["map_information"]["grid_points"]["latitude"])
latmax = parse(Float64,info["output"]["map_information"]["max"]["latitude"])
latmin = parse(Float64,info["output"]["map_information"]["min"]["latitude"])
lonmax = parse(Float64,info["output"]["map_information"]["max"]["longitude"])
lonmin = parse(Float64,info["output"]["map_information"]["min"]["longitude"])
PGAlon = range(lonmin,lonmax,length=Nlon)
PGAlat = range(latmin,latmax,length=Nlat)