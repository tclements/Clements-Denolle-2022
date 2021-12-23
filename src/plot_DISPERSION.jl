using CSV 
using DataFrames
using Glob  
using Plots 

# read Vp / Vs / Rho data 
df = DataFrame(CSV.File(joinpath(@__DIR__,"../data/DISPERSION/disp"),delim=' ',ignorerepeated=true,select=[3,7,8,9], header=0))
model = df[:,[:Column3, :Column7, :Column8, :Column9]] ./ 1000 
rename!(model,["DEPTH","VP", "VS", "RHO"])

# read elastic rayleigh wave data 
kernels = glob("elastic*",joinpath(@__DIR__,"../data/DISPERSION/"))
N = length(kernels)

# plot Vp Vs depth 
p1 = plot(
    model[!,:VP],
    model[!,:DEPTH],
    label="\$V_p\$",
    color=:black,
    linewidth=2,
    alpha=0.75,
    yaxis=:flip,
    xlabel="km/s or g / cm^3",
    ylabel="Depth [km]",
    ylims=(0,4),
    legend=:top
)
plot!(
    p1, 
    model[!,:VS],
    model[!,:DEPTH],
    label="\$V_s\$",
    color=:red,
    linewidth=2,
    alpha=0.75,
    yaxis=:flip,
    ylims=(0,4),
)
plot!(
    p1, 
    model[!,:RHO],
    model[!,:DEPTH],
    label="\$\\rho\$",
    color=:blue,
    yaxis=:flip,
    linewidth=2,
    alpha=0.75,
    ylims=(0,4),
)

# plot dc/db -> change in phase velocity wrt shear wave velocity at each layer 
p2 = plot()
for ii in N:-1:1 
    kernel = kernels[ii]
    period = parse(Float64,replace(split(kernel, "_")[end],".txt"=>""))
    freq = 1 / period 
    df = DataFrame(CSV.File(kernel, delim=' ', ignorerepeated=true, header=2))
    df[!,:DEPTH] = cumsum(df[!,:THICK]) .- df[1,:THICK]
    ind = findall(df[!,:DEPTH] .<= 4)
    df = df[ind,:]
    maxdepth = round(df[argmax(df[!,"dc/db"]),:DEPTH],digits=2)
    plot!(
        p2,
        df[:,"dc/db"],
        df[:,:DEPTH],
        linewidth=2,
        alpha=0.75,
        yaxis=:flip, 
        ylim=(0,4),
        xlabel="\$\\frac{\\partial c}{\\partial \\beta}\$",
        ylabel="Depth [km]",
        label="$(freq) Hz $(maxdepth) km",
        legend=:right,
    )
end
l = @layout [a b]
plot(p1,p2,layout=l,dpi=500)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CILJR-DISPERSION.png"))
