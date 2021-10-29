using Arrow 
using DataFrames 
using Dates 
using Glob
using Plots
using Statistics

function compheat(arrows,day,FIGDIR)
    station = join(split(basename(arrows[1]),'.')[1:2],".")
    println("Printing station $station $(now())")
    # load all components
    EN = Arrow.Table(arrows[1]) |> Arrow.columntable |> DataFrame
    EZ = Arrow.Table(arrows[2]) |> Arrow.columntable |> DataFrame
    NZ = Arrow.Table(arrows[3]) |> Arrow.columntable |> DataFrame

    # get unqiue rows 
    unique!(EN)
    unique!(EZ)
    unique!(NZ)

    # get comp dv/v 
    ALL = Arrow.Table("/media/FOUR/data/DVV-$day-DAY-COMP/2.0-4.0/$(station).arrow") |> DataFrame

    dates = intersect(EN[!,:DATE],EZ[!,:DATE],NZ[!,:DATE])
    if length(dates) == 0 
        return nothing 
    end
    EN = EN[[d in dates for d in EN[!,:DATE]],:]
    EZ = EZ[[d in dates for d in EZ[!,:DATE]],:]
    NZ = NZ[[d in dates for d in NZ[!,:DATE]],:]

    # form matrix of dv/v 
    A = zeros(size(NZ,1),6)
    A[:,1] = EN[!,:DVVPOS]
    A[:,2] = EN[!,:DVVNEG]
    A[:,3] = EZ[!,:DVVPOS]
    A[:,4] = EZ[!,:DVVNEG]
    A[:,5] = NZ[!,:DVVPOS]
    A[:,6] = NZ[!,:DVVNEG]

    # get correlation 
    DVVcor = cor(A)

    # annotations for each comps 
    comps = ["EN POS", "EN NEG", "EZ POS", "EZ NEG", "NZ POS", "NZ NEG"]
    anns = Array{Tuple{Real,Real,Any}}(undef,36)
    for ii = 1:6
        for jj = 1:6
            anns[(ii-1)*6+jj] = (ii-0.75,jj-0.25,text("$(round(DVVcor[ii,jj],digits=2))",8,:white,:top,:left))
        end
    end

    # plot comps + dv/v 
    l = @layout [a; b]
    p1 = scatter(
        ALL[!,:DATE],
        ALL[!,:DVV],
        marker_z=ALL[!,:CC],
        alpha=ALL[!,:CC] ./ 10,
        label = "",
        ylabel="dv/v [%]",
        dpi=500,
        size=(800,400)
    )
    cmap = cgrad(:diverging_bkr_55_10_c35_n256,scale=(-1,1))
    p2 = heatmap(comps,comps,DVVcor,seriescolor=cmap,clims=(-1,1))
    annotate!(anns)
    title!("$station")
    plot(p1, p2, layout=l)
    savepath = joinpath(FIGDIR,"$station.png")
    savefig(savepath)
    return nothing 
end

day = 10 
CORRS = glob("*.arrow*","/media/FOUR/data/DVV-$day-DAY/2.0-4.0/")
FIGDIR = "/media/FOUR/data/FIGURES-DVV-$day-DAY-COMP/"
if !isdir(FIGDIR)
    mkpath(FIGDIR)
end

for ii in 1:3:length(CORRS)-1
    compheat(CORRS[ii:ii+2],day,FIGDIR)
end

# scatter(EN[!,:DATE],EN[!,:DVVPOS],marker_z=EN[!,:CCPOS],alpha=EN[!,:CCPOS] ./ 2, label = "EN POS",ylims=(-3,3))
# scatter(EN[!,:DATE],EN[!,:DVVNEG],marker_z=EN[!,:CCNEG],alpha=EN[!,:CCNEG] ./ 2, label = "EN NEG",ylims=(-3,3))
# scatter(EZ[!,:DATE],EZ[!,:DVVPOS],marker_z=EZ[!,:CCPOS],alpha=EZ[!,:CCPOS] ./ 2, label = "EZ POS",ylims=(-3,3))
# scatter(EZ[!,:DATE],EZ[!,:DVVNEG],marker_z=EZ[!,:CCNEG],alpha=EZ[!,:CCNEG] ./ 2, label = "EZ NEG",ylims=(-3,3))
# scatter(NZ[!,:DATE],NZ[!,:DVVPOS],marker_z=NZ[!,:CCPOS],alpha=NZ[!,:CCPOS] ./ 2, label = "NZ POS",ylims=(-3,3))
# scatter(NZ[!,:DATE],NZ[!,:DVVNEG],marker_z=NZ[!,:CCNEG],alpha=NZ[!,:CCNEG] ./ 2, label = "NZ NEG",ylims=(-3,3))