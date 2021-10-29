using Arrow 
using DataFrames 
using Glob
using Plots 
using Statistics

import GR
GR.inline("png")

# get all dv/v files by station 
OUTDIR = "/media/FOUR/data/COMP-NC-MWCS/2.0-4.0/"
if !isdir(OUTDIR)
    mkpath(OUTDIR)
end
FIGDIR = "/media/FOUR/data/FIGURES-NC-COMP-MWCS/2.0-4.0/"
if !isdir(FIGDIR)
    mkpath(FIGDIR)
end
files = glob("*","/media/FOUR/data/NC-MWCS-DVV/2.0-4.0/")
NCdvv = DataFrame(FILE=files)
NCdvv[:,:STA] = [join(split(basename(f),".")[1:2],".") for f in files]
NCdvv[:,:COMP] = [join(replace.(replace.(split(basename(f),".")[4:4:8],"H"=>""),"B"=>""),"") for f in files]
stas = groupby(NCdvv,:STA)
for sta in stas 
    sta = DataFrame(sta)
    sort!(sta,:COMP)
    EN = Arrow.Table(sta[1,:FILE]) |> DataFrame 
    EZ = Arrow.Table(sta[2,:FILE]) |> DataFrame 
    NZ = Arrow.Table(sta[3,:FILE]) |> DataFrame 

    # get common dates 
    # get intersect of all dates 
    dates = intersect(EN[!,:DATE],EZ[!,:DATE],NZ[!,:DATE])
    EN = EN[[d in dates for d in EN[!,:DATE]],:]
    EZ = EZ[[d in dates for d in EZ[!,:DATE]],:]
    NZ = NZ[[d in dates for d in NZ[!,:DATE]],:]

    # get CC and DV/V
    CCEN = clamp.(1 .- EN[:,:ERRBOTH] * 1e4,0,1)  
    CCEZ = clamp.(1 .- EZ[:,:ERRBOTH] * 1e4,0,1)  
    CCNZ = clamp.(1 .- NZ[:,:ERRBOTH] * 1e4,0,1) 
    CC = CCEN .^ 2 .+ CCEZ .^ 2 .+ CCNZ .^ 2
    DVV = CCEN .^ 2 .* EN[:,:DVVBOTH] .+ CCEZ .^ 2 .* EZ[:,:DVVBOTH] .+ CCNZ .^ 2 .* NZ[:,:DVVBOTH]
    DVV ./= CC 
    CC = (CCEN .^ 3 .+ CCEZ .^ 3 .+ CCNZ .^ 3) ./ CC 

    ALL = DataFrame(DATE=EN[:,:DATE],DVV=DVV,CC=CC)

    # form matrix of dv/v 
    A = zeros(size(NZ,1),3)
    A[:,1] = EN[!,:DVVBOTH]
    A[:,2] = EZ[!,:DVVBOTH]
    A[:,3] = NZ[!,:DVVBOTH]

    # get correlation 
    DVVcor = cor(A)

    # annotations for each comps 
    comps = ["EN","EZ","NZ"]
    anns = Array{Tuple{Real,Real,Any}}(undef,9)
    for ii = 1:3
        for jj = 1:3
            anns[(ii-1)*3+jj] = (ii-0.75,jj-0.25,text("$(round(DVVcor[ii,jj],digits=2))",8,:white,:top,:left))
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
        c=:linear_ternary_red_0_50_c52_n256,
    )
    cmap = cgrad(:diverging_bkr_55_10_c35_n256,scale=(-1,1))
    p2 = heatmap(comps,comps,DVVcor,seriescolor=cmap,clims=(-1,1))
    annotate!(anns)
    title!("$(sta[1,:STA])")
    plot(p1, p2, layout=l)
    figpath = joinpath(FIGDIR,"$(sta[1,:STA]).png")
    savefig(figpath)

    # write MWCS comp 
    savepath = joinpath(OUTDIR,"$(sta[1,:STA]).arrow")
    Arrow.write(savepath,ALL)
end
