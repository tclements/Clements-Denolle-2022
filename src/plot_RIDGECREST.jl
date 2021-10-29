using Arrow 
using CSV
using DataFrames
using Dates 
using Plots 

# read dv/v for relevant stations
CLC = Arrow.Table("/media/FOUR/data/DVV-10-DAY-COMP/2.0-4.0/CI.CLC.arrow") |> DataFrame
JRC2 = Arrow.Table("/media/FOUR/data/DVV-10-DAY-COMP/2.0-4.0/CI.JRC2.arrow") |> DataFrame
LRL = Arrow.Table("/media/FOUR/data/DVV-10-DAY-COMP/2.0-4.0/CI.LRL.arrow") |> DataFrame
WCS2 = Arrow.Table("/media/FOUR/data/DVV-10-DAY-COMP/2.0-4.0/CI.WCS2.arrow") |> DataFrame
WRC2 = Arrow.Table("/media/FOUR/data/DVV-10-DAY-COMP/2.0-4.0/CI.WRC2.arrow") |> DataFrame
mindate = Date(2017,1,1)
maxdate = Date(2021,2,18)
cmap = cgrad(:inferno, (0,1),rev=true)

l = @layout [a; b; c]
p1 = scatter(
    CLC[!,:DATE],
    CLC[!,:DVV],  
    alpha = CLC[!,:CC] ./ 5, 
    seriescolor=cmap,
    marker_z=CLC[!,:CC],
    label="CLC",
    xlims = (mindate,maxdate),
    colorbar=false,
    xtick=false,
    ylabel="dv/v [%]",
    legend=:outerright,
)
vline!(
    p1,
    [Date(2019,7,4)],
    linewidth=2,
    c=:black,
    linestyle=:dash,
    label="",
)
p2 = scatter(
    WCS2[!,:DATE],
    WCS2[!,:DVV], 
    alpha = WCS2[!,:CC] ./ 5, 
    seriescolor=cmap,
    marker_z=WCS2[!,:CC],
    label="WCS2",
    xlims = (mindate,maxdate),
    colorbar=false,
    xtick=false,
    ylabel="dv/v [%]",
    legend=:outerright,
)
vline!(
    p2,
    [Date(2019,7,4)],
    linewidth=2,
    c=:black,
    linestyle=:dash,
    label="",
)
p3 = scatter(
    WRC2[!,:DATE],
    WRC2[!,:DVV], 
    alpha = WRC2[!,:CC] ./ 5, 
    seriescolor=cmap,
    marker_z=WRC2[!,:CC],
    label="WRC2",
    xlims = (mindate,maxdate),
    colorbar=false,
    ylabel="dv/v [%]",
    legend=:outerright,
)
vline!(
    p3,
    [Date(2019,7,4)],
    linewidth=2,
    c=:black,
    linestyle=:dash,
    label="",
)
plot(
    p1, 
    p2, 
    p3, 
    layout = l, 
    xlims = (mindate,maxdate),
    sharex=true,
    dpi=500,
)
savefig("/media/FOUR/data/FINAL-FIGURES/RIDGECREST-DVV.svg")
savefig("/media/FOUR/data/FINAL-FIGURES/RIDGECREST-DVV.png")
# to-do 
# 1 move legends 
# 2 sharex 
# 3 add line at earthquake 
# 4 show decorrelation better 