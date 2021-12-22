using CSV 
using DataFrames
using Dates 
using Plots 
using Plots.Measures

# load wells 
MW1 = CSV.File(joinpath(@__DIR__,"../data/TRC-MW1.csv"),dateformat="Y/m/d") |> DataFrame
MW1[:,:GWL] ./= 3.28

MW16 = CSV.File(joinpath(@__DIR__,"../data/TRC-MW16D.csv"),dateformat="Y/m/d") |> DataFrame
MW16[:,:GWL] ./= 3.28

LWCD = CSV.File(joinpath(@__DIR__,"../data/LWCD-lebec.csv"),dateformat="Y/m/d") |> DataFrame
LWCD[:,:GWL] ./= 3.28

PW90 = CSV.File(joinpath(@__DIR__,"../data/TRC-PW90.csv"),dateformat="Y/m/d") |> DataFrame
PW90[:,:GWL] ./= 3.28

# load Castac Lake groundwater 
LAKE = CSV.File(joinpath(@__DIR__,"../data/CASTAC-LAKE-ELEVATION.csv"),dateformat="Y/m/d") |> DataFrame
LAKE[:,:ELEVATION] ./= 3.28 

# plot monitoring wells 
plot(
    MW1[:,:DATE],
    MW1[:,:GWL],
    label="Monitoring Well 1",
    legend=:topright,
    size=(800,400),
    color=:lightskyblue,
    linewidth=2.5,
    alpha=0.85,
    linestyle=:dash,
    ylabel="Groundwater Elevation [m]",
    left_margin=1cm,
    dpi=250,
)
plot!(
    MW16[:,:DATE],
    MW16[:,:GWL],
    label="Monitoring Well 16D",
    color=:lightskyblue,
    linewidth=2.5,
    alpha=0.85,
)

# plot pumping wells 
plot!(
    PW90[:,:DATE],
    PW90[:,:GWL],
    label="Pumping Well 90",
    color=:indianred,
    linewidth=2.5,
    alpha=0.85,
)
plot!(
    LWCD[:,:DATE],
    LWCD[:,:GWL],
    label="Lebec Pumping Well",
    color=:indianred,
    linewidth=2.5,
    alpha=0.85,
    linestyle=:dash,
)

# plot lake elevation 
plot!(
    LAKE[:,:DATE],
    LAKE[:,:ELEVATION], 
    label="Castac Lake Elevation",
    color=:black,
    linewidth=2.5,
    alpha=0.85,
)
savefig(joinpath(@__DIR__,"../data/FINAL-FIGURES/CASTAC-GWL"))