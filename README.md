# Clements-Denolle-2022

This collection of codes are available to reproduce the work in:

Clements T, Denolle MA. "The Seismic Signature of California's Earthquakes,
Droughts, and Floods", eartharxiv, https://doi.org/10.31223/X5RQ0G.

Please cite the paper if you use the codes.

## Files to reproduce methods and results 

1. src/01-stream-auto-corr.jl: Autocorrelation of seismic noise from California on AWS
2. src/02-single-station-dvv.jl: dv/v of single-station autocorrelations on AWS
3. src/03-dvv-comp.jl: combine EN, EZ, NZ componentsof dv/v into single single time-series
4. src/04-fit-thermo-hydro-models.jl: Fit hydro-elastic and thermo-elastic models to dv/v time series 

## Files to create manuscript figures: 
1. src/05-plot-figure-1.jl: Plot map of California seismic stations 
2. src/06-plot-figure-2-3.jl: Plot single-station autocorrelation and noise spectrum of CI.LJR
3. src/06-plot-figure-2-3.jl: Plot CI.LJR dv/v time-series 
3a. src/07-plot-figure-3-map.jl: Plot map of CI.LJR using GMT.jl 
4. src/08-plot-figure-4.jl: Plot ratio of hydro-elastic to thermo-elastic component in dv/v
5. src/09-plot-figure-5.jl: Plot dv/v and GRACE for 2005 and 2011-2016
6. src/10-plot-figure-6.jl: Plot dv/v near Salton Sea 
7. src/11-plot-figure-7.jl: Plot dv/v relaxation time-series after earthquakes in California
8. src/12-plot-figure-8.jl: Plot dv/v and PGV for Ridgecrest earthquake

## Files for getting/producing data  
1. src/download_zenodo.jl: Downloads data from Zenodo for reproducing figures. Run this before script before running figure scripts. 
2. src/PRISMgetscript.jl: Downloads the PRISM dataset
3. src/bil2netcdf.jl: Converts PRISM .bil files to single netCDF file 
4. src/LJR-psd.jl: Calculates daily power spectral density for station CI.LJR 
5. src/ridgecrest-pgv.jl: Calculates peak ground velcity (PGV) for Ridgecrest earthquake
6. src/make-dispersion.py: Calculate dispersion using Community Velocity Model and 
    Computer Programs in Seismology

## Files to produce supplementary figures 
1. src/plot_CASTAC.jl: Plots groundwater levels in the Castac valley
2. src/plot_DISPERSION.jl: Plots surface wave sensitivity kernels for station CI.LJR 
3. src/plot_LA-DVV.jl: Plots dv/v time series for stations in the Los Angeles Basin  
4. src/plot_NJQ.jl: Plots dv/v and hydrology time series for station CI.NJQ 
5. src/plot_NJQ-map.jl: Plots station location for station CI.NJQ
6. src/plot_RIDGECREST.jl: Plots dv/v for stations that show drop in velocity after Ridgecrest earthquake 
7. src/plot_RXH.jl: Plot dv/v time series for station CI.RXH
velocity after Ridgecrest earthquake 
8. src/plot_SAL.jl: Plot dv/v time series for station CI.SAL
velocity after Ridgecrest earthquake 
9. src/plot_WES.jl: Plot dv/v time series for station CI.WES
