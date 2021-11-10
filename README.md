# Clements-Denolle-2022
Code for replication of Clements &amp; Denolle, 2022 

Files to reproduce

1. src/stream-auto-corr.jl: Autocorrelation of seismic noise from California on AWS
2. src/single-station-dvv.jl: dv/v of single-station autocorrelations on AWS
3. src/dvv-comp.jl: combine EN, EZ, NZ componentsof dv/v into single single time-series
4. src/fit-hydro-models.jl: Fit hydro-elastic and thermo-elastic models to dv/v time series 

Figures: 
1. src/plot-figure-1.jl: Plot map of California seismic stations 
2. src/plot-figure-2-3.jl: Plot single-station autocorrelation and noise spectrum of CI.LJR
3. src/plot-figure-2-3.jl: Plot CI.LJR dv/v time-series 
3a. src/plot-figure-3-map.jl: Plot map of CI.LJR using GMT.jl 
4. src/plot-figure-4.jl: Plot ratio of hydro-elastic to thermo-elastic component in dv/v
5. src/plot-figure-5.jl: Plot dv/v and GRACE for 2005 and 2011-2016
6. src/plot-figure-6.jl: Plot dv/v near Salton Sea 
7. src/plot-figure-7.jl: Plot dv/v relaxation time-series after earthquakes in California
8. src/plot-figure-8.jl: Plot dv/v and PGV for Ridgecrest earthquake