using Dates
using Glob  
using NetCDF 

files = glob("*","/media/FOUR/data/GRACE/")
sort!(files)
N = length(files)
t = Array{String}(undef,N)
gws_inst = zeros(Float32,464,224,N)
rtzsm_inst = similar(gws_inst)
sfsm_inst = similar(gws_inst)
for ii = 1:N
    println(ii,now())
    t[ii] = ncgetatt(files[ii],"time","begin_date")
    gws_inst[:,:,ii] .= ncread(files[ii],"gws_inst")[:,:,1]
    rtzsm_inst[:,:,ii] .= ncread(files[ii],"rtzsm_inst")[:,:,1]
    sfsm_inst[:,:,ii] .= ncread(files[ii],"sfsm_inst")[:,:,1]
end
lon = ncread(files[1],"lon")
lat = ncread(files[1],"lat")

tim = [Dates.datetime2unix(DateTime(Date(d,"yyyymmdd"))) for d in t]
attribs = Dict("units" => "%", "data_min" => -999., "data_max" =>100.)
lonatts = Dict("longname" => "Longitude",
          "units"    => "degrees east")
latatts = Dict("longname" => "Latitude",
          "units"    => "degrees north")
timatts = Dict("longname" => "Time",
          "units"    => "Seconds since 1970")

# convert to netcdf file
filename = "/media/FOUR/data/GRACE.nc"
nccreate(filename,"gws_inst","lat",lat,latatts,"lon",lon,lonatts,"t",tim,timatts,atts=attribs)
ncwrite(permutedims(gws_inst,[2,1,3]),filename,"gws_inst")
nccreate(filename,"rtzsm_inst","lat",lat,latatts,"lon",lon,lonatts,"t",tim,timatts,atts=attribs)
ncwrite(permutedims(rtzsm_inst,[2,1,3]),filename,"rtzsm_inst")
nccreate(filename,"sfsm_inst","lat",lat,latatts,"lon",lon,lonatts,"t",tim,timatts,atts=attribs)
ncwrite(permutedims(sfsm_inst,[2,1,3]),filename,"sfsm_inst")
