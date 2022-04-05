using ArchGDAL
using Arrow
using DataFrames
using Dates
using Glob
using FFTW
using SeisIO 

function rot2d(x::AbstractVector,y::AbstractVector,θ::Real)
    return x .* cos(θ) .- y .* sin(θ), x .* sin(θ) .+ y .* cos(θ)
end

function findmaxtheta(S::SeisData,n::Real)
    θs = 0: π / 2 / n : π / 2 - (π / 2 / n)
    amps = zeros(n)
    for ii in 1:n 
        x1, y1 = rot2d(S.x[1],S.x[2],θs[ii])
        amps[ii] = maximum([abs.(x1);abs.(y1)])
    end
    return θs,amps
end

# read Baja data 
S = read_data("mseed",joinpath(@__DIR__,"../data/14607652.ms"))

# move channel names into dataframe
nslc = string.(hcat(split.(S.name,".")...))
df = DataFrame(
    NAME=S.name,
    NET=nslc[1,:],
    STA=nslc[2,:],
    LOC=nslc[3,:],
    CHAN=nslc[4,:],
    BAND=String.(SubString.(nslc[4,:],1,1)),
    INST=String.(SubString.(nslc[4,:],2,2)),
    COMP=String.(SubString.(nslc[4,:],3,3))
)

# filter for HHE and HHN components
filter!(x ->x.BAND == "H",df) 
filter!(x ->x.INST == "N",df)
filter!(x ->x.COMP != "1",df)
filter!(x ->x.COMP != "2",df)
filter!(x ->x.COMP != "3",df)
filter!(x ->x.COMP != "Z",df)
filter!(x ->x.NET != "CE",df)
filter!(x ->x.NET != "PG",df)
filter!(x ->x.NET != "WR",df)
filter!(x ->x.NET != "NP",df)
filter!(x ->x.STA != "NBS",df) # NBS has corrupted data 
SHH = S[findall(in(df.NAME),S.name)]

# get stationXML for non CI network stations 
# directory for instrument response 
XMLDIR = joinpath(@__DIR__,"../data/XML/")
XMLfiles = glob("*",XMLDIR)

for ii in 1:size(df,1)
    # create file path 
    xmlpath = joinpath(XMLDIR,"$(df[ii,:NET])_$(df[ii,:STA]).xml")

    # check if already exists 
    if in(xmlpath,XMLfiles)
        continue
    end

    # request from IRIS 
    if df[ii,:NET] == "NC"
        src = "NCEDC"
    else
        src = "IRIS"
    end
    try 
        FDSNsta(df[ii,:NAME],xf=xmlpath,s="2010-01-01",t="2021-01-01",src=src)
    catch 
    end
end

# minimal preprocessing 
demean!(SHH)
ungap!(SHH)
detrend!(SHH)
taper!(SHH)
filtfilt!(SHH,rt="Highpass",fh=0.1)

# add location, gain, response and units to each SeisChannel 
for ii in 1:SHH.n 
    # get start, end times 
    s = Date(u2d(SeisIO.starttime(SHH.t[ii],SHH.fs[ii]) .* SeisIO.μs)) |> string
    t = Date(u2d(SeisIO.endtime(SHH.t[ii],SHH.fs[ii]) .* SeisIO.μs)) + Day(1) |> string 

    # load instrument response
    net, sta, loc, chan = split(SHH.name[ii], ".")
    XMLfile = joinpath(XMLDIR,"$(net)_$sta.xml")   
    println("Reading $XMLfile") 
    R = SeisIO.read_sxml(XMLfile,s,t,false,false,2)

    # add to seischannel 
    SHH.loc[ii] = R.loc[1]
    SHH.gain[ii] = R.gain[1]
    SHH.resp[ii] = R.resp[1]
    SHH.units[ii] = R.units[1]
end

# remove instrument response 
# SHHR = remove_resp(SHH)
SHHR = deepcopy(SHH)

# correct for gain and convert to m/s
for ii in 1:SHHR.n 
    SHHR.x[ii] ./= SHHR.gain[ii] 
    SHHR.gain[ii] = 1.0
end

# sort stations 
ind = sortperm(SHHR.id)
SHHR = SHHR[ind]

# cut channels to have same start and end time 
Spga = Array{SeisData}(undef,SHHR.n ÷ 2)
for ii in 1:SHHR.n ÷ 2  
    SEN = deepcopy(SHHR[(ii-1) * 2 + 1:(ii-1) * 2 + 2])
    # fix timing errors -> channels off by one sample 
    μt = round(Int,1 / SEN.fs[1] * 1e6)
    off1 = SEN.t[1][1,2] % μt
    if off1 >= μt / 2 
        SEN.t[1][1,2] += μt - off1
    else
        SEN.t[1][1,2] -= off1 
    end
    off2 = SEN.t[2][1,2] % μt
    if off2 >= μt / 2 
        SEN.t[2][1,2] += μt - off2
    else
        SEN.t[2][1,2] -= off2 
    end
    sync!(SEN,s="last",t="first")
    Spga[ii] = SEN
end

# find PGV in different frequency bands 
freqmin = 2.0 .^ (0:2)
freqmax = 2.0 .^ (1:3)
PGA = DataFrame()
for ii in 1:length(freqmin)

    # iterate over stations 
    for jj in 1:length(Spga)
        SA = deepcopy(Spga[jj])
        net, sta, loc, chan = split(SA.name[1], ".")
        println("Calculating PGV $net.$sta $(freqmin[ii])-$(freqmax[ii]) Hz $(now())")
        
        # filter in pass band 
        filtfilt!(SA,rt="Bandpass",fl=freqmin[ii],fh=freqmax[ii])

        # throw out unrealistic data due to sensor glitch 
        if maximum(abs.([SA.x[1]; SA.x[2]])) > 1000 # > 1000 cm/s^2 
            continue 
        end

        # get PGV and angle 
        SV = convert_seis(SA)
        angles, PGVs = findmaxtheta(SV,100)
        indV = argmax(PGVs)

        # get PGA and angle 
        angles, PGAs = findmaxtheta(SA,100)
        indA = argmax(PGAs)

        append!(PGA,DataFrame(
            NET=net,
            STA=sta,
            LON=SV.loc[1].lon,
            LAT=SV.loc[1].lat,
            PGV=PGVs[indV] .* 100, # convert to cm/s 
            PGA=PGAs[indA] .* 100, # convert to cm/s^2
            FREQMIN=freqmin[ii],
            FREQMAX=freqmax[ii],
            ),
        )
    end
end

# add VS30 values to PGA
# load VS30 data from USGS https://earthquake.usgs.gov/data/vs30/
ds = ArchGDAL.readraster(joinpath(@__DIR__,"../data/VS30/global_vs30.tif"))
geotransform = ArchGDAL.getgeotransform(ds)
vs30lon = range(geotransform[1],step=geotransform[2],length=ds.size[1])
vs30lat = range(geotransform[4],step=geotransform[6],length=ds.size[2])

# subset to area with PGA measurements 
lonind = findall(minimum(PGA[:,:LON]) - 0.5 .< vs30lon .< maximum(PGA[:,:LON]) + 0.5)
latind = findall(minimum(PGA[:,:LAT]) - 0.5 .< vs30lat .< maximum(PGA[:,:LAT]) + 0.5)
vs30lon = vs30lon[lonind]
vs30lat = vs30lat[latind]
dataset = ds[lonind,latind,1]

# set oceans values to zero 
dataset[dataset .== 600.0] .= 0.0

# find nearest PGA values 
VS30 = zeros(size(PGA,1))
for ii in 1:length(VS30)
    lonind = argmin(abs.(PGA[ii,:LON] .- vs30lon))
    latind = argmin(abs.(PGA[ii,:LAT] .- vs30lat))
    VS30[ii] = dataset[lonind,latind,1]
end

# get shear strain 
PGA[:,:VS30] = VS30
PGA[:,:STRAIN] = PGA[:,:PGV] ./ PGA[:,:VS30]

# write to disk in arrow format 
Arrow.write(joinpath(@__DIR__,"../data/baja-PGV.arrow"),PGA)