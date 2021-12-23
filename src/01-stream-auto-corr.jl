# this script: 
# 1. streams files from the scedc-pds public bucket in s3 to a compute instance on ec2 
# 2. computes the single-station autocorrelation using SeisNoise.jl
# 3. saves single-station autocorrelations in a bucket in s3

using Distributed
addprocs()

@everywhere begin

using AWS
using AWSS3
using Dates
using DataFrames
using SCEDC
using SCEDCCorr
using SeisIO
using SeisNoise
using Serialization
using Statistics 

function daygrouper(files)
    N = length(files)

    # get net, sta, chan , location, date 
    nslcd = replace.(basename.(files),".ms"=>"")
    net = Array{String}(undef,N)
    sta = Array{String}(undef,N)
    loc = Array{String}(undef,N)
    chan = Array{String}(undef,N)
    date = Array{DateTime}(undef,N)
    for ii = 1:N
        net[ii] = nslcd[ii][1:2]
        sta[ii] = replace(nslcd[ii][3:7],"_"=>"")
        chan[ii] = nslcd[ii][8:10]
        loc[ii] = replace(nslcd[ii][11:13],"_"=>"")
        date[ii] = yyyyjjj2date(nslcd[ii][14:end])
    end

    # create dataframe 
    df = DataFrame(
        "PATH"=>files,
        "NET"=>net,
        "STA"=>sta,
        "LOC"=>loc,
        "CHAN"=>chan,
        "DATE"=>date,
    )

    # group channels by net, sta, loc, channel, date
    groups = groupby(df,[:NET,:STA,:LOC,:DATE])

    outfiles = []
    for g in groups
        # check that each station has BHE, BHN, BHZ components
        if g[!,:CHAN] == ["BHE","BHN","BHZ"]
            push!(outfiles,string.(g[!,:PATH]))
        end
    end
    return outfiles
end

# from parallelism.jl 
function retry_check(delay_state, err)
    # Below each condition to retry is listed along with an explanation about why
    # retrying should/might work.
    should_retry = (
        # Worker death is normally stocastic, if not then doesn't matter how many
        # retries as it will rapidly kill all workers
        err isa ProcessExitedException ||
        # If we are in the middle of fetching data and the process is killed we could
        # get an ArgumentError saying that the stream was closed or unusable.
        # So same as above.
        err isa ArgumentError && occursin("stream is closed or unusable", err.msg) ||
        # In general IO errors can be transient and related to network blips
        err isa Base.IOError
    )
    return should_retry
end

function combine_corr(CS)
    Ncols = length(CS)
    Nrows = size(CS[1][1].corr,1)
    C1data = zeros(eltype(CS[1][1].corr),Nrows,Ncols)
    C2data = similar(C1data)
    C3data = similar(C1data)
    t = zeros(eltype(CS[1][1].t),Ncols)
    for ii = 1:length(CS)
        C1data[:,ii] .= CS[ii][1].corr[:]
        C2data[:,ii] .= CS[ii][2].corr[:]
        C3data[:,ii] .= CS[ii][3].corr[:]
        t[ii] = CS[ii][1].t[1]
        println("Added day $(u2d(t[ii])) $ii of $(length(CS)) $(now())")
    end
    
    # allocate new CorrData
    C1 = deepcopy(CS[1][1])
    C2 = deepcopy(CS[1][2])
    C3 = deepcopy(CS[1][3])

    # update data 
    C1.corr = C1data
    C2.corr = C2data
    C3.corr = C3data

    # update time 
    C1.t = t 
    C2.t = t 
    C3.t = t

    return C1,C2,C3
end

cc_len = 1800
cc_step = 450
fs = 40.
freqmin = 0.5
freqmax = 19.
responsefreq = 0.4
maxlag = 32.
upload_bucket = "bh-auto-corr"
end 

# get station list 
aws = AWSConfig(region="us-west-2")
stations = s3_get(aws,upload_bucket,"stations")
stations = String(stations)
stations = split(stations,"\n")
for ii = 1:length(stations)
    station = stations[ii]
    net = station[1:2]
    sta = station[3:end]

    # get files 
    athenabucket = "scedc-athena-queries"
    query = "seedchan LIKE 'BH%' and net = '$net' and sta = '$sta' and location = '--'"
    filelist = athenaquery(athenabucket,query)
    sort!(filelist)

    # get metadata for each file 
    filesizes = pmap(
        x->
        parse(
            Int,
            s3_get_meta(
                aws,
                "scedc-pds",
                x
            )["Content-Length"]
        ),
        filelist
    )

    # filter based on file sizes
    med = median(filesizes) 
    ind = findall(filesizes .> med / 10)
    files = filelist[ind]

    # group by day 
    sort!(files)
    files = unique(files)
    outfiles = daygrouper(files)

    # get the resp 
    XMLfile = "FDSNstationXML/CI/$(net)_$sta.xml"
    stream = s3_get(aws,"scedc-pds",XMLfile,raw=true)
    xsta = String(stream);
    RESP = SeisIO.FDSN_sta_xml(xsta, false, "0001-01-01T00:00:00", "9999-12-31T12:59:59", 0)
    ind = findall(occursin.("BH",RESP.id))
    RESP = RESP[ind]
    update_resp_t!(RESP)

    # auto correlate 
    CS = pmap(
        x -> stream_autocorr(
            x,
            fs,
            cc_len,
            cc_step,
            freqmin,
            freqmax,
            maxlag,
            RESP,
            responsefreq=responsefreq,
        ), 
        outfiles,
        retry_check=retry_check,
        retry_delays=ExponentialBackOff(n=3),
    )

    # combine auto-correlations 
    C1,C2,C3 = combine_corr(CS)

    # upload to s3 
    for C in [C1,C2,C3]
        path = "CORR/$(C.name)"
        io = IOBuffer()
        serialize(io,C)
        s3_put(upload_bucket,path,io.data)
        close(io)
    end
end
# # read from s3 
# s = IOBuffer(s3_get(aws,upload_bucket,"CORR/$(C1.name)"))
# seekstart(s)
# C = deserialize(s)