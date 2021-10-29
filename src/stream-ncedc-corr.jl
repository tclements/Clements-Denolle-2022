using Distributed
addprocs()

@everywhere begin

using Arrow 
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

function daygrouper(df::DataFrame)
    groups = groupby(df,[:NET,:STA,:LOC,:DATE])
    outfiles = []
    for g in groups
        push!(outfiles,string.(g[!,:KEY]))
    end
    return outfiles 
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
data_bucket = "seisbasin"
end 

# get station list 
aws = AWSConfig(region="us-west-2")
hquake = AWSConfig(region="us-west-2",creds=AWSCredentials(profile="hquake"))
df = s3_get(aws,"bh-auto-corr","ncedc.arrow") |> Arrow.Table |> DataFrame 
for netsta in groupby(df,[:NET,:STA])
    net = netsta[1,:NET]
    sta = netsta[1,:STA]

    # get files 
    ind = sortperm(netsta,:KEY)
    netsta = netsta[ind,:]

    # get metadata for each file 
    filesizes = pmap(
        x->
        parse(
            Int,
            s3_get_meta(
                hquake,
                "seisbasin",
                x
            )["Content-Length"]
        ),
        netsta[:,:KEY],
        on_error = e -> (isa(e,AWS.AWSExceptions.AWSException) ? 0 : rethrow()),
    )

    # filter based on file sizes
    ind = findall(.!iszero.(filesizes))
    netsta = netsta[ind,:]
    filesizes = filesizes[ind]
    med = median(filesizes) 
    ind = findall(filesizes .> med / 10)
    netsta = netsta[ind,:]

    if size(netsta,1) == 0 
        continue
    end

    # group by day 
    ind = sortperm(netsta[:,:KEY])
    netsta = netsta[ind,:]
    files = unique(netsta[:,:KEY])

    # group channels by net, sta, loc, channel, date
    outfiles = daygrouper(netsta)
    @eval @everywhere outfiles=$outfiles

    # get the resp 
    if net in ["NN","TA","NP","AZ","G","II","IM","YU","YN","BC","SN"]
        src = "IRIS"
    elseif net in ["NC","BK"]
        src = "NCEDC"
    end

    RESP = SeisData()
    try
        RESP = FDSNsta("$net.$sta..HH?,$net.$sta..BH?",src=src,s="2000-01-01",t="2021-03-30")
    catch 
        RESP = FDSNsta("$net.$sta..HH?,$net.$sta..BH?",src="IRIS",s="2000-01-01",t="2021-03-30")
    end
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
            aws_config=hquake,
            bucket=data_bucket,
            filetype="seisio",
            resp=false,
        ), 
        outfiles,
        retry_check=retry_check,
        retry_delays=ExponentialBackOff(n=3),
    )

    # combine auto-correlations 
    C1,C2,C3 = combine_corr(CS)

    # upload to s3 
    for C in [C1,C2,C3]
        path = "NCCORR/$(C.name)"
        io = IOBuffer()
        serialize(io,C)
        s3_put(aws,upload_bucket,path,io.data)
        close(io)
    end
end
# # read from s3 
# s = IOBuffer(s3_get(aws,upload_bucket,"CORR/$(C1.name)"))
# seekstart(s)
# C = deserialize(s)