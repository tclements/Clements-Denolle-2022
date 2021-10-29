using Distributed
addprocs()

@everywhere using AWS: @service 
@everywhere @service S3
@everywhere begin 
    using AWS, AWSS3, Arrow, DataFrames, Statistics, Dates, Plots 

    function s3_list_all_keys(aws::AWSConfig,bucket::String,prefix::String)
        req = S3.list_objects_v2(bucket,Dict("prefix"=>prefix),aws_config=aws)
        allkeys = [req["Contents"][ii]["Key"] for ii = 1:length(req["Contents"])]

        while req["IsTruncated"] == "true"
            req = S3.list_objects_v2(
                bucket,
                Dict(
                    "prefix"=>prefix,
                    "continuation-token"=>req["NextContinuationToken"]
                ),
                aws_config=aws,
            )
            append!(allkeys,[req["Contents"][ii]["Key"] for ii = 1:length(req["Contents"])])
        end
        return allkeys
    end


    function get_dvv(aws::AWSConfig, bucket::String, stakeys)
        EN = Arrow.Table(s3_get(aws,bucket,stakeys[1])) |> DataFrame
        EZ = Arrow.Table(s3_get(aws,bucket,stakeys[2])) |> DataFrame
        NZ = Arrow.Table(s3_get(aws,bucket,stakeys[3])) |> DataFrame

            # get intersect of all dates 
        dates = intersect(EN[!,:DATE],EZ[!,:DATE],NZ[!,:DATE])
        EN = EN[[d in dates for d in EN[!,:DATE]],:]
        EZ = EZ[[d in dates for d in EZ[!,:DATE]],:]
        NZ = NZ[[d in dates for d in NZ[!,:DATE]],:]
        return EN, EZ, NZ 
    end

    function weighteddvv(EN::DataFrame, EZ::DataFrame, NZ::DataFrame)
        CC = EN[!,:CCPOS] .^ 2 .+ EZ[!,:CCPOS] .^ 2 .+ NZ[!,:CCPOS] .^ 2 .+
                EN[!,:CCNEG] .^ 2 .+ EZ[!,:CCNEG] .^ 2 .+ NZ[!,:CCNEG] .^ 2

        DVV = (EN[!,:CCPOS] .^ 2 .* EN[!,:DVVPOS]) .+ (EZ[!,:CCPOS] .^ 2 .*
            EZ[!,:DVVPOS]) .+ (NZ[!,:CCPOS] .^ 2 .* NZ[!,:DVVPOS]) .+
            (EN[!,:CCNEG] .^ 2 .* EN[!,:DVVNEG]) .+ (EZ[!,:CCNEG] .^ 2 .*
                    EZ[!,:DVVNEG]) .+ (NZ[!,:CCNEG] .^ 2 .* NZ[!,:DVVNEG])

        DVV ./= CC

        CC = (EN[!,:CCPOS] .^ 3 .+ EZ[!,:CCPOS] .^ 3 .+ NZ[!,:CCPOS] .^ 3 .+
            EN[!,:CCNEG] .^ 3 .+ EZ[!,:CCNEG] .^ 3 .+ NZ[!,:CCNEG] .^ 3) ./ CC

        return DVV, CC 
    end

    function plot_dvv(
        aws::AWSConfig,
        bucket::String,
        stakeys::AbstractArray,
        freqmin::Real,
        freqmax::Real,
        FIGURES::String,
    )
        freqsta = join(split(basename(stakeys[1]),".")[1:2],".")
        OUTDIR = joinpath(FIGURES,"$freqmin-$freqmax")
        if !isdir(OUTDIR)
            mkpath(OUTDIR)
        end
        println("$freqsta $freqmin-$freqmax $(now())")

        # get dvv
        EN, EZ, NZ = get_dvv(aws, bucket, stakeys)

        # compute CC from Hobiger, 2014
        DVV, CC = weighteddvv(EN,EZ,NZ)

        # scatter plot DVV
        scatter(
            EN[!,:DATE],
            DVV,
            alpha=max.(0,CC ./ 2),
            ylims=(-3.5,3.5),
            ylabel="dv/v (%)",
            label="",
            title="$freqsta $freqmin-$freqmax Hz",
            size=(1000,400),
            marker_z=CC,
            seriescolor=:Reds_9,
            show=false
        )
        savepath = joinpath(OUTDIR,"$freqsta.png")
        savefig(savepath)
        return nothing
    end

    # load dv/v data
    aws = AWSConfig(region="us-west-2")
    upload_bucket = "bh-auto-corr"
    FIGURES = "/media/FOUR/data/FIGURES/"
    allkeys = s3_list_all_keys(aws,upload_bucket,"DVV/")
    freqranges = unique([dirname(a) for a in allkeys])
    threshold = 0.5 # minimum cc value 
end

# loop through all frequencies 
for ii in 1:length(freqranges)

    # get files for each frequency range
    freqrange = freqranges[ii]
    freqmin, freqmax = parse.(Float64,(split(basename(freqrange),"-")))
    freqind = findall(contains.(allkeys,freqrange))
    freqkeys = allkeys[freqind]

    # get unique stations 
    freqstas = unique([join(split(basename(f),".")[1:2],".") for f in freqkeys])
    stakeys = [sort(freqkeys[ii:ii+2]) for ii = 1:3:length(freqkeys)-1]

    pmap(x -> plot_dvv(aws,upload_bucket,x,freqmin,freqmax,FIGURES),stakeys)
end
