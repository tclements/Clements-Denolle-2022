using DataFrames
using Dates 
using Glob 
using SCEDCCorr
using SeisIO
using SeisNoise

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

# PARAMETERS
ROOT = expanduser("/media/FOUR/data")
CORROUT = joinpath(ROOT,"CORR")
XMLDIR = joinpath(ROOT,"XML")
if !isdir(CORROUT)
    mkpath(CORROUT)
end
cc_len = 1800
cc_step = 450
fs = 40.
freqmin = 0.5
freqmax = 19.
responsefreq = 0.4
maxlag = 20.

# FINDING PAIRS
files = glob("continuous_waveforms/*/*/*",ROOT)
files = [f for f in files if occursin("ADO",f)]
files = size_check(files)
outfiles = daygrouper(files)

# 30 minute windows 
EN30, EZ30, NZ30 = sc_all(
    outfiles[1],
    fs,
    cc_len,
    cc_step,
    freqmin,
    freqmax,
    maxlag,
    CORROUT,
    XMLDIR,
    responsefreq=responsefreq, 
)

# 15 minute windows 
EN15, EZ15, NZ15 = sc_all(
    outfiles[1],
    fs,
    900,
    225,
    freqmin,
    freqmax,
    maxlag,
    CORROUT,
    XMLDIR,
    responsefreq=responsefreq, 
)

# 7.5 minute windows 
EN7, EZ7, NZ7 = sc_all(
    outfiles[1],
    fs,
    450,
    112.5,
    freqmin,
    freqmax,
    maxlag,
    CORROUT,
    XMLDIR,
    responsefreq=responsefreq, 
)

# filter autocorrelations 2-4 Hz 
EN30clean = clean_up(EN30,2.,4.)
EN15clean = clean_up(EN15,2.,4.)
EN7clean = clean_up(EN7,2.,4.)

# stack using either linear or robust stack 
EN30lin = SeisNoise.stack(EN30clean)
EN30rob = SeisNoise.stack(EN30clean,stacktype=robuststack)
EN15lin = SeisNoise.stack(EN15clean)
EN15rob = SeisNoise.stack(EN15clean,stacktype=robuststack)
EN7lin = SeisNoise.stack(EN7clean)
EN7rob = SeisNoise.stack(EN7clean,stacktype=robuststack)

# plot 30 minute robust vs linear  
lags = -EN15.maxlag:1/EN15.fs:EN15.maxlag
plot(
    lags,
    [abs_max(EN30lin.corr[:]), abs_max(EN30rob.corr[:])],
    label=["30-min Linear" "30-min Robust"], 
    xlabel = "Lag [s]",
)

# zoom in around zero 
ind = findall(abs.(lags) .< 2)
plot(
    lags[ind],
    [abs_max(EN30lin.corr[ind,:]), abs_max(EN30rob.corr[ind,:])],
    label=["30-min Linear" "30-min Robust"], 
    xlabel = "Lag [s]",
)

# plot 15 minute robust vs linear  
plot(
    lags,
    [abs_max(EN15lin.corr[:]), abs_max(EN15rob.corr[:])],
    label=["15-min Linear" "15-min Robust"], 
    xlabel = "Lag [s]",
)

# zoom in around zero 
plot(
    lags[ind],
    [abs_max(EN15lin.corr[ind,:]), abs_max(EN15rob.corr[ind,:])],
    label=["15-min Linear" "15-min Robust"], 
    xlabel = "Lag [s]",
)

# plot 7 minute robust vs linear  
plot(
    lags,
    [abs_max(EN7lin.corr[:]), abs_max(EN7rob.corr[:])],
    label=["7-min Linear" "7-min Robust"], 
    xlabel = "Lag [s]",
)

# zoom in around zero 
plot(
    lags[ind],
    [abs_max(EN7lin.corr[ind,:]), abs_max(EN7rob.corr[ind,:])],
    label=["7-min Linear" "7-min Robust"], 
    xlabel = "Lag [s]",
)

# do one month of auto-correlations 
N = length(outfiles)
CSlin = Array{Array{CorrData,1}}(undef,N)
CSrob = Array{Array{CorrData,1}}(undef,N)
for ii = 1:N
    println("Correlating $(basename(outfiles[ii][1])) $(now())")
    C1,C2,C3 = sc_all(
        outfiles[ii],
        fs,
        900,
        225,
        freqmin,
        freqmax,
        maxlag,
        CORROUT,
        XMLDIR,
        responsefreq=responsefreq, 
    ) 
    ENlin, EZlin, NZlin = map(SeisNoise.stack,[C1,C2,C3])
    CSlin[ii] = [ENlin, EZlin, NZlin]
    ENrob, EZrob, NZrob = map(robuststack,[C1,C2,C3])
    CSrob[ii] = [ENrob, EZrob, NZrob]
end

# need to filter then linear / robust stack 
ENlin, EZlin, ENlin = combine_corr(CSlin)
ENrob, EZrob, ENrob = combine_corr(CSrob)

