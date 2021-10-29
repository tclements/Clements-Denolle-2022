using SeisNoise, SeisIO, DataFrames, Glob, SCEDCCorr, ProgressBars

function daygrouper(files)
    N = length(files)
    ii = 1 
    outfiles = []
    while ii < N - 2
        if dirname(files[ii]) == dirname(files[ii + 1]) == dirname(files[ii+2])
            push!(outfiles,files[ii:ii+2])
            ii += 3
        else
            ii += 1
        end
    end
    return outfiles
end

# PARAMETERS
ROOT = expanduser("/media/FOUR/data")
CORROUT = joinpath(ROOT,"BH-CORR")
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
files = glob("continuous_waveforms/*/*/*BH*",ROOT)
files = size_check(files)
outfiles = daygrouper(files)

# run through each pair 
N = length(outfiles)
for ii in ProgressBar(1:N)
    sc_all(outfiles[ii],
	       fs,
	       cc_len,
	       cc_step,
	       freqmin,
	       freqmax,
	       maxlag,
	       CORROUT,
           XMLDIR,
           responsefreq=responsefreq
    )
end

