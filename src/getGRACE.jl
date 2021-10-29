using ProgressBars
# download GRACE data 
base_url = "https://nasagrace.unl.edu/GRACE/NASApublication/data/"

# get file names 
files = readlines("/media/FOUR/data/grace-list.txt")

# create output directory 
OUTDIR = "/media/FOUR/data/GRACE"
if !isdir(OUTDIR)
    mkpath(OUTDIR)
end

# download files 
for ii in ProgressBar(1:length(files)) 
    infile = joinpath(base_url,files[ii])
    outfile = joinpath(OUTDIR,files[ii])
    download(infile,outfile)
    sleep(0.1)
end

