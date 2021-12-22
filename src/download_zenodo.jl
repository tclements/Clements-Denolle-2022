# this script downloads the replication data from zenodo 
# downloading from zenodo may take a while depending on internet connection....
using InfoZIP

# location of data.zip file 
zenodo_url = "https://zenodo.org/record/5794562/files/data.zip"

# save data.zip at ~/Clements-Denolle-2022/data.zip 
savefile = joinpath(dirname(@__FILE__), "..", basename(zenodo_url))

# download data.zip 
download(zenodo_url, savefile)

# unzip data.zip into ~/Clements-Denolle-2022/data/
datadir = joinpath(dirname(@__FILE__), "..")
InfoZIP.unzip(savefile, datadir)