using Dates

# download precipiation
base_url = "http://services.nacse.org/prism/data/public/4km/ppt"
start = Date(1985,1,1)
stop = Date(now())
date_range = start:Day(1):stop

for ii in eachindex(date_range)
    daypath = joinpath(base_url,Dates.format(date_range[ii],"YYYYmmdd"))
    println("Downloading day $(date_range[ii])")
    command = `wget --content-disposition $daypath`
    run(command)
    sleep(0.1)
end

# download mean temperature
base_url = "http://services.nacse.org/prism/data/public/4km/tmean"
start = Date(1985,1,1)
stop = Date(now())
date_range = start:Day(1):stop

for ii in eachindex(date_range)
    daypath = joinpath(base_url,Dates.format(date_range[ii],"YYYYmmdd"))
    println("Downloading day $(date_range[ii])")
    command = `wget --content-disposition $daypath`
    run(command)
    sleep(0.1)
end
