#using ncdf files. prob not needed

require(ncdf4)
nc_data <- nc_open('../../Data/Complete_TAVG_EqualArea.nc')
# Save the print(nc) dump to a text file
{
  sink('../../Data/Complete_TAVG_EqualArea.txt')
  print(nc_data)
  sink()
}

lon <- ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude", verbose = F)
time<-ncvar_get(nc_data, "time")
temp.array <- ncvar_get(nc_data, "temperature")
dim(temp.array)
timelastyear<-time[(3240-12):3240]
time[3240]
one<-mean(temp.array[,(3240-12):3240])
plot(lat,one)




nc_data <- nc_open('../../Data/Raw_TAVG_EqualArea.nc')
# Save the print(nc) dump to a text file
{
  sink('../../Data/Raw_TAVG_EqualArea.txt')
  print(nc_data)
  sink()
}

lon <- ncvar_get(nc_data, "longitude")
lat <- ncvar_get(nc_data, "latitude", verbose = F)
time<-ncvar_get(nc_data, "time")
temp.array <- ncvar_get(nc_data, "temperature")
dim(temp.array)
timelastyear<-time[(3240-12):3240]
time[3240]
one<-temp.array[,3000]
plot(lat,one)
