#http://www-das.uwyo.edu/~geerts/cwx/notes/chap16/geo_clim.html 
#sketchy sources but vaguely reasonable results 
#london -51.5 degrees N 
#not using anymore
latitudes<-seq(-90,90,length.out = 100)
#directions<-c(rep("N",length(latitudes)/2),rep("S",length(latitudes)/2))
latitude_temp<-data.frame(matrix(nrow=100,ncol=5))
colnames(latitude_temp)<-c("latitude","low","high","temp","peak")
latitude_temp$latitude<-latitudes
#latitude_temp$direction<-directions
latitude_temp$low<-NA
latitude_temp$high<-NA
latitude_temp$temp<-NA
latitude_temp$peak<-NA

for (i in 1:nrow(latitude_temp)){
  #direction<-latitude_temp$direction[i]
  degrees<-latitude_temp$latitude[i]
  if (degrees>0){
    lat_difference<-degrees-16
    if (lat_difference<0){
      lat_difference<-0
     }
  temp<-27-0.86*lat_difference
  peak<-197
  }
  if (degrees<0){
    lat_difference<-abs(degrees)-20
   if (lat_difference<0){
      lat_difference<-0
     }
  temp<-27-0.63*lat_difference
  peak<-16
  }

range<-0.4*abs(degrees)
high<-temp+range/2
low<-temp-range/2
latitude_temp$low[i]<-low
latitude_temp$high[i]<-high
latitude_temp$temp[i]<-temp
latitude_temp$peak[i]<-peak


}
plot(latitude_temp$latitude,latitude_temp$temp)
plot(latitude_temp$latitude,abs(latitude_temp$high-latitude_temp$low))

write.csv(latitude_temp,"../../Data/simulated_latitude_temp.csv")


read.csv("../../Data/simulated_latitude_temp.csv")
