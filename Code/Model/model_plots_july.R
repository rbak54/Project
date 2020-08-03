require("deSolve")
require("ggplot2")

long_data<-read.csv("../../Data/deyle_edited_all.csv")
short_data<-read.csv("../../Data/deyle_edited_summ.csv")


long_data$day_year<-(long_data$week-1)*7+3.5
long_data$day_year<-(long_data$week-1)*7+3.5


short_data$day_peak<-(short_data$peak-1)*7+3.5
short_data$day_peak<-(short_data$peak-1)*7+3.5
#require("tidyr")
library("tidyverse")


#for accessing min temp
#also want to compare flu and temp
#need to edit belpw
data_influenza<-read.csv("../../Data/deyle_data.txt")
data_influenza_order<group_by(data_influenza,variable,country,year) %>% mutate(week=row_number())
data_influenza_order_max<- data_influenza_order %>% group_by(variable,country,year) %>% mutate(no=max(week))
data_influenza_full_years<-data_influenza_order_max[which(data_influenza_order_max$no>51),]
temp<-data_influenza_full_years[which(data_influenza_full_years$country==location),]
temp<-temp[which(temp$variable=="flu"),]
temp$index_week<-seq(1,nrow(temp),by=1)
week_1<-temp[which(temp$week==1),]
temp<-temp[which(temp$value!="NaN"),]
#this would work if before remove NAN
#also need to index
flu_long<-long_data[which(long_data$variable=="flu"),]
temp_long<-long_data[which(long_data$variable=="T"),]
#long_data_less<-long_data[c("country",variable,week)]



just_peak_day<-short_data[,c("country","variable","day_peak")]
wider_just_peak_day<-spread(just_peak_day,variable, day_peak)
head(wider_just_peak_day)
ggplot(data=wider_just_peak_day,aes(x=T,y=flu))+geom_point()+theme_bw()


just_high<-short_data[,c("country","variable","high")]
wider_high<-spread(just_high,variable, high)
ggplot(data=wider_high,aes(x=T,y=flu))+geom_point()+theme_bw()
#averaged this is pretty meaningless
just_low<-short_data[,c("country","variable","low")]
wider_low<-spread(just_low,variable, low)
ggplot(data=wider_low,aes(x=T,y=flu))+geom_point()+theme_bw()




head(long_data)
long_data_summ_peak_value<-long_data %>% group_by(year,country,variable) %>% summarise(year_max=max(value))
long_data_summ_peak_day<-long_data %>% group_by(year,country,variable) %>% summarise(year_peak_day=day_year[which.max(value)])
wider_peak_value<-spread(long_data_summ_peak_value,variable,year_max)
wider_peak_day<-spread(long_data_summ_peak_day,variable,year_peak_day)
wider_peak_day
ggplot(data=wider_peak_value,aes(T,flu,color=country))+geom_point()+theme(legend.position = "none") 
ggplot(data=wider_peak_day,aes(T,flu,color=country))+geom_point()+theme(legend.position = "none") 



for (location in unique(wider_peak_day$country)){
  png(paste0("../../Results/Plots/country_peaks/country_peaks_",location,".png"))
  par(mfrow=c(2,1))
  temp<-wider_peak_day[which(wider_peak_day$country==location),]
  plot(temp$T,temp$flu,xlab="day with peak temperature",ylab="day with peak flu",pch=16,cex=0.5)
  temp<-wider_peak_value[which(wider_peak_value$country==location),]
  plot(temp$T,temp$flu,xlab="peak temperature",ylab="peak flu",pch=16,cex=0.5)
  graphics.off()
}


#should be using min temperature 


#getting data into nice form
data_wider<-read.csv("../../Data/deyle_edited_wide.csv")
data_wider<-data_wider[which(data_wider$T!="NAN"),]
data_wider<-data_wider[which(data_wider$flu!="NAN"),]
data_wider
data_wider_means<-data_wider %>% group_by(country,week) %>% summarise(meanflu=mean(flu),meantemp=mean(T))
data_wider_means$day<-(data_wider_means$week-1)*7+3.5
data_wider_means
#for day
data_wider_means_summ<-data_wider_means %>% group_by(country) %>% summarise(minflu=min(meanflu),maxflu=max(meanflu),peakflu=day[which.max(meanflu)],troughflu=day[which.min(meanflu)],minT=min(meantemp),maxT=max(meantemp),peakT=day[which.max(meantemp)], troughT=day[which.min(meantemp)])
# for week
data_wider_means_summ<-data_wider_means %>% group_by(country) %>% summarise(minflu=min(meanflu),maxflu=max(meanflu),peakflu=week[which.max(meanflu)],troughflu=week[which.min(meanflu)],minT=min(meantemp),maxT=max(meantemp),peakT=week[which.max(meantemp)], troughT=week[which.min(meantemp)])

data_wider_means_summ
#using mean over all years for each country. and the min/maxes of these values. i.e the dayweek with the highest
plot(data_wider_means_summ$troughT,data_wider_means_summ$peakflu)# peak flu time and low temp time related
plot(data_wider_means_summ$peakT,data_wider_means_summ$peakflu)  # negative correlation between high temp time and high flu time
plot(data_wider_means_summ$maxT,data_wider_means_summ$maxflu) #when max  temperature is higher,less flu
plot(data_wider_means_summ$minT,data_wider_means_summ$maxflu) #when min temperature is higher,less flu
plot(data_wider_means_summ$maxT,data_wider_means_summ$minflu) #
plot(data_wider_means_summ$minT,data_wider_means_summ$minflu) #


a<-lm(data_wider_means_summ$peakflu~data_wider_means_summ$troughT)  # increasing minimum T 
summary(a)


ggplot(data=data_wider_means,aes(meantemp,meanflu,color=country))+geom_point()+theme(legend.position ="none")
for (location in unique(data_wider_means$country)){
  temp<-data_wider_means[which(data_wider_means$country==location),]#
  png(paste0("../../Results/Plots/country_week_means/",location,".png"))
  plot(temp$meantemp,temp$meanflu)
  graphics.off()
}

#not as clear but still good
for (location in unique(data_wider$country)){
  temp<-data_wider[which(data_wider$country==location),]#
  png(paste0("../../Results/Plots/country_week/",location,".png"))
  plot(temp$T,temp$flu)
  graphics.off()
}


data_wider_summ<-data_wider %>% group_by(country, year) %>% summarise(minflu=min(flu),maxflu=max(flu),peakflu=week[which.max(flu)],troughflu=week[which.min(flu)],minT=min(T),maxT=max(T),peakT=week[which.max(T)], troughT=week[which.min(T)])
data_wider_summ
for (location in unique(data_wider_summ$country)){
  temp<-data_wider_summ[which(data_wider_summ$country==location),]
  png(paste0("../../Results/Plots/country_year_summ/",location,".png"))
  par(mfrow=c(3,2))
  plot(temp$minT,temp$maxflu)
  plot(temp$maxT,temp$maxflu)
  plot(temp$troughT,temp$peakflu)
  plot(temp$peakT,temp$peakflu)
  plot(temp$troughT,temp$maxflu)
  
  graphics.off()
}
