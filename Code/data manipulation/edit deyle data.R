#5th july
#modifying delye data- tidy, average over multiple years, and obtain min , max and week where peak occurs


data_influenza<-read.csv("../../Data/deyle_data.txt")
require("tidyverse")
require("ggplot2")
#alreasdy per capita

#data_influenza<-data_influenza[which(data_influenza$country=="Luxembourg"),]
data_influenza_order<-dplyr::group_by(data_influenza,variable,country,year) %>% mutate(week=row_number())
data_influenza_order_max<- data_influenza_order %>% group_by(variable,country,year) %>% mutate(no=max(week))
data_influenza_full_years<-data_influenza_order_max[which(data_influenza_order_max$no>51),]
data_edited_wide<-data_influenza_full_years %>% spread(variable,value)
write.csv(data_edited_wide,"../../Data/deyle_edited_wide.csv")


data_influenza_full_value<-data_influenza_full_years[which(data_influenza_full_years$value!="NaN"),]
data_influenza_full_value_mean<-data_influenza_full_value %>% group_by(variable,country,week) %>% mutate(avg=mean(value))
data_influenza_full_value_mean[which(data_influenza_full_value_mean$variable=="T"),"avg"]<-(data_influenza_full_value_mean[which(data_influenza_full_value_mean$variable=="T"),"avg"]-32)*(5/9)

shorter_all<-data_influenza_full_value_mean
#shorter_all[which(shorter_all$variable=="T"),"avg"]<-(shorter_all[which(shorter_all$variable=="T"),"avg"]-32)*(5/9)
shorter_all_summ<-shorter_all %>% group_by(country,variable) %>%  summarise(low=min(avg), high=max(avg), peak=week[which.max(avg)] )

#now need to find mean flu for each week of the year in each country
flu<-data_influenza_full_value_mean[which(data_influenza_full_value_mean$variable=="flu"),]

data_flu_summ<-flu %>% group_by(variable,country) %>% summarise(min_value=min(avg),max_value=max(avg),max_week=which.max(avg))

temp<-data_influenza_full_value_mean[which(data_influenza_full_value_mean$variable=="T"),]
#data_temp_summ<-temp %>% group_by(variable,country) %>% summarise(min_value=min(avg),max_value=max(avg),max_week=which.max(avg))

ah<-data_influenza_full_value_mean[which(data_influenza_full_value_mean$variable=="AH"),]
#data_ah_summ<-ah %>% group_by(variable,country) %>% summarise(min_value=min(avg),max_value=max(avg),max_week=which.max(avg))

prec<-data_influenza_full_value_mean[which(data_influenza_full_value_mean$variable=="PRCP"),]
#data_prec_summ<-prec %>% group_by(variable,country) %>% summarise(min_value=min(avg),max_value=max(avg),max_week=which.max(avg))


#T, PRCP, AH 
ggplot(data=flu,aes(x=week,y=avg,group=country,col=country)) +geom_point() + geom_line() + theme(legend.position = "none") 
ggplot(data=temp,aes(x=week,y=avg,group=country,col=country)) +geom_point() + geom_line() +theme_bw() + theme(legend.position = "none")
ggplot(data=prec,aes(x=week,y=avg,group=country,col=country)) +geom_point() + geom_line() + theme(legend.position = "none") 
ggplot(data=ah,aes(x=week,y=avg,group=country,col=country)) +geom_point() + geom_line() + theme(legend.position = "none") 


 write.csv(shorter_all,"../../Data/deyle_edited_all.csv")
 write.csv(shorter_all_summ,"../../Data/deyle_edited_summ.csv")
 
 temp_summ<-shorter_all_summ[which(shorter_all_summ$variable=="T"),]
flu_summ<-shorter_all_summ[which(shorter_all_summ$variable=="flu"),]
prec_summ<-shorter_all_summ[which(shorter_all_summ$variable=="PRCP"),]
ah_summ<-shorter_all_summ[which(shorter_all_summ$variable=="AH"),]

for (location in unique(data_influenza_full_years$country)){
  temp<-data_influenza_full_years[which(data_influenza_full_years$country==location),]
  temp<-temp[which(temp$variable=="flu"),]
  temp$index_week<-seq(1,nrow(temp),by=1)
  week_1<-temp[which(temp$week==1),]
  temp<-temp[which(temp$value!="NaN"),]
  png(paste0("../../Results/Plots/country_flu/",location,".png"))
  plot(temp$index_week,temp$value,cex=0.5,ylim=c(0,0.01))
  abline(v=unique(week_1$index_week),col="red")
  graphics.off()
}
#most yearly but some really messy




  
