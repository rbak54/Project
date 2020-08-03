
#specifically when summary not given by model. but unlikely to be necessary fgiven can do this in integration now
########################################################################
#reload integration results and add further info. make more summary tibbles.
out<-read.csv("../../Results/model_results.csv")
out<-as.data.frame(out)
data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ.csv")
out[,"country"]<-data_wider_means_summ[out[,"country"],"country"]
out[,"year"]<-ceiling((out[,"time"])/365)
out[,"day"]<-(out[,"time"]-(out[,"year"]-1)*365)
out[,"week"]<-ceiling((out[,"day"])/7)
model_wide<-out
model_means<-model_wide %>% group_by(country,week,mismatch) %>% summarise(meanI=mean(I),meanR0=mean(R0),meantemp=mean(temperature))
model_means_summ<-model_means %>% group_by(country,mismatch) %>% summarise(minI=min(meanI),maxI=max(meanI),meanI=mean(meanI),peakI=week[which.max(meanI)],troughI=week[which.min(meanI)],minR0=min(meanR0),meanR0=mean(meanR0),maxR0=max(meanR0),peakR0=week[which.max(meanR0)],troughR0=week[which.min(meanR0)],minT=min(meantemp),maxT=max(meantemp),peakT=week[which.max(meantemp)], troughT=week[which.min(meantemp)])
model_summ_year<- model_wide %>% group_by(country,mismatch,year) %>% summarise(minI=min(I),maxI=max(I),peakI=week[which.max(I)],meanI=mean(I),troughI=week[which.min(I)],minR0=min(R0),maxR0=max(R0),meanR0=mean(R0),peakR0=week[which.max(R0)],troughR0=week[which.min(R0)],minT=min(temperature),maxT=max(temperature),peakT=week[which.max(temperature)], troughT=week[which.min(temperature)])
write.csv(model_means,"../../Results/model_means.csv")
write.csv(model_means_summ,"../../Results/model_means_summ.csv")
write.csv(model_summ_year,"../../Results/model_summ_year.csv")



###########

