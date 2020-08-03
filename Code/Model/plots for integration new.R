require(ggplot2)

data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ_POP.csv")
data_wider_summ<-read.csv("../../Data/data_wider_summ_POP.csv")
data_wider_means<-read.csv("../../Data/data_wider_means_POP.csv")
require(tidyverse)
sims=80
parms=list(climate_label="Temperature")
correlation_df<-read.csv(paste0("../../Results/fromfunction/cors/",sims,parms[["climate_label"]],"correlation_dataframe.csv"))
std <- function(x) sd(x)/sqrt(length(x))
correlation_df$Mismatch<-as.factor(correlation_df$mismatch)

bests<-correlation_df %>% group_by(country) %>% summarise(best=Mismatch[which.max(corsI)],.groups="keep")
correlation_df_means<-as_tibble(correlation_df) %>% group_by(Mismatch) %>% mutate(means=mean(corsI),errors=std(corsI))

png(paste0("../../Results/Plots/boxandwhisker",sims,parms[["climate_label"]],".png"))
print(ggplot(data=correlation_df, aes(x= Mismatch,y= corsI)) +geom_boxplot()  +theme_bw()+xlab("Mismatch") +ylab("Mean Correlation"))
graphics.off()

png(paste0("../../Results/Plots/erros",sims,parms[["climate_label"]],".png"))
print(ggplot(data=correlation_df_means, aes(x= Mismatch,y= means)) +geom_point()+theme_bw()+
        geom_errorbar(aes(ymin=means+errors,ymax=means-errors))+xlab("Mismatch") +ylab("Mean Correlation"))
graphics.off()

correlation_df_means<-as_tibble(correlation_df) %>% group_by(Mismatch,lat) %>% summarise(means=mean(corsI),errors=std(corsI),.groups="keep")
png(paste0("../../Results/Plots/lats",sims,parms[["climate_label"]],".png"))
print(ggplot(data=correlation_df_means, aes(x= abs(lat), col=Mismatch,y= means)) +scale_fill_discrete(name = "Dose")+geom_point()+theme_bw()+
        geom_errorbar(aes(ymin=means+errors,ymax=means-errors)) +xlab("Absolute Value of Latitude") +ylab("Mean Correlation") )
graphics.off()
model_means<-read.csv(paste0("../../Results/fromfunction/",sims,parms[["climate_label"]],".csv"))


m<-model_means %>% group_by(country,mismatch) %>% summarise(meanmeanR0=mean(meanR0),meanmeanI=mean(meanI),maxmeanR0=max(meanR0),maxmeanI=max(meanI),errors=std(meanR0),.groups="keep")
png(paste0("../../Results/Plots/mismatchr0",sims,parms[["climate_label"]],".png"))
print(ggplot(data=m,aes(mismatch,meanmeanR0,group=country, colour = as.factor(country)))+geom_point()+geom_line()+
  theme_bw()+theme(legend.position = "none")+ylab("Mean R0")+geom_errorbar(aes(ymin=meanmeanR0+errors,ymax=meanmeanR0-errors),width=0.02))
graphics.off()
png(paste0("../../Results/Plots/mismatchr0",sims,parms[["climate_label"]],".png"))
print(ggplot(data=m,aes(mismatch,meanmeanR0,group=country, colour = as.factor(country)))+geom_point()+geom_line()+
        theme_bw()+theme(legend.position = "none")+ylab("Mean R0")+geom_errorbar(aes(ymin=meanmeanR0+errors,ymax=meanmeanR0-errors),width=0.02))
graphics.off()
m<-model_means %>% group_by(country,mismatch,week) %>% summarise(meanmeanR0=mean(meanR0),meanmeanI=mean(meanI),maxmeanR0=max(meanR0),maxmeanI=max(meanI),errorsI=std(meanI),errorsR=std(meanR0),.groups="keep")
model_means<-m
##########################################################
#R0 over time for model ,  meanI over time for model, flu over time for data
#ggplot
for (i in unique(m$country)){
  temp<-m[which(m$country==i),]
  temp2<-data_wider_means[which(data_wider_means$country==unique(data_wider_means$country)[i]),]
#print(ggplot(data=temp,aes(week,meanR0,color=as.factor(mismatch)),group=as.factor(mismatch))+geom_point()+geom_line())
   png(paste0("../../Results/Plots/R0_time/",i,"R0.png"))
 print(ggplot(data=temp,aes(week,meanmeanR0,color=as.factor(mismatch)),group=as.factor(mismatch))+
    geom_errorbar(aes(ymin=meanmeanR0+errorsR,ymax=meanmeanR0-errorsR))+geom_point()+geom_line()+theme_bw()+
    ylab("mean R0"))
 graphics.off()
 
  png(paste0("../../Results/Plots/R0_time/",i,"I.png"))
   print(ggplot(data=temp,aes(week,meanmeanI,color=as.factor(mismatch)),group=as.factor(mismatch))+
    geom_errorbar(aes(ymin=meanmeanI+errorsI,ymax=meanmeanI-errorsI))+geom_point()+geom_line()+theme_bw()+
    ylab("mean I"))
   graphics.off()
   
   png(paste0("../../Results/Plots/R0_time/",i,"flu.png"))
    print(ggplot(data = temp2,aes(week,meanflu))+geom_point() +geom_line()+theme_bw()+ylab("mean per capita flu"))
    graphics.off()
    
  #BASER
    for (i in unique(model_means$country)){
      png(paste0("../../Results/Plots/R0_time/",i,".png"))
      par(mfrow=c(3,1))
      temp<-model_means[which(model_means$country==i),]
      plot(x=temp$week,y=temp$meanmeanR0,"n",xlab="week",ylab="mean R0")
      cols<-rainbow(length(unique(temp$mismatch)))
      
      for (j in 1:length(unique(temp$mismatch))){
        temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
        lines(y=temptemp$meanmeanR0,x=temptemp$week,col=cols[j],pch=16) 
      }
      legend("topright",legend=c(unique(temp$mismatch)),col=c(cols),pch=16)
      plot(x=temp$week,y=temp$meanmeanI,"n",xlab="week",ylab="mean I")
      for (j in 1:length(unique(temp$mismatch))){
        temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
        lines(y=temptemp$meanmeanI,x=temptemp$week,col=cols[j],pch=16) 
      }
      legend("topright",legend=c(unique(temp$mismatch)),col=c(cols),pch=16)
      temp<-data_wider_means[which(data_wider_means$country==unique(data_wider_means$country)[i]),]
      plot(temp$week,temp$meanflu,type="s",xlab="week",ylab="mean flu cases per capita")
      graphics.off()
    }
  # par(mfrow=c(3,1))
  # plot(x=temp$week,y=temp$meanR0,"n")
  # cols<-rainbow(length(unique(temp$mismatch)))
  # 
  # for (j in 1:length(unique(temp$mismatch))){
  #   temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
  #   lines(y=temptemp$meanR0,x=temptemp$week,col=cols[j],pch=16) 
  # }
  # legend("topright",legend=c(unique(temp$mismatch)),col=c(cols),pch=16)
  # plot(x=temp$week,y=temp$meanI,"n")
  # for (j in 1:length(unique(temp$mismatch))){
  #   temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
  #   lines(y=temptemp$meanI,x=temptemp$week,col=cols[j],pch=16) 
  # }
  # legend("topright",legend=c(unique(temp$mismatch)),col=c(cols),pch=16)
  # temp<-data_wider_means[which(data_wider_means$country==unique(data_wider_means$country)[i]),]
  # plot(temp$week,temp$meanflu,type="s")
  # graphics.off()
}
##############