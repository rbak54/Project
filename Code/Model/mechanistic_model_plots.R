#data and functions 


#source("modelling_functions.R")
require("tidyverse")
require("ggplot2")
data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ.csv")
data_wider_summ<-read.csv("../../Data/data_wider_summ.csv")
data_wider_means<-read.csv("../../Data/data_wider_means.csv")

#reload summary results
model_means<-read.csv("../../Results/model_means.csv")
model_means_summ<-read.csv("../../Results/model_means_summ.csv")
model_summ_year<-read.csv("../../Results/model_summ_year.csv")

##########################################################################
#plots

  #all model means- average of models for each timepoint 
ggplot(data=model_means, aes(meantemp,meanI,col=mismatch))+geom_point() #too complex to be useful
ggplot(data=model_means,aes(meantemp,meanI,color=country))+geom_point()+theme(legend.position ="none")
ggplot(data=model_means,aes(meantemp,meanR0,color=country))+geom_point()+theme(legend.position ="none")
ggplot(data=model_means,aes(meantemp,meanR0,color=mismatch,group=mismatch))+geom_point() 

  #summary -DATA
plot(data_wider_means_summ$troughT,data_wider_means_summ$peakflu)# peak flu time and low temp time related
plot(data_wider_means_summ$peakT,data_wider_means_summ$peakflu)  # negative correlation between high temp time and high flu time- looks almost lke 30 week lag
plot(data_wider_means_summ$maxT,data_wider_means_summ$maxflu) #when max  temperature is higher,less flu
plot(data_wider_means_summ$maxT,data_wider_means_summ$minflu) #?
plot(data_wider_means_summ$minT,data_wider_means_summ$maxflu) #when min temperature is higher,less flu
plot(data_wider_means_summ$minT,data_wider_means_summ$minflu)#?

  #summary -models -I
ggplot(data=model_means_summ, aes(troughT,peakI,col=mismatch))+geom_point()# peak flu time and low temp time related.  mismatch determines lag. mismatch od 0 looks most like data
ggplot(data=model_means_summ, aes(peakT,peakI,col=mismatch))+geom_point() #oposite to expected. seems to suggest that mismatch determines extent of lag but peak temp does map lag
ggplot(data=model_means_summ, aes(maxT,maxI,col=mismatch))+geom_point()# for some mismatches, increase temp increases I, for some it decreases 
ggplot(data=model_means_summ, aes(maxT,minI,col=mismatch))+geom_point() #max/min used not all that insightful.  country means more insigntful
ggplot(data=model_means_summ, aes(minT,maxI,col=mismatch))+geom_point()
ggplot(data=model_means_summ, aes(minT,minI,col=mismatch))+geom_point()

 #summary- models- R0
ggplot(data=model_means_summ, aes(troughT,peakR0,col=mismatch))+geom_point()# peak flu time and low temp time related
ggplot(data=model_means_summ, aes(peakT,peakR0,col=mismatch))+geom_point()
ggplot(data=model_means_summ, aes(maxT,maxR0,col=mismatch))+geom_point() #these show importance of mismatch
ggplot(data=model_means_summ, aes(maxT,minR0,col=mismatch))+geom_point()# but also always decrease when look at R0, above confounded by seir dynamics
ggplot(data=model_means_summ, aes(minT,maxR0,col=mismatch))+geom_point()
ggplot(data=model_means_summ, aes(minT,minR0,col=mismatch))+geom_point()

#too complex to be v helpful
ggplot(data=data_wider_means,aes(meantemp,meanflu,color=country))+geom_point()+theme(legend.position ="none")
ggplot(data=model_means,aes(week,meanR0,col=mismatch))+geom_point()
ggplot(data=data_wider,aes(y=flu,x=T))+geom_point()

ggplot(data=model_means,aes(meantemp,meanI,col=as.factor(mismatch)))+geom_point()
#plots of temp and I/R0. these are extremely hard to find match with mismatch. bu
#something happens around 15degrees for 0 mismatch
#why so different to data. seems fairly likely thaty my constact rate thing isn't so good. or data not great or missing something key EG humidity
par(mfrow=c(3,2))
for(i in unique(model_means$mismatch)){
  temporary<-model_means[which(model_means$mismatch==i),]
  plot(temporary$meantemp,temporary$meanI,xlim=c(-30,40),ylim=c(0.05,0.35))
}
for(i in unique(model_means$mismatch)){
  temporary<-model_means[which(model_means$mismatch==i),]
  plot(temporary$meantemp,temporary$meanR0)
}
plot(data_wider_means$meantemp,data_wider_means$meanflu,ylim =c(0.01,0.1))
#mismatch and R0/incidence
model_means_summ$mismatch<-as.factor(model_means_summ$mismatch)
plot(model_means_summ$mismatch,model_means_summ$meanI) #LOW OR HIGH MORE DISEASE
plot(model_means_summ$mismatch,model_means_summ$meanR0)# LOW HIGHER r0 , INTERMEDIATE LOWER
plot(model_means_summ$mismatch,model_means_summ$maxI) #HIGHEST PEAK SIZE AT INTERMEDIATE, LOWER PEAK SIZE FOR HIGH MISMATCH
plot(model_means_summ$mismatch,model_means_summ$maxR0)#LOWEST MAX r0 AT INTERMEDIATE??
model_summ_year$mismatch<-as.factor(model_summ_year$mismatch)
plot(model_summ_year$mismatch,model_summ_year$meanI)
plot(model_summ_year$mismatch,model_summ_year$meanR0)
plot(model_summ_year$mismatch,model_summ_year$maxI) 
plot(model_summ_year$mismatch,model_summ_year$maxR0)#THIS IS SUBSTATIONALLY DIFFERENT TO ABOVE
#when you look at I intermediate mismatch has the least disease
#when you look at R0 , full mismatch lower R0 then 0 mismatch but lowest when R0 is 0.5
#

#lowest mismatch
#lower mismatches have worse disease for most of these plots 
#plots for individual countries
############################################################################
### flu incidence at different temperatures using weekly summaries. - R0, I (different mismatch) and data

# below shows huge effect of mismatch. but also that a low mismatch pattern is most common
for (i in unique(model_means$country)){
  png(paste0("../../Results/Plots/country_week_means/",i,".png"))
  par(mfrow=c(2,2))
  temp<-model_means[which(model_means$country==i),]
  plot(temp$meantemp,temp$meanI,"n")
  cols<-rainbow(length(unique(temp$mismatch)))
  for (j in 1:length(unique(temp$mismatch))){
    temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
    lines(y=temptemp$meanI,x=temptemp$meantemp,col=cols[j],pch=16,type="b") 
    
}
  legend("topright",legend=c(unique(temp$mismatch)),col=c(cols),pch=16)
  plot(x=temp$meantemp,temp$meanR0,"n")
  cols<-rainbow(length(unique(temp$mismatch)))
  for (j in 1:length(unique(temp$mismatch))){
    temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
    lines(y=temptemp$meanR0,x=temptemp$meantemp,col=cols[j],pch=16,type="b") 
  }
  temp<-data_wider_means[which(data_wider_means$country==i),]
  plot(temp$meantemp,temp$meanflu)
  
  graphics.off()
}
# mismatch big but not consistent effect on time of peak. low mismatch often means trough of T close to peak of I


################################################################
#week at which highest flu vs week at which lowest t

for (i in unique(model_summ_year$country)){
  png(paste0("../../Results/Plots/country_week_means_summ/",i,".png"))
  par(mfrow=c(2,2))
  
  temp<-model_summ_year[which(model_summ_year$country==i),]
  cols<-rainbow(length(unique(temp$mismatch)))
  
  plot(temp$troughT,temp$peakI,"n")
  for (j in 1:length(unique(temp$mismatch))){
    temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
    points(y=temptemp$peakI,x=jitter(temptemp$troughT,3),col=cols[j],pch=16) 
  }
  legend("topright",legend=c(unique(temp$mismatch)),col=c(cols),pch=16)
  plot(x=temp$troughT,temp$peakR0,"n")
  for (j in 1:length(unique(temp$mismatch))){
    temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
    points(y=temptemp$peakR0,x=jitter(temptemp$troughT,3),col=cols[j],pch=16) 
  }
  legend("topright",legend=c(unique(temp$mismatch)),col=c(cols),pch=16)

  temp<-data_wider_summ[which(data_wider_summ$country==i),]
  plot(temp$troughT,temp$peakflu)

  graphics.off()
}
#added jitter to see roughly how many points. but all same x value 
############################################################
#model_means_summ$mismatch<-as.factor(model_means_summ$mismatch)


##########################################################
#temp and flu for all data
for (location in unique(data_wider$country)){
  temp<-data_wider[which(data_wider$country==location),]#
  png(paste0("../../Results/Plots/country_week/",location,".png"))
  plot(temp$T,temp$flu)
  graphics.off()
}
#########################################################

#not all that helpful. summary stat graphs for each year for each country 
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
##########################################################

model_means<-read.csv("../../Results/model_results_summary.csv")
##########################################################
#R0 over time for model ,  meanI over time for model, flu over time for data
for (i in unique(model_means$country)){
  png(paste0("../../Results/Plots/R0_time/",i,".png"))
  par(mfrow=c(3,1))
  temp<-model_means[which(model_means$country==i),]
  plot(x=temp$week,y=temp$meanR0,"n",xlab="week",ylab="mean R0")
  cols<-rainbow(length(unique(temp$mismatch)))
  
  for (j in 1:length(unique(temp$mismatch))){
    temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
    lines(y=temptemp$meanR0,x=temptemp$week,col=cols[j],pch=16) 
  }
  legend("topright",legend=c(unique(temp$mismatch)),col=c(cols),pch=16)
  plot(x=temp$week,y=temp$meanI,"n",xlab="week",ylab="mean I")
  for (j in 1:length(unique(temp$mismatch))){
    temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
    lines(y=temptemp$meanI,x=temptemp$week,col=cols[j],pch=16) 
  }
  legend("topright",legend=c(unique(temp$mismatch)),col=c(cols),pch=16)
  temp<-data_wider_means[which(data_wider_means$country==unique(data_wider_means$country)[i]),]
  plot(temp$week,temp$meanflu,type="s",xlab="week",ylab="mean flu cases per capita")
  graphics.off()
}

###############################################################
