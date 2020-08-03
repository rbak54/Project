source("modelling_functions.R")
require("tidyverse")
require("ggplot2")
##################################
#additional functions
non_neg<-function(x){
  for (i in 1:length(x)){
    if (x[i]<0){
      x[i]=0
    }
  }
  return(x)
}

find_temp<-function(time){
  ### assuming time is dot to dot
  upperpos<-Position(function(x) x >time, country_subset$day_tot)
  days<-country_subset$day_tot[c(upperpos-1,upperpos)]
  temps<-country_subset$T[c(upperpos-1,upperpos)]#print(day)
  gradient<-(temps[2]-temps[1])/(days[2]-days[1])
  return(temps[1]+gradient*(time-days[1]))
  #print("CLimate=")
}


SEIR_model_mod <- function(t=time, values=start, parms) {
  #function where c changes with the climate and humidity at each timepoint
  #initial values
  S <- values[1]
  E <- values[2]
  I <- values[3]
  R <- values[4]
  N <- values[5]
  #parameters
  
  for (parameter in names(parms)){
    assign(parameter,parms[[parameter]])
  }
  
  for (variable in names(Climate_Variables)){
    assign(variable,Climate_Variables[[variable]])
  }
  
  Climate<-find_temp(t)
  cr_value=cr_climate(Climate =  Climate,Max_Coordinates_cr = c(Max_Climate_cr,Max_cr), range_C=range_C)
  q_value= q_climate(q0=q0,g=g, Climate=Climate)
  beta_value = beta(c_r = cr_value,q =q_value ,d = d,h =h ,epsilon =epsilon )    

  #seir model
  dS = nu * N - beta_value * I * S / N - mu * S + f * R
  dE = beta_value * S * I / N - (sigma + mu) * E
  dI = sigma * E - (mu + gamma) * I * (1/(1-p))
  dR = gamma * I - mu * R - f * R
  dN = nu * N - (p/(1-p)) *(mu + gamma) * I - mu * (S + E + I + R)
  
  list(c(dS, dE, dI, dR, dN))
}

################################################

###############################################################
#open data
data_wider<-read.csv("../../Data/deyle_edited_wide.csv")
################################################################

################################################################
#data tidying and reformatting 
#remove missing flu and temp values
data_wider<-data_wider[which(data_wider$T!="NAN"),]
data_wider<-data_wider[which(data_wider$flu!="NAN"),]
data_wider$day_of_year<-(data_wider$week-1)*7+3.5
data_wider$day_tot<-(data_wider$year-min(data_wider$year)-1)*365 +data_wider$day_of_year
#convert to celcius
data_wider$T<-(data_wider$T-32)*(5/9)
#find weekly mean flu and temp over the years
data_wider_means<-data_wider %>% group_by(country,week) %>% summarise(meanflu=mean(flu),meantemp=mean(T))
#convert week to day
data_wider_means$day<-(data_wider_means$week-1)*7+3.5

# find summary for each country- week and values for maxima and mimima
data_wider_means_summ<-data_wider_means %>% group_by(country) %>% summarise(minflu=min(meanflu),maxflu=max(meanflu),peakflu=week[which.max(meanflu)],troughflu=week[which.min(meanflu)],minT=min(meantemp),maxT=max(meantemp),peakT=week[which.max(meantemp)], troughT=week[which.min(meantemp)])
# convert to data frame
data_wider_means_summ<-as.data.frame(data_wider_means_summ)
#summary for each year for each country
data_wider_summ<-data_wider %>% group_by(country, year) %>% summarise(minflu=min(flu),maxflu=max(flu),peakflu=week[which.max(flu)],troughflu=week[which.min(flu)],minT=min(T),maxT=max(T),peakT=week[which.max(T)], troughT=week[which.min(T)])
#################################################################
#if above edits break -see other



#################################################################
#integration


extra_cols<-7
out=matrix(ncol=6+extra_cols)  
parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = 1, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="Temperature",
              g=0.085,q0=-9.079,Climate_Variables=NA)
#parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.005, gamma =0.25,f=0.1,
   #           N = 1, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="Temperature",
  #            g=0.085,q0=-9.079,Climate_Variables=NA)
#parms = list( mu = 0,sigma = 0.68 ,p = 0.00, gamma =0.25,f=0.1,
 #             N = 1, nu =0, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="Temperature",
  #            g=0.085,q0=-9.079,Climate_Variables=NA)
start = c(S = (1-1e-3)*parms[["N"]],
          E = 0.00*parms[["N"]],
          I =(1e-3)*parms[["N"]],
          R = 0*parms[["N"]],
          N  = parms[["N"]])
#remember to change h and epsilon because sensitivity analysis


for (location_index in 1:nrow(data_wider_means_summ)){
  location<-data_wider_means_summ[location_index,"country"]
  print(location)
  peak_contact_seq<-seq(data_wider_means_summ[location_index,"minT"],data_wider_means_summ[location_index,"maxT"],length.out=3)
  country_subset<-data_wider[which(data_wider$country==location),]
  country_subset$day_tot<-country_subset$day_of_year+(country_subset$year-min(country_subset$year,na.rm = T))*365
  Climate_Variables_Min=list(Max_Climate_cr=NA, range_C=c(data_wider_means_summ[location_index,"minT"],data_wider_means_summ[location_index,"maxT"]))
  for (i in peak_contact_seq){
    #difference between temperature where contact rate is highest and lowest temperature in range (i.e virus does best survivasl or virus does best contact)
    #mismatch= 0 is when contact rate is highest at low temp 
    #parms[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakT"]*7,range_C=c(data_wider_means_summ[location_index,"minT"],data_wider_means_summ[location_index,"maxT"]),Max_Climate_cr=i)
    Climate_Variables_Min[["Max_Climate_cr"]]=i
    parms[["Climate_Variables"]]=Climate_Variables_Min
    #time<-country_subset$day_tot
    time<-seq(min(country_subset$day_tot),max(country_subset$day_tot),by=1)
    mismatch=abs((i-data_wider_means_summ[location_index,"minT"])/(data_wider_means_summ[location_index,"maxT"]-data_wider_means_summ[location_index,"minT"]))
    temp = ode(   y=start,  parms =  parms ,  times=time, func=SEIR_model_mod)
    #temp<-temp[1:(nrow(temp)-1),]
    temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minT"],data_wider_means_summ[location_index,"maxT"],data_wider_means_summ[location_index,"peakT"],mismatch,NA,NA)),nrow=nrow(temp),ncol=extra_cols,byrow=T)
    colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","temperature","R0")
    #added bc breaking. temp fix!!!
   #  for (j in 1:nrow(temp_extra)){
    #    temp_extra[j,"temperature"]<-find_temp(temp[j,"time"])
     # }
    temp_extra[,"temperature"]<-sapply(temp[,"time"],find_temp)
    temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[,"temperature"],parms=parms, Climate_Variables_Temp=Climate_Variables_Min, max_R0_Req=F)
    png(paste0("../../Results/Plots/model_series_temp/",data_wider_means_summ[location_index,"country"],gsub("\\.","",mismatch),".png"))
    plottime(temp)
    graphics.off()
    temp<-cbind(temp,temp_extra)
   # out<-rbind(temp_mean,out)
    out<-rbind(out,temp)
  }
}

out<-out[-1,]
write.csv(out,"../../Results/model_results_temporal.csv")
########################################################################

########################################################################
###############################################################
#open data
data_wider<-read.csv("../../Data/deyle_edited_wide.csv")
################################################################

################################################################
#data tidying and reformatting 
#remove missing flu and temp values
data_wider<-data_wider[which(data_wider$T!="NAN"),]
data_wider<-data_wider[which(data_wider$flu!="NAN"),]
#convert to celcius
data_wider$T<-(data_wider$T-32)*(5/9)
#find weekly mean flu and temp over the years
data_wider_means<-data_wider %>% group_by(country,week) %>% summarise(meanflu=mean(flu),meantemp=mean(T))
#convert week to day
data_wider_means$day<-(data_wider_means$week-1)*7+3.5

# find summary for each country- week and values for maxima and mimima
data_wider_means_summ<-data_wider_means %>% group_by(country) %>% summarise(minflu=min(meanflu),maxflu=max(meanflu),peakflu=week[which.max(meanflu)],troughflu=week[which.min(meanflu)],minT=min(meantemp),maxT=max(meantemp),peakT=week[which.max(meantemp)], troughT=week[which.min(meantemp)])

# convert to data frame
data_wider_means_summ<-as.data.frame(data_wider_means_summ)
#summary for each year for each country
data_wider_summ<-data_wider %>% group_by(country, year) %>% summarise(minflu=min(flu),maxflu=max(flu),peakflu=week[which.max(flu)],troughflu=week[which.min(flu)],minT=min(T),maxT=max(T),peakT=week[which.max(T)], troughT=week[which.min(T)])


#backup if above breaks-2
#################################################################
###############################################################

#reload integration results and add further info. make more summary tibbles.
out<-read.csv("../../Results/model_results_temporal.csv")
out<-as.data.frame(out)
out[,"country"]<-data_wider_means_summ[out[,"country"],"country"]
out[,"year"]<-ceiling((out[,"time"])/365)
out[,"day"]<-(out[,"time"]-(out[,"year"]-1)*365)
out[,"week"]<-ceiling((out[,"day"])/7)

model_wide<-out
model_wide$mismatch<-round(model_wide$mismatch,3)
model_means<-model_wide %>% group_by(country,week,mismatch) %>% summarise(meanI=mean(I,na.rm = T),meanR0=mean(R0,na.rm = T),meantemp=mean(temperature,na.rm = T))
model_means_summ<-model_means %>% group_by(country,mismatch) %>% summarise(minI=min(meanI),maxI=max(meanI),peakI=week[which.max(meanI)],troughI=week[which.min(meanI)],minR0=min(meanR0),maxR0=max(meanR0),peakR0=week[which.max(meanR0)],troughR0=week[which.min(meanR0)],minT=min(meantemp),maxT=max(meantemp),peakT=week[which.max(meantemp)], troughT=week[which.min(meantemp)])
model_summ_year<- model_wide %>% group_by(country,mismatch,year) %>% summarise(minI=min(I),maxI=max(I),peakI=week[which.max(I)],troughI=week[which.min(I)],minR0=min(R0),maxR0=max(R0),peakR0=week[which.max(R0)],troughR0=week[which.min(R0)],minT=min(temperature),maxT=max(temperature),peakT=week[which.max(temperature)], troughT=week[which.min(temperature)])
write.csv(model_means,"../../Results/model_means_temporal.csv")
write.csv(model_means_summ,"../../Results/model_means_summ_temporal.csv")
write.csv(model_summ_year,"../../Results/model_summ_year_temporal.csv")

#reload summary results
model_means<-read.csv("../../Results/model_means_temporal.csv")
model_means_summ<-read.csv("../../Results/model_means_summ_temporal.csv")
model_summ_year<-read.csv("../../Results/model_summ_year_temporal.csv")


#####################################################################

#####################################################################
#plots C+P!!!-see 3

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


#plots for individual countries
############################################################################
### flu incidence at different temperatures using weekly summaries. - R0, I (different mismatch) and data

# below shows huge effect of mismatch. but also that a low mismatch pattern is most common
for (i in unique(model_means$country)){
  png(paste0("../../Results/Plots/country_week_means_T/",i,".png"))
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
  png(paste0("../../Results/Plots/country_week_means_summ_T/",i,".png"))
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
  png(paste0("../../Results/Plots/country_week_T/",location,".png"))
  plot(temp$T,temp$flu)
  graphics.off()
}
#########################################################

#not all that helpful. summary stat graphs for each year for each country 
for (location in unique(data_wider_summ$country)){
  temp<-data_wider_summ[which(data_wider_summ$country==location),]
  png(paste0("../../Results/Plots/country_year_summ_T/",location,".png"))
  par(mfrow=c(3,2))
  plot(temp$minT,temp$maxflu)
  plot(temp$maxT,temp$maxflu)
  plot(temp$troughT,temp$peakflu)
  plot(temp$peakT,temp$peakflu)
  plot(temp$troughT,temp$maxflu)
  
  graphics.off()
}
##########################################################


##########################################################
#R0 over time for model ,  meanI over time for model, flu over time for data
for (i in unique(model_means$country)){
  png(paste0("../../Results/Plots/R0_time_T/",i,".png"))
  par(mfrow=c(3,1))
  temp<-model_means[which(model_means$country==i),]
  plot(x=temp$week,y=temp$meanR0,"n")
  for (j in 1:length(unique(temp$mismatch))){
    temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
    points(y=temptemp$meanR0,x=temptemp$week,col=cols[j],pch=16) 
  }
  legend("topright",legend=c(unique(temp$mismatch)),col=c(cols),pch=16)
  plot(x=temp$week,y=temp$meanI,"n")
  for (j in 1:length(unique(temp$mismatch))){
    temptemp<-temp[which(temp$mismatch==(unique(temp$mismatch))[j]),]
    points(y=temptemp$meanI,x=temptemp$week,col=cols[j],pch=16) 
  }
  legend("topright",legend=c(unique(temp$mismatch)),col=c(cols),pch=16)
  temp<-data_wider_means[which(data_wider_means$country==i),]
  plot(temp$week,temp$meanflu)
  graphics.off()
}
###############################################################

###############################################################
#correlation between model and data for each mismatch for each location

correlation_df<-as.data.frame(matrix(nrow=nrow(model_means_summ),ncol=3))
colnames(correlation_df)=c("Country","Mismatch","Correlation")
correlation_df[,"Country"]<-model_means_summ[,"country"]
correlation_df[,"Mismatch"]<-model_means_summ[,"mismatch"]
for (i in 1:nrow(correlation_df)){
  png(paste0("../../Results/Plots/correlation_plots_T/",i,".png"))
  
  country<-correlation_df[i,"Country"]
  mismatch<-correlation_df[i,"Mismatch"]
  
  data_sub<-data_wider_means[which(data_wider_means$country==country),]
  model_sub<-model_means[which(model_means$country==country),]
  model_sub<-model_sub[which(model_sub$mismatch==mismatch),]
  model_sub<-model_sub[c(1:nrow(data_sub)),]
  plot(model_sub$meanI,data_sub$meanflu)
  graphics.off()
  correlation_df[i,"Correlation"]<-cor.test(model_sub$meanI,data_sub$meanflu)[4]
}

cor2<-correlation_df%>% group_by(Mismatch) %>% summarise(mean(Correlation))
plot(cor2$Mismatch,cor2$`mean(Correlation)`)
correlation_df$Mismatch<-as.factor(correlation_df$Mismatch)
plot(correlation_df$Mismatch,correlation_df$Correlation)

#################################################################

