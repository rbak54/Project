#data and functions 


source("modelling_functions.R")
require("tidyverse")
require("ggplot2")
#source("latincube.R")
Y<-read.csv("../../Data/latincube.csv")
source("data_sorting.R")
source("gazeteer.R")
data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ.csv")
data_wider_summ<-read.csv("../../Data/data_wider_summ.csv")

#################################################################
#integration initialise
extra_cols<-7
#out=matrix(ncol=5+extra_cols)  

time = seq(1,365*20, by=2)
parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = 1, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="RH",
              g=0.0209,q0=-21.98,Climate_Variables=NA)
start = c(S = (1-1e-5)*parms[["N"]],
          E = 0.00*parms[["N"]], 
          I =(1e-5)*parms[["N"]],
          R = 0*parms[["N"]])
################################################################

###################################################################
#integration loop
# a=Sys.time()
# 
# for (location_index in 1:nrow(data_wider_means_summ)){
#   peak_contact_seq<-seq(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],length.out=3)
#   for (i in peak_contact_seq){
#     #difference between relative_humidity where contact rate is highest and lowest relative_humidity in range (i.e virus does best survivasl or virus does best contact)
#     #mismatch= 0 is when contact rate is highest at low temp 
#     parms[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakRH"]*7,range_C=c(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"]),Max_Climate_cr=i)
#     mismatch=abs((i-parms[["Climate_Variables"]][["range_C"]][1])/(parms[["Climate_Variables"]][["range_C"]][2]-parms[["Climate_Variables"]][["range_C"]][1]))
#     #test time<-as.vector(read.csv("../../Results/time.csv"))[,2]
#     mismatch<-round(mismatch,2)
#     temp = ode(    y = start,    time = time,    func = SEIR_model,    parms = parms )
#     temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],data_wider_means_summ[location_index,"peakRH"],mismatch,NA,NA)),nrow=nrow(temp),ncol=extra_cols,byrow=T)
#     colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","relative_humidity","R0")
#     temp_extra[,"relative_humidity"]<-Climate_Time_Function(time = time[1:nrow(temp)],min=parms[["Climate_Variables"]][["range_C"]][1],max=parms[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms[["Climate_Variables"]][["time_at_peak"]] )
#     temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(temp)),"relative_humidity"],parms=parms, Climate_Variables_Temp=parms[["Climate_Variables"]], max_R0_Req=F)
#     png(paste0("../../Results/Plots/model_series/",data_wider_means_summ[location_index,"country"],gsub("\\.","",mismatch),".png"))
#     plottime(temp)
#     graphics.off()
#     temp<-cbind(temp,temp_extra)
#     out<-rbind(out,temp)
#   }
# }
# out<-out[-1,]
# Sys.time()-a
# write.csv(out,"../../Results/model_results_H.csv")
############################################################################################################################
#integration summary
a=Sys.time()
extra_cols<-8
model_means=matrix(ncol =6)
names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH")
for (location_index in 1:nrow(data_wider_means_summ)){
  peak_contact_seq<-seq(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],length.out=5)
  for (i in peak_contact_seq){
    #difference between relative_humidity where contact rate is highest and lowest relative_humidity in range (i.e virus does best survivasl or virus does best contact)
    #mismatch= 0 is when contact rate is highest at low temp 
    parms[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakRH"]*7,range_C=c(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"]),Max_Climate_cr=i)
    mismatch=abs((i-parms[["Climate_Variables"]][["range_C"]][1])/(parms[["Climate_Variables"]][["range_C"]][2]-parms[["Climate_Variables"]][["range_C"]][1]))
    #test time<-as.vector(read.csv("../../Results/time.csv"))[,2]
    mismatch<-round(mismatch,2)
    temp = ode(    y = start,    time = time,    func = SEIR_model,    parms = parms )
    #png(paste0("../../Results/Plots/model_series/",data_wider_means_summ[location_index,"country"],gsub("\\.","",mismatch),".png"))
    #plottime(temp)
    #graphics.off()
    temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],data_wider_means_summ[location_index,"peakRH"],mismatch,NA,NA,NA)),nrow=nrow(temp),ncol=extra_cols,byrow=T)
    colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","relative_humidity","R0","week")
    temp_extra[,"relative_humidity"]<-Climate_Time_Function(time = time[1:nrow(temp)],min=parms[["Climate_Variables"]][["range_C"]][1],max=parms[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms[["Climate_Variables"]][["time_at_peak"]] )
    temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(temp)),"relative_humidity"],parms=parms, Climate_Variables_Temp=parms[["Climate_Variables"]], max_R0_Req=F)
    year<-ceiling((temp[,"time"])/365)
    day<-(temp[,"time"]-(year-1)*365)
    temp_extra[,"week"]<-ceiling(day/7)
    temp_extra[,"country"]<-location_index
    temp<-cbind(temp,temp_extra)
    temp<-as_tibble(temp)
    model_means_temp<- temp %>% group_by(country,week,mismatch) %>% summarise(meanI=mean(I),meanR0=mean(R0),meanRH=mean(relative_humidity),.groups="drop")
    #    names(model_means_RH)<-c("country","week","mismatch","meanI","meanR0","meanRH")
    #    names(model_means_RH)<-c("country","week","mismatch","meanI","meanR0","meanRH")
    
    model_means<-bind_rows(model_means,model_means_temp)
    #names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH")
    #names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
    

  }
}
model_means<-model_means[-1,]
write.csv(model_means,"../../Results/model_results_summary_H.csv")
#Sys.time()-a
#this way is quicker for more data
#####################################################################################################################
#####################################################################################################################
#integration combinations

parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = 1, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="RH",
              g=0.0209,q0=-21.98,Climate_Variables=NA)
parms_temp[["sigma"]]<-NA
parms_temp[["h"]]<-NA
parms_temp[["mu"]]<-NA
parms_temp[["f"]]<-NA
start = c(S = (1-1e-5)*parms_temp[["N"]],
          E = 0.00*parms_temp[["N"]], 
          I =(1e-5)*parms_temp[["N"]],
          R = 0*parms_temp[["N"]])
Y<-read.csv("../../Data/latincube.csv")
extra_cols<-9
model_means=matrix(ncol =7)
names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH","combination")
a<-Sys.time()
for (combination in 1:nrow(Y)){
  parms_temp[["sigma"]]<-Y$sigma[combination]
  parms_temp[["h"]]<-Y$h[combination]
  parms_temp[["mu"]]<-Y$mu[combination]
  parms_temp[["f"]]<-Y$f[combination]
  for (location_index in 1:nrow(data_wider_means_summ)){
    peak_contact_seq<-seq(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],length.out=5)
    for (i in peak_contact_seq){
      #difference between relative_humidity where contact rate is highest and lowest relative_humidity in range (i.e virus does best survivasl or virus does best contact)
      #mismatch= 0 is when contact rate is highest at low temp 
      parms_temp[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakRH"]*7,range_C=c(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"]),Max_Climate_cr=i)
      mismatch=abs((i-parms_temp[["Climate_Variables"]][["range_C"]][1])/(parms_temp[["Climate_Variables"]][["range_C"]][2]-parms_temp[["Climate_Variables"]][["range_C"]][1]))
      #test time<-as.vector(read.csv("../../Results/time.csv"))[,2]
      mismatch<-round(mismatch,2)
      temp = ode(    y = start,    time = time,    func = SEIR_model,    parms = parms_temp ,)
      temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],data_wider_means_summ[location_index,"peakRH"],mismatch,NA,NA,NA,combination)),nrow=nrow(temp),ncol=extra_cols,byrow=T)
      colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","relative_humidity","R0","week","combination")
      temp_extra[,"relative_humidity"]<-Climate_Time_Function(time = time[1:nrow(temp)],min=parms_temp[["Climate_Variables"]][["range_C"]][1],max=parms_temp[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms_temp[["Climate_Variables"]][["time_at_peak"]] )
      temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(temp)),"relative_humidity"],parms=parms_temp, Climate_Variables_Temp=parms_temp[["Climate_Variables"]], max_R0_Req=F)
      year<-ceiling((temp[,"time"])/365)
      day<-(temp[,"time"]-(year-1)*365)
      temp_extra[,"week"]<-ceiling(day/7)
      temp_extra[,"country"]<-location_index
      temp<-cbind(temp,temp_extra)
      temp<-as_tibble(temp)
      model_means_temp<- temp %>% group_by(country,week,mismatch,combination) %>% summarise(meanI=mean(I),meanR0=mean(R0),meanRH=mean(relative_humidity),.groups="drop")
      model_means<-bind_rows(model_means,model_means_temp)
    }
  }
}

model_means<-model_means[-1,]
Sys.time()-a
write.csv(model_means,"../../Results/model_results_combinations_summary_H.csv")
#15 minutes for 10 comninations






#
#################################################################
#integration  N changing
population_deyle<-read.csv("../../Data/population/API_SP.POP.TOTL_DS2_en_csv_v2_1217749/populations_sel.csv")
data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ_POP.csv")
time = seq(1,365*20, by=1)
parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
              N = NA, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="relative_humidity",
              g=0.085,q0=-9.079,Climate_Variables=NA)
extra_cols<-7
out=matrix(ncol=5+extra_cols)  

data_wider_means_summ<-data_wider_means_summ[which(!is.na(population_deyle$V2)),]
population_deyle<-population_deyle[which(!is.na(population_deyle$V2)),]
# a=Sys.time()
# #integration loop
# for (location_index in 1:nrow(data_wider_means_summ)){
#   #currently 2019 but could change
#   parms[["N"]]=population_deyle[location_index,3]
#   start = c(S = (1-1e-4)*parms[["N"]],
#             E = 0.00*parms[["N"]], 
#             I =(1e-4)*parms[["N"]],
#             R = 0*parms[["N"]])
#   peak_contact_seq<-seq(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],length.out=3)
#   for (i in peak_contact_seq){
#     #difference between relative_humidity where contact rate is highest and lowest relative_humidity in range (i.e virus does best survivasl or virus does best contact)
#     #mismatch= 0 is when contact rate is highest at low RH 
#     parms[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakRH"]*7,range_C=c(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"]),Max_Climate_cr=i)
#     mismatch=abs((i-parms[["Climate_Variables"]][["range_C"]][1])/(parms[["Climate_Variables"]][["range_C"]][2]-parms[["Climate_Variables"]][["range_C"]][1]))
#     #test time<-as.vector(read.csv("../../Results/time.csv"))[,2]
#     mismatch<-round(mismatch,2)
#     temp = ode(    y = start,    time = time,    func = SEIR_model,    parms = parms )
#     temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],data_wider_means_summ[location_index,"peakRH"],mismatch,NA,NA)),nrow=nrow(temp),ncol=7,byrow=T)
#     colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","relative_humidity","R0")
#     temp_extra[,"relative_humidity"]<-Climate_Time_Function(time = time[1:nrow(temp)],min=parms[["Climate_Variables"]][["range_C"]][1],max=parms[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms[["Climate_Variables"]][["time_at_peak"]] )
#     temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(temp)),"relative_humidity"],parms=parms, Climate_Variables_Temp=parms[["Climate_Variables"]], max_R0_Req=F)
#     png(paste0("../../Results/Plots/model_series/",data_wider_means_summ[location_index,"country"],gsub("\\.","",mismatch),".png"))
#     plottime(temp)
#     graphics.off()
#     temp<-cbind(temp,temp_extra)
#     out<-rbind(out,temp)
#   }
# }
# out<-out[-1,]
# Sys.time()-a
# write.csv(out,"../../Results/model_results_population.csv")
a=Sys.time()
extra_cols<-8
model_means=matrix(ncol =6)
names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH")
for (location_index in 1:nrow(data_wider_means_summ)){
  parms[["N"]]=population_deyle[location_index,3]
  start = c(S = (1-1e-4)*parms[["N"]],
            E = 0.00*parms[["N"]], 
            I =(1e-4)*parms[["N"]],
            R = 0*parms[["N"]])
  peak_contact_seq<-seq(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],length.out=5)
  for (i in peak_contact_seq){
    #difference between relative_humidity where contact rate is highest and lowest relative_humidity in range (i.e virus does best survivasl or virus does best contact)
    #mismatch= 0 is when contact rate is highest at low RH 
    parms[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakRH"]*7,range_C=c(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"]),Max_Climate_cr=i)
    mismatch=abs((i-parms[["Climate_Variables"]][["range_C"]][1])/(parms[["Climate_Variables"]][["range_C"]][2]-parms[["Climate_Variables"]][["range_C"]][1]))
    #test time<-as.vector(read.csv("../../Results/time.csv"))[,2]
    mismatch<-round(mismatch,2)
    temp = ode(    y = start,    time = time,    func = SEIR_model,    parms = parms )
    temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],data_wider_means_summ[location_index,"peakRH"],mismatch,NA,NA,NA)),nrow=nrow(temp),ncol=extra_cols,byrow=T)
    colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","relative_humidity","R0","week")
    temp_extra[,"relative_humidity"]<-Climate_Time_Function(time = time[1:nrow(temp)],min=parms[["Climate_Variables"]][["range_C"]][1],max=parms[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms[["Climate_Variables"]][["time_at_peak"]] )
    temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(temp)),"relative_humidity"],parms=parms, Climate_Variables_Temp=parms[["Climate_Variables"]], max_R0_Req=F)
    year<-ceiling((temp[,"time"])/365)
    day<-(temp[,"time"]-(year-1)*365)
    temp_extra[,"week"]<-ceiling(day/7)
    temp_extra[,"country"]<-location_index
    temp<-cbind(temp,temp_extra)
    temp<-as_tibble(temp)
    model_means_temp<- temp %>% group_by(country,week,mismatch) %>% summarise(meanI=mean(I),meanR0=mean(R0),meanRH=mean(relative_humidity),.groups="drop")
    #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
    #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
    
    model_means<-bind_rows(model_means,model_means_temp)
    #names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH")
    #names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
    
    #png(paste0("../../Results/Plots/model_series/",data_wider_means_summ[location_index,"country"],gsub("\\.","",mismatch),".png"))
    #plottime(temp)
    #graphics.off()
  }
}
model_means<-model_means[-1,]
write.csv(model_means,"../../Results/model_results_population_summary_H.csv")
Sys.time()-a
#this way is quicker for more data
########################################################################################################
#integrations N changing and combinations -modelmeans

population_deyle<-read.csv("../../Data/population/API_SP.POP.TOTL_DS2_en_csv_v2_1217749/populations_sel.csv")
data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ_POP.csv")
Y<-read.csv("../../Data/latincube.csv")

time = seq(1,365*20, by=1)
extra_cols<-9
model_means=matrix(ncol =7)
names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH","combination")

parms_temp = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
                   N = 1, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="relative_humidity",
                   g=0.085,q0=-9.079,Climate_Variables=NA)
parms_temp[["sigma"]]<-NA
parms_temp[["h"]]<-NA
parms_temp[["mu"]]<-NA
parms_temp[["f"]]<-NA
parms_temp[["N"]]<-NA
a<-Sys.time()



for (combination in 1:nrow(Y)){
  parms_temp[["sigma"]]<-Y$sigma[combination]
  parms_temp[["h"]]<-Y$h[combination]
  parms_temp[["mu"]]<-Y$mu[combination]
  parms_temp[["f"]]<-Y$f[combination]
  for (location_index in 1:nrow(data_wider_means_summ)){
    parms_temp[["N"]]=population_deyle[location_index,3]
    start = c(S = (1-1e-4)*parms_temp[["N"]],
              E = 0.00*parms_temp[["N"]], 
              I =(1e-4)*parms_temp[["N"]],
              R = 0*parms_temp[["N"]])
    peak_contact_seq<-seq(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],length.out=5)
    for (i in peak_contact_seq){
      #difference between relative_humidity where contact rate is highest and lowest relative_humidity in range (i.e virus does best survivasl or virus does best contact)
      #mismatch= 0 is when contact rate is highest at low RH 
      
      parms_temp[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakRH"]*7,range_C=c(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"]),Max_Climate_cr=i)
      mismatch=abs((i-parms_temp[["Climate_Variables"]][["range_C"]][1])/(parms_temp[["Climate_Variables"]][["range_C"]][2]-parms_temp[["Climate_Variables"]][["range_C"]][1]))
      #test time<-as.vector(read.csv("../../Results/time.csv"))[,2]
      mismatch<-round(mismatch,2)
      temp = ode(    y = start,    time = time,    func = SEIR_model,    parms = parms_temp)
      temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],data_wider_means_summ[location_index,"peakRH"],mismatch,NA,NA,NA,combination)),nrow=nrow(temp),ncol=extra_cols,byrow=T)
      colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","relative_humidity","R0","week","combination")
      temp_extra[,"relative_humidity"]<-Climate_Time_Function(time = time[1:nrow(temp)],min=parms_temp[["Climate_Variables"]][["range_C"]][1],max=parms_temp[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms_temp[["Climate_Variables"]][["time_at_peak"]] )
      temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(temp)),"relative_humidity"],parms=parms_temp, Climate_Variables_Temp=parms_temp[["Climate_Variables"]], max_R0_Req=F)
      year<-ceiling((temp[,"time"])/365)
      day<-(temp[,"time"]-(year-1)*365)
      temp_extra[,"week"]<-ceiling(day/7)
      temp_extra[,"country"]<-location_index
      temp<-cbind(temp,temp_extra)
      temp<-as_tibble(temp)
      model_means_temp<- temp %>% group_by(country,week,mismatch,combination) %>% summarise(meanI=mean(I),meanR0=mean(R0),meanRH=mean(relative_humidity),.groups="drop")
      #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
      #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
      
      model_means<-bind_rows(model_means,model_means_temp)
      #names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH")
      #names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
      
      #png(paste0("../../Results/Plots/model_series/",data_wider_means_summ[location_index,"country"],gsub("\\.","",mismatch),".png"))
      #plottime(temp)
      #graphics.off()
    }
  }
}
model_means<-model_means[-1,]





write.csv(model_means,"../../Results/model_results_population_combination_summary_H.csv")



Sys.time()-a
