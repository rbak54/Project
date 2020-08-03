#data and functions 
source("modelling_functions.R")
  require("tidyverse")
require("ggplot2")
source("latincube.R")
source("data_sorting.R")
source("gazeteer.R")

data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ_POP.csv")
data_wider_summ<-read.csv("../../Data/data_wider_summ_POP.csv")
data_wider_means<-read.csv("../../Data/data_wider_means_POP.csv")
#integration  N changing
population_deyle<-read.csv("../../Data/population/API_SP.POP.TOTL_DS2_en_csv_v2_1217749/populations_sel.csv")
latlong=read.csv("../../Data/latlong/latlong_sel_short.csv")
pop<-read.csv("../../Data/population/API_SP.POP.TOTL_DS2_en_csv_v2_1217749/populations_sel.csv")




integration_general<-function(parms,sims,time){
  Y<-make_lhs(n=sims,parms = parms)
  write.csv(Y,"../../Data/latincube.csv")
  extra_cols<-9
  model_means=matrix(ncol =7)
  names(model_means)<-c("country","week","mismatch","meanI","meanR0","meantemp","combination")
  parms_temp<-parms
  parms_temp[["sigma"]]<-NA
  parms_temp[["h"]]<-NA
  parms_temp[["mu"]]<-NA
  parms_temp[["f"]]<-NA
  parms_temp[["N"]]<-NA

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
          
      if(parms[["climate_label"]]=="Temperature"){
          peak_contact_seq<-seq(data_wider_means_summ[location_index,"minT"],data_wider_means_summ[location_index,"maxT"],length.out=5)
          for (i in peak_contact_seq){
          #difference between temperature where contact rate is highest and lowest temperature in range (i.e virus does best survivasl or virus does best contact)
          #mismatch= 0 is when contact rate is highest at low temp 
          parms_temp[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakT"]*7,range_C=c(data_wider_means_summ[location_index,"minT"],data_wider_means_summ[location_index,"maxT"]),Max_Climate_cr=i)
           mismatch=(i-parms_temp[["Climate_Variables"]][["range_C"]][1])/(parms_temp[["Climate_Variables"]][["range_C"]][2]-parms_temp[["Climate_Variables"]][["range_C"]][1])
      #test time<-as.vector(read.csv("../../Results/time.csv"))[,2]
          mismatch<-round(mismatch,2)
          temp = ode(    y = start,    time = time,    func = SEIR_model,    parms = parms_temp)
           temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minT"],data_wider_means_summ[location_index,"maxT"],data_wider_means_summ[location_index,"peakT"],mismatch,NA,NA,NA,combination)),nrow=nrow(temp),ncol=extra_cols,byrow=T)
           colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","temperature","R0","week","combination")
           temp_extra[,"temperature"]<-Climate_Time_Function(time = time[1:nrow(temp)],min=parms_temp[["Climate_Variables"]][["range_C"]][1],max=parms_temp[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms_temp[["Climate_Variables"]][["time_at_peak"]] )
           temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(temp)),"temperature"],parms=parms_temp, Climate_Variables_Temp=parms_temp[["Climate_Variables"]], max_R0_Req=F)
           year<-ceiling((temp[,"time"])/365)
           day<-(temp[,"time"]-(year-1)*365)
            temp_extra[,"week"]<-ceiling(day/7)
            temp_extra[,"country"]<-location_index
           temp<-cbind(temp,temp_extra)
            temp<-as_tibble(temp)
           #model_means_temp<- temp %>% group_by(country,week,mismatch,combination) %>% summarise(meanI=mean(I),meanR0=mean(R0),meantemp=mean(temperature),.groups="drop")
           model_means_temp<- temp %>% group_by(country,week,mismatch,combination) %>% summarise(meanI=mean(I/(S+E+I+R)),meanR0=mean(R0),meantemp=mean(temperature),.groups="drop")
           
                 #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meantemp")
      #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meantemp")
      
            model_means<-bind_rows(model_means,model_means_temp)
      #names(model_means)<-c("country","week","mismatch","meanI","meanR0","meantemp")
      #names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meantemp")
      
      #png(paste0("../../Results/Plots/model_series/",data_wider_means_summ[location_index,"country"],gsub("\\.","",mismatch),".png"))
      #plottime(temp)
      #graphics.off()
          }
      }
    if(parms[["climate_label"]]=="RH"){
        peak_contact_seq<-seq(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],length.out=5)
        for (i in peak_contact_seq){
          #difference between relative_humidity where contact rate is highest and lowest relative_humidity in range (i.e virus does best survivasl or virus does best contact)
          #mismatch= 0 is when contact rate is highest at low RH 
          parms_temp[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakRH"]*7,range_C=c(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"]),Max_Climate_cr=i)
          mismatch=(i-parms_temp[["Climate_Variables"]][["range_C"]][1])/(parms_temp[["Climate_Variables"]][["range_C"]][2]-parms_temp[["Climate_Variables"]][["range_C"]][1])
          #test time<-as.vector(read.csv("../../Results/time.csv"))[,2]
          mismatch<-round(mismatch,2)
          RH = ode(    y = start,    time = time,    func = SEIR_model,    parms = parms_temp)
          temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],data_wider_means_summ[location_index,"peakRH"],mismatch,NA,NA,NA,combination)),nrow=nrow(RH),ncol=extra_cols,byrow=T)
          colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","relative_humidity","R0","week","combination")
          temp_extra[,"relative_humidity"]<-Climate_Time_Function(time = time[1:nrow(RH)],min=parms_temp[["Climate_Variables"]][["range_C"]][1],max=parms_temp[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms_temp[["Climate_Variables"]][["time_at_peak"]] )
          temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(RH)),"relative_humidity"],parms=parms_temp, Climate_Variables_Temp=parms_temp[["Climate_Variables"]], max_R0_Req=F)
          year<-ceiling((RH[,"time"])/365)
          day<-(RH[,"time"]-(year-1)*365)
          temp_extra[,"week"]<-ceiling(day/7)
          temp_extra[,"country"]<-location_index
          RH<-cbind(RH,temp_extra)
          RH<-as_tibble(RH)
          model_means_temp<- RH %>% group_by(country,week,mismatch,combination) %>% summarise(meanI=mean(I),meanR0=mean(R0),meanRH=mean(relative_humidity),.groups="drop")
          #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
          #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
          
          model_means<-bind_rows(model_means,model_means_temp)
          #names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH")
          #names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
          
          #png(paste0("../../Results/Plots/model_series/",data_wider_means_summ[location_index,"country"],gsub("\\.","",mismatch),".png"))
          #plottime(RH)
          #graphics.off()
        }
      }
  }
  }
model_means<-model_means[-1,]
return(model_means)
}


correlations<-function(mean,country){
  data_sub<-data_wider_means[which(data_wider_means$country==country),"meanflu"]
  return(unname(cor.test(mean[1:52],data_sub[1:52])[4][[1]]))
  
}

#addidng latitude and pop info
correlation_function<-function(model_means,parms){
  model_means$country<-data_wider_means_summ[model_means$country,"country"]
  correlation_df<- as_tibble(model_means)  %>% group_by(country,combination,mismatch) 
#  correlation_df<-na.omit(correlation_df$)
  if (parms[["climate_label"]]=="Temperature"){
    correlation_df <- correlation_df %>%  summarise(corsI=correlations(meanI,unique(country)),corsR=correlations(meanR0,unique(country)),maxs=max(meantemp),mins=min(meantemp),means=mean(meantemp),.groups="keep")
  }
  if (parms[["climate_label"]]=="RH"){
    correlation_df <- correlation_df %>%  summarise(corsI=correlations(meanI,unique(country)),corsR=correlations(meanR0,unique(country)),maxs=max(meanRH),mins=min(meanRH),means=mean(meanRH),.groups="keep")
  }
  #correlation_df$mismatch<-as.factor(correlation_df$mismatch)
  #for (i in unique(correlation_df$combination)){
  #  png(paste0("../../Results/Plots/comboplots/",i,".png"))
  #  toplot<-correlation_df[which(correlation_df$combination==i),]
  #  plot(toplot$mismatch,toplot$cors)
  #  graphics.off()
  #}
  #addidng latitude and pop info
  
  matrix_extra<-as.data.frame(matrix(nrow=nrow(correlation_df),ncol=3))
  colnames(matrix_extra)<-c("lat","long","pop")
  for (i in 1:nrow(matrix_extra)){
    matrix_extra[i,c(1:2)]<-latlong[which(latlong$V1==correlation_df$country[i]),c(3:4)]
    matrix_extra[i,3]<-pop[which(pop$V1==correlation_df$country[i]),"V2"]
    
  }
  correlation_df<-as.data.frame(correlation_df)
  correlation_df<-cbind(correlation_df,matrix_extra)
  #correlation_df$mismatch<-as.factor(correlation_df$mismatch)
  return(correlation_df)
}

std <- function(x) sd(x)/sqrt(length(x))

parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
                   N = NA, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="Temperature",
                   g=0.085,q0=-9.079,Climate_Variables=NA)
#parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
#           N = 1, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="RH",
 #          g=0.0209,q0=-21.98,Climate_Variables=NA)
#sims_range<-c(1, 5,10,20,40,80,160)
#sims_range<-c(5,10,20,40,80,160)

sims_range<-c(1)

for (sss in sims_range){

sims=sss
time = seq(1,365*10, by=1)

model_means<-integration_general(parms,sims,time)
write.csv(model_means,paste0("../../Results/fromfunction/",sims,parms[["climate_label"]],".csv"))
correlation_df<-correlation_function(model_means,parms)
write.csv(correlation_df,paste0("../../Results/fromfunction/cors/",sims,parms[["climate_label"]],"correlation_dataframe.csv"))

bests<-correlation_df %>% group_by(country) %>% summarise(best=mismatch[which.max(corsI)],.groups="keep")
correlation_df_means<-as_tibble(correlation_df) %>% group_by(mismatch) %>% mutate(means=mean(corsI),errors=std(corsI))

png(paste0("../../Results/Plots/boxandwhisker",sims,parms[["climate_label"]],".png"))
print(ggplot(data=correlation_df, aes(x= as.factor(mismatch),y= corsI)) +geom_boxplot())
graphics.off()

png(paste0("../../Results/Plots/erros",sims,parms[["climate_label"]],".png"))
print(ggplot(data=correlation_df_means, aes(x= as.factor(mismatch),y= means)) +geom_point()+
  geom_errorbar(aes(ymin=means+errors,ymax=means-errors)))
graphics.off()

correlation_df_means<-as_tibble(correlation_df) %>% group_by(mismatch,lat) %>% summarise(means=mean(corsI),errors=std(corsI),.groups="keep")
png(paste0("../../Results/Plots/lats",sims,parms[["climate_label"]],".png"))
print(ggplot(data=correlation_df_means, aes(x= abs(lat), col=as.factor(mismatch),y= means)) +geom_point()+
  geom_errorbar(aes(ymin=means+errors,ymax=means-errors)))
graphics.off()

correlation_df$Mismatch<-as.factor(correlation_df$mismatch)

bests<-correlation_df %>% group_by(country) %>% summarise(best=Mismatch[which.max(corsI)],.groups="keep")
correlation_df_means<-as_tibble(correlation_df) %>% group_by(Mismatch) %>% mutate(means=mean(corsI),errors=std(corsI))

pdf(paste0("../../Results/Plots/boxandwhisker",sims,parms[["climate_label"]],".pdf"))
print(ggplot(data=correlation_df, aes(x= Mismatch,y= corsI)) +geom_boxplot()  +theme_bw()+xlab("Mismatch") +ylab("Mean Correlation"))
graphics.off()

pdf(paste0("../../Results/Plots/erros",sims,parms[["climate_label"]],".pdf"))
print(ggplot(data=correlation_df_means, aes(x= Mismatch,y= means)) +geom_point()+theme_bw()+
        geom_errorbar(aes(ymin=means+errors,ymax=means-errors))+xlab("Mismatch") +ylab("Mean Correlation"))
graphics.off()

correlation_df_means<-as_tibble(correlation_df) %>% group_by(Mismatch,lat) %>% summarise(means=mean(corsI),errors=std(corsI),.groups="keep")
pdf(paste0("../../Results/Plots/lats",sims,parms[["climate_label"]],".pdf"))
print(ggplot(data=correlation_df_means, aes(x= abs(lat), col=Mismatch,y= means)) +scale_fill_discrete(name = "Dose")+geom_point()+theme_bw()+
        geom_errorbar(aes(ymin=means+errors,ymax=means-errors)) +xlab("Absolute Value of Latitude") +ylab("Mean Correlation") )
graphics.off()

}



#lots some plots here :(
#maxs-mins vs means
#also temp, min, max time at max

#also looked at interactions between these and mismatch.nothing much 

)


# 
# parms = list( mu = 2.06e-5,sigma = 0.68 ,p = 0.001, gamma =0.25,f=0.1,
#               N = 1, nu = 5.07e-5, h=0.25 / 24 ,epsilon= 0.05, d=4/24,Max_cr=29.97,climate_label="RH",
#               g=0.0209,q0=-21.98,Climate_Variables=NA)
# 
# 
# plot(correlation_df$lat,correlation_df$corsI)
# plot(correlation_df$mismatch,correlation_df$corsI)
# plot(log(correlation_df$pop),correlation_df$corsI)
# 
# 
# # think about whether best is useful 
# # #########################################
# # #plots for H
# # 
# # model_means<-read.csv("../../Results/model_results_population_combination_summary_H_10.csv")
# # 
# # model_means$country<-data_wider_means_summ[model_means$country,"country"]
# # correlation_df<- as_tibble(model_means)  %>% group_by(country,combination,mismatch) 
# # correlation_df<-na.omit(correlation_df)
# # #correlation_df <- correlation_df %>%  summarise(cors=correlations(meanR0,unique(country)),.groups="keep")
# #correlation_df <- correlation_df %>%  summarise(cors=correlations(meanR0,unique(country)),maxs=max(meantemp),mins=min(meantemp),means=mean(meantemp),.groups="keep")
# correlation_df <- correlation_df %>%  summarise(corsI=correlations(meanI,unique(country)),corsR=correlations(meanR0,unique(country)),maxs=max(meanRH),mins=min(meanRH),means=mean(meanRH),.groups="keep")
# #correlation_df$mismatch<-as.factor(correlation_df$mismatch)
# #for (i in unique(correlation_df$combination)){
# #  png(paste0("../../Results/Plots/comboplots/",i,".png"))
# #  toplot<-correlation_df[which(correlation_df$combination==i),]
# # plot(toplot$mismatch,toplot$cors)
# #  graphics.off()
# #}
# 
# matrix_extra<-as.data.frame(matrix(nrow=nrow(correlation_df),ncol=3))
# colnames(matrix_extra)<-c("lat","long","pop")
# for (i in 1:nrow(matrix_extra)){
#   matrix_extra[i,c(1:2)]<-latlong[which(latlong$V1==correlation_df$country[i]),c(3:4)]
#   matrix_extra[i,3]<-pop[which(pop$V1==correlation_df$country[i]),"V2"]
#   
# }
# 
# correlation_df<-as.data.frame(correlation_df)
# correlation_df<-cbind(correlation_df,matrix_extra)
# correlation_df$mismatch<-as.factor(correlation_df$mismatch)
# 
# 
# write.csv(correlation_df,"../../Results/correlation_dataframe_H.csv")
# 
# plot(correlation_df$lat,correlation_df$cors)
# plot(correlation_df$mismatch,correlation_df$cors)
# plot(log(correlation_df$pop),correlation_df$cors
# 
# ###############
# #humidity
# data_wider_means_summ<-data_wider_means_summ[which(!is.na(population_deyle$V2)),]
# population_deyle<-population_deyle[which(!is.na(population_deyle$V2)),]
# parms<-parms_main
# a=Sys.time()
# extra_cols<-8
# model_means=matrix(ncol =6)
# names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH")
# for (location_index in 1:nrow(data_wider_means_summ)){
#   parms[["N"]]=population_deyle[location_index,3]
#   start = c(S = (1-1e-4)*parms[["N"]],
#             E = 0.00*parms[["N"]], 
#             I =(1e-4)*parms[["N"]],
#             R = 0*parms[["N"]])
#   peak_contact_seq<-seq(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],length.out=5)
#   for (i in peak_contact_seq){
#     #difference between RH where contact rate is highest and lowest RH in range (i.e virus does best survivasl or virus does best contact)
#     #mismatch= 0 is when contact rate is highest at low RH 
#     parms[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakRH"]*7,range_C=c(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"]),Max_Climate_cr=i)
#     mismatch=(i-parms[["Climate_Variables"]][["range_C"]][1])/(parms[["Climate_Variables"]][["range_C"]][2]-parms[["Climate_Variables"]][["range_C"]][1])
#     #test time<-as.vector(read.csv("../../Results/time.csv"))[,2]
#     mismatch<-round(mismatch,2)
#     RH = ode(    y = start,    time = time,    func = SEIR_model,    parms = parms )
#     temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],data_wider_means_summ[location_index,"peakRH"],mismatch,NA,NA,NA)),nrow=nrow(RH),ncol=extra_cols,byrow=T)
#     colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","relative_humidity","R0","week")
#     temp_extra[,"relative_humidity"]<-Climate_Time_Function(time = time[1:nrow(RH)],min=parms[["Climate_Variables"]][["range_C"]][1],max=parms[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms[["Climate_Variables"]][["time_at_peak"]] )
#     temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(RH)),"relative_humidity"],parms=parms, Climate_Variables_Temp=parms[["Climate_Variables"]], max_R0_Req=F)
#     year<-ceiling((RH[,"time"])/365)
#     day<-(RH[,"time"]-(year-1)*365)
#     temp_extra[,"week"]<-ceiling(day/7)
#     temp_extra[,"country"]<-location_index
#     RH<-cbind(RH,temp_extra)
#     RH<-as_tibble(RH)
#     model_means_temp<- RH %>% group_by(country,week,mismatch) %>% summarise(meanI=mean(I),meanR0=mean(R0),meanRH=mean(relative_humidity),.groups="drop")
#     #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
#     #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
#     
#     model_means<-bind_rows(model_means,model_means_temp)
#     #names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH")
#     #names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
#     
#     #png(paste0("../../Results/Plots/model_series/",data_wider_means_summ[location_index,"country"],gsub("\\.","",mismatch),".png"))
#     #plottime(RH)
#     #graphics.off()
#   }
# }
# model_means<-model_means[-1,]
# write.csv(model_means,"../../Results/model_results_population_summary_H.csv")
# 
# 
# #integrations N changing and combinations -modelmeans
# 
# Y<-read.csv("../../Data/latincube.csv")
# 
# extra_cols<-9
# model_means=matrix(ncol =7)
# names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH","combination")
# parms_temp<-parms_main
# parms_temp[["sigma"]]<-NA
# parms_temp[["h"]]<-NA
# parms_temp[["mu"]]<-NA
# parms_temp[["f"]]<-NA
# parms_temp[["N"]]<-NA
# a<-Sys.time()
# 
# for (combination in 1:nrow(Y)){
#   parms_temp[["sigma"]]<-Y$sigma[combination]
#   parms_temp[["h"]]<-Y$h[combination]
#   parms_temp[["mu"]]<-Y$mu[combination]
#   parms_temp[["f"]]<-Y$f[combination]
#   for (location_index in 1:nrow(data_wider_means_summ)){
#     parms_temp[["N"]]=population_deyle[location_index,3]
#     start = c(S = (1-1e-4)*parms_temp[["N"]],
#               E = 0.00*parms_temp[["N"]], 
#               I =(1e-4)*parms_temp[["N"]],
#               R = 0*parms_temp[["N"]])
#     peak_contact_seq<-seq(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],length.out=5)
#     for (i in peak_contact_seq){
#       #difference between relative_humidity where contact rate is highest and lowest relative_humidity in range (i.e virus does best survivasl or virus does best contact)
#       #mismatch= 0 is when contact rate is highest at low RH 
#       parms_temp[["Climate_Variables"]]= list(time_at_peak=data_wider_means_summ[location_index,"peakRH"]*7,range_C=c(data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"]),Max_Climate_cr=i)
#       mismatch=(i-parms_temp[["Climate_Variables"]][["range_C"]][1])/(parms_temp[["Climate_Variables"]][["range_C"]][2]-parms_temp[["Climate_Variables"]][["range_C"]][1])
#       #test time<-as.vector(read.csv("../../Results/time.csv"))[,2]
#       mismatch<-round(mismatch,2)
#       RH = ode(    y = start,    time = time,    func = SEIR_model,    parms = parms_temp)
#       temp_extra<-matrix(rep(c(location_index,data_wider_means_summ[location_index,"minRH"],data_wider_means_summ[location_index,"maxRH"],data_wider_means_summ[location_index,"peakRH"],mismatch,NA,NA,NA,combination)),nrow=nrow(RH),ncol=extra_cols,byrow=T)
#       colnames(temp_extra)<-c("country","lower","upper","peak_week","mismatch","relative_humidity","R0","week","combination")
#       temp_extra[,"relative_humidity"]<-Climate_Time_Function(time = time[1:nrow(RH)],min=parms_temp[["Climate_Variables"]][["range_C"]][1],max=parms_temp[["Climate_Variables"]][["range_C"]][2],time_at_peak =parms_temp[["Climate_Variables"]][["time_at_peak"]] )
#       temp_extra[,"R0"]<-find_R0_function(Climate=temp_extra[c(1:nrow(RH)),"relative_humidity"],parms=parms_temp, Climate_Variables_Temp=parms_temp[["Climate_Variables"]], max_R0_Req=F)
#       year<-ceiling((RH[,"time"])/365)
#       day<-(RH[,"time"]-(year-1)*365)
#       temp_extra[,"week"]<-ceiling(day/7)
#       temp_extra[,"country"]<-location_index
#       RH<-cbind(RH,temp_extra)
#       RH<-as_tibble(RH)
#       model_means_temp<- RH %>% group_by(country,week,mismatch,combination) %>% summarise(meanI=mean(I),meanR0=mean(R0),meanRH=mean(relative_humidity),.groups="drop")
#       #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
#       #    names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
#       
#       model_means<-bind_rows(model_means,model_means_temp)
#       #names(model_means)<-c("country","week","mismatch","meanI","meanR0","meanRH")
#       #names(model_means_temp)<-c("country","week","mismatch","meanI","meanR0","meanRH")
#       
#       #png(paste0("../../Results/Plots/model_series/",data_wider_means_summ[location_index,"country"],gsub("\\.","",mismatch),".png"))
#       #plottime(RH)
#       #graphics.off()
#     }
#   }
# }
# model_means<-model_means[-1,]
# write.csv(model_means,"../../Results/model_results_population_combination_summary_H.csv")
# 
# Sys.time()-a


#BUT WHAT ABOUT COVID
#parms = list( mu = 2.06e-5,sigma = 1/5 ,p = 0.01, gamma =1/21,f=,
    #          N = NA, nu = 5.07e-5, h=0.25/24 ,epsilon= , d=4/24,Max_cr=29.97,climate_label="Temperature",
    #          g=,q0=,Climate_Variables=NA)

