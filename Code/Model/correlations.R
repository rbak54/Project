require("tidyverse")
require("ggplot2")
data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ_POP.csv")
data_wider_summ<-read.csv("../../Data/data_wider_summ_POP.csv")
data_wider_means<-read.csv("../../Data/data_wider_means_POP.csv")
###############################################################
#correlation between model and data for each mismatch for each location
#model_means<-read.csv("../../Results/model_results_summary.csv")
#for just one set of parameters

model_means<-read.csv("../../Results/model_results_population_summary.csv")


model_means[,"country"]<-data_wider_means_summ[model_means[,"country"],"country"]

model_means_summ<-model_means %>% group_by(country,mismatch) %>% summarise()
correlation_df<-as.data.frame(matrix(nrow=nrow(model_means_summ),ncol=3))
colnames(correlation_df)=c("Country","Mismatch","Correlation")
correlation_df[,c(1:2)]<-model_means_summ
for (i in 1:nrow(correlation_df)){
  #png(paste0("../../Results/Plots/correlation_plots/",i,".png"))
  country<-correlation_df[i,"Country"]
  mismatch<-correlation_df[i,"Mismatch"]
  data_sub<-data_wider_means[which(data_wider_means$country==country),]
  model_sub<-model_means[which(model_means$country==country),]
  model_sub<-model_sub[which(model_sub$mismatch==mismatch),]
  model_sub<-model_sub[c(1:nrow(data_sub)),]
  #plot(model_sub$meanI,data_sub$meanflu)
  #graphics.off()
  correlation_df[i,"Correlation"]<-cor.test(model_sub$meanR0,data_sub$meanflu)[4]
}

cor2<-correlation_df%>% group_by(Mismatch) %>% summarise(mean(Correlation))
plot(cor2$Mismatch,cor2$`mean(Correlation)`) #median better?
correlation_df$Mismatch<-as.factor(correlation_df$Mismatch)
plot(correlation_df$Mismatch,correlation_df$Correlation)
#
bests<-correlation_df %>% group_by(Country) %>% summarise(best=Mismatch[which.max(Correlation)])
ggplot(data=bests,aes(best))+geom_histogram(stat="count")
#how to quantify this
cor.test(correlation_df$Correlation,as.numeric(as.character(correlation_df$Mismatch)))
#cor.test(as.numeric(as.character(medians$Mismatch)),medians$medians)
#some sort of linear or non linear mpde;
write.csv(correlation_df,"../../Results/correlation_dataframe_1.csv")


model_means<-read.csv("../../Results/model_results_population_summary_H.csv")
model_means[,"country"]<-data_wider_means_summ[model_means[,"country"],"country"]

model_means_summ<-model_means %>% group_by(country,mismatch) %>% summarise()
correlation_df<-as.data.frame(matrix(nrow=nrow(model_means_summ),ncol=3))
colnames(correlation_df)=c("Country","Mismatch","Correlation")
correlation_df[,c(1:2)]<-model_means_summ
for (i in 1:nrow(correlation_df)){
  #png(paste0("../../Results/Plots/correlation_plots/",i,".png"))
  country<-correlation_df[i,"Country"]
  mismatch<-correlation_df[i,"Mismatch"]
  data_sub<-data_wider_means[which(data_wider_means$country==country),]
  model_sub<-model_means[which(model_means$country==country),]
  model_sub<-model_sub[which(model_sub$mismatch==mismatch),]
  model_sub<-model_sub[c(1:nrow(data_sub)),]
  #plot(model_sub$meanI,data_sub$meanflu)
  #graphics.off()
  correlation_df[i,"Correlation"]<-cor.test(model_sub$meanR0,data_sub$meanflu)[4]
}
write.csv(correlation_df,"../../Results/correlation_dataframe_H_1.csv")




#################################################################
#correlation using model means-varying combinations
# #ignore pop
# data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ.csv")
# data_wider_summ<-read.csv("../../Data/data_wider_summ.csv")
# data_wider_means<-read.csv("../../Data/data_wider_means.csv")

# model_means<-read.csv("../../Results/model_results_combinations_summary.csv")
# #OR
# 
# model_means<-read.csv("../../Results/model_results_combinations_summary_H.csv")
# #OR
# 
# 
# correlations<-function(mean,country){
#   data_sub<-data_wider_means[which(data_wider_means$country==country),"meanflu"]
#   return(unname(cor.test(mean[1:52],data_sub[1:52])[4][[1]]))
#   
# }
# 
# 
# model_means$country<-data_wider_means_summ[model_means$country,"country"]
# correlation_df<- as_tibble(model_means)  %>% group_by(country,combination,mismatch) 
# correlation_df<-na.omit(correlation_df)
# #correlation_df <- correlation_df %>%  summarise(cors=correlations(meanI,unique(country)),.groups="keep")
# correlation_df <- correlation_df %>%  summarise(cors=correlations(meanI,unique(country)),range=max(meantemp)-min(meantemp),.groups="keep")
# correlation_df <- correlation_df %>%  summarise(cors=correlations(meanI,unique(country)),range=max(meanRH)-min(meanRH),.groups="keep")
# correlation_df$mismatch<-as.factor(correlation_df$mismatch)
# for (i in unique(correlation_df$combination)){
#   png(paste0("../../Results/Plots/comboplots/",i,".png"))
# toplot<-correlation_df[which(correlation_df$combination==i),]
#   plot(toplot$mismatch,toplot$cors)
#   graphics.off()
# }
# #someting big effect- seeing if mu
# 
# plot(correlation_df$mismatch,correlation_df$cors)
# #
# correlation_df
# bests<-correlation_df %>% group_by(country) %>% summarise(best=mismatch[which.max(cors)])
# ggplot(data=bests,aes(best))+geom_histogram(stat="count")
# #hoq to quantify this
# #cor.test(as.numeric(as.character(medians$Mismatch)),medians$medians)
# #some sort of linear or non linear mpde;
# #obtaining latlong
# latlong=read.csv("../../Data/latlong/latlong_sel.csv")
# matrix_latlong<-as.data.frame(matrix(nrow=nrow(correlation_df),ncol=2))
# colnames(matrix_latlong)<-c("lat","long")
# for (i in 1:nrow(matrix_latlong)){
#   matrix_latlong[i,]<-latlong[which(latlong$V1==correlation_df$country[i]),c(3:4)]
# }
# 
# correlation_df<-as.data.frame(correlation_df)
# correlation_df<-cbind(correlation_df,matrix_latlong)
# correlation_df$mismatch<-as.factor(correlation_df$mismatch)
# plot(correlation_df$lat,correlation_df$cors)
# plot(correlation_df$mismatch,correlation_df$cors)
# #########################




#THIS IS THE ONE
#plots for population 

data_wider_means_summ<-read.csv("../../Data/data_wider_means_summ_POP.csv")
data_wider_summ<-read.csv("../../Data/data_wider_summ_POP.csv")
data_wider_means<-read.csv("../../Data/data_wider_means_POP.csv")

correlations<-function(mean,country){
  data_sub<-data_wider_means[which(data_wider_means$country==country),"meanflu"]
  return(unname(cor.test(mean[1:52],data_sub[1:52])[4][[1]]))
  
}

#addidng latitude and pop info
latlong=read.csv("../../Data/latlong/latlong_sel_short.csv")
pop<-read.csv("../../Data/population/API_SP.POP.TOTL_DS2_en_csv_v2_1217749/populations_sel.csv")


model_means<-read.csv("../../Results/model_results_population_summary.csv")

model_means<-read.csv("../../Results/model_results_population_combination_summary.csv")

model_means$country<-data_wider_means_summ[model_means$country,"country"]
correlation_df<- as_tibble(model_means)  %>% group_by(country,combination,mismatch) 
correlation_df<-na.omit(correlation_df)
correlation_df <- correlation_df %>%  summarise(corsI=correlations(meanI,unique(country)),corsR=correlations(meanR0,unique(country)),maxs=max(meantemp),mins=min(meantemp),means=mean(meantemp),.groups="keep")

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
correlation_df$mismatch<-as.factor(correlation_df$mismatch)

write.csv(correlation_df,"../../Results/correlation_dataframe.csv")

plot(correlation_df$lat,correlation_df$cors)
plot(correlation_df$mismatch,correlation_df$cors)
plot(log(correlation_df$pop),correlation_df$cors)



#########################################
#plots for H

model_means<-read.csv("../../Results/model_results_population_combination_summary_H_10.csv")

model_means$country<-data_wider_means_summ[model_means$country,"country"]
correlation_df<- as_tibble(model_means)  %>% group_by(country,combination,mismatch) 
correlation_df<-na.omit(correlation_df)
#correlation_df <- correlation_df %>%  summarise(cors=correlations(meanR0,unique(country)),.groups="keep")
#correlation_df <- correlation_df %>%  summarise(cors=correlations(meanR0,unique(country)),maxs=max(meantemp),mins=min(meantemp),means=mean(meantemp),.groups="keep")
correlation_df <- correlation_df %>%  summarise(corsI=correlations(meanI,unique(country)),corsR=correlations(meanR0,unique(country)),maxs=max(meanRH),mins=min(meanRH),means=mean(meanRH),.groups="keep")
#correlation_df$mismatch<-as.factor(correlation_df$mismatch)
#for (i in unique(correlation_df$combination)){
#  png(paste0("../../Results/Plots/comboplots/",i,".png"))
#  toplot<-correlation_df[which(correlation_df$combination==i),]
 # plot(toplot$mismatch,toplot$cors)
#  graphics.off()
#}

matrix_extra<-as.data.frame(matrix(nrow=nrow(correlation_df),ncol=3))
colnames(matrix_extra)<-c("lat","long","pop")
for (i in 1:nrow(matrix_extra)){
  matrix_extra[i,c(1:2)]<-latlong[which(latlong$V1==correlation_df$country[i]),c(3:4)]
  matrix_extra[i,3]<-pop[which(pop$V1==correlation_df$country[i]),"V2"]
  
}

correlation_df<-as.data.frame(correlation_df)
correlation_df<-cbind(correlation_df,matrix_extra)
correlation_df$mismatch<-as.factor(correlation_df$mismatch)


write.csv(correlation_df,"../../Results/correlation_dataframe_H.csv")

plot(correlation_df$lat,correlation_df$cors)
plot(correlation_df$mismatch,correlation_df$cors)
plot(log(correlation_df$pop),correlation_df$cors)

# model_means<-read.csv("../../Results/model_results_population_combination_summary_H_10.csv")
# 
# 
# 
# 
# model_means$country<-data_wider_means_summ[model_means$country,"country"]
# correlation_df<- as_tibble(model_means)  %>% group_by(country,combination,mismatch) 
# correlation_df<-na.omit(correlation_df)
# #correlation_df <- correlation_df %>%  summarise(cors=correlations(meanR0,unique(country)),.groups="keep")
# correlation_df <- correlation_df %>%  summarise(cors=correlations(meanR0,unique(country)),maxs=max(meantemp),mins=min(meantemp),means=mean(meantemp),.groups="keep")
# correlation_df <- correlation_df %>%  summarise(cors=correlations(meanR0,unique(country)),maxs=max(meanRH),mins=min(meanRH),means=mean(meanRH),.groups="keep")
# correlation_df$mismatch<-as.factor(correlation_df$mismatch)
# for (i in unique(correlation_df$combination)){
#   png(paste0("../../Results/Plots/comboplots/",i,".png"))
#   toplot<-correlation_df[which(correlation_df$combination==i),]
#   plot(toplot$mismatch,toplot$cors)
#   graphics.off()
# }
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
# plot(correlation_df$lat,correlation_df$cors)
# plot(correlation_df$mismatch,correlation_df$cors)
# plot(log(correlation_df$pop),correlation_df$cors)
# 
# write.csv(correlation_df,"../../Results/correlation_dataframe_R0.csv")
# write.csv(correlation_df,"../../Results/correlation_dataframe_H.csv")
# 
# 
# write.csv(correlation_df,"../../Results/correlation_dataframe_R0_H.csv")


