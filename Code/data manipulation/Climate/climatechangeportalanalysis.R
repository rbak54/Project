#from climate change portal
#deyla prob better
#195 countries as simple as possible                                  
require(dplyr)
require(tidyr)
data_climate<-read.csv("../../../Data/tas_1991_2016_AFG_ALB_DZA_AND_AGO_ATG_ARG_ARM_AUS_AUT_AZE_BHS_BHR_BGD_BRB_BLR_BEL_BLZ_BEN_BTN_BOL_BIH_BWA_BRA_BRN_BGR_BFA_BDI_KHM_CMR_CAN_CPV_CAF_TCD_CHL_CHN_COL_COM_COD_COG_CRI_CIV_HRV_CUB_CYP_CZE_DNK_DJI_DMA_DOM_ECU_EGY_SLV.csv")
data_climate<-tidyr::separate(data_climate,col=Statistics,into=c("not_sure","Month","Average"),sep="\\ ")

days<-  c(31,28,31,30,31,30,31,31,30,31,30,31)  
day_of_year<-matrix(c(c(1:12),days,rep(NA,12)),nrow=12,ncol=4,byrow = F)
for (i in 1:12){
  day_of_year[i,3]<-sum(day_of_year[0:(i-1),2])+ceiling((day_of_year[i,2]/2))
}

length(unique(data_climate$Country))*12
data_climate_means<-data_climate %>% group_by(Country,Month) %>% summarise(means=mean(Temperature))


data_climate_min_max<- data_climate_means %>% group_by(Country) %>% summarise(low=min(means),max_temp_month=Month[which.max(means)], high=max(means))  
data_climate_min_max$peak<-  day_of_year[match(data_climate_min_max$max_temp_month,month.abb),3]
data_climate_min_max<-as.data.frame(data_climate_min_max)

#data_climate_min_max$max_day<-data_climate_min_max$max_temp_month
#ata_climate_min_max$max_day<-match(data_climate_min_max$max_day)

write.csv(data_climate_min_max,"../../../Data/climatechangeportalresults")




