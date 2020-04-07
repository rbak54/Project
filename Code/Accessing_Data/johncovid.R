
datak<-read.csv("../Data/Covid/JOHN HOPKINS/COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv")
require(ggplot2)
datak[1,]
ggplot(data=datak,aes())
datak[order(datak["Country.Region"]),1:2]
datak[which(datak$Country.Region=="US"),]
#[order(datak["Country.Region"]),1:2]

data<-as.data.frame(matrix(ncol=ncol(datak)))
colnames(data)<-colnames(datak)
for(i in unique(datak$Country.Region)){
  sub<-datak[which(datak$Country.Region==i),]
  if (nrow(sub)>5){
   data<-rbind(data,sub) 
  }
}
