dataUK<-read.csv("../../Sandbox/Covid/uk/covid-19-uk-data/data/covid-19-cases-uk.csv")
#PROB MOST USEFUL
dataUK<-read.csv("../Data/Covid/uk/covid-19-uk-data/data/covid-19-totals-uk.csv")
dataUK<-read.csv("../Data/Covid/uk/covid-19-uk-data/data/covid-19-indicators-uk.csv")
#ABOVE HAS DEATHS TESTS AND CASES FOR DIFFERENT UK COUNTRIES

dataloc<-dataUK[which(dataUK$Area=="Coventry"),]
dataloc<-dataUK[which(dataUK$Area=="Warwickshire"),]
dataloc<-dataUK[which(dataUK$Area=="Birmingham"),]
dataloc<-dataUK[which(dataUK$Area=="Blackburn with Darwen"),]
dataloc<-dataUK[which(dataUK$Area=="West Berkshire"),]
dataloc<-dataUK[which(dataUK$Area=="Nottingham"),]
dataloc<-dataUK[which(dataUK$Area=="Greater Glasgow and Clyde"),]
head(dataloc)
dataloc[,6]<-c(1:nrow(dataloc))
dataloc[,7]<-NA
for (i in 2:nrow(dataloc)){
  dataloc[i,7]<-dataloc[i,5]-dataloc[i-1,5]
}
plot(dataloc$V6,dataloc$TotalCases)
plot(dataloc$V6,log(dataloc$TotalCases))
plot(dataloc$V6,dataloc$V7)
plot(dataloc$V6,log(dataloc$TotalCases))

sort(unique(dataUK$Area))
