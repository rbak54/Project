#italy
dataI<-read.csv("../Data/Covid/Italy/COVID-19/dati-regioni/dpc-covid19-ita-regioni-latest.csv")
dataI<-read.csv("../Data/Covid/Italy/COVID-19/dati-province/dpc-covid19-ita-province-latest.csv")
#more detail is province
dataI<-read.csv("../Data/Covid/Italy/COVID-19/dati-province/dpc-covid19-ita-province.csv")
##this has enough data for time series
head(dataI)
tail(dataI)
