dataA<-read.csv("../Data/Covid/usa/covid-tracking-data/data/states_current.csv")
#only state level
head(dataA)
#(domr fo aboout website- poss more info there)
dataA<-read.csv("../Data/Covid/usa/covid-tracking-data/data/states_daily_4pm_et.csv")
#not sure whether this or above better


dataA<-read.csv("../Data/Covid/usa/covid-tracking-data/data/states_info.csv")
#info aboout website- poss more info there)


US<-datak[which(datak$Country.Region=="US"),]
unique(US$Province.State)

       