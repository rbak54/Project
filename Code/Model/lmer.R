require(lme4)
C<-read.csv("../../Results/correlation_dataframe_I.csv")
C<-read.csv("../../Results/correlation_dataframe_I_H.csv")

C<-read.csv("../../Results/correlation_dataframe_R0.csv")
C<-read.csv("../../Results/correlation_dataframe_R0_H.csv")

#seems to be same for both
c.lmer<-lmer(cors~(1|country)+mismatch+lat+long+pop+(1|combination),data=C)
summary(c.lmer)
c.lm<-c.lmer<-lm(cors~country+mismatch+lat+long+pop+combination,data=C)
summary(c.lm)
c.lm<-c.lmer<-lm(cors~mismatch+lat+long+pop+combination,data=C)
summary(c.lm)

plot(C$combination,C$cors)
plot(C$lat,C$cors)
#these graphs are really interesting!!
#absolute lat seems to determine mismatch
#>30 0 is best
#<30 chaos. something else effecting it. WHAT IS IT???
C$mismatch<-as.factor(C$mismatch)

ggplot(data=C,aes(lat,cors,col=mismatch))+geom_point()
ggplot(data=C,aes(abs(lat),cors,col=mismatch))+geom_point()
#from making graphs this doesn't seem to be mediated by anything in my dataset

C$mismatch<-as.numeric(as.character(C$mismatch))
low<-C[which(C$mismatch<0.5),]
ggplot(data=low,aes(lat,cors))+geom_point()
med<-C[which(C$mismatch==0.5),]
ggplot(data=med,aes(lat,cors))+geom_point()
high<-C[which(C$mismatch>0.5),]
ggplot(data=high,aes(lat,cors))+geom_point()




plot(C$long,C$cors)
C$mismatch<-as.factor(C$mismatch)
plot(C$mismatch,C$cors)

#latitudes -20-0-20 v all over the place correlations
#outside of that seem to be split into 3 distinct groups
#-those distinct groups are prob miusmatch, see 1 country eg 
#low, mid and high, but for different countries mismatch may mean different things

D<-C[which(C$country=="United Kingdom"),]
plot(D$combination,D$cors)
D$mismatch<-as.factor(D$mismatch)
ggplot(data=D,aes(combination,cors,col=mismatch))+geom_point()
plot(C$mismatch,C$cors)



