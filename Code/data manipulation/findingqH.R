#Harper







data<-read.csv("../../Data/converted harper 1961.csv",header =F)
# in this dataset nil becomes 0, trace or n.d becomes NA,and the midpoint of the range is used
#where there is a range I have used the midpoint
data<-as.data.frame(data)
data

viability<-as.vector(as.matrix((data[4:10,])))
formatted_data<-as.data.frame(matrix(nrow=length(viability),ncol = 11))
formatted_data[,5]<-viability
formatted_data[,4]<-rep(c(1/60/60/24,1/12/24,1/2/24,1/24,4/24,6/24,23/24),length(viability)/7)
formatted_data[,6]<-rep(1:11,each=7)

formatted_data[,1]<-rep(as.vector(as.matrix((data[1,]))),each=7)
formatted_data[,2]<-rep(as.vector(as.matrix((data[2,]))),each=7)
formatted_data[,3]<-rep(as.vector(as.matrix((data[3,]))),each=7)
colnames(formatted_data)<-c("Temp","Hum","Rep","Time","Viability","Experiment","q1","v0","q2","aic1","aic2")
formatted_data

formatted_data$AH<-c(rep(0.00193,7),rep(0.00410,7),rep(0.00659,7),rep(0.00414,7),rep(0.00691,7),rep(0.01007,7),rep(0.01283,7),rep(0.01599,7),rep(0.00677,7),rep(0.01693,7),rep(0.02742,7))


require("minpack.lm")

times<-seq(0,1,length.out = 100)
par(mfrow=c(3,4))
for (i in (unique(formatted_data$Experiment))){
  data_subset<-formatted_data[which(formatted_data[,6]==i),]
  #this is model 2 
  model<-nlsLM(Viability~100*exp(q*Time),data=data_subset,start=list(q=-10))
  #plot(times,100*exp(times*summary(model)$parameters[1,1]),"s",ylim=c(0,100))
  #points(data_subset$Time,data_subset$Viability,col="blue",cex=1)
  formatted_data[which(formatted_data[,6]==i),11]<-AIC(model) 
  
  formatted_data[which(formatted_data[,6]==i),9]<-summary(model)$parameters[1,1]
  #this is model one
  model<-nlsLM(Viability~v0*exp(q*Time),data=data_subset,start=list(v0=100,q=-5),lower = c(v0=0  ,q=-300),upper = c(v0=100  ,q=000))
  #
  #without limits results look better? not sure which makes more biological sense
  
  formatted_data[which(formatted_data[,6]==i),7]<-summary(model)$parameters[2,1]
  formatted_data[which(formatted_data[,6]==i),8]<-summary(model)$parameters[1,1] 
  formatted_data[which(formatted_data[,6]==i),10]<-AIC(model)
  plot(times,summary(model)$parameters[1,1]*exp(times*summary(model)$parameters[2,1]),"s",ylim=c(0,100))
  points(data_subset$Time,data_subset$Viability,col="blue",cex=1)
}
formatted_data
shrunk_data<-formatted_data[seq(1,77,by=7),c(1,2,3,6,7,8,9,10,11,12)]
#https://planetcalc.com/2167/
shrunk_data
#first model best for all but 1 eg so just use first model?
par(mfrow=c(2,3))
plot(shrunk_data$Temp,shrunk_data$q1)
plot(shrunk_data$Temp,shrunk_data$q2)
plot(shrunk_data$Temp,shrunk_data$v0)
plot(shrunk_data$Hum,shrunk_data$q1)
plot(shrunk_data$Hum,shrunk_data$q2)
plot(shrunk_data$Hum,shrunk_data$v0)

plot(shrunk_data$AH,shrunk_data$q1)
plot(shrunk_data$AH,shrunk_data$q2)

par(mfrow=c(2,2))
plot(shrunk_data$Temp,shrunk_data$q1)

plot(shrunk_data$Hum,shrunk_data$q1)


shrunk_data_repeats<-as.data.frame(matrix(NA,ncol=ncol(shrunk_data),nrow = sum(shrunk_data$Rep)))
colnames(shrunk_data_repeats)=colnames(shrunk_data)
index<-0
for (i in 1:nrow(shrunk_data)) {
  shrunk_data_repeats[(index+1):(index+shrunk_data$Rep[i]) , ] <-shrunk_data[i,]
  index<-index+shrunk_data$Rep[i]
}

Hums<-seq(0,100,length.out = 1000)
Hum_Model<-nlsLM(q1~q0*exp(Hum*g),data=shrunk_data_repeats[,c(2,5)],start = list(g=0,q0=0))
g<-summary(Hum_Model)$coefficients[1,1]
q0<-summary(Hum_Model)$coefficients[2,1]

plot(Hums,q0*exp(g*Hums),"s")
points(shrunk_data_repeats$Hum,shrunk_data_repeats$q1)
shrunk_data_repeats

#result of all this is that equation for how b responds to temp is 
#(b0 * exp (g * Temps))
q0
#-21.98
g
#0.0209
###AH
plot(shrunk_data_repeats$AH,shrunk_data_repeats$q1)

AH<-seq(min(formatted_data$AH),max(formatted_data$AH),length.out = 1000)
AH_Model<-nlsLM(q1~q0*exp(AH*g),data=shrunk_data_repeats[,c(5,10)],start = list(g=0,q0=0))
g<-summary(AH_Model)$coefficients[1,1]
q0<-summary(AH_Model)$coefficients[2,1]

g 
#73.00948
q0
#-30.30238

#plot(shrunk_data_repeats$AH,shrunk_data_repeats$q1)
plot(AH,q0*exp(g*AH),"s")
points(shrunk_data_repeats$AH,shrunk_data_repeats$q1)
shrunk_data_repeats
