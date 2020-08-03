#Harper







  data<-read.csv("../../Data/converted harper 1961.csv",header =F)
# in this dataset nil becomes 0, trace or n.d becomes NA,and the midpoint of the range is used
#where there is a range I have used the midpoint
data<-as.data.frame(data)


viability<-as.vector(as.matrix((data[4:10,])))
formatted_data<-as.data.frame(matrix(nrow=length(viability),ncol = 11))
formatted_data[,5]<-viability
formatted_data[,4]<-rep(c(1/60/60/24,1/12/24,1/2/24,1/24,4/24,6/24,23/24),length(viability)/7)
formatted_data[,6]<-rep(1:11,each=7)

formatted_data[,1]<-rep(as.vector(as.matrix((data[1,]))),each=7)
formatted_data[,2]<-rep(as.vector(as.matrix((data[2,]))),each=7)
formatted_data[,3]<-rep(as.vector(as.matrix((data[3,]))),each=7)
colnames(formatted_data)<-c("Temp","Hum","Rep","Time","Viability","Experiment","q1","v0","q2","aic1","aic2")



require("minpack.lm")

times<-seq(0,1,length.out = 10000)
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
 plot(times,summary(model)$parameters[1,1]*exp(times*summary(model)$parameters[2,1]),type="s",ylim=c(0,100),ylab="Viability")
 points(data_subset$Time,data_subset$Viability,col="blue",cex=1)
 }

seq(1,77,by=7)
shrunk_data<-formatted_data[seq(1,77,by=7),c(1,2,3,6,7,8,9,10,11)]


shrunk_data
#first model best for all but 1 eg so just use first model?
par(mfrow=c(2,3))
plot(shrunk_data$Temp,shrunk_data$q1)
plot(shrunk_data$Temp,shrunk_data$q2)
plot(shrunk_data$Temp,shrunk_data$v0)
plot(shrunk_data$Hum,shrunk_data$q1)
plot(shrunk_data$Hum,shrunk_data$q2)
plot(shrunk_data$Hum,shrunk_data$v0)



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

Temps<-seq(-10,40,length.out = 1000)
Temp_Model<-nlsLM(q1~q0*exp(Temp*g),data=shrunk_data_repeats[,c(1,5)],start = list(g=0,q0=0))
g<-summary(Temp_Model)$coefficients[1,1]
q0<-summary(Temp_Model)$coefficients[2,1]
  
plot(Temps,q0*exp(g*Temps),"s")
points(shrunk_data_repeats$Temp,shrunk_data_repeats$q1)
shrunk_data_repeats

#result of all this is that equation for how b responds to temp is 
  #(b0 * exp (g * Temps))
q0
#0.572
g
#0.164
  
# or if apply limits
#b0=9.079
#g=0.085  

#dublineau
#watch glass ; suypplementary
data<-read.csv("../../Data/dublineaub.csv",header =T)
#if <1.67 have replaced with then NA
data<-as.data.frame(data)


viability<-as.vector(as.matrix((data[,4:9])))
formatted_data<-as.data.frame(matrix(nrow=length(viability),ncol = 11))
formatted_data[,5]<-viability
formatted_data[,4]<-rep(c(0.5/24,0,1,2,3,8),length(viability)/6)
formatted_data[,6]<-rep(1:12,each=6)

formatted_data[,1]<-rep(as.vector(as.matrix((data[,1]))),each=6)
formatted_data[,2]<-rep(as.vector(as.matrix((data[,2]))),each=6)
formatted_data[,3]<-rep(as.vector(as.matrix((data[,3]))),each=6)
colnames(formatted_data)<-c("strain","Temp","dtheoretical","Time","tcid","Experiment","b1","v0","b2","aic1","aic2")


for (i in 1:nrow(formatted_data)){
  
  if (formatted_data[i,4]==0){
    
    if (formatted_data[i,2]==4){
      formatted_data[i,4]<-((17+(40/60))/24)
    }
    if (formatted_data[i,2]==25){
      formatted_data[i,4]<-((7+(45/60))/24)
    }
    if (formatted_data[i,2]==35){
      formatted_data[i,4]<-((5+(15/60))/24)
    }
  }
}

formatted_data$tcid<-10^(formatted_data$tcid)
formatted_data$dtheoretical<-10^(formatted_data$dtheoretical)

 require("minpack.lm")

times<-seq(0,10,length.out = 100)
par(mfrow=c(3,4))
for (i in (unique(formatted_data$Experiment))){
  data_subset<-formatted_data[which(formatted_data[,6]==i),]
  theoretical<-data_subset[1,3]
  model<-nlsLM(tcid~theoretical*exp(-b*Time),data=data_subset,start=list(b=00))
  plot(times,theoretical*exp(-times*summary(model)$parameters[1,1]),"s")
  points(data_subset$Time,data_subset$tcid,col="blue",cex=1)
  formatted_data[which(formatted_data[,6]==i),10]<-AIC(model) 
  
  formatted_data[which(formatted_data[,6]==i),7]<-summary(model)$parameters[1,1]
  
  model<-nlsLM(tcid~t0*exp(-b*Time),data=data_subset,start=list(t0=theoretical,b=0))
  #,lower=c(0,0),upper = c(100,1000)) 
  #without limits results look better? not sure which makes more biological sense
  
  formatted_data[which(formatted_data[,6]==i),9]<-summary(model)$parameters[2,1]
  formatted_data[which(formatted_data[,6]==i),8]<-summary(model)$parameters[1,1] 
  formatted_data[which(formatted_data[,6]==i),11]<-AIC(model)
  #plot(times,summary(model)$parameters[1,1]*exp(-times*summary(model)$parameters[2,1]),"s")
  #points(data_subset$Time,data_subset$tcid,col="blue",cex=1)
}

seq(1,77,by=7)
shrunk_data<-formatted_data[seq(1,72,by=6),c(1,2,3,4,6,7,8,9,10,11)]

shrunk_data
#first modeformatted_datal best for all but 1 eg so just use first model?
par(mfrow=c(2,3))
plot(shrunk_data$Temp,shrunk_data$b1)
plot(shrunk_data$Temp,shrunk_data$b2)
plot(shrunk_data$Temp,shrunk_data$v0)


par(mfrow=c(2,2))
plot(shrunk_data$Temp,shrunk_data$b1)

#wrong?


Temps<-seq(-10,40,length.out = 1000)
g<-summary(nlsLM(b2~b0*exp(Temp*g),data=shrunk_data[,c(2,8)],start = list(g=0,b0=0)))$coefficients[1,1]
b0<-summary(nlsLM(b2~b0*exp(Temp*g),data=shrunk_data[,c(2,8)],start = list(g=0,b0=0)))$coefficients[2,1]

plot(Temps,b0*exp(g*Temps),"s")
points(shrunk_data$Temp,shrunk_data$b2)



