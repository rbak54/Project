source("../Model/modelling_functions.R")


Temperature<-seq(-10,40,length.out = 10000)
Contact_Rate_0<-cr_normal(Temperature,m = c(-10,30),stype=1,c = c(-10,40))
Contact_Rate_05<-cr_normal(Temperature,m = c(15,30),stype=1,c = c(-10,40))
Contact_Rate_1<-cr_normal(Temperature,m = c(40,30),stype=0,c = c(-10,40))
cols<-rainbow(n=3)
png("mismatch.png")
par(mar = c(5, 5, 3, 5))
plot(Temperature,-9.079*exp(0.085*Temperature),"s",xlab="Temperature",ylab="Growth Rate")
par(new=TRUE)
plot(Temperature,Contact_Rate_0,"l",xaxt="n",yaxt="n",ylab="",xlab="",col=cols[1],lty=1)
par(new=TRUE)

plot(Temperature,Contact_Rate_05,"l",xaxt="n",yaxt="n",ylab="",xlab="",col=cols[2],lty=1)
par(new=TRUE)

plot(Temperature,Contact_Rate_1,"l",xaxt="n",yaxt="n",ylab="",xlab="",col=cols[3],lty=1)

axis(side=4)
mtext("Contact Rate", side = 4, line = 3)
legend("bottom", c("Growth Rate", "Mismatch= 0","Mismatch=  0.5","Mismatch= 1"),
       col = c("black", cols), lty = c(1, 1,1,1))
#legend("bottom", c("Growth Rate", "Mismatch= 0"),
#       col = c("black", cols[1]), lty = c(1, 1))ssd
graphics.off()
cols
