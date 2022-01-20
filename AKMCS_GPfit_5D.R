rm(list = ls())
graphics.off()
library(GPfit)
library(lhs)
library(fields) # Make grid plots
setwd("D:/Penn state/paper2/Sensitivity")
Reliability<-function (X) {
  X[ ,1]*X[ ,2]+2*X[ ,1]*X[ ,3]-X[ ,4]*X[ ,5]+0.5
}

d = 5
set.seed(1)

#direct Sobol
n <- 25000
X1 <- data.frame(matrix(runif(d*n), nrow = n))
y1 <- Reliability(X1)
plot(density(y1),lwd=2,xlab="y",main="")
length(which(y1<0))/25000

X2 <- data.frame(matrix(runif(d*n), nrow = n))
First_order<-matrix(NA,nrow=1000,ncol=5)
First_order_up<-matrix(NA,nrow=1000,ncol=5)
First_order_low<-matrix(NA,nrow=1000,ncol=5)
Second_order_dum<-rep(NA,1000)
Second_order_dumup<-rep(NA,1000)
Second_order_dumlow<-rep(NA,1000)

Second_order<-matrix(NA,nrow=1000,ncol=10)
Second_order_up<-matrix(NA,nrow=1000,ncol=10)
Second_order_low<-matrix(NA,nrow=1000,ncol=10)

for (i in 1:1000){
  N = n/1000*i
  S<-sensitivity::sobol(model = Reliability, X1 = X1[1:N, ], X2 = X2[1:N, ],
                        order=2, nboot = 100)
  First_order[i, ] <- S$S$original[1:5]
  First_order_up[i, ] <- S$S$`max. c.i.`[1:5]
  First_order_low[i, ] <- S$S$`min. c.i.`[1:5]
  Second_order[i, ] <- S$S$original[6:15]
  Second_order_up[i, ] <- S$S$`max. c.i.`[6:15]
  Second_order_low[i, ] <- S$S$`min. c.i.`[6:15]
  print(i)
}

save(First_order,file="Data/5D/First_order")
save(First_order_up,file="Data/5D/First_order_up")
save(First_order_low,file="Data/5D/First_order_low")

save(Second_order,file="Data/5D/Second_order")
save(Second_order_up,file="Data/5D/Second_order_up")
save(Second_order_low,file="Data/5D/Second_order_low")


library(randomcoloR)
k <- 5
Col <- distinctColorPalette(k)
pie(rep(1, k), col=Col)
Col<-c("Red","Green","Blue","Orange","Cyan")

n=25000*4
X=seq(n/1000,n,by=n/1000)
par(mar=c(4,5.1,1.6,2.1))
plot(X,First_order[ ,1]*100,type="l",lwd=2,xlab="Sample size",ylab="Percentage of explained variance (%)",
     cex.lab=1.5,cex.axis=1.5,ylim=c(-30,100),xlim=c(100,n),log="x",col="red",xaxt="n")
lines(X,First_order[ ,1]*100,lwd=2,col=Col[1])
lines(X,First_order[ ,2]*100,lwd=2,col=Col[2])
lines(X,First_order[ ,3]*100,lwd=2,col=Col[3])
lines(X,First_order[ ,4]*100,lwd=2,col=Col[4])
lines(X,First_order[ ,5]*100,lwd=2,col=Col[5])
lines(X,First_order_up[ ,1]*100,lty=2,lwd=1,col=Col[1])
lines(X,First_order_up[ ,2]*100,lty=2,lwd=1,col=Col[2])
lines(X,First_order_up[ ,3]*100,lty=2,lwd=1,col=Col[3])
lines(X,First_order_up[ ,4]*100,lty=2,lwd=1,col=Col[4])
lines(X,First_order_up[ ,5]*100,lty=2,lwd=1,col=Col[5])
lines(X,First_order_low[ ,1]*100,lty=2,lwd=1,col=Col[1])
lines(X,First_order_low[ ,2]*100,lty=2,lwd=1,col=Col[2])
lines(X,First_order_low[ ,3]*100,lty=2,lwd=1,col=Col[3])
lines(X,First_order_low[ ,4]*100,lty=2,lwd=1,col=Col[4])
lines(X,First_order_low[ ,5]*100,lty=2,lwd=1,col=Col[5])
axis(1, at = c(10,100,1000,10000,100000),cex.axis=1.5)
legend(40000,110,lwd=c(2,2,2,2,2),lty=c(1,1,1,1,1),col = Col,
       legend = c("X1","X2","X3","X4","X5"),cex=1,bty="n")


Col <- rainbow(10)
n=25000*4
X=seq(n/1000,n,by=n/1000)
par(mar=c(4,5.1,1.6,2.1))
plot(X,Second_order[ ,1]*100,type="l",lwd=2,xlab="Sample size",ylab="Percentage of explained variance (%)",
     cex.lab=1.5,cex.axis=1.5,ylim=c(-10,40),xlim=c(100,n),log="x",col="red",xaxt="n")
lines(X,Second_order[ ,1]*100,lwd=2,col=Col[1])
lines(X,Second_order[ ,2]*100,lwd=2,col=Col[2])
lines(X,Second_order[ ,3]*100,lwd=2,col=Col[3])
lines(X,Second_order[ ,4]*100,lwd=2,col=Col[4])
lines(X,Second_order[ ,5]*100,lwd=2,col=Col[5])
lines(X,Second_order[ ,6]*100,lwd=2,col=Col[6])
lines(X,Second_order[ ,7]*100,lwd=2,col=Col[7])
lines(X,Second_order[ ,8]*100,lwd=2,col=Col[8])
lines(X,Second_order[ ,9]*100,lwd=2,col=Col[9])
lines(X,Second_order[ ,10]*100,lwd=2,col=Col[10])
axis(1, at = c(10,100,1000,10000,100000),cex.axis=1.5)
legend(7000,40,lwd=c(2,2,2,2,2),lty=c(1,1,1,1,1),col = Col[1:5],
       legend = c("X1&X2","X1&X3","X1&X4","X1&X5","X2&X3"),cex=1,bty="n")
legend(30000,40,lwd=c(2,2,2,2,2),lty=c(1,1,1,1,1),col = Col[6:10],
       legend = c("X2&X4","X2&X5","X3&X4","X3&X5","X4&X5"),cex=1,bty="n")


#Kriging
Kriging <- function (X){
  a<-predict(GPmodel,X)
  a$Y_hat
}
n = 100
x = maximinLHS(n, d)
y = Reliability(x)
GPmodel = GP_fit(x, y)
levels<-seq(0,1,length=100)
Grid<-expand.grid(levels,levels)
Data<-data.frame(X1=Grid$Var1,X2=Grid$Var2,X3=0.5,X4=0.5,X5=0.5)
Y_fitted<-Kriging(Data)
Y_true<-Reliability(Data)
colfunc <- colorRampPalette(c("white", "black"))
cols <- colfunc(101)
par(mar=c(4,5.1,1.6,2.1))
image.plot(levels,levels, matrix(Y_true,ncol=100),xaxt="n",yaxt="n", xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
axis(1, at = c(0,1),cex.axis=1.5)
axis(2, at = c(0,1),cex.axis=1.5)
par(mar=c(4,5.1,1.6,2.1))
image.plot(levels,levels, matrix(Y_fitted-Y_true,ncol=100),xaxt="n",yaxt="n", xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
axis(1, at = c(0,1),cex.axis=1.5)
axis(2, at = c(0,1),cex.axis=1.5)

Kriging_First_order<-matrix(NA,nrow=25,ncol=5)
Kriging_First_order_up<-matrix(NA,nrow=25,ncol=5)
Kriging_First_order_low<-matrix(NA,nrow=25,ncol=5)
Kriging_Second_order<-matrix(NA,nrow=25,ncol=10)
Kriging_Second_order_up<-matrix(NA,nrow=25,ncol=10)
Kriging_Second_order_low<-matrix(NA,nrow=25,ncol=10)

load("Data/5D/Kriging_First_order")
load("Data/5D/Kriging_First_order_up")
load("Data/5D/Kriging_First_order_low")
load("Data/5D/Kriging_Second_order")
load("Data/5D/Kriging_Second_order_up")
load("Data/5D/Kriging_Second_order_low")

for (i in 1:25){
  n = 10*i
  x = maximinLHS(n, d)
  y = Reliability(x)
  GPmodel = GP_fit(x, y)
  
  S<-sensitivity::sobol(model = Kriging, X1 = X1, X2 = X2,
                        order=2, nboot = 100,conf=0.9)
  Kriging_First_order[i, ] <- S$S$original[1:5]
  Kriging_First_order_up[i, ] <- S$S$`max. c.i.`[1:5]
  Kriging_First_order_low[i, ] <- S$S$`min. c.i.`[1:5]
  Kriging_Second_order[i, ] <- S$S$original[6:15]
  Kriging_Second_order_up[i, ] <- S$S$`max. c.i.`[6:15]
  Kriging_Second_order_low[i, ] <- S$S$`min. c.i.`[6:15]
  save(Kriging_First_order,file="Data/5D/Kriging_First_order")
  save(Kriging_First_order_up,file="Data/5D/Kriging_First_order_up")
  save(Kriging_First_order_low,file="Data/5D/Kriging_First_order_low")
  save(Kriging_Second_order,file="Data/5D/Kriging_Second_order")
  save(Kriging_Second_order_up,file="Data/5D/Kriging_Second_order_up")
  save(Kriging_Second_order_low,file="Data/5D/Kriging_Second_order_low")
  print(i)
}


S1<-sensitivity::sobol(model = Kriging, X1 = X1[1:1000, ], X2 = X2[1:1000, ],
                      order=2, nboot = 100)
S2<-sensitivity::sobol(model = Kriging, X1 = X1[1:1000, ], X2 = X2[1:1000, ],
                       order=2, nboot = 100,conf=0.95)
