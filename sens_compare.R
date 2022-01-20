# Sensitivity comparison of multiple methods
# remove all existing variables and plots
rm(list = ls())
graphics.off()
# load the required libraries
library(lhs) #Generate Latin Hypercube samples
library(fields) # Make grid plots
library(sensitivity)
library(randtoolbox) #Generate Sobol sequence (Quasi MC)
library(gstat)
library(sp)

Reliability<-function (X) {
  floor(X[ ,1]*X[ ,2]/0.75)
  #X[ ,1]*X[ ,2]
}

Kriging <- function (X){
  coordinates(X) <- ~ X1+X2
  KrigSurface <- predict(TheGStat, model=FittedModel, newdata=X)
  mean <- KrigSurface$Rel.pred
}

#N <- seq(1000,20000,by=1000)
N_K <- seq(25,25*50,by=25)
S1_Kriging <- S2_Kriging <- rep(NA,length(N_K))
S1_Kriging_up <- S2_Kriging_up <- rep(NA,length(N_K))
S1_Kriging_low <- S2_Kriging_low <- rep(NA,length(N_K))
set.seed(42)
for (i in 30: length(N_K)){
  #Kriging
  n <- N_K[i]
  x <- data.frame(randomLHS(n,2))
  y <- Reliability(x)
  Data <- data.frame(X1=x$X1,X2=x$X2,y=y)
  coordinates(Data)=~X1+X2
  load("FittedModel")
  TheGStat <- gstat(id="Rel", formula=y ~ 1, data=Data,model=FittedModel )
  
  #Randomly divide into two samples
  X <- data.frame(randomLHS(100000,2))
  ind <- sample(c(rep(TRUE,50000), rep(FALSE,50000)), replace=FALSE)
  X1 <- X[ind, ]
  X2 <- X[!ind, ]
  S<-sensitivity::sobol(model = Kriging, X1 = X1, X2 = X2,
                        order=2, nboot = 100)
  S1_Kriging[i] <- S$S$original[1]
  S2_Kriging[i] <- S$S$original[2]
  S3_Kriging[i] <- S$S$original[3]
  S1_Kriging_up[i]<-S$S$`max. c.i.`[1]
  S2_Kriging_up[i]<-S$S$`max. c.i.`[2]
  S3_Kriging_up[i] <- S$S$`max. c.i.`[3]
  S1_Kriging_low[i]<-S$S$`min. c.i.`[1]
  S2_Kriging_low[i]<-S$S$`min. c.i.`[2]
  S3_Kriging_low[i] <- S$S$`min. c.i.`[3]
}
S3_Kriging<-1-S1_Kriging-S2_Kriging
S3_Kriging_low<-1-S1_Kriging_up-S2_Kriging_up
S3_Kriging_up<-1-S1_Kriging_low-S2_Kriging_low

save(S1_Kriging,file="Data/S1_Kriging")
save(S2_Kriging,file="Data/S2_Kriging")
save(S1_Kriging_up,file="Data/S1_Kriging_up")
save(S2_Kriging_up,file="Data/S2_Kriging_up")
save(S1_Kriging_low,file="Data/S1_Kriging_low")
save(S2_Kriging_low,file="Data/S2_Kriging_low")

load("Data/S1")
load("Data/S2")
load("Data/S3")
load("Data/S1_up")
load("Data/S2_up")
load("Data/S3_up")
load("Data/S1_low")
load("Data/S2_low")
load("Data/S3_low")
load("Data/S1_Kriging")
load("Data/S2_Kriging")
#load("S3_Kriging")
load("Data/S1_Kriging_up")
load("Data/S2_Kriging_up")
#load("S3_Kriging_up")
load("Data/S1_Kriging_low")
load("Data/S2_Kriging_low")
#load("S3_Kriging_low")
load("N_K")
n=10^5
X=seq(n/100,n,by=n/100)
par(mar=c(4,5.1,1.6,2.1))
plot(X,S1*100,type="l",lwd=2,xlab="Sample size",ylab="Percentage of explained variance (%)",
     cex.lab=1.5,cex.axis=1.5,ylim=c(-10,120),xlim=c(50,n),log="x",col="red",xaxt="n")
polygon(c(X,rev(X)),c(S1_up*100,rev(S1_low*100)),col=adjustcolor("rosybrown1",alpha.f = 0.5), border = NA)
polygon(c(X,rev(X)),c(S2_up*100,rev(S2_low*100)),col=adjustcolor("lightblue1",alpha.f = 0.5), border = NA)
polygon(c(X,rev(X)),c(S3_up*100,rev(S3_low*100)),col=adjustcolor("darkseagreen1",alpha.f = 0.5), border = NA)
polygon(c(N_K,rev(N_K)),c(S1_Kriging_up*100,rev(S1_Kriging_low*100)),col=adjustcolor("rosybrown1",alpha.f = 0.5), border = NA)
polygon(c(N_K,rev(N_K)),c(S2_Kriging_up*100,rev(S2_Kriging_low*100)),col=adjustcolor("lightblue1",alpha.f = 0.5), border = NA)
polygon(c(N_K,rev(N_K)),c(S3_Kriging_up*100,rev(S3_Kriging_low*100)),col=adjustcolor("darkseagreen1",alpha.f = 0.5), border = NA)
lines(X,S1*100,lwd=2,col="red")
lines(X,S2*100,lwd=2,col="blue")
lines(X,S3*100,lwd=2,col="green")

lines(N_K,S1_Kriging*100,lwd=2,lty=5,col="red")
lines(N_K,S2_Kriging*100,lwd=2,lty=5,col="blue")
lines(N_K,S3_Kriging*100,lwd=2,lty=5,col="green")
axis(1, at = c(10,100,1000,10000,100000,1000000),cex.axis=1.5)
abline(h=14.5,lty=2) # theoretical value of S1 and S2
legend("topleft",lwd=c(2,2,2),lty=c(1,1,1),col = c("red","blue","green"),
       legend = c("X1","X2","Interaction"),cex=0.8,bty="n")
legend(150,125,lwd=c(2,2,2,1),lty=c(5,5,5,2),col = c("red","blue","green","black"),
       legend = c("X1 by Kriging","X2 by Kriging","Interaction by Kriging","Theoretical value"),cex=0.8,bty="n")
legend(44,105,fill=c("rosybrown1","lightblue1","darkseagreen1"),border = c("black","black","black"),
       legend = c("X1 95% CI","X2 95% CI","Interaction 95% CI"),cex=0.8,bty="n")


#plot(N,S1_Sobol,log="x",type="l",xlab="sample size",ylab="% of total variance",cex.lab=1.2,cex.axis=1.2)
plot(N_K,S1_Kriging*100,type="l",col="red",lwd=2,ylim=c(13,30),xlab="sample size",ylab="% of total variance",cex.lab=1.2,cex.axis=1.2)
lines(N_K,S2_Kriging*100,col="green",lwd=2)
legend("topright",lty=c(1,1),lwd=c(2,2),col=c("red","green"),legend=c("X1","X2"),cex=1.5,bty="n")