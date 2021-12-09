#Perform AKMCS algorithm by gstat package
#2D example problem
#Author: Haochen Ye, hxy46@psu.edu
rm(list = ls())
graphics.off()
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
thres <- 0.5
S1_Kriging <- S2_Kriging <- rep(NA,100)
total_size <- rep (NA,100)

Kriging <- function (X){
  coordinates(X) <- ~ X1+X2
  KrigSurface <- predict(TheGStat, model=FittedModel, newdata=X)
  mean <- KrigSurface$Rel.pred
}

set.seed(42)
#Repeat for many times to get the distribution of sensitivity indices
total_run <- 1
for (k in 1:total_run){
level <- seq(0, 1, length=100)
N <- 10000
X_ALL <- expand.grid(level,level)
Grid <- data.frame(X1=X_ALL$Var1,X2=X_ALL$Var2)
n <- 30
n_init <- n
indx <- sample(N,n)
x <- Grid[indx, ]
X <- Grid[-indx, ]
y <- Reliability(x)
Data <- data.frame(X1=x$X1,X2=x$X2,y=y)
coordinates(Data)=~X1+X2
TheVariogram=variogram(y~1, data=Data)
load("FittedModel")

TheGStat <- gstat(id="Rel", formula=y ~ 1, data=Data,model=FittedModel )
XX <- X
coordinates(XX) <- ~ X1+X2
KrigSurface <- predict(TheGStat, model=FittedModel, newdata=XX)
#gstat.cv(TheGStat,nfold=10)
mean <- KrigSurface$Rel.pred
std <- sqrt(KrigSurface$Rel.var)
P_fail <- (sum(mean>thres)+sum(y>thres))/N
U <- abs(mean-thres)/std
if (is.na(P_fail))
  P_fail <- 0
if (is.na(U[1]))
  U <- 0
print(min(U))
while ((min(U)<=0.05) | P_fail==0){
  n <- n+1
  i <- which(U==min(U))
  if (length(i)>1){
    i <- sample(i,1)
  }
  y <- append(y,Reliability(X[i, ]))
  x <- rbind.data.frame(x,X[i, ])
  X <- X[-i, ]
  XX <- X
  coordinates(XX) <- ~ X1+X2
  Data <- data.frame(X1=x$X1,X2=x$X2,y=y)
  coordinates(Data)=~X1+X2
  TheVariogram=variogram(y~1, data=Data)
  #FittedModel <- fit.variogram(TheVariogram, vgm(psill=0.015, model="Exp", nugget=0.05, range=0.1))
  #FittedModel$psill[1] <- 0.01
  TheGStat <- gstat(id="Rel", formula=y ~ 1, data=Data,model=FittedModel )
  KrigSurface <- predict(TheGStat, model=FittedModel, newdata=XX)
  mean <- KrigSurface$Rel.pred
  std <- sqrt(KrigSurface$Rel.var)
  P_fail <- (sum(mean>thres)+sum(y>thres))/N
  if (is.na(P_fail))
    P_fail <- 0
  U <- abs(mean-thres)/std
  if (is.na(U[1]))
    U <- 0
  print(c(n,min(U),P_fail))
}

Grid <- data.frame(X1=X_ALL$Var1,X2=X_ALL$Var2)
y=Reliability(Grid)
y=matrix(y,ncol=100)
coordinates(Grid) <- ~ X1+X2
KrigSurface <- predict(TheGStat, model=FittedModel, newdata=Grid)

colfunc <- colorRampPalette(c("white", "red"))
cols <- colfunc(101)
par(mar=c(4,5.1,1.6,2.1))
image.plot(level,level, matrix(KrigSurface$Rel.pred,ncol=100),xaxt="n",yaxt="n", xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
#image.plot(level,level, abs(Y-y),xaxt="n",yaxt="n", xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)

Y <- KrigSurface$Rel.pred
points(x$X1[1:n_init],x$X2[1:n_init],pch=21,bg=cols[ceiling((y[1:n_init])*100)+1])
points(x$X1[(n_init+1):n],x$X2[(n_init+1):n],pch=21,bg=cols[ceiling((y[(n_init+1):n])*100)+1])
lines(level,0.75/level,lwd=2)

  #Randomly divide into two samples
  X <- data.frame(randomLHS(100000,2))
  ind <- sample(c(rep(TRUE,50000), rep(FALSE,50000)), replace=FALSE)
  X1 <- X[ind, ]
  X2 <- X[!ind, ]
  S<-sensitivity::sobol(model = Kriging, X1 = X1, X2 = X2,
                        order=2, nboot = 100)
  S1_Kriging[k] <- S$S$original[1]
  S2_Kriging[k] <- S$S$original[2]
  total_size[k] <- n
}

h1=hist(S1_Kriging,breaks = seq(0,0.4, length.out = 11))
h2=hist(S2_Kriging,breaks = seq(0,0.4, length.out = 11))

plot( h1, col=rgb(0,0,1,1/4), xlim=c(0.05,0.4),main="",
      xlab="Sobol first-order sensitivity index")  # first histogram
plot( h2, col=rgb(1,0,0,1/4), xlim=c(0.05,0.4), add=T)  # second
legend("topright",fill=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),legend=c("X1","X2"),bty="n")