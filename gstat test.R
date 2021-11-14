#gstat package

rm(list = ls())
graphics.off()
library(lhs)
library(gstat)
library(sp)
library(fields)
Reliability<-function (X) {
  #floor(X[ ,1]*X[ ,2]/0.75)
  X[ ,1]*X[ ,2]
}
thres <- 0.75
level <- seq(0, 1, length=100)
N <- 10000
X_ALL <- expand.grid(level,level)
Grid <- data.frame(X1=X_ALL$Var1,X2=X_ALL$Var2)
n <- 20
n_init <- n
indx <- sample(N,n)
x <- Grid[indx, ]
X <- Grid[-indx, ]
y <- Reliability(x)
Data <- data.frame(X1=x$X1,X2=x$X2,y=y)
coordinates(Data)=~X1+X2
TheVariogram=variogram(y~1, data=Data)
FittedModel <- fit.variogram(TheVariogram, vgm("Exp"))
TheGStat <- gstat(id="Rel", formula=y ~ 1, data=Data,model=FittedModel )
XX <- X
coordinates(XX) <- ~ X1+X2
KrigSurface <- predict(TheGStat, model=FittedModel, newdata=XX)
mean <- KrigSurface$Rel.pred
std <- sqrt(KrigSurface$Rel.var)
P_fail <- (sum(mean>thres)+sum(y>thres))/N
U <- abs(mean-thres)/std
if (is.na(P_fail))
  P_fail <- 0
if (is.na(U[1]))
  U <- 0
while ((min(U)<=0.03) | P_fail==0){
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
  FittedModel <- fit.variogram(TheVariogram, vgm("Exp"))
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
coordinates(Grid) <- ~ X1+X2
KrigSurface <- predict(TheGStat, model=FittedModel, newdata=Grid)

colfunc <- colorRampPalette(c("white", "red"))
cols <- colfunc(100)
image.plot(level,level, matrix(KrigSurface$Rel.pred,ncol=100),xaxt="n",yaxt="n", xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
points(x$X1[1:n_init],x$X2[1:n_init],pch=20,cex=2)
points(x$X1[(n_init+1):n],x$X2[(n_init+1):n],pch=20,cex=2,col="blue")
lines(level,0.75/level,lwd=2)