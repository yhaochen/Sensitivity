# This script shows an example of adaptive learning reliability method combining Kriging and Monte
#   Carlo simulation (AK-MCS)
# Author: Haochen Ye (hxy46@psu.edu)
#------------------
# The Kriging model follows the example in Gaussian_process.R
# References: Echard et al, 2011; Gramacy, 2020

# remove all existing variables and plots
rm(list = ls())
graphics.off()
# load the required libraries
library(plgp)
library(mvtnorm)  # Multivariate normal distribution library
library(lhs) #Generate Latin Hypercube samples
library(fields) # Make grid plots
library(mistral)
# Adding a very small number to aviod ill-conditioned matrix calculation that may happen 
#   for square matrix sometimes. (Neal, 1998)
g <- 10^(-13)#sqrt(.Machine$double.eps)
# Problem definition
# Note in this example we modify the output as continuous, success is positive, fail is negative
#   limit state is zero. This is a requirement of AK-MCS algorithm.

# Kriging part: 
# Choose a small subsample from a larger sample and fit a Kriging model
# Initial sample size is recommended as (n+1)(n+2)/2 = 3*4/2 = 6
Reliability<-function (X) {
  floor(X[ ,1]*X[ ,2]/0.75)
  #X[ ,1]*X[ ,2]
}
thres <- 0.5
# 100 total samples, start with 3 samples and fit the rest 97 by Kriging
N <- 5000
CV <- rep(NA, 10)
for (k in 1:1){
  n=100
  n_init <- n
  X_all <- data.frame(randomLHS(N,2))
  Y_all <- Reliability(X_all)
  # Fit a variogram to the selected samples, and then adopt the fitted variogram as the
  #   covariance function.
  # x is the samples we use, X is the rest samples
  x <- X_all[c(1:n), ]
  X <- X_all[c((n+1):N), ]
  Y_all <- Y_all[c((n+1):N)]
  y <- Reliability(x)
  d <- distance(x)
  theta <- 1
  Sigma <- exp(-d/theta) + diag(g, nrow(x))
  DX <- distance(X)
  SX <- exp(-DX/theta) + diag(g, ncol(DX))
  D <- distance(X, x)
  S <- exp(-D/theta)
  # Insert into the linear predictor
  Si <- solve(Sigma)
  mup <- S %*% Si %*% y
  Sigmap <- SX - S %*% Si %*% t(S)
  
  P_fail <- (sum(mup>=thres)+sum(y>=thres))/N
  # Adaptive learning part:
  # Now we build the Kriging model, the next step is to add an additional point that minimize the
  #   learning function U(x)
  U <- abs(thres-mup)/diag(Sigmap)
  while ((min(U)<=-1) | (P_fail==0)){
    n <- n+1
    i <- which(U==min(U))
    if (length(i)>1){
      i <- sample(i,1)
    }
    x <- rbind.data.frame(x,X[i, ])
    y <- Reliability(x)
    d <- distance(x)
    Sigma <- exp(-d/theta) + diag(g, nrow(x))
    # Update remaining samples
    X <- X[-i, ]
    Y_all <- Y_all[-i]
    DX <- distance(X)
    SX <- exp(-DX/theta) + diag(g, ncol(DX))
    # Update distances matrices
    D <- distance(X, x)
    S <- exp(-D/theta)
    Si <- solve(Sigma)
    mup <- S %*% Si %*% y
    Sigmap <- SX - S %*% Si %*% t(S)
    Diag <- diag(Sigmap)
    U <- abs(thres-mup)/abs(diag(Sigmap))
    P_fail <- (sum(mup>=thres)+sum(y>=thres))/N
    print(c(n,min(U)))
  }
  CV[k] <- mean(sqrt(var(Y_all-mup)))
}
level <- seq(0, 1, length=50)
X_ALL <- expand.grid(level,level)
DX_ALL <- distance(X_ALL)
SX_ALL <- exp(-DX_ALL/theta) + diag(g, ncol(DX_ALL))
D_ALL <- distance(X_ALL, x)
S_ALL <- exp(-D_ALL/theta)
# Insert into the linear predictor
mup_ALL <- S_ALL %*% Si %*% y
Sigma_ALL <- SX_ALL - S_ALL %*% Si %*% t(S_ALL)
DIAG_ALL <- diag(Sigma_ALL)
U_ALL <- abs(thres-mup_ALL)/abs(diag(Sigma_ALL))
colfunc <- colorRampPalette(c("white", "red"))
cols <- colfunc(100)
par(mfrow=c(2,2))
par(mar=c(2.1,4.4,1.6,3))
image.plot(level,level, matrix(mup_ALL,ncol=50),xaxt="n",yaxt="n", xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
axis(1, at=0:1, labels=c(0,1),cex.axis=1.5)
axis(2, at=0:1, labels=c(0,1),cex.axis=1.5)
points(x$X1[1:n_init],x$X2[1:n_init],pch=20,cex=2)
points(x$X1[(n_init+1):N],x$X2[(n_init+1):N],pch=20,cex=2,col="blue")
lines(level,0.75/level,lwd=2)
par(mar=c(2.1,4.4,1.6,4.5))
image.plot(level,level, matrix(DIAG_ALL,ncol=50),xaxt="n",yaxt="n", xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
axis(1, at=0:1, labels=c(0,1),cex.axis=1.5)
axis(2, at=0:1, labels=c(0,1),cex.axis=1.5)
points(x$X1[1:n_init],x$X2[1:n_init],pch=20,cex=2)
points(x$X1[(n_init+1):N],x$X2[(n_init+1):N],pch=20,cex=2,col="blue")
lines(level,0.75/level,lwd=2)
par(mar=c(2.1,4.4,1.6,4.5))
image.plot(level,level, matrix(1/U_ALL,ncol=50),xaxt="n",yaxt="n", xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
axis(1, at=0:1, labels=c(0,1),cex.axis=1.5)
axis(2, at=0:1, labels=c(0,1),cex.axis=1.5)
points(x$X1[1:n_init],x$X2[1:n_init],pch=20,cex=2)
points(x$X1[(n_init+1):N],x$X2[(n_init+1):N],pch=20,cex=2,col="blue")
lines(level,0.75/level,lwd=2)




rm(list = ls())
graphics.off()
# load the required libraries
library(plgp)
library(mvtnorm)  # Multivariate normal distribution library
library(lhs) #Generate Latin Hypercube samples
library(fields) # Make grid plots
library(mistral)
g <- 10^(-13)#sqrt(.Machine$double.eps)
Reliability<-function (X) {
  floor(X[ ,1]*X[ ,2]/0.75)
  #X[ ,1]*X[ ,2]
}

n=100
x <- data.frame(randomLHS(n,2))
y <- Reliability(x)
d <- distance(x)
theta <- 1
Sigma <- exp(-d/theta) + diag(g, nrow(x))

level <- seq(0, 1, length=50)
X <- expand.grid(level,level)

D <- distance(X, x)
S <- exp(-D/theta)
Si <- solve(Sigma)
mup <- S %*% Si %*% y

colfunc <- colorRampPalette(c("white", "red"))
cols <- colfunc(100)
image.plot(level,level, matrix(mup,ncol=50),xaxt="n",yaxt="n", xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
axis(1, at=0:1, labels=c(0,1),cex.axis=1.5)
axis(2, at=0:1, labels=c(0,1),cex.axis=1.5)
points(x$X1,x$X2,pch=20,cex=2)
