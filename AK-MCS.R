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
library(DiceKriging)
# Adding a very small number to aviod ill-conditioned matrix calculation that may happen 
#   for square matrix sometimes. (Neal, 1998)
eps <- sqrt(.Machine$double.eps)
# Problem definition
# Note in this example we modify the output as continuous, success is positive, fail is negative
#   limit state is zero. This is a requirement of AK-MCS algorithm.

# Kriging part: 
# Choose a small subsample from a larger sample and fit a Kriging model
# Initial sample size is recommended as (n+1)(n+2)/2 = 3*4/2 = 6
Reliability<-function (X) {
  0.75-X[ ,1]*X[ ,2]
}
# 100 total samples, start with 3 samples and fit the rest 97 by Kriging
N <- 15000
n <- 3
n_init <- n
X_all <- data.frame(randomLHS(N,2))
# x is the samples we use, X is the rest samples
x <- X_all[c(1:n), ]
X <- X_all[c((n+1):N), ]
#y <- branin(x)$X2-10
y <- Reliability(x)
d <- distance(x)
Sigma <- exp(-d) + diag(eps, nrow(x))
DX <- distance(X)
SX <- exp(-DX) + diag(eps, ncol(DX))
D <- distance(X, x)
S <- exp(-D)
# Insert into the linear predictor
Si <- solve(Sigma)
mup <- S %*% Si %*% y
Sigmap <- SX - S %*% Si %*% t(S)
P_fail <- (sum(mup<0)+sum(y<0))/N
COV <- sqrt((1-P_fail)/(P_fail*N))
# Adaptive learning part:
# Now we build the Kriging model, the next step is to add an additional point that minimize the
#   learning function U(x)
U <- abs(mup)/diag(Sigmap)
while ((min(abs(U))<=2) | (COV>0.05)){
  n <- n+1
  i <- which(U==min(U))
  x <- rbind.data.frame(x,X[i, ])
  y <- Reliability(x)
  #y <- branin(x)$X2-10
  d <- distance(x)
  Sigma <- exp(-d) + diag(eps, nrow(x))
  # Update remaining samples
  X <- X[-i, ]
  DX <- distance(X)
  SX <- exp(-DX) + diag(eps, ncol(DX))
  # Update distances matrices
  D <- distance(X, x)
  S <- exp(-D)
  Si <- solve(Sigma)
  mup <- S %*% Si %*% y
  Sigmap <- SX - S %*% Si %*% t(S)
  Diag <- diag(Sigmap)
  U <- abs(mup)/diag(Sigmap)
  P_fail <- (sum(mup<0)+sum(y<0))/N
  COV <- sqrt((1-P_fail)/(P_fail*N))
  print(c(n,min(U),COV))
}

level <- seq(0, 1, length=50)
X_ALL <- expand.grid(level,level)
DX_ALL <- distance(X_ALL)
SX_ALL <- exp(-DX_ALL) + diag(eps, ncol(DX_ALL))
D_ALL <- distance(X_ALL, x)
S_ALL <- exp(-D_ALL)
# Insert into the linear predictor
mup_ALL <- S_ALL %*% Si %*% y
cols <- rev(heat.colors(128))
image.plot(level,level, matrix(mup_ALL,ncol=50), xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
points(x$X1[1:n_init],x$X2[1:n_init],pch=20,cex=2)
points(x$X1[(n_init+1):n],x$X2[(n_init+1):n],pch=20,cex=2,col="blue")
lines(level,0.75/level,lwd=2)
