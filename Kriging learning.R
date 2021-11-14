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
library(gstat)
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
  floor(X[ ,1]*X[ ,2]/0.75)
}
# 100 total samples, start with 3 samples and fit the rest 97 by Kriging
N <- 10000
n <- 20
n_init <- n
X_all <- data.frame(randomLHS(N,2))
# Fit a variogram to the selected samples, and then adopt the fitted variogram as the
#   covariance function.

# Two functions used to optimize hyper-parameters of the variogram
# theta: decay rate;  g: nugget
gradnl <- function(par, D, Y)
{
  ## extract parameters
  theta <- par[1]
  g <- par[2]
  
  ## calculate covariance quantities from data and parameters
  n <- length(Y)
  K <- exp(-D/theta) + diag(g, n)
  Ki <- solve(K)
  dotK <- K*D/theta^2
  KiY <- Ki %*% Y
  
  ## theta component
  dlltheta <- (n/2) * t(KiY) %*% dotK %*% KiY / (t(Y) %*% KiY) - 
    (1/2)*sum(diag(Ki %*% dotK))
  
  ## g component
  dllg <- (n/2) * t(KiY) %*% KiY / (t(Y) %*% KiY) - (1/2)*sum(diag(Ki))
  
  ## combine the components into a gradient vector
  return(-c(dlltheta, dllg))
}

nl <- function(par, D, Y) 
{
  theta <- par[1]                                       ## change 1
  g <- par[2]
  n <- length(Y)
  K <- exp(-D/theta) + diag(g, n)                       ## change 2
  Ki <- solve(K)
  ldetK <- determinant(K, logarithm=TRUE)$modulus
  ll <- - (n/2)*log(t(Y) %*% Ki %*% Y) - (1/2)*ldetK
  return(-ll)
}

# x is the samples we use, X is the rest samples
x <- X_all[c(1:n), ]
X <- X_all[c((n+1):N), ]
y <- Reliability(x)
d <- distance(x)
output <- optim(c(0.1, 0.1*var(y)),nl, gradnl, method="L-BFGS-B", 
                lower=eps, upper=c(10, var(y)), D=d, Y=y)
theta <- output$par[1]
g <- output$par[2]
Sigma <- exp(-d/theta) + diag(g, nrow(x))
DX <- distance(X)
SX <- exp(-DX/theta) + diag(g, ncol(DX))
D <- distance(X, x)
S <- exp(-D/theta)
# Insert into the linear predictor
Si <- solve(Sigma)
mup <- S %*% Si %*% y
tau <- drop(t(y) %*% Si %*% y / n)
Sigmap <-tau*(SX - S %*% Si %*% t(S)) 

P_fail <- (sum(mup>0.5)+sum(y>0.5))/N
# Adaptive learning part:
# Now we build the Kriging model, the next step is to add an additional point that minimize the
#   learning function U(x)
U <- abs(0.5-mup)/diag(Sigmap)
while ((min(U)<=2) | P_fail==0){
  n <- n+1
  i <- which(U==min(U))
  if (length(i)>1){
    i <- sample(i,1)
  }
  x <- rbind.data.frame(x,X[i, ])
  y <- Reliability(x)
  d <- distance(x)
  output <- optim(c(0.1, 0.1*var(y)),nl, gradnl, method="L-BFGS-B", 
                  lower=eps, upper=c(10, var(y)), D=d, Y=y)
  theta <- output$par[1]
  g <- output$par[2]
  Sigma <- exp(-d/theta) + diag(g, nrow(x))
  # Update remaining samples
  X <- X[-i, ]
  DX <- distance(X)
  SX <- exp(-DX/theta) + diag(g, ncol(DX))
  # Update distances matrices
  D <- distance(X, x)
  S <- exp(-D/theta)
  Si <- solve(Sigma)
  mup <- S %*% Si %*% y
  tau <- drop(t(y) %*% Si %*% y / n)
  Sigmap <-tau*(SX - S %*% Si %*% t(S)) 
  Diag <- diag(Sigmap)
  U <- abs(0.5-mup)/diag(Sigmap)
  P_fail <- (sum(mup>0.5)+sum(y>0.5))/N
  print(c(n,min(U)))
}

level <- seq(0, 1, length=50)
X_ALL <- expand.grid(level,level)
DX_ALL <- distance(X_ALL)
SX_ALL <- exp(-DX_ALL/theta) + diag(g, ncol(DX_ALL))
D_ALL <- distance(X_ALL, x)
S_ALL <- exp(-D_ALL/theta)
# Insert into the linear predictor
mup_ALL <- S_ALL %*% Si %*% y
cols <- rev(heat.colors(128))
image.plot(level,level, matrix(mup_ALL,ncol=50), xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
points(x$X1[1:n_init],x$X2[1:n_init],pch=20,cex=2)
points(x$X1[(n_init+1):n],x$X2[(n_init+1):n],pch=20,cex=2,col="blue")
lines(level,0.75/level,lwd=2)