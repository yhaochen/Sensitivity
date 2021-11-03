# This code illustrates a Gaussian process emulator model example
# Author: Haochen Ye (hxy46@psu.edu)
#------------------

# remove all existing variables and plots
rm(list = ls())
graphics.off()
# load the required libraries
library(plgp)
library(mvtnorm)  # Multivariate normal distribution library
library(lhs) #Generate Latin Hypercube samples
library(fields) # Make grid plots

# A simple 1-D Gaussian process example:
# input samples
n <- 100 
X <- matrix(seq(0, 10, length=n), ncol=1)

# Distance between each input pair, used to define covariance matrix
D <- distance(X)
# Adding a very small number to aviod ill-conditioned matrix calculation that may happen 
#   for square matrix sometimes. (Neal, 1998)
eps <- sqrt(.Machine$double.eps)
# Assume the covariance decays exponentially with respect to distance
Sigma <- exp(-D) + diag(eps, n)
image(Sigma)
# Generate random output values by mean and sigma (assume mean is 1 everywhere) 
# Now Y is a random smooth function under a Gaussian process prior, but still a finite
#     dimension realization (n = 100).
Y <- mvtnorm::rmvnorm(1, rep(1,100), sigma=Sigma)
plot(X, Y, type="l",cex.axis=1.5,cex.lab=1.5)
# Y can be very different because we are not making restrictions to any sample point:
Y <- mvtnorm::rmvnorm(n=5,rep(1,100), sigma=Sigma)
matplot(X, t(Y),lwd=2, type="l", ylab="Y",cex.axis=1.5,cex.lab=1.5)

# Now assume we have some sample points and want to use a Gaussian process to predict the 
#   conditional distribution of unobserved outputs.
# 8 samples with equal space between every neighbour 
n <- 8
x <- matrix(seq(0, 10, length=n), ncol=1)
# Assigned sample outputs
y <- c(0,0.4,1.4,0.8,-0.2,1,1.5,2.3)
plot(x,y,pch=20,cex.axis=1.5,cex.lab=1.5)
# Sample distances and covariance
d <- distance(x) 
Sigma <- exp(-d) + diag(eps, ncol(d))
# Unobserved data points
X <- matrix(seq(0, 10, length=100), ncol=1)
# The entire prediction process is actually a linear predictor.
# To predict the conditional distribution of outputs, we further need (1) the covariance of 
#     unobserved data (SX); (2) the covariance between samples and unobserved data (S).
DX <- distance(X)
SX <- exp(-DX) + diag(eps, ncol(DX))
D <- distance(X, x)
S <- exp(-D)
# Calculate the inverse of observed covariance matrix, and then apply the linear predictor.
# The added eps ensures the solve() function works properly.
Si <- solve(Sigma)
mup <- S %*% Si %*% y
Sigmap <- SX - S %*% Si %*% t(S)

# Now we can generate outputs by the calculated mu and covariance
Y <- mvtnorm::rmvnorm(1, mup, Sigmap)
plot(x,y,pch=20,ylim=c(-1,3),cex.axis=1.5,cex.lab=1.5)
lines(X,Y)
# This is an example of possible outputs conditioned on samples, we add the mean values, which 
#     are the best estimates.
lines(X,mup,col="blue")
# We can also add 95% confidence intervals:
up <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
low <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))
lines(X,up,col="red",lty=2)
lines(X,low,col="red",lty=2)
legend("topleft",col=c("black","black","blue","red"),lty=c(NA,1,1,2),pch=c(20,NA,NA,NA),bty="n",
       legend=c("Sample","Example of a function fit","Best estimate","95% uncertainty range"))

# Then we apply it to our 2-d reliability problem
Reliability<-function (X) {
  floor(X[ ,1]*X[ ,2]/0.75)
}
# First sample 2-d inputs, calculate corresponding output, distance, covariance
n <- 100
x <- data.frame(randomLHS(n,2))
y <- Reliability(x)
d <- distance(x)
Sigma <- exp(-d) + diag(eps, nrow(x))
# Then select some evenly distributed unobserved data points and calculate covariance
level <- seq(0, 1, length=50)
X <- expand.grid(level,level)
DX <- distance(X)
SX <- exp(-DX) + diag(eps, ncol(DX))
D <- distance(X, x)
S <- exp(-D)
# Insert into the linear predictor
Si <- solve(Sigma)
mup <- S %*% Si %*% y
Sigmap <- SX - S %*% Si %*% t(S)
# We can visualize the mean value of these unobserved data
cols <- rev(heat.colors(128))
image.plot(level,level, matrix(mup,ncol=50), xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
points(x$X1,x$X2,pch=20)
lines(level,0.75/level,col="red",lwd=2)
# The standard deviation plot
image.plot(level,level, matrix(sqrt(diag(Sigmap)),ncol=50), xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
points(x$X1,x$X2,pch=20)
lines(level,0.75/level,col="red",lwd=2)
