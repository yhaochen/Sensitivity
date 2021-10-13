#An example of Sobol and Morris method applied to a simple reliability problem

rm(list = ls())
graphics.off()
library(sensitivity)
library(lhs)
library(randtoolbox)
#Problem definition: return yes if the product of two numbers >=0.75
Reliability<-function (X) {
  floor(X[ ,1]*X[ ,2]/0.75)
}

#Plot
x <- seq(0.75,1,by=0.01)
plot(x,0.75/x,type="l",lwd=2,xlim=c(0,1),ylim=c(0,1),xlab="x1",ylab="x2")
polygon(c(0.75,1,rev(x)),c(1,1,rev(0.75/x)),col="red",border = NA)
#Sobol analysis

#Base sample size (note it needs two samples of equal size in "sensitivity" package)
set.seed(1)
n <- 50000

#Random sampling
X1 <- data.frame(matrix(runif(2*n), nrow = n))
X2 <- data.frame(matrix(runif(2*n), nrow = n))

#Quasi MC -- Sobol sequence
X <- data.frame(randtoolbox::sobol(n = 2*n, dim = 2))
#Randomly separate into two samples
ind <- sample(c(rep(TRUE,n), rep(FALSE,n)), replace=FALSE)
X1 <- X[ind, ]
X2 <- X[!ind, ]

#Latin Hypercube sampling
X1 <- data.frame(randomLHS(n,2))
X2 <- data.frame(randomLHS(n,2))

#Sobol index with different base sample size
S1 <- S2 <- rep(NA,100)
for (i in 1:100){
  S<-sensitivity::sobol(model = Reliability, X1 = X1[1:(n/100*i), ], X2 = X2[1:(n/100*i), ],  nboot = 100)
  S1[i]<-S$S$original[1]
  S2[i]<-S$S$original[2]
}

#Plot sensitivity index
par(mar=c(4,5.1,1.6,2.1))
plot(seq(n/50,2*n,by=n/50),S1,type="l",lwd=2,xlab="Sample size",ylab="First-order sensitivity index",
     cex.lab=2,cex.axis=2,ylim=c(0,0.4),col="red")
lines(seq(n/50,2*n,by=n/50),S2,lwd=2,col="green")
abline(h=0.147,lty=2)
legend("topright",lwd=c(2,2),col = c("red","green"),legend = c("X1","X2"),cex=2,bty="n")
title(main = "Random MC",cex=2)

#Morris method (number of model runs = r*(factors+1))
M <- morris(model = Reliability, factors = 2, r = 100,
            design = list(type = "oat", levels = 200, grid.jump = 100))
print(M)
plot(M)

