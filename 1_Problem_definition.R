# Scenario definition (dimension, model, single evaluation time, LHS samples)
rm(list = ls())
graphics.off()
library(lhs)
setwd("/storage/work/h/hxy46/Sensitivity")
set.seed(1)
d = 2 # Dimension of the model (number of input parameter)

# The model is defined as a sum (altogether d terms) of multiple products of random two parameters
# para generates the parameters in the model
# coef generates the coefficient of each product

# Example: when d=2, para is (1,1,2,1), coef is (0.2, 0.5), then the model will be:
#   y = 0.2x1x1+0.5x1x2

para <- matrix(runif(d*2), nrow = 2)
para <- ceiling(para*d)
coef <- runif(d)

Reliability<-function (X) {
  S<-0
  for (i in 1:d){
    dim1<-para[1,i]
    dim2<-para[2,i]
    a<-coef[i]*X[ ,dim1]*X[ ,dim2]
    S<-S+a
  }
  S
}

# Two LHS samples, each with a sample size N. These two samples can be used for Sobol 
#    analysis. A subset of the first sample is also used for Kriging and AKMCS methods.
N <- 50000
X1 <- data.frame(matrix(runif(d*N), nrow = N))
X2 <- data.frame(matrix(runif(d*N), nrow = N))
y<- Reliability(X1)

# Define a success/failure threshold as 2% of lower tail 
q <- 0.02
thres<-quantile(y,q)

# These are the indices of samples at "failure" region
true_neg_indx<-which(y<thres)

