#42D problem test
rm(list = ls())
graphics.off()
library(GPfit)
library(lhs)
library(sensitivity)
setwd("/storage/work/h/hxy46/Sensitivity")

set.seed(1)
para <- matrix(runif(20*2), nrow = 2)
coef <- runif(20)
para <- ceiling(para*42)
Reliability<-function (X) {
  S<-0
  for (i in 1:20){
    dim1<-para[1,i]
    dim2<-para[2,i]
    a<-coef[i]*X[ ,dim1]*X[ ,dim2]
    S<-S+a
  }
  S-1
}


d = 42

#direct Sobol
n <- 150000
X1 <- data.frame(matrix(runif(d*n), nrow = n))
y<-Reliability(X1)
X2 <- data.frame(matrix(runif(d*n), nrow = n))

#load("Data/42D/First_order")
#load("Data/42D/First_order_up")
#load("Data/42D/First_order_low")

Sobol_size=c(seq(1000,25000,by=1000),seq(30000,150000,by=5000))
Sobol_time<-rep(NA,length(Sobol_size))
Sobol_fail_prob<-rep(NA,length(Sobol_size)) # failure probability
start.time <- Sys.time()

First_order<-matrix(NA,nrow=50,ncol=d)
First_order_up<-matrix(NA,nrow=50,ncol=d)
First_order_low<-matrix(NA,nrow=50,ncol=d)
# 
# Second_order<-matrix(NA,nrow=2000,ncol=10)
# Second_order_up<-matrix(NA,nrow=2000,ncol=10)
# Second_order_low<-matrix(NA,nrow=2000,ncol=10)

for (i in 1:50){
  N = Sobol_size[i]
  S<-sensitivity::sobol(model = Reliability, X1 = X1[1:N, ], X2 = X2[1:N, ],
                        order=2, nboot = 100)
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  Sobol_time[i]<-as.numeric(time.taken)
  First_order[i, ] <- S$S$original[1:d]
  First_order_up[i, ] <- S$S$`max. c.i.`[1:d]
  First_order_low[i, ] <- S$S$`min. c.i.`[1:d]
  Sobol_fail_prob[i]<- length(which(y[1:N]<0))/N*100
  save(Sobol_time,file = "Data/42D/Sobol/Sobol_time")
  save(Sobol_fail_prob,file = "Data/42D/Sobol/Sobol_fail_prob")
  save(First_order,file="Data/42D/Sobol/Sobol_First_order")
  save(First_order_up,file="Data/42D/Sobol/Sobol_First_order_up")
  save(First_order_low,file="Data/42D/Sobol/Sobol_First_order_low")
  print(i)
}

