rm(list = ls())
graphics.off()
library(sensitivity)
library(lhs)
setwd("/storage/work/h/hxy46/Sensitivity")
Reliability<-function (X) {
  X[ ,1]*X[ ,2]+2*X[ ,1]*X[ ,3]-X[ ,4]*X[ ,5]+0.5
}

d = 6 #add a sixth dummy variable

#direct Sobol
n <- 100000
X1 <- data.frame(matrix(runif(d*n), nrow = n))
y<- Reliability(X1)
X2 <- data.frame(matrix(runif(d*n), nrow = n))
Sobol_size<-seq(n/2000,n,by=n/2000)
Sobol_First_order<-matrix(NA,nrow=length(Sobol_size),ncol=d)
Sobol_First_order_up<-matrix(NA,nrow=length(Sobol_size),ncol=d)
Sobol_First_order_low<-matrix(NA,nrow=length(Sobol_size),ncol=d)
Sobol_time<-rep(NA,length(Sobol_size))
Sobol_fail_prob<-rep(NA,length(Sobol_size)) # failure probability
for (i in 1:2000){
  start.time <- Sys.time()
  N = n/2000*i
  S<-sensitivity::sobol(model = Reliability, X1 = X1[1:N, ], X2 = X2[1:N, ],
                        order=2, nboot = 100)
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  Sobol_time[i]<-as.numeric(time.taken)
  Sobol_First_order[i, ] <- S$S$original[1:d]
  Sobol_First_order_up[i, ] <- S$S$`max. c.i.`[1:d]
  Sobol_First_order_low[i, ] <- S$S$`min. c.i.`[1:d]
  Sobol_fail_prob[i]<- length(which(y[1:N]<0))/N*100
  save(Sobol_time,file = "Data/5D/Sobol/Sobol_time")
  save(Sobol_fail_prob,file = "Data/5D/Sobol/Sobol_fail_prob")
  save(Sobol_First_order,file="Data/5D/Sobol/Sobol_First_order")
  save(Sobol_First_order_up,file="Data/5D/Sobol/Sobol_First_order_up")
  save(Sobol_First_order_low,file="Data/5D/Sobol/Sobol_First_order_low")
  print(i)
}

