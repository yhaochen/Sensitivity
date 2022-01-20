rm(list = ls())
graphics.off()
library(GPfit)
library(lhs)
setwd("/storage/work/h/hxy46/Sensitivity")
Reliability<-function (X) {
  X[ ,1]*X[ ,2]+2*X[ ,1]*X[ ,3]-X[ ,4]*X[ ,5]+0.5
}

d = 6 #add a sixth dummy variable

#direct Sobol
N <- 50000
X1 <- data.frame(matrix(runif(d*N), nrow = N))
X2 <- data.frame(matrix(runif(d*N), nrow = N))
y<- Reliability(X1)
true_neg_indx<-which(y<0)
#Kriging
Kriging <- function (X){
  a<-predict(GPmodel,X)
  a$Y_hat
}

Kriging_size<-c(10:100)
Step<-length(Kriging_size)
Kriging_First_order<-matrix(NA,nrow=Step,ncol=d)
Kriging_First_order_up<-matrix(NA,nrow=Step,ncol=d)
Kriging_First_order_low<-matrix(NA,nrow=Step,ncol=d)
Kriging_time<-rep(NA,Step)
Kriging_fail_prob<-rep(NA,Step) # failure probability
Kriging_wrong_pos<-rep(NA,Step) # % of wrong success (should be failure)
Kriging_wrong_neg<-rep(NA,Step) # % of wrong failure (should be success)
for (i in 1:Step){
  start.time <- Sys.time()
  n = 9+i
  x = maximinLHS(n, d)
  y = Reliability(x)
  GPmodel = GP_fit(x, y)
  yfit<-Kriging(X1)
  Kriging_fail_prob[i]<-length(which(yfit<0))/N*100
  Kriging_wrong_pos[i]<-length(which(yfit[true_neg_indx]>0))/length(true_neg_indx)*100
  Kriging_wrong_neg[i]<-length(which(yfit[-true_neg_indx]<0))/(N-length(true_neg_indx))*100
  S<-sensitivity::sobol(model = Kriging, X1 = X1, X2 = X2,
                        order=2, nboot = 100)
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  Kriging_time[i]<-as.numeric(time.taken)
  Kriging_First_order[i, ] <- S$S$original[1:d]
  Kriging_First_order_up[i, ] <- S$S$`max. c.i.`[1:d]
  Kriging_First_order_low[i, ] <- S$S$`min. c.i.`[1:d]
  save(Kriging_time,file="Data/5D/Kriging/Kriging_time")
  save(Kriging_fail_prob,file="Data/5D/Kriging/Kriging_fail_prob")
  save(Kriging_wrong_neg,file = "Data/5D/Kriging/Kriging_wrong_neg")
  save(Kriging_wrong_pos,file = "Data/5D/Kriging/Kriging_wrong_pos")
  save(Kriging_First_order,file="Data/5D/Kriging/Kriging_First_order")
  save(Kriging_First_order_up,file="Data/5D/Kriging/Kriging_First_order_up")
  save(Kriging_First_order_low,file="Data/5D/Kriging/Kriging_First_order_low")
  print(i)
}