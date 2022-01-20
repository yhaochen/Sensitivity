rm(list = ls())
graphics.off()
library(GPfit)
library(lhs)
library(sensitivity)
setwd("/storage/work/h/hxy46/Sensitivity")
Reliability<-function (X) {
  X[ ,1]*X[ ,2]+2*X[ ,1]*X[ ,3]-X[ ,4]*X[ ,5]+0.5
}

AKMCS <- function (X){
  a<-predict(GPmodel,X)
  a$Y_hat
}

d = 6
set.seed(1)


#AKMCS
n=1000
n_init<-9
x_all <- maximinLHS(n, d)
indx <- sample(n,n_init)
x <- x_all[indx, ]
x_rest <- x_all[-indx, ]
y <- Reliability(x)
GPmodel = GP_fit(x, y)
a<-predict(GPmodel,x_rest)
mean<-a$Y_hat
std<-sqrt(a$MSE)
U<-mean/std

AKMCS_size<-c(10:100)
Step<-length(AKMCS_size)
AKMCS_First_order<-matrix(NA,nrow=Step,ncol=d)
AKMCS_First_order_up<-matrix(NA,nrow=Step,ncol=d)
AKMCS_First_order_low<-matrix(NA,nrow=Step,ncol=d)
AKMCS_time<-rep(NA,Step)
AKMCS_fail_prob<-rep(NA,Step) # failure probability
AKMCS_wrong_pos<-rep(NA,Step) # % of wrong success (should be failure)
AKMCS_wrong_neg<-rep(NA,Step) # % of wrong failure (should be success)

N <- 50000
X1 <- data.frame(matrix(runif(d*N), nrow = N))
X2 <- data.frame(matrix(runif(d*N), nrow = N))

y1 <- Reliability(X1)
true_neg_indx<-which(y1<0)


for (i in 1:Step){
  start.time <- Sys.time()
  k<-which(U==min(U))
  if (length(k)>1){
    k <- sample(k,1)
  }
  x_add<-matrix(x_rest[k, ],ncol=d)
  x_rest<-x_rest[-k, ]
  y_add<-Reliability(x_add)
  y <- append(y,y_add)
  x <- rbind(x,x_add)
  GPmodel = GP_fit(x, y)
  a<-predict(GPmodel,x_rest)
  mean<-a$Y_hat
  b<-predict(GPmodel,X1)
  yfit<-b$Y_hat
  AKMCS_fail_prob[i]<-length(which(yfit<0))/N*100
  AKMCS_wrong_pos[i]<-length(which(yfit[true_neg_indx]>0))/length(true_neg_indx)*100
  AKMCS_wrong_neg[i]<-length(which(yfit[-true_neg_indx]<0))/(N-length(true_neg_indx))*100
  std<-sqrt(a$MSE)
  U<-mean/std
  S<-sensitivity::sobol(model = AKMCS, X1 = X1, X2 = X2,
                        order=2, nboot = 100)
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  AKMCS_time[i]<-as.numeric(time.taken)
  AKMCS_First_order[i, ] <- S$S$original[1:d]
  AKMCS_First_order_up[i, ] <- S$S$`max. c.i.`[1:d]
  AKMCS_First_order_low[i, ] <- S$S$`min. c.i.`[1:d]
  save(AKMCS_time,file="Data/5D/AKMCS/AKMCS_time")
  save(AKMCS_fail_prob,file="Data/5D/AKMCS/AKMCS_fail_prob")
  save(AKMCS_wrong_neg,file = "Data/5D/AKMCS/AKMCS_wrong_neg")
  save(AKMCS_wrong_pos,file = "Data/5D/AKMCS/AKMCS_wrong_pos")
  save(AKMCS_First_order,file="Data/5D/AKMCS/AKMCS_First_order")
  save(AKMCS_First_order_up,file="Data/5D/AKMCS/AKMCS_First_order_up")
  save(AKMCS_First_order_low,file="Data/5D/AKMCS/AKMCS_First_order_low")
  print(i)
}
