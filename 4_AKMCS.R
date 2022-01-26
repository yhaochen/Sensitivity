# AKMCS method
rm(list = ls())
graphics.off()
setwd("/storage/work/h/hxy46/Sensitivity")
library(sensitivity)
library(lhs)
library(GPfit)
# Load problem definition
source("1_Problem_definition.R")
set.seed(1)
# Create a folder for d dimension test scenario if not yet created
folder<-paste("/storage/work/h/hxy46/Sensitivity/Data/",d,"D/AKMCS",sep="")
dir.create(file.path(folder), showWarnings = FALSE)

AKMCS_step<-91
AKMCS_size<-c(10:100)
save(AKMCS_size,file = paste(folder,"/AKMCS_size",sep=""))
AKMCS_First_order<-matrix(NA,nrow=AKMCS_step,ncol=d)
# Record two times, one for AKMCS emulator, one for Sobol based on this emulator
# Note that for AKMCS, time for emulator part should be cumulative
AKMCS_Krigingtime<-rep(NA,AKMCS_step)
AKMCS_Soboltime<-rep(NA,AKMCS_step)
AKMCS_fail_prob<-rep(NA,AKMCS_step) # estimated failure probability

# For emulator based method, also record their wrong classification 
AKMCS_wrong_pos<-rep(NA,AKMCS_step) # % of wrong success (should be failure)
AKMCS_wrong_neg<-rep(NA,AKMCS_step) # % of wrong failure (should be success)

#Kriging function, used in Sobol analysis based on Kriging emulator
Kriging <- function (X){
  a<-predict(GPmodel,X)
  a$Y_hat
}

# AKMCS initial state
n_init<-9
indx <- sample(N,n_init)
x <- X1[indx, ]
x_rest <- X1[-indx, ]
y <- Reliability(x)
GPmodel = GP_fit(x, y)
a<-predict(GPmodel,x_rest)
mean<-a$Y_hat
std<-sqrt(a$MSE)
U<-(mean-thres)/std

for (i in 1:AKMCS_step){
  start.time <- Sys.time()
  k<-which(U==min(U))
  if (length(k)>1){
    k <- sample(k,1)
  }
  x_add<-x_rest[k, ]
  x_rest<-x_rest[-k, ]
  y_add<-Reliability(x_add)
  y <- append(y,y_add)
  x <- rbind(x,x_add)
  GPmodel = GP_fit(x, y)
  a<-predict(GPmodel,x_rest)
  mean<-a$Y_hat
  b<-predict(GPmodel,X1)
  yfit<-b$Y_hat
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  AKMCS_Krigingtime[i]<-as.numeric(time.taken)
  AKMCS_fail_prob[i]<-length(which(yfit<0))/N*100
  AKMCS_wrong_pos[i]<-length(which(yfit[true_neg_indx]>0))/length(true_neg_indx)*100
  AKMCS_wrong_neg[i]<-length(which(yfit[-true_neg_indx]<0))/(N-length(true_neg_indx))*100
  std<-sqrt(a$MSE)
  U<-(mean-thres)/std
  start.time <- Sys.time()
  S<-sensitivity::sobol(model = Kriging, X1 = X1, X2 = X2,
                        order=1, nboot = 100)
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  AKMCS_Soboltime[i]<-as.numeric(time.taken)
  
  AKMCS_First_order[i, ] <- S$S$original[1:d]
  save(AKMCS_Krigingtime,file = paste(folder,"/AKMCS_Krigingtime",sep=""))
  save(AKMCS_Soboltime,file = paste(folder,"/AKMCS_Soboltime",sep=""))
  save(AKMCS_fail_prob,file = paste(folder,"/AKMCS_fail_prob",sep=""))
  save(AKMCS_wrong_neg,file = paste(folder,"/AKMCS_wrong_neg",sep=""))
  save(AKMCS_wrong_pos,file = paste(folder,"/AKMCS_wrong_pos",sep=""))
  save(AKMCS_First_order,file = paste(folder,"/AKMCS_First_order",sep=""))
  print(i)
}
