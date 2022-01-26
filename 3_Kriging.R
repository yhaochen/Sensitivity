# Kriging emulator
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
folder<-paste("/storage/work/h/hxy46/Sensitivity/Data/",d,"D/Kriging",sep="")
dir.create(file.path(folder), showWarnings = FALSE)

Kriging_step<-91
Kriging_size<-c(10:100)
save(Kriging_size,file = paste(folder,"/Kriging_size",sep=""))
Kriging_First_order<-matrix(NA,nrow=Kriging_step,ncol=d)
# Record two times, one for Kriging emulator, one for Sobol based on this emulator
Kriging_time<-rep(NA,Kriging_step)
Kriging_Soboltime<-rep(NA,Kriging_step)
Kriging_fail_prob<-rep(NA,Kriging_step) # estimated failure probability

# For emulator based method, also record their wrong classification 
Kriging_wrong_pos<-rep(NA,Kriging_step) # % of wrong success (should be failure)
Kriging_wrong_neg<-rep(NA,Kriging_step) # % of wrong failure (should be success)

#Kriging function, used in Sobol analysis based on Kriging emulator
Kriging <- function (X){
  a<-predict(GPmodel,X)
  a$Y_hat
}

for (i in 1:Kriging_step){
  n <- Kriging_size[i]
  x <- X1[c(1:n), ]
  y <- Reliability(x)
  
  start.time <- Sys.time()
  GPmodel <- GP_fit(x, y)
  yfit<-Kriging(X1)
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  Kriging_time[i]<-as.numeric(time.taken)
  
  Kriging_fail_prob[i]<-length(which(yfit<0))/N*100
  Kriging_wrong_pos[i]<-length(which(yfit[true_neg_indx]>0))/length(true_neg_indx)*100
  Kriging_wrong_neg[i]<-length(which(yfit[-true_neg_indx]<0))/(N-length(true_neg_indx))*100
  
  start.time <- Sys.time()
  S<-sensitivity::sobol(model = Kriging, X1 = X1, X2 = X2,
                        order=1, nboot = 100)
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  Kriging_Soboltime[i]<-as.numeric(time.taken)
  
  Kriging_First_order[i, ] <- S$S$original[1:d]
  save(Kriging_time,file = paste(folder,"/Kriging_time",sep=""))
  save(Kriging_Soboltime,file = paste(folder,"/Kriging_Soboltime",sep=""))
  save(Kriging_fail_prob,file = paste(folder,"/Kriging_fail_prob",sep=""))
  save(Kriging_wrong_neg,file = paste(folder,"/Kriging_wrong_neg",sep=""))
  save(Kriging_wrong_pos,file = paste(folder,"/Kriging_wrong_pos",sep=""))
  save(Kriging_First_order,file = paste(folder,"/Kriging_First_order",sep=""))
  print(i)
}
