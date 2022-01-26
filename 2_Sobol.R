# Sobol analysis
rm(list = ls())
graphics.off()
setwd("/storage/work/h/hxy46/Sensitivity")
library(sensitivity)
library(lhs)
# Load problem definition
source("1_Problem_definition.R")
set.seed(1)
# Create a folder for d dimension test scenario if not yet created
folder<-paste("/storage/work/h/hxy46/Sensitivity/Data/",d,"D/Sobol",sep="")
dir.create(file.path(folder), showWarnings = FALSE)

# Define total steps for Sobol, then record computational time, first-order Sobol sensitivity index
#     estimated failure probability for each step (each sample size)
Sobol_step<-2000
Sobol_size<-seq(N/Sobol_step,N,by=N/Sobol_step)
save(Sobol_size,file = paste(folder,"/Sobol_size",sep=""))
Sobol_First_order<-matrix(NA,nrow=Sobol_step,ncol=d)
Sobol_time<-rep(NA,Sobol_step)
LHS_fail_prob<-rep(NA,Sobol_step) # estimated failure probability by LHS
for (i in 1:Sobol_step){
  start.time <- Sys.time()
  N <- Sobol_size[i]
  S<-sensitivity::sobol(model = Reliability, X1 = X1[1:N, ], X2 = X2[1:N, ],
                        order=1, nboot = 100)
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  Sobol_time[i]<-as.numeric(time.taken)
  Sobol_First_order[i, ] <- S$S$original[1:d]
  LHS_fail_prob[i]<- length(which(y[1:N]<thres))/N*100
  save(Sobol_time,file = paste(folder,"/Sobol_time",sep=""))
  save(LHS_fail_prob,file = paste(folder,"/LHS_fail_prob",sep=""))
  save(Sobol_First_order,file= paste(folder,"/Sobol_First_order",sep=""))
  print(i)
}

