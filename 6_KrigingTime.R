# This script records the computational time of Kriging under different dimensions and sample sizes
# Record the computational time of each method under different dimensions (model evaluation time is very short)
rm(list = ls())
graphics.off()
setwd("/storage/work/h/hxy46/Sensitivity")
library(sensitivity)
library(lhs)
library(GPfit)
set.seed(1)
# Create a folder for time record if not yet created
folder<-"/storage/work/h/hxy46/Sensitivity/Data/KrigingTime"
dir.create(file.path(folder), showWarnings = FALSE)

#Kriging function, used in Sobol analysis based on Kriging emulator
Kriging <- function (X){
  a<-predict(GPmodel,X)
  a$Y_hat
}
N<-50000
# Record computational time under different dimensions and different samples
Dim <- c(2:15)
Sample <- seq(10,100,by=10)
Kriging_time<-matrix(NA,nrow = length(Dim),ncol = length(Sample))
for (j in 1:length(Dim)) {
  d <- Dim[j]
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
  X1 <- data.frame(matrix(runif(d*N), nrow = N))
  
  # Kriging time: assume we have altogether 50000 points need to be fitted
  for (k in 1:length(Sample)){
    x <- X1[c(1:Sample[k]), ]
    y <- Reliability(x)
    start.time <- Sys.time()
    GPmodel <- GP_fit(x, y)
    yfit<-Kriging(X1)
    end.time <- Sys.time()
    time.taken <- difftime(end.time,start.time,units = "secs")
    Kriging_time[j,k]<-as.numeric(time.taken)
    save(Kriging_time,file = paste(folder,"/Kriging_time",sep = ""))
  }
} 