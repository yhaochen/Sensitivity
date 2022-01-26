# Record the computational time of each method under different dimensions (model evaluation time is very short)
rm(list = ls())
graphics.off()
setwd("/storage/work/h/hxy46/Sensitivity")
library(sensitivity)
library(lhs)
library(GPfit)
set.seed(1)
# Create a folder for time record if not yet created
folder<-"/storage/work/h/hxy46/Sensitivity/Data/Time"
dir.create(file.path(folder), showWarnings = FALSE)
N <- 50000

#Kriging function, used in Sobol analysis based on Kriging emulator
Kriging <- function (X){
  a<-predict(GPmodel,X)
  a$Y_hat
}

# Record computational time under different dimensions
step<-30
Sobol_time<-rep(NA,step)
Kriging_time<-rep(NA,step)
for (d in 2: step) {
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
  X2 <- data.frame(matrix(runif(d*N), nrow = N))
  y<- Reliability(X1)
  # Define a success/failure threshold as 2% of lower tail 
  q <- 0.02
  thres<-quantile(y,q)
  
  # Sobol time:
  start.time <- Sys.time()
  S<-sensitivity::sobol(model = Reliability, X1 = X1[1:N, ], X2 = X2[1:N, ],
                        order=1, nboot = 100)
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  Sobol_time[d]<-as.numeric(time.taken)
  save(Sobol_time,file = paste(folder,"/Sobol_time",sep=""))
  # Kriging time: Assume we take 50 initial samples and fit X1 (50000 samples)
  x <- X1[c(1:50), ]
  y <- Reliability(x)
  start.time <- Sys.time()
  GPmodel <- GP_fit(x, y)
  yfit<-Kriging(X1)
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  Kriging_time[d]<-as.numeric(time.taken)
  save(Kriging_time,file = paste(folder,"/Krging_time",sep=""))
} 