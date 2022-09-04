# Sobol based on the Kriging emulator
# Note: this script can take an extremely long running time. This is because high-dimensional test problems take
#       a long time to build the Kriging emulator, and we need to search for the required sample size 
#       (which means we need to repeat the emulation process). It would be more efficient to separate high
#       dimension tests by modifying the vector "D".

# For code check, keep only the first two terms in "D" (avoid high-dimensional runs) 
#       so that the total time is acceptable. (see the block with dashed lines below)

# Remove all existing environment and plots
rm(list = ls())
graphics.off()

# Set a working directory, please set it to your own working folder when testing
dir <- "/storage/work/h/hxy46/Sensitivity"
setwd(dir = dir)

# Load the required packages
library(GPfit)
library(lhs)
library(sensobol)

# Set a random seed
set.seed(1)

# The dimensions we consider
D <- c(2,5,10,15,20,30)

#-------------------------- edit the code in this block if you want to rerun the code for a check
# A switch for code check, set the number to "1" if for code check, otherwise set to "0" or other numbers
Check_switch <- 1
if (Check_switch == 1){
  D <- c(2,5)
}
#-----------------------------------------

#Kriging function, used when calling the emulator in the Sobol analysis
Kriging <- function (X){
  a <- predict(GPmodel,X)
  a$Y_hat
}

# Define the test model in each dimension, build the Kriging emulator and perform the Sobol analysis
for (k in 1:length(D)){
  # model dimension
  d <- D[k]
  # model
  Reliability<-function (X) {
    S = 0
    for (i in 1:d){
      if (i%%3==1){
        a = 1
        b = X[i]
      } else if (i%%3==2){
        a = 10
        b = X[i]^2*X[i-1]
      } else {
        a = 100
        b = X[i]
      }
      S = S + a*b
    }
    return(S)
  }
  # Begin with 10*dimension as a base sample size
  Kriging_size <- 10*d
  
  # The training samples for emulator quality control
  # Select 20,000 samples by Latin Hypercube Sampling (LHS)
  x_test <- randomLHS(20000,d)
  
  # Folder for d dimension test scenario
  folder <- paste("Example_Data/",d,"D/Kriging",sep="")
  dir.create(folder, recursive = TRUE)
  
  # A while loop includes the stopping criterion
  while (1>0){
    # Begin with the base sample size (also sampled by LHS)
    X_GP <- randomLHS(Kriging_size,d)
    # Evaluate the model to get their true outputs
    Y_GP <- apply(X_GP,1,Reliability)
    
    # Record the Kriging model evaluation time
    start.time <- Sys.time()
    GPmodel <- GP_fit(X_GP,Y_GP)
    end.time <- Sys.time()
    
    # Emulator convergence check
    
    # First get the standard error of all the training samples
    a <- predict(GPmodel,x_test)
    std <- sqrt(a$MSE)
    
    # If all the standard errors are less than or equal to 1, then convergence is reached, else add sample size
    if (max(std)<=1){
      break
    } else {
      Kriging_size <- Kriging_size + d
    }
  }
  
  # Time for building the emulator
  T_Kriging<-difftime(end.time,start.time,units = "secs")
  
  # Save this time and the required sample size
  save(T_Kriging,file = paste(folder,"/T_Kriging",sep=""))
  save(Kriging_size,file = paste(folder,"/Kriging_size",sep=""))
  
  # Then perform the sensitivity analysis
  # The sample size in sensitivity analysis is directly estimated by the required sample size of the original
  #     model for convenience.
  load(paste("/Example_Data/Sobol_convergencesize",sep = ""))
  N <- Sobol_convergesize[k]/(d+2+d*(d-1)/2)
  
  # Record the time for sensitivity analysis
  start.time<-Sys.time()
  mat <- sobol_matrices(N = N, params = as.character(c(1:d)), order = "second")
  Y_S <- Kriging(mat)
  
  # The sensitivity indices
  S_Kriging <- sobol_indices(Y=Y_S,N=N,params = as.character(c(1:d)),
                     boot=TRUE,R=100,order="second")
  end.time <- Sys.time()
  
  # Save the results
  T_KrigingSobol <- difftime(end.time,start.time,units = "secs")
  save(T_KrigingSobol,file = paste(folder,"/T_KrigingSobol",sep=""))
  save(S_Kriging,file=paste(folder,"/S_Kriging",sep=""))
}
