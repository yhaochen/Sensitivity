# Sobol based on the Adaptive Kriging combined with Monte Carlo Sampling (AKMCS) method
# Note: This script also takes an extremely long time when dealing with high-dimensional
#       models. Again change the dimension vector D for code replication check. 

# Remove all existing environment and plots
rm(list = ls())
graphics.off()

# Set a working directory, please set it to your own working folder when testing
dir <- "/storage/work/h/hxy46/Sensitivity"
setwd(dir = dir)

# Load the required packages
library(lhs)
library(GPfit)
library(sensobol)

# Set a random seed
set.seed(1)

# The dimensions we consider
D <- c(2,5,10,15,20,30)

#-------------------------- edit the code in this block if you want to rerun the code for a check
# A switch for code check, set the number "1" if for code check, otherwise set to "0" or other numbers
Check_switch <- 1
if (Check_switch == 1){
  D <- c(2,5)
}
#-----------------------------------------

# Define the test model in each dimension, apply AKMCS and perform the Sobol analysis
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

  #Kriging function
  Kriging <- function (X){
    a <- predict(GPmodel,X)
    a$Y_hat
  }

  folder<-paste("./Example_Data/",d,"D/AKMCS",sep="")
  if (!dir.exists(folder)){
    dir.create(folder, recursive = TRUE)
  }
    
  # Start recording the time from AKMCS initial state
  # AKMCS also begins with 20,000 training samples
  start.time <- Sys.time()
  candidate_size <- 20000
  X <- randomLHS(candidate_size,d)
  
  # Save these training samples
  save(X,file = paste(folder,"/initial_sample",sep=""))
  
  # Begin with 12 random samples from these training samples
  n_init <- 12
  indx <- sample(candidate_size,n_init)
  
  # Update the used samples and remaining samples
  x <- X[indx, ]
  x_rest <- X[-indx, ]
  
  # Evaluate model outputs and fit a Kriging model
  y <- apply(x,1,Reliability)
  GPmodel <- GP_fit(x, y)
  a <- predict(GPmodel,x_rest)
  
  # U is the learning function, which is simply the standard error here
  U <- sqrt(a$MSE)
  
  # Stopping criterion: all the remaining samples have standard errors larger than 1
  # If the criterion is not reached, pick the next sample adaptively based on the learning function
  while (1>0){
    # Find which sample has the largest standard error
    m <- which(U==max(U))
    if (length(m)>1){
      m <- sample(m,1)
    }
    
    # Add that sample and update the remaining samples
    x_add <- x_rest[m, ]
    x_rest <- x_rest[-m, ]
    
    # Evaluate the output of that sample and update
    y_add <- Reliability(x_add)
    y <- append(y,y_add)
    x <- rbind(x,x_add)
    
    # Fit the Kriging model again
    GPmodel <- GP_fit(x, y)
    a <- predict(GPmodel,x_rest)
    
    # Get the learning function again
    U <- sqrt(a$MSE)
    
    # End the loop if the stopping criterion is fulfilled
    if (max(U)<=1){
      break
    }
  }
  
  # Then record the total used time
  end.time <- Sys.time()
  time.taken <- difftime(end.time,start.time,units = "secs")
  T_AKMCS <- as.numeric(time.taken)
  save(T_AKMCS,file = paste(folder,"/T_AKMCS",sep=""))
  
  # Next perform the sensitivity analysis, and again directly get the convergence size of standard Sobol
  load(paste("./Example_Data/Sobol_convergencesize",sep = ""))
  N <- floor(Sobol_convergesize[k]/(d+2+d*(d-1)/2))
  
  # Time for sensitivity analysis
  start.time <- Sys.time()
  mat <- sobol_matrices(N = N, params = as.character(c(1:d)), order = "second")
  Y_S <- Kriging(mat)
  S_AKMCS <- sobol_indices(Y=Y_S,N=N,params = as.character(c(1:d)),
                           boot=TRUE,R=100,order="second")
  end.time<-Sys.time()
  T_AKMCSSobol<-difftime(end.time,start.time,units = "secs")
  
  # Save the time and sensitivity indices
  save(T_AKMCSSobol,file = paste(folder,"/T_AKMCSSobol",sep=""))
  save(S_AKMCS,file=paste(folder,"/S_AKMCS",sep=""))
}