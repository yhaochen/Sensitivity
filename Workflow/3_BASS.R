# Sobol based on the Bayesian Adaptive Spline Surface (BASS) method
# Note: the Sobol analysis in this script directly uses the function in "BASS" package instead of sensobol package

# This script also takes a time when dealing with high-dimensional models (but not extremely long).
#       Again change the dimension vector D for code replication check.
# Remove all existing environment and plots
rm(list = ls())
graphics.off()

# Set a working directory, please set it to your own working folder when testing
dir <- "/storage/work/h/hxy46/Sensitivity"
setwd(dir = dir)

# Load the required packages
library(BASS)
library(lhs)

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

# Define the model in each dimension and apply BASS method
for (k in 1:length(D)){
  # Model dimension
  d=D[k]
  S_total <- rep(NA,d)

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
  
  # Folder for d dimension test scenario
  folder <- paste("./Example_Data/",d,"D/BASS",sep="")
  if (!dir.exists(folder)){
    dir.create(file.path(folder), showWarnings = FALSE)
  }
    
  # Use 20,000 LHS training data points to test emulator quality
  x_test <- randomLHS(20000,d)
  
  # Similar to Kriging method, we begin with 10 times model dimension samples
  while (1>0) {
    sample_size <- 10*d
    X <- randomLHS(sample_size, d)
    Y <- apply(X, 1, Reliability)
    
    # Use a MCMC size of 500,000, burn-in period of 100,000, record the output every 1,000 steps
    mcmc_size <- 500000
    # Record the time of BASS emulation
    start.time <- Sys.time()
    mod <- bass(X, Y, nmcmc = mcmc_size, nburn = 100000, thin = 1000) # fit BASS model
    end.time <- Sys.time()
    T_BASS <- difftime(end.time,start.time, units = "secs")
    
    # Check emulator quanlity by the same stopping criterion: stop sampling if all the training data's standard
    #     errors are less than one
    P <- predict(mod, x_test, verbose = TRUE)
    V <- apply(P, 2, var)
    
    # If still need to take more samples, add the sample size by d
    if (max(V)>1){
      sample_size <- sample_size + d
    }
    # otherwise perform sensitivity analysis and record the time
    else{
      start.time <- Sys.time()
      S_BASS <- sobol(mod, verbose = FALSE)
      end.time <- Sys.time()
      T_BASSSobol <- difftime(end.time,start.time,units = "secs")
      
      # Because we are not using sensobol package here, we need an additional convergence check for the Sobol
      #     analysis. Usually if the emulator's quality is good enough, the sensitivity analysis results should
      #     be stable. We add a convergence check here to make sure the results are reliable. If the results are 
      #     not converged, add the BASS sample size by d.
      for (j in 1:d){
        S_total[j] <- quantile(S_BASS$T[ ,j],probs = 0.975) - quantile(S_BASS$T[ ,j],probs = 0.025)
      }
      
      # If converged, save the emulation time, sensitivity analysis time, sensitivity indices and the sample size
      if (all(S_total <= 0.05)){
        save(T_BASS, file = paste(folder, "/T_BASS", sep=""))
        save(T_BASSSobol,file = paste(folder,"/T_BASSSobol",sep=""))
        save(S_BASS,file = paste(folder,"/S_BASS",sep=""))
        save(sample_size,file = paste(folder,"/BASS_size",sep=""))
        break
      } else {
        sample_size <- sample_size + d
      }
    }
  }
}
