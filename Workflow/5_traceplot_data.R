# This script prepares the required data of Figure 2 in the paper

# Remove all existing environment and plots
rm(list = ls())
graphics.off()

# Set a working directory, please set it to your own working folder when testing
dir <- "/storage/work/h/hxy46/Sensitivity"
setwd(dir = dir)

# Load the required packages
library(sensobol)
library(lhs)
library(GPfit)
library(BASS)

# Plot the results in the example of 5D problem
d<-5
folder<-paste("Example_Data/",d,"D",sep="")

# Test model in 5D
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

#Kriging function, used in Kriging methods
Kriging <- function (X){
  a<-predict(GPmodel,X)
  a$Y_hat
}

# Converged size of the 5D test model
load("Data/Sobol_convergencesize")
N_K <- Sobol_convergesize[2]/(d+2+d*(d-1)/2)

# Sizes used for Sobol trace plot
# Standard Sobol's results do not vary by random seeds, other methods' results do. 
Size_S <- c(seq(40,1000,by=40),seq(1000,100000,by=1000))
# Upper and lower bound of 95% CI of the total-order indices
T_S_high <-rep(NA, length(Size_S))
T_S_low <- rep(NA, length(Size_S))
T_S <- rep(NA, length(Size_S))

# Evaluate the results of the standard Sobol method
for (i in 1:length(Size_S)){
  N <- floor(Size_S[i]/(d+2+d*(d-1)/2))
  mat <- sobol_matrices(N = N, params = as.character(c(1:d)), order = "second")
  X <- split(t(mat), rep(1:dim(mat)[1], each = d))
  Y <- sapply(X, Reliability)
  S <- sobol_indices(Y=Y,N=N,params = as.character(c(1:d)),
                     boot=TRUE,R=100,order="second")
  T_S_high[i] <- S$results$high.ci[8]
  T_S_low[i] <- S$results$low.ci[8]
  
  T_S[i] <- S$results$original[8]
}
save(T_S,file = paste(folder,"/Traceplot/T_S",sep=""))
save(T_S_high,file = paste(folder,"/Traceplot/T_S_high",sep=""))
save(T_S_low,file = paste(folder,"/Traceplot/T_S_low",sep=""))

# Sizes used for Kriging traceplot
# 5 different seeds for the emulation methods
Size_K <- seq(20,80,by=2)
# Total-order indices
T_K <- matrix(NA,nrow=5,ncol=length(Size_K))

# Sizes used for AKMCS traceplot
Size_A <- seq(13,50,by=1)
T_A <- matrix(NA,nrow=5,ncol=length(Size_A))

# Sizes used for BASS traceplot
Size_B <- seq(20,100,by=5)
T_B <- matrix(NA,nrow=5,ncol=length(Size_B))

# Test for 5 seeds for the emulation methods
for (seed in 1:5){
  set.seed(seed)
  
  # Kriging
  for (i in 1:length(Size_K)){
    # Take the corresponding sample size, fit a Kriging model and perform Sobol analysis
    X_GP <- randomLHS(Size_K[i],d)
    Y_GP <- apply(X_GP,1,Reliability)
    GPmodel <- GP_fit(X_GP,Y_GP)
    mat <- sobol_matrices(N = N_K, params = as.character(c(1:d)), order = "second")
    Y_K <- Kriging(mat)
    S <- sobol_indices(Y=Y_K,N=N_K,params = as.character(c(1:d)),
                               boot=TRUE,R=100,order="second")
    T_K[seed,i] <- S$results$original[8]
  }
  # Save the results of the best estimates
  save(T_K,file = paste(folder,"/Traceplot/T_K",sep=""))
  
  # AKMCS records the results after each step
  # 20,000 training points
  candidate_size <- 20000
  X <- randomLHS(candidate_size,d)
  # Begin with 12 initial samples
  n_init<-12
  indx <- sample(candidate_size,n_init)
  x <- X[indx, ]
  x_rest <- X[-indx, ]
  # Evaluate their outputs, fit a Kriging model
  y <- apply(x,1,Reliability)
  GPmodel <- GP_fit(x, y)
  a<-predict(GPmodel,x_rest)
  # Learning function
  U<-sqrt(a$MSE)
  # Take the next sample adaptively, repeat these steps
  for (i in 1:length(Size_A)){
    k<-which(U==max(U))
    if (length(k)>1){
      k <- sample(k,1)
    }
    x_add<-x_rest[k, ]
    x_rest<-x_rest[-k, ]
    y_add<-Reliability(x_add)
    y <- append(y,y_add)
    x <- rbind(x,x_add)
    GPmodel <- GP_fit(x, y)
    a<-predict(GPmodel,x_rest)
    U<-sqrt(a$MSE)
    
    # Save the sensitivity analysis results
    mat <- sobol_matrices(N = N_K, params = as.character(c(1:d)), order = "second")
    Y_A <- Kriging(mat)
    S <- sobol_indices(Y=Y_A,N=N_K,params = as.character(c(1:d)),
                               boot=TRUE,R=100,order="second")
    T_A[seed,i] <- S$results$original[8]
  }
  save(T_A,file = paste(folder,"/Traceplot/T_A",sep=""))
  
  # BASS method simply takes different initial sample sizes
  for (i in 1:length(Size_A)){
    X <- randomLHS(Size_B[i], d)
    Y <- apply(X, 1, Reliability)
    mcmc_size <- 500000
    mod <- bass(X, Y, nmcmc = mcmc_size, nburn = 100000, thin = 1000) # fit BASS model
    S_BASS <- sobol(mod, verbose = FALSE)
    # Take the mean of the ensemble as the best estimate
    T_B[seed,i] <- mean(S_BASS$T[ ,3])
  }
  save(T_B,file = paste(folder,"/Traceplot/T_B",sep=""))
}
