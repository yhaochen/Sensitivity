# Sensitivity comparison
# remove all existing variables and plots
rm(list = ls())
graphics.off()
# load the required libraries
library(lhs) #Generate Latin Hypercube samples
library(fields) # Make grid plots
library(sensitivity)
library(randtoolbox) #Generate Sobol sequence (Quasi MC)
library(gstat)
library(sp)

Reliability<-function (X) {
  floor(X[ ,1]*X[ ,2]/0.75)
}

Kriging <- function (X){
  coordinates(X) <- ~ X1+X2
  KrigSurface <- predict(TheGStat, model=FittedModel, newdata=X)
  mean <- KrigSurface$Rel.pred
}

N <- c(5000,10000,25000,50000,100000,250000,500000)
N_K <- c(10,25,50,100,250,500,1000)
S1_Sobol <- S2_Sobol <- S1_Morris <- S2_Morris <- S1_Kriging <- S2_Kriging <- rep(NA,length(N))

for (i in 1: length(N)){
  #Sobol
  n <- floor(N[i]/4)
  X <- data.frame(randomLHS(n,2))
  #Randomly divide into two samples
  ind <- sample(c(rep(TRUE,n/2), rep(FALSE,n/2)), replace=FALSE)
  X1 <- X[ind, ]
  X2 <- X[!ind, ]
  S<-sensitivity::sobol(model = Reliability, X1 = X1, X2 = X2,
                        order=2, nboot = 10)
  S1_Sobol[i]<-S$S$original[1]
  S2_Sobol[i]<-S$S$original[2]
  
  #Morris
  n <- floor(N[i]/3)
  M <- morris(model = Reliability, factors = 2, r = n,
              design = list(type = "oat", levels = 1000, grid.jump = 100))
  S1_Morris[i] <- apply(M$ee, 2, sd)[1]
  S2_Morris[i] <- apply(M$ee, 2, sd)[2]
  
  #Kriging
  n <- N_K[i]
  x <- data.frame(randomLHS(n,2))
  y <- Reliability(x)
  Data <- data.frame(X1=x$X1,X2=x$X2,y=y)
  coordinates(Data)=~X1+X2
  TheVariogram=variogram(y~1, data=Data)
  TheVariogramModel <- vgm(psill=0.025, model="Exp", nugget=0.001, range=0.5)
  FittedModel <- fit.variogram(TheVariogram, model=TheVariogramModel)
  TheGStat <- gstat(id="Rel", formula=y ~ 1, data=Data,model=FittedModel )
  
  #Randomly divide into two samples
  X <- data.frame(randomLHS(100000,2))
  ind <- sample(c(rep(TRUE,50000), rep(FALSE,50000)), replace=FALSE)
  X1 <- X[ind, ]
  X2 <- X[!ind, ]
  S<-sensitivity::sobol(model = Kriging, X1 = X1, X2 = X2,
                        order=2, nboot = 100)
  S1_Kriging[i] <- S$S$original[1]
  S2_Kriging[i] <- S$S$original[2]
}

level <- seq(0, 1, length=50)
X_ALL <- expand.grid(level,level)
Y_ALL <- Kriging(X_ALL)
colfunc <- colorRampPalette(c("white", "red"))
cols <- colfunc(100)
image.plot(level,level, matrix(Y_ALL,ncol=50),xaxt="n",yaxt="n", xlab="x1", ylab="x2", col=cols,cex.lab=1.5,cex.axis=1.5)
