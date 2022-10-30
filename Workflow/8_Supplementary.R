# This script makes the supplementary figures
# Note: this script requires the full data. 

# Remove all existing environment and plots
rm(list = ls())
graphics.off()

# Set a working directory, please set it to your own working folder when testing
dir <- commandArgs(trailingOnly=TRUE)
setwd(dir = dir)

# Load the required packages
library(plot.matrix)
library(RColorBrewer)

# Tested dimension, method names, and evaluation time
tested_D_num <- c(2,5,10,15,20,30)
tested_D <- c("2D","5D","10D","15D","20D","30D")
tested_M <- c("Kriging","AKMCS","BASS")
tested_eval_time <- c(10^(-6),10^(-5),0.0001,0.001,0.01,0.1,1,10,60,3600,3600*6,3600*12,3600*24)

# Total time of each method & each scenario
Time_Sobol <- matrix(NA,nrow = length(tested_D),ncol = length(tested_eval_time))
Time_Kriging <- matrix(NA,nrow = length(tested_D),ncol = length(tested_eval_time))
Time_AKMCS <- matrix(NA,nrow = length(tested_D),ncol = length(tested_eval_time))
Time_BASS <- matrix(NA,nrow = length(tested_D),ncol = length(tested_eval_time))

# Load all related results for each test scenario
for (i in 1:length(tested_D)) {
  folder <- paste("./Data/",tested_D[i],sep="")
  
  # General:
  load(paste(folder,"/avg_eval_time",sep=""))
  
  # Sobol:
  load(paste(folder,"/Sobol/T_Sobol",sep=""))
  load(paste(folder,"/Sobol/S_Sobol",sep=""))
  Sobol_convergesize <- S$C
  
  # Kriging:
  load(paste(folder,"/Kriging/T_Kriging",sep=""))
  load(paste(folder,"/Kriging/T_KrigingSobol",sep=""))
  load(paste(folder,"/Kriging/Kriging_size",sep=""))
  
  # AKMCS:
  load(paste(folder,"/AKMCS/AKMCS_time",sep=""))
  load(paste(folder,"/AKMCS/x",sep=""))
  AKMCS_size <- dim(x)[1]
  load(paste(folder,"/AKMCS/T_AKMCSSobol",sep=""))
  
  # BASS:
  load(paste(folder,"/BASS/T_BASS",sep=""))
  load(paste(folder,"/BASS/T_BASSSobol",sep=""))
  load(paste(folder,"/BASS/BASS_size",sep=""))
  BASS_size <- sample_size
  
  # Calculation of the computational time in each scenario
  # Sensitivity analysis time + model evaluation adjusted time + emulation time
  for (j in 1:length(tested_eval_time)) {
    Time_Sobol[i,j] <- T_Sobol + (tested_eval_time[j]-avg_eval_time)*Sobol_convergesize
    Time_Kriging[i,j] <- T_Kriging + (tested_eval_time[j]-avg_eval_time)*Kriging_size+ T_KrigingSobol
    Time_AKMCS[i,j] <- T_AKMCS + (tested_eval_time[j]-avg_eval_time)*AKMCS_size + T_AKMCSSobol
    Time_BASS[i,j] <- T_BASS + (tested_eval_time[j]-avg_eval_time)*BASS_size + T_BASSSobol
  }
}

# An adjustment for standard Sobol's time: for fast test models, it is likely the time becomes negative due
#     to random noises in the time recording.
Time_Sobol[Time_Sobol < 0] <- Time_Sobol[1,1]

PercentTime_Kriging <- Time_Sobol/Time_Kriging
PercentTime_AKMCS <- Time_Sobol/Time_AKMCS
PercentTime_BASS <- Time_Sobol/Time_BASS

# Label of evaluation time
eval_time_lab <- c("1us","10us","0.1ms","1ms","10ms","0.1s","1s","10s","1min","1h","6h","12h","1d")

# Create a folder to save figures
folder <- "./Figures"
if (!dir.exists(folder)){
  dir.create(folder, recursive = TRUE)
}

#--------------------------------------
# Errors in sensitivity index of each method

# Sizes used for standard Sobol traceplot
Size_S <- c(seq(40,1000,by=40),seq(1000,100000,by=1000))
# Sizes used for Kriging traceplot
Size_K <- seq(2,80,by=2)
# Sizes used for AKMCS traceplot
Size_A <- seq(5,50,by=1)
# Sizes used for BASS traceplot
Size_B <- seq(2,100,by=2)

# Load the saved sensitivity indices from script 5
load("./Data/5D/Traceplot/T_S")
load("./Data/5D/Traceplot/T_S_high")
load("./Data/5D/Traceplot/T_S_low")

load("./Data/5D/Traceplot/T_K")
load("./Data/5D/Traceplot/T_B")
load("./Data/5D/Traceplot/T_A")

# Load the convergence size of the 5D test model
load("./Data/Sobol_convergencesize")
C_S <- Sobol_convergesize[2]
load("./Data/5D/Kriging/Kriging_size")
C_K <- Kriging_size
load("./Data/5D/AKMCS/x")
C_A <- dim(x)[1]
load("./Data/5D/BASS/BASS_size")
C_B <- sample_size

# Error in sensitivity index 
# For standard Sobol, take the last result (when converged) as the "true" result
# For emulation-based models, take the mean of 5 random seeds as the result under each sample size.
E_K <- colMeans(T_K) - T_S[length(T_S)]
E_A <- colMeans(T_A) - T_S[length(T_S)]
E_B <- colMeans(T_B) - T_S[length(T_S)]
pdf(file = "./Figures/Figure_S1.pdf",width = 11,height = 7)
par(mar=c(5,5,2.6,3))
plot(Size_K,E_K,type="l",col="purple",lwd=1.5,xlab="Sample size",
     ylab = "Error in sensitivity index",xlim = c(min(c(Size_A,Size_B,Size_K)),max(c(Size_A,Size_B,Size_K))),
     ylim = c(min(c(E_A,E_B,E_K)),max(c(E_A,E_B,E_K))),cex.lab=1.5,cex.axis=1.5)
lines(Size_B,E_B,col="red",lwd=1.5)
lines(Size_A,E_A,col="blue",lwd=1.5)
abline(h = 0, lty=2,lwd=0.5)
legend("bottomright",lty=c(1,1,1),lwd=c(1.5,1.5,1.5),col=c("purple","blue","red"),bty="n",
       legend = c("Kriging","AKMCS","BASS"),cex=1.5)
dev.off()



#-------------------------------------
#Grid plot of Kriging gain
Mat = floor(log10(PercentTime_Kriging))
Cols<-brewer.pal(n = 11, name = "RdYlGn")
rownames(Mat)<-tested_D_num
colnames(Mat)<-eval_time_lab
pdf(file = "./Figures/Figure_S2.pdf",width = 14,height = 7)
par(mar=c(5,5,2.6,7))
plot(Mat[nrow(Mat):1, ],breaks = c(-7:4),col=Cols,xlab="Time of single run",
     ylab="Number of input parameters",digits = 1,fmt.cell='%.0f',max.col=100,
     cex.lab=1.5,cex.axis=1,main="",fmt.key = "%+.0f")
mtext('Order of magnitude', side=3, line=1, at=15)
mtext('Speed gain of Kriging',side=3, line=-0.5, at=15)
dev.off()

#-------------------------------------
#Grid plot of AKMCS gain
Mat = floor(log10(PercentTime_AKMCS))
Cols<-brewer.pal(n = 11, name = "RdYlGn")
rownames(Mat)<-tested_D_num
colnames(Mat)<-eval_time_lab
pdf(file = "./Figures/Figure_S3.pdf",width = 14,height = 7)
par(mar=c(5,5,2.6,7))
plot(Mat[nrow(Mat):1, ],breaks = c(-8:4),col=c("orangered4",Cols),xlab="Time of single run",
     ylab="Number of input parameters",digits = 1,fmt.cell='%.0f',max.col=100,
     cex.lab=1.5,cex.axis=1,main="",fmt.key = "%+.0f")
mtext('Order of magnitude', side=3, line=1, at=15)
mtext('Speed gain of AKMCS',side=3, line=-0.5, at=15)
dev.off()

#-------------------------------------
#Grid plot of BASS gain
Mat = floor(log10(PercentTime_BASS))
Cols<-brewer.pal(n = 8, name = "RdYlGn")
rownames(Mat)<-tested_D_num
colnames(Mat)<-eval_time_lab
pdf(file = "./Figures/Figure_S4.pdf",width = 14,height = 7)
par(mar=c(5,5,2.6,7))
plot(Mat[nrow(Mat):1, ],breaks = c(-4:4),col=Cols,xlab="Time of single run",
     ylab="Number of input parameters",digits = 1,fmt.cell='%.0f',max.col=100,
     cex.lab=1.5,cex.axis=1,main="",fmt.key = "%+.0f")
mtext('Order of magnitude', side=3, line=1, at=15)
mtext('Speed gain of BASS',side=3, line=-0.5, at=15)
dev.off()

#----------------------------------------------
# 
# Speed gain compared with brute force Sobol method (Figure S4)
Mat <- Time_Sobol
textMat <- matrix(NA,nrow=nrow(Mat),ncol=ncol(Mat))
for (i in 1:length(tested_D)) {
  for (j in 1:length(tested_eval_time)){
    Mat[i,j] <- min(Time_Sobol[i,j],Time_Kriging[i,j],
                    Time_AKMCS[i,j],Time_BASS[i,j],na.rm=TRUE)
    if (Mat[i,j]==Time_Sobol[i,j]){
      textMat[i,j] <- "S"
    } 
    if (Mat[i,j]==Time_Kriging[i,j]){
      textMat[i,j] <- "K"
    } 
    if (Mat[i,j]==Time_AKMCS[i,j]){
      textMat[i,j] <- "A"
    } 
    if (Mat[i,j]==Time_BASS[i,j]){
      textMat[i,j] <- "B"
    }
  }
}
Mat <- Time_Sobol/Mat
Mat <- floor(log10(Mat))
rownames(Mat)<-tested_D_num
colnames(Mat)<-eval_time_lab
Cols<-brewer.pal(n = 4, name = "Reds")
pdf(file = "./Figures/Figure_S5.pdf",width = 14,height = 7)
par(mar=c(5,5,2.6,6))
time<-plot(Mat[nrow(Mat):1, ],breaks=c(0:4),digits=1,fmt.cell='%.0f',
           xlab="Time of single run",ylab="Number of input parameters",col=Cols,
           text.cell=list(pos=3, cex=1),max.col=100,
           cex.lab=1.5,cex.axis=1,main="",fmt.key = "%.0f")
for (i in 1:nrow(Mat)) {
  for (j in 1:ncol(Mat)) {
    args<-time$cell.text[[nrow(Mat)+1-i,j]]
    args$labels <- textMat[i,j]
    args$cex    <- 1.5
    args$pos    <- 1
    do.call(text, args)
  }
}
mtext("Orders of",side = 3, at = 15)
mtext("magnitude",side = 3, at = 15, line = -1)
dev.off()

#-------------------
# speed gain of best method compared with the second best method (Figure S5)
Mat<-Time_Sobol
textMat <- matrix(NA,nrow=nrow(Mat),ncol=ncol(Mat))
for (i in 1:nrow(Mat)){
  for (j in 1:ncol(Mat)){
    vec <- c(Time_Sobol[i,j],Time_Kriging[i,j],Time_AKMCS[i,j],Time_BASS[i,j])
    Text <- c("S","K","A","B")
    vec <- vec[!is.na(vec)]
    r <- order(vec)
    Mat[i,j]<-vec[r[2]]/vec[r[1]]
    textMat[i,j]<-paste(Text[r[1]],"(",Text[r[2]],")",sep="")
  }
}
Mat <- floor(log10(Mat))
rownames(Mat)<-tested_D_num
colnames(Mat)<-eval_time_lab
Cols<-brewer.pal(n = 4, name = "Reds")
pdf(file = "./Figures/Figure_S6.pdf",width = 14,height = 7)
par(mar=c(5,5,2.6,5.5))
Fast<-plot(Mat[nrow(Mat):1, ],breaks = c(0:4),col=Cols,digits = 1,fmt.cell='%.0f',
     xlab="Time of single run",ylab="Number of input parameters",
     text.cell=list(pos=3, cex=1),cex.lab=1.5,cex.axis=1,main="",fmt.key = "%.0f")
for (i in 1:nrow(Mat)) {
  for (j in 1:ncol(Mat)) {
    args<-Fast$cell.text[[nrow(Mat)+1-i,j]]
    args$labels <- textMat[i,j]
    args$cex    <- 1.5
    args$pos    <- 1
    do.call(text, args)
  }
}
mtext('Order of magnitude', side=3, line=1, at=15)
mtext(expression(paste("(",plain("T")[plain("2nd")]," / ",plain("T")[plain("1st")],")")),
      side=3, line=-0.5, at=15)
dev.off()

#-----------------------
