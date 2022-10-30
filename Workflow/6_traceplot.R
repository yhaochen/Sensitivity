# This script plots the traceplot (Figure 2 of the paper)

# Remove all existing environment and plots
rm(list = ls())
graphics.off()

# Set a working directory, please set it to your own working folder when testing
dir <- commandArgs(trailingOnly=TRUE)
setwd(dir = dir)

# Load the required package
library(RColorBrewer)

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

# Create a folder to save figures
folder <- "./Figures"
if (!dir.exists(folder)){
  dir.create(folder, recursive = TRUE)
}

# 4 panels of trace plots + 1 line showing convergence locations
pdf(file = paste("./Figures/Figure_2.pdf",sep=""),width = 18,height = 12)
layout(matrix(c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE))
par(mar=c(4,5,6,2.6))
plot(Size_S,T_S,type="l",col="seagreen",xlab="Sample size",ylab="Sensitivity",ylim=c(0.8,1.2),
     cex.axis=1.7,cex.lab=1.7)
lines(Size_S,T_S_high,lty = 2, col = "seagreen")
lines(Size_S,T_S_low,lty = 2, col = "seagreen")
legend("topright",lty = c(1,2), col = c("seagreen","seagreen"), 
       legend = c("Standard Sobol","95% CI"), bty = "n", cex = 2)
mtext("a",side = 3, line = 1, at = 0, cex = 1.5)
arrows(C_S,0.85,C_S,0.95,length = 0.1, col = "seagreen")

par(mar=c(4,5,6,2.6))
plot(Size_K,T_K[1, ],type="l",col="purple",ylim = c(0.8,1.2),
     xlab="Sample size",ylab="Sensitivity",cex.axis=1.7,cex.lab=1.7)
for (i in 2:5){
  lines(Size_K,T_K[i, ],col="purple")
}
legend("topright",lty = 1, col = "purple", legend = "Kriging", bty = "n", cex = 2)
mtext("b",side = 3, line = 1, at = min(Size_K), cex = 1.5)
arrows(C_K,0.85,C_K,0.95,length = 0.1,col="purple")
abline(h = T_S[125], lty = 2, lwd = 0.5)

par(mar=c(4,5,2.6,2.6))
plot(Size_B,T_B[1, ],type="l",col="blue",ylim = c(0.8,1.2),
     xlab="Sample size",ylab="Sensitivity",cex.axis=1.7,cex.lab=1.7)
for (i in 2:5){
  lines(Size_B,T_B[i, ],col="blue")
}
legend("topright",lty = 1, col = "blue", legend = "BASS", bty = "n", cex = 2)
mtext("c",side = 3, line = 1, at = min(Size_B), cex = 1.5)
arrows(C_B,0.85,C_B,0.95,length = 0.1,col="blue")
abline(h = T_S[125], lty = 2, lwd = 0.5)

par(mar=c(4,5,2.6,2.6))
plot(Size_A,T_A[1, ],type="l",col="red",ylim = c(0.8,1.2),
     xlab="Sample size",ylab="Sensitivity",cex.axis=1.7,cex.lab=1.7)
for (i in 2:5){
  lines(Size_A,T_A[i, ],col="red")
}
legend("topright",lty = 1, col = "red", legend = "AKMCS", bty = "n", cex = 2)
mtext("d",side = 3, line = 1, at = min(Size_A), cex = 1.5)
arrows(C_A,0.85,C_A,0.95,length = 0.1,col="red")
abline(h = T_S[125], lty = 2, lwd = 0.5)

plot(0,0,type = "n", xaxt = "n", yaxt = "n", bty="n", xlab = "", ylab="",
     xlim=c(10, 100000), ylim=c(0, 0.7),log="x")
axis(1, at = c(10,100,1000,10000,100000),labels = c(10,100,1000,10000,100000), pos = 0.5,
     cex.axis = 2, cex.lab = 2)
points(C_S, 0.5, col = "seagreen", pch = 20, cex = 2)
points(C_K, 0.5, col = "purple", pch = 20, cex = 2)
points(C_B, 0.5, col = "blue", pch = 20, cex = 2)
points(C_A, 0.5, col = "red", pch = 20, cex = 2)
text(C_S, 0.6, labels = "Sobol", col = "seagreen", cex = 2)
text(C_K, 0.6, labels = "Kriging", col = "purple", cex = 2)
text(C_B, 0.4, labels = "BASS", col = "blue", cex = 2)
text(C_A, 0.6, labels = "AKMCS", col = "red", cex = 2)
text(1000, 0.2, labels = "Required sample size for convergence", col = "black", cex = 2)
dev.off()