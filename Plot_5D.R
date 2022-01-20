#Result summary
rm(list = ls())
graphics.off()
library(fields) # Make grid plots

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}

#Load direct sobol results:
load("Data/5D/Sobol/Sobol_First_order")
load("Data/5D/Sobol/Sobol_First_order_up")
load("Data/5D/Sobol/Sobol_First_order_low")
load("Data/5D/Sobol/Sobol_time")
load("Data/5D/Sobol/Sobol_fail_prob")
n<-400000
Sobol_size<-seq(n/2000,n,by=n/2000)

#Load direct Kriging results:
load("Data/5D/Kriging/Kriging_First_order")
load("Data/5D/Kriging/Kriging_First_order_up")
load("Data/5D/Kriging/Kriging_First_order_low")
load("Data/5D/Kriging/Kriging_fail_prob")
load("Data/5D/Kriging/Kriging_time")
load("Data/5D/Kriging/Kriging_wrong_neg")
load("Data/5D/Kriging/Kriging_wrong_pos")
for (i in 2:51){
  Kriging_time[i]<-Kriging_time[i]-Kriging_time[i-1]
}
Kriging_size<-c(10:100)

#Load AKMCS results:
load("Data/5D/AKMCS/AKMCS_First_order")
load("Data/5D/AKMCS/AKMCS_First_order_up")
load("Data/5D/AKMCS/AKMCS_First_order_low")
load("Data/5D/AKMCS/AKMCS_fail_prob")
load("Data/5D/AKMCS/AKMCS_time")
load("Data/5D/AKMCS/AKMCS_wrong_neg")
load("Data/5D/AKMCS/AKMCS_wrong_pos")
AKMCS_size<-c(10:100)

# The palette with black:
Col <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Finish loading part
#---------------------------------------------
#Plot part


#Sensitivity indices by different methods
#Only sobol
par(mar=c(4,5,1.6,2))
plot(sobol_size, First_order[ ,1],type="l",ylim=c(-0.1,0.8),col=Col[1],lwd=2,xlab="Sample size",
     ylab="Sobol index",xaxt="n",cex.axis=1.5,cex.lab=1.5,log = "x")
axis(1, at = c(10,100,1000,10000,100000),cex.axis=1.5)
lines(sobol_size, First_order[ ,2],col=Col[2],lwd=2)
lines(sobol_size, First_order[ ,3],col=Col[3],lwd=2)
lines(sobol_size, First_order[ ,4],col=Col[4],lwd=2)
lines(sobol_size, First_order[ ,5],col=Col[5],lwd=2)
lines(sobol_size, First_order[ ,6],lwd=2)
lines(sobol_size, First_order_up[ ,1],col=Col[1],lty=2,lwd=0.2)
lines(sobol_size, First_order_up[ ,2],col=Col[2],lty=2,lwd=0.2)
lines(sobol_size, First_order_up[ ,3],col=Col[3],lty=2,lwd=0.2)
lines(sobol_size, First_order_up[ ,4],col=Col[4],lty=2,lwd=0.2)
lines(sobol_size, First_order_up[ ,5],col=Col[5],lty=2,lwd=0.2)
lines(sobol_size, First_order_up[ ,6],lty=2,lwd=1)
lines(sobol_size, First_order_low[ ,1],col=Col[1],lty=2,lwd=0.2)
lines(sobol_size, First_order_low[ ,2],col=Col[2],lty=2,lwd=0.2)
lines(sobol_size, First_order_low[ ,3],col=Col[3],lty=2,lwd=0.2)
lines(sobol_size, First_order_low[ ,4],col=Col[4],lty=2,lwd=0.2)
lines(sobol_size, First_order_low[ ,5],col=Col[5],lty=2,lwd=0.2)
lines(sobol_size, First_order_low[ ,6],lty=2,lwd=1)

conv<-matrix(NA,nrow=length(sobol_size)-9,ncol=5)
for (i in 1:(length(sobol_size)-9)){
  for (j in 1:5)
  conv[i,j] <- abs(sum(First_order[c(i:(i+4)),j])-sum(First_order[c((i+5):(i+9)),j]))
}
conv_sobol<-max(which(!(conv[ ,1:5] < rep(0.02,5)))) %% (length(sobol_size)-9) +1
abline(v=sobol_size[conv_sobol],lty=3)

#Kriging
par(mar=c(4,5,1.6,2))
plot(Kriging1_size, Kriging1_First_order[ ,1],type="l",ylim=c(-0.1,0.8),col=Col[1],lwd=2,xlab="Sample size",
     ylab="Kriging1 index",cex.axis=1.5,cex.lab=1.5)
lines(Kriging1_size, Kriging1_First_order[ ,2],col=Col[2],lwd=2)
lines(Kriging1_size, Kriging1_First_order[ ,3],col=Col[3],lwd=2)
lines(Kriging1_size, Kriging1_First_order[ ,4],col=Col[4],lwd=2)
lines(Kriging1_size, Kriging1_First_order[ ,5],col=Col[5],lwd=2)
lines(Kriging1_size, Kriging1_First_order_up[ ,1],col=Col[1],lty=2,lwd=0.2)
lines(Kriging1_size, Kriging1_First_order_up[ ,2],col=Col[2],lty=2,lwd=0.2)
lines(Kriging1_size, Kriging1_First_order_up[ ,3],col=Col[3],lty=2,lwd=0.2)
lines(Kriging1_size, Kriging1_First_order_up[ ,4],col=Col[4],lty=2,lwd=0.2)
lines(Kriging1_size, Kriging1_First_order_up[ ,5],col=Col[5],lty=2,lwd=0.2)
lines(Kriging1_size, Kriging1_First_order_low[ ,1],col=Col[1],lty=2,lwd=0.2)
lines(Kriging1_size, Kriging1_First_order_low[ ,2],col=Col[2],lty=2,lwd=0.2)
lines(Kriging1_size, Kriging1_First_order_low[ ,3],col=Col[3],lty=2,lwd=0.2)
lines(Kriging1_size, Kriging1_First_order_low[ ,4],col=Col[4],lty=2,lwd=0.2)
lines(Kriging1_size, Kriging1_First_order_low[ ,5],col=Col[5],lty=2,lwd=0.2)
conv<-matrix(NA,nrow=length(Kriging1_size)-9,ncol=5)
for (i in 1:(length(Kriging1_size)-9)){
  for (j in 1:5)
    conv[i,j] <- abs(sum(Kriging1_First_order[c(i:(i+4)),j])-sum(Kriging1_First_order[c((i+5):(i+9)),j]))
}
conv_Kriging<-max(which(!(conv[ ,1:5] < rep(0.02,5)))) %% (length(Kriging1_size)-9) +1
abline(v=Kriging1_size[conv_Kriging],lty=3)



#AKMCS
par(mar=c(4,5,1.6,2))
plot(AKMCS_size, AKMCS_First_order[ ,1],type="l",ylim=c(-0.1,0.8),col=Col[1],lwd=2,xlab="Sample size",
     ylab="AKMCS index",cex.axis=1.5,cex.lab=1.5)
lines(AKMCS_size, AKMCS_First_order[ ,2],col=Col[2],lwd=2)
lines(AKMCS_size, AKMCS_First_order[ ,3],col=Col[3],lwd=2)
lines(AKMCS_size, AKMCS_First_order[ ,4],col=Col[4],lwd=2)
lines(AKMCS_size, AKMCS_First_order[ ,5],col=Col[5],lwd=2)
lines(AKMCS_size, AKMCS_First_order_up[ ,1],col=Col[1],lty=2,lwd=0.2)
lines(AKMCS_size, AKMCS_First_order_up[ ,2],col=Col[2],lty=2,lwd=0.2)
lines(AKMCS_size, AKMCS_First_order_up[ ,3],col=Col[3],lty=2,lwd=0.2)
lines(AKMCS_size, AKMCS_First_order_up[ ,4],col=Col[4],lty=2,lwd=0.2)
lines(AKMCS_size, AKMCS_First_order_up[ ,5],col=Col[5],lty=2,lwd=0.2)
lines(AKMCS_size, AKMCS_First_order_low[ ,1],col=Col[1],lty=2,lwd=0.2)
lines(AKMCS_size, AKMCS_First_order_low[ ,2],col=Col[2],lty=2,lwd=0.2)
lines(AKMCS_size, AKMCS_First_order_low[ ,3],col=Col[3],lty=2,lwd=0.2)
lines(AKMCS_size, AKMCS_First_order_low[ ,4],col=Col[4],lty=2,lwd=0.2)
lines(AKMCS_size, AKMCS_First_order_low[ ,5],col=Col[5],lty=2,lwd=0.2)

conv<-matrix(NA,nrow=length(AKMCS_size)-9,ncol=5)
for (i in 1:(length(AKMCS_size)-9)){
  for (j in 1:5)
    conv[i,j] <- abs(sum(AKMCS_First_order[c(i:(i+4)),j])-sum(AKMCS_First_order[c((i+5):(i+9)),j]))
}
conv_AKMCS<-max(which(!(conv[ ,1:5] < rep(0.02,5)))) %% (length(AKMCS_size)-9) +1
abline(v=AKMCS_size[conv_AKMCS],lty=3)


#All together:
polygonindx<-seq(25,1995,by=5)
par(mar=c(4,5,1.6,2))
plot(sobol_size, First_order[ ,1],type="l",ylim=c(-0.3,0.9),xlim=c(10,100000),col=Col[1],lwd=2,xlab="Sample size",
     ylab="Sobol index",xaxt="n",cex.axis=1.5,cex.lab=1.5,log = "x")
axis(1, at = c(10,100,1000,10000,100000),cex.axis=1.5)
polygon(c(sobol_size[polygonindx],rev(sobol_size[polygonindx])), c(First_order_up[polygonindx,1],rev(First_order_low[polygonindx,1])),
        border="NA", col = t_col(Col[1],90))
lines(sobol_size, First_order[ ,2],col=Col[2],lwd=2)
polygon(c(sobol_size[polygonindx],rev(sobol_size[polygonindx])), c(First_order_up[polygonindx,2],rev(First_order_low[polygonindx,2])),
        border="NA", col = t_col(Col[2],90))
lines(sobol_size, First_order[ ,3],col=Col[3],lwd=2)
polygon(c(sobol_size[polygonindx],rev(sobol_size[polygonindx])), c(First_order_up[polygonindx,3],rev(First_order_low[polygonindx,3])),
        border="NA", col = t_col(Col[3],90))
lines(sobol_size, First_order[ ,4],col=Col[4],lwd=2)
polygon(c(sobol_size[polygonindx],rev(sobol_size[polygonindx])), c(First_order_up[polygonindx,4],rev(First_order_low[polygonindx,4])),
        border="NA", col = t_col(Col[4],90))
lines(sobol_size, First_order[ ,5],col=Col[5],lwd=2)
polygon(c(sobol_size[polygonindx],rev(sobol_size[polygonindx])), c(First_order_up[polygonindx,5],rev(First_order_low[polygonindx,5])),
        border="NA", col = t_col(Col[5],90))
lines(sobol_size, First_order[ ,6],lwd=2)
polygon(c(sobol_size[polygonindx],rev(sobol_size[polygonindx])), c(First_order_up[polygonindx,6],rev(First_order_low[polygonindx,6])),
        border="NA", col = t_col("grey",50))
abline(v=sobol_size[conv_sobol],lty=3)


lines(Kriging1_size, Kriging1_First_order[ ,1],col=Col[1],lwd=2,lty=3)
polygon(c(Kriging1_size,rev(Kriging1_size)),c(Kriging1_First_order_up[ ,1],rev(Kriging1_First_order_low[ ,1])),
        border = "NA",col = t_col(Col[1],90))
lines(Kriging1_size, Kriging1_First_order[ ,2],col=Col[2],lwd=2,lty=3)
polygon(c(Kriging1_size,rev(Kriging1_size)),c(Kriging1_First_order_up[ ,2],rev(Kriging1_First_order_low[ ,2])),
        border = "NA",col = t_col(Col[2],90))
lines(Kriging1_size, Kriging1_First_order[ ,3],col=Col[3],lwd=2,lty=3)
polygon(c(Kriging1_size,rev(Kriging1_size)),c(Kriging1_First_order_up[ ,3],rev(Kriging1_First_order_low[ ,3])),
        border = "NA",col = t_col(Col[3],90))
lines(Kriging1_size, Kriging1_First_order[ ,4],col=Col[4],lwd=2,lty=3)
polygon(c(Kriging1_size,rev(Kriging1_size)),c(Kriging1_First_order_up[ ,4],rev(Kriging1_First_order_low[ ,4])),
        border = "NA",col = t_col(Col[4],90))
lines(Kriging1_size, Kriging1_First_order[ ,5],col=Col[5],lwd=2,lty=3)
polygon(c(Kriging1_size,rev(Kriging1_size)),c(Kriging1_First_order_up[ ,5],rev(Kriging1_First_order_low[ ,5])),
        border = "NA",col = t_col(Col[5],90))
abline(v=Kriging1_size[conv_Kriging],lty=3)
legend(50,0.95,lwd=rep(2,6),col=c(Col,"black"),legend = c("X1 (Sobol)","X2 (Sobol)","X3 (Sobol)","X4 (Sobol)","X5 (Sobol)","Dummy variable"),bty="n",cex=1.5)
legend(300,0.95,lwd=rep(2,5),col=c(Col),lty=rep(3,5),legend = c("X1 (Kriging)","X2 (Kriging)","X3 (Kriging)","X4 (Kriging)","X5 (Kriging)"),bty="n",cex=1.5)
legend(300,0.71,fill=t_col(Col[1],90),legend="95% CI by 100 bootstraps",bty="n",cex=1.5)



# in one column
#polygonindx<-seq(25,995,by=5)
Sobol_time<-Sobol_time+Sobol_size*60
Kriging_time<-Kriging_time+Kriging_size*60
AKMCS_time<-AKMCS_time+AKMCS_size*60


par(mfrow=c(3, 1))
plot(Sobol_size, Sobol_First_order[ ,1],type="l",ylim=c(0,1),xlim=c(200,400000),col=Col[1],lwd=2,xlab="Sample size",
     ylab="Sobol index",xaxt="n",cex.axis=1.5,cex.lab=1.5,log = "x")
axis(1, at = c(10,100,1000,10000,100000),cex.axis=1.5)
axis(3, at = Sobol_size[seq(10,1960,by=50)], labels = round(Sobol_time[seq(10,1960,by=50)]),cex.axis=1.5)
mtext("Sobol computational time (s)",3,line=2.5,at=7500,cex = 1)
lines(Sobol_size, Sobol_First_order[ ,2],col=Col[2],lwd=2)
lines(Sobol_size, Sobol_First_order[ ,3],col=Col[3],lwd=2)
lines(Sobol_size, Sobol_First_order[ ,4],col=Col[4],lwd=2)
lines(Sobol_size, Sobol_First_order[ ,5],col=Col[5],lwd=2)
lines(Sobol_size, Sobol_First_order[ ,6],lwd=2)
legend(30000,1,lwd=rep(2,3),col=c(Col[1:3]),legend = c("X1","X2","X3"),bty="n",cex=1.5)
legend(80000,1,lwd=rep(2,3),col=c(Col[4:5],"black"),legend = c("X4","X5","Dummy variable"),bty="n",cex=1.5)
#legend("topright",fill=t_col(Col[1],90),legend="95% CI by 100 bootstraps",bty="n",cex=1.5)
#abline(v=Sobol_size[conv_sobol],lty=3)

plot(Kriging_size, Kriging_First_order[ ,1],type="l",col=Col[1],xlim=c(10,60),ylim=c(0,1),lwd=2,lty=1,xlab="Sample size",
     ylab="Sobol index",xaxt="n",cex.axis=1.5,cex.lab=1.5)
axis(1, at = c(10,20,30,40,50,60,70,80,90,100),cex.axis=1.5)
axis(3, at = Kriging_size, labels = round(Kriging_time),cex.axis=1.5)
mtext("Kriging computational time (s)",3,line=2.5,at=35,cex = 1)
lines(Kriging_size, Kriging_First_order[ ,2],col=Col[2],lwd=2,lty=1)
lines(Kriging_size, Kriging_First_order[ ,3],col=Col[3],lwd=2,lty=1)
lines(Kriging_size, Kriging_First_order[ ,4],col=Col[4],lwd=2,lty=1)
lines(Kriging_size, Kriging_First_order[ ,5],col=Col[5],lwd=2,lty=1)
lines(Kriging_size, Kriging_First_order[ ,6],lwd=2)

#abline(v=Kriging_size[conv_Kriging],lty=3)


plot(AKMCS_size, AKMCS_First_order[ ,1],type="l",col=Col[1],xlim=c(10,60),ylim=c(0,1),lwd=2,lty=1,xlab="Sample size",
     ylab="Sobol index",xaxt="n",cex.axis=1.5,cex.lab=1.5)
axis(1, at = c(10,20,30,40,50,60,70,80,90,100),cex.axis=1.5)
axis(3, at = AKMCS_size, labels = round(AKMCS_time),cex.axis=1.5)
mtext("AKMCS computational time (s)",3,line=2.5,at=35,cex = 1)
lines(AKMCS_size, AKMCS_First_order[ ,2],col=Col[2],lwd=2,lty=1)
lines(AKMCS_size, AKMCS_First_order[ ,3],col=Col[3],lwd=2,lty=1)
lines(AKMCS_size, AKMCS_First_order[ ,4],col=Col[4],lwd=2,lty=1)
lines(AKMCS_size, AKMCS_First_order[ ,5],col=Col[5],lwd=2,lty=1)
lines(AKMCS_size, AKMCS_First_order[ ,6],lwd=2)

#abline(v=AKMCS_size[conv_AKMCS],lty=3)



# Compare the failure probability versus computational time:
par(mar=c(4,5,5,2))
plot(Kriging_time,Kriging_fail_prob,type="l",lwd=2,lty=1,xlab="Computational time (s)",xlim=c(0,15000),ylim=c(0,7.5),
     ylab="Estimated failure probability (%)",cex.axis=1.5,cex.lab=1.5,col=Col[2])
lines(AKMCS_time,AKMCS_fail_prob,col=Col[3])
axis(3, at = AKMCS_time, labels = AKMCS_size,cex.axis=1.5,col.ticks = Col[3],col.axis=Col[3])
mtext("AKMCS sample size",3,line=2.5,at=7500,cex = 1.5,col=Col[3])

axis(3, at = Kriging_time, labels = Kriging_size,cex.axis=1.5,col=Col[2],col.ticks = Col[2],
     col.axis=Col[2],pos=7)
mtext("Kriging sample size",3,line=-2.5,at=12500,cex = 1.5,col=Col[2])


lines(Sobol_time,Sobol_fail_prob,col=Col[1])
axis(3, at = Sobol_time[seq(1,1951,by=50)], labels = Sobol_size[seq(1,1951,by=50)],
     cex.axis=1.5,col=Col[1],col.axis=Col[1],pos=6.2)
mtext("Brute force LHS sample size",3,line=-4.5,at=8500,cex = 1.5,col = Col[1])
legend(11000,6.3,lty = rep(1,3),col = Col[1:3],legend = c("LHS Sobol","Kriging","AKMCS"),bty="n",cex = 1.5)



# Compare the wrong negative/positive rate
par(mar=c(4,5,5,2))
par(mfrow=c(2, 1))

plot(AKMCS_time,AKMCS_wrong_neg,type="l",lwd=2,lty=1,xlab="Computational time (s)",ylim=c(0,5),
     ylab="False failure rate (%)",cex.axis=1.5,cex.lab=1.5,col=Col[3])
lines(Kriging_time,Kriging_wrong_neg,col=Col[2])
axis(3, at = Kriging_time, labels = Kriging_size,cex.axis=1.5,col.axis=Col[2],col.ticks = Col[2])
axis(3, at = AKMCS_time, labels = AKMCS_size,cex.axis=1.5,col.axis=Col[3],col.ticks = Col[3],pos=4.3)
mtext("Kriging sample size",3,line=2.5,at=7500,cex = 1.5,col=Col[2])
mtext("AKMCS sample size",3,line=-4.5,at=7500,cex = 1.5,col=Col[3])

legend(12000,4,lty = rep(1,2),col = Col[2:3],legend = c("Kriging","AKMCS"),bty="n",cex = 1.5)

plot(AKMCS_time,AKMCS_wrong_pos,type="l",lwd=2,lty=1,xlab="Computational time (s)",ylim=c(0,115),
     ylab="False success rate (%)",cex.axis=1.5,cex.lab=1.5,col=Col[3])
lines(Kriging_time,Kriging_wrong_pos,col=Col[2])
axis(3, at = Kriging_time, labels = Kriging_size,cex.axis=1.5,col.axis=Col[2],col.ticks = Col[2])
axis(3, at = AKMCS_time, labels = AKMCS_size,cex.axis=1.5,col.axis=Col[3],col.ticks = Col[3],pos=100)
mtext("Kriging sample size",3,line=2.5,at=7500,cex = 1.5,col=Col[2])
mtext("AKMCS sample size",3,line=-4.5,at=7500,cex = 1.5,col=Col[3])

legend(12000,80,lty = rep(1,2),col = Col[2:3],legend = c("Kriging","AKMCS"),bty="n",cex = 1.5)

#Horse race plot
d=5
Eva_time=0.01 #(s)

layout(mat = matrix(c(1,2,3),nrow = 3,ncol=1),heights = c(1,1,5))
par(mar=c(4,5,5,2))
plot(c(Eva_time),c(0.5),xlim=c(0.01,86400),ylim=c(0,1),axes = FALSE,
     xlab = "",ylab = "", log="x",pch=20,cex=2) 
axis(1, at=c(0.01, 0.1, 1, 10, 60, 600, 3600, 86400),cex.label=1.5,
     labels = c("0.01s","0.1s","1s","10s","1min","10min","1h","1d"))
mtext("Computational time per evaluation",3,line = 2)

par(mar=c(4,5,5,2))
plot(c(d),c(0.5),xlim=c(0,50),ylim=c(0,1),axes = FALSE,
     xlab = "",ylab = "", pch=20,cex=2) 
axis(1, at=c(1:50))
mtext("Model dimension (number of parameters)",3,line = 2)

plot(Sobol_size,Sobol_First_order[ ,1],type = "l",xlim=c(10,1000000),ylim=c(0,1),log="x",col="blue",lwd=2,
     lty=1,xlab="Sample size",ylab="Sobol index",xaxt="n",cex.axis=1.5,cex.lab=1.5)
axis(1,at=c(10,100,1000,10000,100000,1000000),cex.axis=1.5)
lines(Kriging_size,Kriging_First_order[ ,1],lwd=2,lty=3,col="blue")
lines(AKMCS_size,AKMCS_First_order[ ,1],lwd=2,lty=5,col="blue")
lines(Sobol_size,Sobol_First_order[ ,3],lwd=2,lty=3,col="red")
lines(Kriging_size,Kriging_First_order[ ,3],lwd=2,lty=3,col="red")
lines(AKMCS_size,AKMCS_First_order[ ,3],lwd=2,lty=5,col="red")
legend("topright",lty = c(1,3,5,1,3,5),col = c("blue","blue","blue","red","red","red"),bty="n",
       lwd=rep(2,6),legend = c("Sobol (X1)","Kriging + Sobol (X1)","AKMCS + Sobol (X1)",
                               "Sobol (X3)","Kriging + Sobol (X3)","AKMCS + Sobol (X3)"),cex=1.5)

plot(Sobol_time,Sobol_First_order[ ,1],pch=20,xlim=c(0.1,100000),ylim=c(0,1),log="x",col="blue",lwd=2,
     lty=1,xlab="Computational time",ylab="Sobol index",xaxt="n",cex.axis=1.5,cex.lab=1.5)
axis(1,at=c(0.01, 0.1, 1, 10, 60, 600, 3600, 86400),
     labels = c("0.01s","0.1s","1s","10s","1min","10min","1h","1d"),cex.axis=1.5)
lines(Kriging_time,Kriging_First_order[ ,1],lwd=2,lty=3,col="blue")
lines(AKMCS_time,AKMCS_First_order[ ,1],lwd=2,lty=5,col="blue")
lines(Sobol_time,Sobol_First_order[ ,3],lwd=2,lty=3,col="red")
lines(Kriging_time,Kriging_First_order[ ,3],lwd=2,lty=3,col="red")
lines(AKMCS_time,AKMCS_First_order[ ,3],lwd=2,lty=5,col="red")
legend("topright",lty = c(1,3,5,1,3,5),col = c("blue","blue","blue","red","red","red"),bty="n",
       lwd=rep(2,6),legend = c("Sobol (X1)","Kriging + Sobol (X1)","AKMCS + Sobol (X1)",
                               "Sobol (X3)","Kriging + Sobol (X3)","AKMCS + Sobol (X3)"),cex=1.5)
