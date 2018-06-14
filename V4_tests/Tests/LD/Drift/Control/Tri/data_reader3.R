setwd("~/Documents/LD_Tests/Drift/Control")

dataA <- read.table("tabA.csv",header=TRUE)
dataB <- read.table("tabB.csv",header=TRUE)
dataC <- read.table("tabD.csv",header=TRUE)

plot(0:59,dataA$mean_r,type="l",col="firebrick2",ylim=c(0,1.2),
     ylab="DL",xlab="t")
points(0:59,dataB$mean_r,type="l",col="mediumorchid")
points(0:59,dataC$mean_r,type="l",col="dodgerblue4")
legend("topright",c("size=100","size=500","size=1000"),
       lty=c(1,1,1),col=c("firebrick2","mediumorchid","dodgerblue4"))
