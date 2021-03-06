setwd("~/Documents/LD_Tests/Drift/Fitness/Ntri")

dataA <- read.table("tabA.csv",header=TRUE)
dataB <- read.table("tabB.csv",header=TRUE)
dataC <- read.table("tabC.csv",header=TRUE)
dataD <- read.table("tabD.csv",header=TRUE)

plot(0:59,dataA$mean_r,type="l",col="firebrick2",ylim=c(0,0.8),
     ylab="DL",xlab="t")
points(0:59,dataB$mean_r,type="l",col="mediumorchid")
points(0:59,dataC$mean_r,type="l",col="dodgerblue4")
points(0:59,dataD$mean_r,type="l",col="forestgreen")
legend("topright",cex=0.5,c("size=100","size=500","size=1000","size=10 000"),
       lty=c(1,1,1,1),col=c("firebrick2","mediumorchid","dodgerblue4",
                            "forestgreen"))
