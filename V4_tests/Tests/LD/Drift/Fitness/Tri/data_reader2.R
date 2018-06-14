setwd("~/Documents/LD_Tests/Drift/Fitness")

dataA <- read.table("tabA.csv",header=TRUE)
dataB <- read.table("tabB.csv",header=TRUE)
dataC <- read.table("tabC.csv",header=TRUE)
dataD <- read.table("tabD.csv",header=TRUE)

plot(0:124,dataA$mean_r,type="l",col="firebrick2",ylim=c(0,1.5),
     ylab="DL",xlab="t")
points(0:199,dataB$mean_r,type="l",col="mediumorchid")
points(0:199,dataC$mean_r,type="l",col="dodgerblue4")
points(0:199,dataD$mean_r,type="l",col="forestgreen")
legend("topright",c("size=100","size=500","size=1000","size=10 000"),
       lty=c(1,1,1,1),col=c("firebrick2","mediumorchid","dodgerblue4",
                            "forestgreen"))
