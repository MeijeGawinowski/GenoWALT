rm(list=ls())
library(pegas)

setwd("~/Documents/LD_Tests/Drift/Control/sizeD/Rep1")

vect_r1 <- c()

for (i in 0:59){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r1 <- abs(Mat[1])
  vect_r1 <- c(vect_r1,val_r1)
}

plot(0:59,vect_r1)

setwd("~/Documents/LD_Tests/Drift/Control/sizeD/Rep2")
vect_r2 <- c()

for (i in 0:59){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r2 <- abs(Mat[1])
  vect_r2 <- c(vect_r2,val_r2)
}

plot(0:59,vect_r2)

setwd("~/Documents/LD_Tests/Drift/Control/sizeD/Rep3")
vect_r3 <- c()

for (i in 0:59){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r3 <- abs(Mat[1])
  vect_r3 <- c(vect_r3,val_r3)
}

plot(0:59,vect_r3)

plot(0:59,vect_r1,type="l",col=1)
points(0:59,vect_r2,type="l",col=2)
points(0:59,vect_r3,type="l",col=3)

data_D <- data.frame(vect_r1,vect_r2,vect_r3)
mean_r <- c()
for (i in 1:dim(data_D)[1]){
  row_data <- as.numeric(data_D[i,])
  m_r <- mean(row_data)
  mean_r <- c(mean_r,m_r)
}

newdata_D <- cbind(data_D,mean_r)
plot(0:59,newdata_D$mean_r)

write.table(newdata_D,"~/Documents/LD_Tests/Drift/Control/tabD.csv")
