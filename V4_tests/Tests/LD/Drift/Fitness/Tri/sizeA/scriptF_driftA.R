rm(list=ls())
library(pegas)

setwd("~/Documents/LD_Tests/Drift/Fitness/sizeA/Rep1")

vect_r1 <- c()

for (i in 0:124){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r1 <- abs(Mat[1])
  vect_r1 <- c(vect_r1,val_r1)
}

plot(0:124,vect_r1)

setwd("~/Documents/LD_Tests/Drift/Fitness/sizeA/Rep2")
vect_r2 <- c()

for (i in 0:124){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r2 <- abs(Mat[1])
  vect_r2 <- c(vect_r2,val_r2)
}

plot(0:124,vect_r2)

setwd("~/Documents/LD_Tests/Drift/Fitness/sizeA/Rep3")
vect_r3 <- c()

for (i in 0:124){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r3 <- abs(Mat[1])
  vect_r3 <- c(vect_r3,val_r3)
}

plot(0:124,vect_r3)

setwd("~/Documents/LD_Tests/DDrift/Fitness/sizeA/Rep4")
vect_r4 <- c()

for (i in 0:124){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r4 <- abs(Mat[1])
  vect_r4 <- c(vect_r4,val_r4)
}

plot(0:124,vect_r4)

setwd("~/Documents/LD_Tests/Drift/Fitness/sizeA/Rep5")
vect_r5 <- c()

for (i in 0:124){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r5 <- abs(Mat[1])
  vect_r5 <- c(vect_r5,val_r5)
}

plot(0:124,vect_r5)

plot(0:124,vect_r1,type="l",col=1,ylim=c(0,1))
points(0:124,vect_r2,type="l",col=2)
points(0:124,vect_r3,type="l",col=3)
points(0:124,vect_r4,type="l",col=4)
points(0:124,vect_r5,type="l",col=5)

data_A <- data.frame(vect_r1,vect_r2,vect_r3,vect_r4,vect_r5)
mean_r <- c()
for (i in 1:dim(data_A)[1]){
  row_data <- as.numeric(data_A[i,])
  m_r <- mean(row_data)
  mean_r <- c(mean_r,m_r)
}

newdata_A <- cbind(data_A,mean_r)
plot(0:124,newdata_A$mean_r)

write.table(newdata_A,"~/Documents/LD_Tests/Drift/Fitness/tabA.csv")
