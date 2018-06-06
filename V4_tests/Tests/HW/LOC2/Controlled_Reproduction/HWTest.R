#######################
####### allo = 0 ######
#######################

rm(list=ls())
library(genepop)
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC2/Controlled_Reproduction/allo0")

pval1 <- c()
SE1 <- c()
WC1 <- c()

pval2 <- c()
SE2 <- c()
WC2 <- c()


for (i in 0:9){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row1 <- strsplit(readLines(res_name)[31]," ")[[1]]
  row2 <- strsplit(readLines(res_name)[32]," ")[[1]]
  # row2 <- strsplit(row," ")[[1]]
  pval1 <- c(pval1,as.numeric(row1[9]))
  SE1 <- c(SE1,as.numeric(row1[11]))
  WC1 <- c(WC1,as.numeric(row1[14]))
  pval2 <- c(pval2,as.numeric(row2[9]))
  SE2 <- c(SE2,as.numeric(row2[11]))
  WC2 <- c(WC2,as.numeric(row2[14]))
}

m_pval <- (pval1+pval2)/2
m_SE <- (SE1+SE2)/2
m_WC <- (WC1+WC2)/2

gen <- 0:9

data <- data.frame(gen,pval1,pval2,m_pval,SE1,SE2,m_SE,WC1,WC2,m_WC)
write.table(data,"HW_allo0.csv")

#########################
####### allo = 0.5 ######
#########################

rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC2/Controlled_Reproduction/allo0.5")
pval1 <- c()
SE1 <- c()
WC1 <- c()

pval2 <- c()
SE2 <- c()
WC2 <- c()


for (i in 0:9){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row1 <- strsplit(readLines(res_name)[31]," ")[[1]]
  row2 <- strsplit(readLines(res_name)[32]," ")[[1]]
  # row2 <- strsplit(row," ")[[1]]
  pval1 <- c(pval1,as.numeric(row1[9]))
  SE1 <- c(SE1,as.numeric(row1[11]))
  WC1 <- c(WC1,as.numeric(row1[14]))
  pval2 <- c(pval2,as.numeric(row2[9]))
  SE2 <- c(SE2,as.numeric(row2[11]))
  WC2 <- c(WC2,as.numeric(row2[14]))
}

m_pval <- (pval1+pval2)/2
m_SE <- (SE1+SE2)/2
m_WC <- (WC1+WC2)/2

gen <- 0:9

data <- data.frame(gen,pval1,pval2,m_pval,SE1,SE2,m_SE,WC1,WC2,m_WC)
write.table(data,"HW_allo0.5.csv")

#######################
####### allo = 1 ######
#######################

rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC2/Controlled_Reproduction/allo1")
pval1 <- c()
SE1 <- c()
WC1 <- c()

pval2 <- c()
SE2 <- c()
WC2 <- c()


for (i in 0:9){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row1 <- strsplit(readLines(res_name)[31]," ")[[1]]
  row2 <- strsplit(readLines(res_name)[32]," ")[[1]]
  # row2 <- strsplit(row," ")[[1]]
  pval1 <- c(pval1,as.numeric(row1[9]))
  SE1 <- c(SE1,as.numeric(row1[11]))
  WC1 <- c(WC1,as.numeric(row1[14]))
  pval2 <- c(pval2,as.numeric(row2[9]))
  SE2 <- c(SE2,as.numeric(row2[11]))
  WC2 <- c(WC2,as.numeric(row2[14]))
}

m_pval <- (pval1+pval2)/2
m_SE <- (SE1+SE2)/2
m_WC <- (WC1+WC2)/2

gen <- 0:9

data <- data.frame(gen,pval1,pval2,m_pval,SE1,SE2,m_SE,WC1,WC2,m_WC)
write.table(data,"HW_allo1.csv")


###################
####### PLOT ######
###################

rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC2/Controlled_Reproduction")

data1 <- read.table("allo0/HW_allo0.csv")
data2 <- read.table("allo0.5/HW_allo0.5.csv")
data3 <- read.table("allo1/HW_allo1.csv")

pdf("LOC2_Control_1000.pdf")
par(mfrow=c(3,3))
plot(data1$gen,data1$pval1,type='o',main="Locus 1, Fitness",ylab="p-value",
     xlab="Generations",col="firebrick2",ylim=c(0,1.3))
points(data1$gen,data2$pval1,type="o",col="dodgerblue4")
points(data1$gen,data3$pval1,type="o",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$pval2,type='o',main="Locus 2, Fitness",ylab="p-value",
     xlab="Generations",col="firebrick2",ylim=c(0,1.3))
points(data1$gen,data2$pval2,type="o",col="dodgerblue4")
points(data1$gen,data3$pval2,type="o",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$m_pval,type='o',main="Mean loci, Fitness",ylab="p-value",
     xlab="Generations",col="firebrick2",ylim=c(0,1.3))
points(data1$gen,data2$m_pval,type="o",col="dodgerblue4")
points(data1$gen,data3$m_pval,type="o",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$SE1,type='o',main="Locus 1, Fitness",ylab="SE",
     xlab="Generations",col="firebrick2",ylim=c(0,0.02))
points(data1$gen,data2$SE1,type="o",col="dodgerblue4")
points(data1$gen,data3$SE1,type="o",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$SE2,type='o',main="Locus 2, Fitness",ylab="SE",
     xlab="Generations",col="firebrick2",ylim=c(0,0.025))
points(data1$gen,data2$SE2,type="o",col="dodgerblue4")
points(data1$gen,data3$SE2,type="o",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$m_SE,type='o',main="Mean loci, Fitness",ylab="SE",
     xlab="Generations",col="firebrick2",ylim=c(0,0.02))
points(data1$gen,data2$m_SE,type="o",col="dodgerblue4")
points(data1$gen,data3$m_SE,type="o",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$WC1,type='o',main="Locus 1, Fitness",ylab="W&C Fis",
     xlab="Generations",col="firebrick2",ylim=c(-0.15,1.3))
points(data1$gen,data2$WC1,type="o",col="dodgerblue4")
points(data1$gen,data3$WC1,type="o",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$WC2,type='o',main="Locus 2, Fitness",ylab="W&C Fis",
     xlab="Generations",col="firebrick2",ylim=c(-0.15,1.3))
points(data1$gen,data2$WC2,type="o",col="dodgerblue4")
points(data1$gen,data3$WC2,type="o",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$m_WC,type='o',main="Mean loci, Fitness",ylab="W&C Fis",
     xlab="Generations",col="firebrick2",ylim=c(-0.15,1.3))
points(data1$gen,data2$m_WC,type="o",col="dodgerblue4")
points(data1$gen,data3$m_WC,type="o",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))
dev.off()

