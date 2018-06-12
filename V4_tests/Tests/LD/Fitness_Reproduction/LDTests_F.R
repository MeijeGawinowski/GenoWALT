#######################
####### allo = 0 ######
#######################


rm(list=ls())
library(pegas)

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Fitness_Reproduction/allo0")

Delta <- c() 

for (i in 0:99){
  name <- paste('pop',i,'.gen',sep="")
  pegas_file <- read.genepop(name)
  loc_file <- as.loci(pegas_file)
  Mat <- LD2(loc_file)$Delta
  Delta <- c(Delta,Mat[2])
}

gen <- 0:99

data <- data.frame(gen,Delta)
write.table(data,"LD_allo0.csv")

#########################
####### allo = 0.5 ######
#########################


rm(list=ls())
library(pegas)

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Fitness_Reproduction/allo0.5")

Delta <- c() 

for (i in 0:99){
  name <- paste('pop',i,'.gen',sep="")
  pegas_file <- read.genepop(name)
  loc_file <- as.loci(pegas_file)
  Mat <- LD2(loc_file)$Delta
  Delta <- c(Delta,Mat[2])
}

gen <- 0:99

data <- data.frame(gen,Delta)
write.table(data,"LD_allo0.5.csv")

#######################
####### allo = 1 ######
#######################


rm(list=ls())
library(pegas)

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Fitness_Reproduction/allo1")

Delta <- c() 

for (i in 0:99){
  name <- paste('pop',i,'.gen',sep="")
  pegas_file <- read.genepop(name)
  loc_file <- as.loci(pegas_file)
  Mat <- LD2(loc_file)$Delta
  Delta <- c(Delta,Mat[2])
}

gen <- 0:99

data <- data.frame(gen,Delta)
write.table(data,"LD_allo1.csv")

###################
####### PLOT ######
###################


rm(list=ls())

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Fitness_Reproduction")

data1 <- read.table("allo0/LD_allo0.csv")
data2 <- read.table("allo0.5/LD_allo0.5.csv")
data3 <- read.table("allo1/LD_allo1.csv")


par(mfrow=c(2,2))
plot(data1$gen,data1$Delta,main="LD, allo=0, size=1000",xlab="Generations",ylab="LD",col="firebrick2")
plot(data1$gen,data2$Delta,main="LD, allo=0.5, size=1000",xlab="Generations",ylab="LD",col="dodgerblue4")
plot(data1$gen,data3$Delta,main="LD, allo=1, size=1000",xlab="Generations",ylab="LD",col="forestgreen")
