rm(list=ls())
library(genepop)

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/LOC2Fitness_Reproduction/allo0/")

pval <- c()
SE <- c()

for (i in 0:9){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name1 <- paste('pop',i,'.txt.DIS',sep="")
  res_name1 <- paste('pop',i,'.txt.TAB',sep="")
  test_LD(name,which='deficit',res_name1)
  write_LD_tables(data_file,res_name2)
  row1 <- strsplit(readLines(res_name1)[15]," ")[[1]]
  pval <- c(pval,as.numeric(row1[17]))
  SE <- c(SE,as.numeric(row1[22]))
}

gen <- 0:9

data <- data.frame(gen,pval,SE)
write.table(data,"LD_allo0.csv")
