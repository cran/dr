source("dr/R/dr.R")
source("luca.R")
dr.compute(x,y,rep(1,100),method="ire",nslices=10)->out1
x1 <- x + .1*matrix(rnorm(100*10),ncol=10)
dr.compute(x1,y,rep(1,100),method="ire",nslices=11)->out1
x2 <- matrix(rnorm(100*10),ncol=10)
dr.compute(x2,y,rep(1,100),method="ire",nslices=11)->out2

ais <- read.table("dr/data/ais.txt",header=TRUE)
a1 <- dr(LBM~.-Sex-Label-Sport,data=ais,method="ire")
a2 <- update(a1,subset=1:100)