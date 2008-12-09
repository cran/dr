### use setwd for your own directory, not mine.
setwd("C://Documents and Settings//sandy//My Documents//dr3//dr//Rnew")

setwd("z:/R/source/dr3/dr/Rnew")
library(alr3)
data(ais)
ais$one <- rep(1,202)
source("dr.R")
source("psir.R")
source("psave.R")
# SIR
summary(s0 <- dr(LBM~log(SSF)+log(Wt)+log(Hg)+log(Ht)+log(WCC)+log(RCC)+
  log(Hc)+log(Ferr),nslices=10,data=ais))
summary(s1 <- update(s0,group=~Sex,nslices=10))
# SAVE
summary(s2 <- update(s0,method="save"))
summary(s3 <- update(s2,group=~Sex, method="save",nslices=10))

# using the library
library(dr)
data(ais)
ais$one <- rep(1,202)
# SIR
model <- LBM~log(SSF)+log(Wt)+log(Hg)+log(Ht)+log(WCC)+log(RCC)+log(Hc)+log(Ferr)
summary(sir1 <- dr(model,nslices=10,data=ais))
summary(sir2 <- update(sir1,group=~Sex))
# SAVE
summary(save1 <- update(sir1,method="save"))
summary(save2 <- update(save1,group=~Sex))
# IRE
(ire1 <- update(sir1,method="ire"))
(ire2 <- update(ire1,group=~Sex))

drop1(sir1,update=FALSE)
drop1(sir2,update=FALSE)
drop1(save1,update=FALSE)
drop1(save1,update=FALSE,test="normal")
drop1(save2,update=FALSE)  # NOT ordered
drop1(ire1,update=FALSE)
drop1(ire2,update=FALSE)

dr.step(save2)
dr.step(sir1,stop=0.05)

ais$g4 <- sample(1:4,202,replace=TRUE)
ais$g2 <- sample(1:2,202,replace=TRUE)
summary(s4 <- update(s3,group=~factor(g4)))
summary(s5 <- update(s3,group=~factor(Sex):factor(g2)))

dr.step(s1)
s6 <- update(s1,~log(Ferr)+log(RCC)+log(SSF)+log(Wt))

drop1(s0,update=FALSE)
drop1(s1,update=FALSE)
drop1(s2,update=FALSE)
drop1(s3,update=FALSE)

dr.coordinate.test(s0,~.-log(WCC),d=3)
dr.coordinate.test(s1,~.-log(WCC),d=3)

dr.coordinate.test(s2,~.-log(WCC),d=3)
dr.coordinate.test(s3,~.-log(WCC),d=3)
