setwd('Z:/R/source/dr3')
source("work/dr.R")
source("work/psir.R")
source("work/psave.R")
data(ais,package="dr")
summary(s0 <- dr(LBM~log(SSF)+log(Wt)+log(Hg)+log(Ht)+log(WCC)+log(RCC)+
  log(Hc)+log(Ferr),data=ais,slice.function=dr.slices.arc,
  chi2approx="wood",numdir=4,method="sir"))
summary(s1 <- update(s0,group=~Sex))
summary(s2 <- dr(LBM~log(SSF)+log(Wt)+log(Hg)+log(Ht)+log(WCC)+log(RCC)+
  log(Hc)+log(Ferr),data=ais,slice.function=dr.slices.arc,group=~Sex,
  chi2approx="wood",numdir=4,method="sir"))
