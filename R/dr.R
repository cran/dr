#####################################################################
#
#     Dimension reduction methods for R and Splus
#     Revised in July, 2004 by Sandy Weisberg and Jorge de la Vega
#     copyright 2001, 2004, Sanford Weisberg
#     sandy@stat.umn.edu
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#####################################################################
#
# Since R and Splus are not identical, different code is sometimes
# needed.  I have tested S-Plus 2000 and S-Plus 6.0 for Linux
# In R, the variable version$language is "R"; in SPlus, it is NULL

#####################################################################
#     dr is the primary function
#####################################################################
dr <- 
function(formula, data, subset, na.action=na.fail, weights, ...)
{  
    mf <- match.call(expand.dots=FALSE)
    mf$... <- NULL # ignore ...
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    mt <- attr(mf,"terms")
    Y <- model.extract(mf,"response")
    X <- model.matrix(mt, mf)
    y.name <- if(is.matrix(Y)) colnames(Y) else
        as.character(attr(mt, "variables")[2])
    int <- match("(Intercept)", dimnames(X)[[2]], nomatch=0)
    if (int > 0) X <- X[, -int, drop=FALSE] # drop the intercept from X
    weights <- mf$"(weights)"               
    if(is.null(weights)){                   # create weights if not present
         weights <- rep(1,dim(X)[1])} else
        {
         if(any(weights<0))stop("Negative weights")
         pos.weights <- weights > 1.e-30
         weights <- weights[pos.weights]
         weights <- dim(X)[1]*weights/sum(weights)
         X <- X[pos.weights,]
         Y <- if(is.matrix(Y)) Y[pos.weights,] else Y[pos.weights]}
    ans <- dr.compute(X,Y,weights=weights,...)
    ans$call <- match.call()
    ans$y.name <- y.name
    ans
}

dr.compute <-
function(x,y,weights,method="sir",...)
{ 
    if (NROW(y) != nrow(x))     #check lengths
        stop("The response and predictors have differing number of observations")
    classname<- if (is.matrix(y))    #Set the class name
        c(paste("m",method,sep=""),method) else method
    genclassname<-"dr"
    sweights <- sqrt(weights)
    qr.z <- qr(scale(apply(x,2,function(a, sweights) a  * sweights, sweights), 
                     center=TRUE, scale=FALSE))  # QR decomp of WEIGHTED x matrix
#initialize the object and then call the fit method
    ans <- list(x=x,y=y,weights=weights,method=method,cases=NROW(y),qr=qr.z)
    class(ans) <-  c(classname,genclassname)   #set the class
    ans <- dr.fit(object=ans,...)   # ... contains args for the method
    class(ans) <-  c(classname,genclassname) # reassign class name
    ans$x <- ans$x[,qr.z$pivot[1:qr.z$rank]] # reorder x and reduce to full rank
    ans #return the object
}
    
#######################################################
#    accessor functions
#######################################################

dr.x <- function(object) {object$x}
dr.wts <- function(object) {object$weights}
dr.qr <- function(object) {object$qr}
dr.Q <- function(object){ qr.Q(dr.qr(object))[,1:object$qr$rank] }
dr.R <- function(object){ qr.R(dr.qr(object))[1:object$qr$rank,1:object$qr$rank]}
dr.z <- function(object) { sqrt(object$cases) * dr.Q(object) }
dr.y.name <- function(object) {object$y.name}

#####################################################################
#     Fitting function
#####################################################################
dr.fit <- function(object,numdir=4,...) UseMethod("dr.fit")

dr.fit.default <-function(object,numdir=4,...){ 
    M <- dr.M(object,...)  # Get the kernel matrix M, method specific
    D <- if(dim(M$M)[1] == dim(M$M)[2]) eigen(M$M) else {
          if(ncol(M$M)==1) eigen(M$M %*% t(M$M)) else eigen(t(M$M) %*% M$M)} 
    or <- rev(order(abs(D$values)))
    evalues <- D$values[or]
    raw.evectors <- D$vectors[,or]  
    evectors <- backsolve(sqrt(object$cases)*dr.R(object),raw.evectors)
    evectors <- if (is.matrix(evectors)) evectors else matrix(evectors,ncol=1)
    evectors <- apply(evectors,2,function(x) x/sqrt(sum(x^2)))
    names <- colnames(dr.x(object))[1:object$qr$rank]
    dimnames(evectors)<-
         list(names, paste("Dir", 1:NCOL(evectors), sep=""))
    aa<-c( object, list(evectors=evectors,evalues=evalues, 
                numdir=min(numdir,dim(evectors)[2],dim(M$M)[1]),
                raw.evectors=raw.evectors), M)
    return(aa)
}

#####################################################################
###
###  dr methods. Each method REQUIRES a dr.M function, and
###  may also have a dr.y function and an function method.
###
#####################################################################

#####################################################################
###  generic functions
#####################################################################

dr.M <- function(object, ...){UseMethod("dr.M")}
dr.y <- function(object) {UseMethod("dr.y")}
dr.y.default <- function(object){object$y}
dr.test <- function(object, nd){  UseMethod("dr.test")}
dr.test.default <-function(object, nd) {NULL}

#####################################################################
#     OLS
#####################################################################

dr.M.ols <- function(object,...) {
 ols <- t(dr.z(object)) %*% (sqrt(dr.wts(object)) * dr.y(object))
 return(list( M= ols %*% t(ols), numdir=1))}
       

#####################################################################
#     Sliced Inverse Regression and multivariate SIR
#####################################################################

dr.M.sir <-function(object,nslices=NULL,slice.info=NULL,...) {
    z <- dr.z(object)
    y <- dr.y(object)
# get slice information    
    h <- if (!is.null(nslices)) nslices else max(8, NCOL(z)+3)
    slices<- if(is.null(slice.info)) dr.slices(y,h) else slice.info
# initialize slice means matrix
    zmeans <- matrix(0,slices$nslices,NCOL(z))
    slice.weight <- rep(0,slices$nslices)  # NOT rep(0,NCOL(z))
    wts <- dr.wts(object)
# compute weighted means within slice 
    wmean <- function (x, wts) { sum(x * wts) / sum (wts) }
    for (j in 1:slices$nslices){
      sel <- slices$slice.indicator==j
      zmeans[j,]<- apply(z[sel,],2,wmean,wts[sel])
      slice.weight[j]<-sum(wts[sel])}
# get M matrix for sir
    M <- t(zmeans) %*% apply(zmeans,2,"*",slice.weight)/ sum(slice.weight)
    return (list (M=M,slice.info=slices))
}

dr.M.msir <-function(...) {dr.M.sir(...)}

dr.test.sir<-function(object,nd) {
#compute the sir test statistic for the first nd directions
    e<-sort(object$evalues)
    p<-length(object$evalues)
    n<-object$cases
    st<-df<-pv<-0
    nt <- min(p,nd)
    for (i in 0:nt-1)
      {st[i+1]<-n*(p-i)*mean(e[seq(1,p-i)])
       df[i+1]<-(p-i)*(object$slice.info$nslices-i-1)
       pv[i+1]<-1-pchisq(st[i+1],df[i+1])
      }
    z<-data.frame(cbind(st,df,pv))
    rr<-paste(0:(nt-1),"D vs >= ",1:nt,"D",sep="")
    dimnames(z)<-list(rr,c("Stat","df","p-value")) 
    z
}


#####################################################################
#     Sliced Average Variance Estimation
#####################################################################

dr.M.save <-function(object,nslices=NULL,slice.info=NULL,...) {
    z <- dr.z(object)
    y <- dr.y(object)
    wts <- dr.wts(object)
# get slice information
    h <- if (!is.null(nslices)) nslices else max(8, ncol(z)+3)
    slices<- if(is.null(slice.info)) dr.slices(y,h) else slice.info
# initialize M
    M <- matrix(0,NCOL(z),NCOL(z))
    ws <-rep(0,slices$nslices)
# Compute weighted within-slice covariance matrices, skipping any slice with
# total weight smaller than 1
    wvar <- function (x, w){
     (if (sum(w) > 1) {(length(w)-1)/(sum(w)-1)} else {0}) *
              var(sweep(x,1,sqrt(w),"*"))}
# Compute M = sum (ws-1)(I-C)^2/sum(ws-1)
    for (j in 1:slices$nslices){
      ind <- slices$slice.indicator==j
      IminusC <- diag(rep(1,NCOL(z))) -  wvar(z[ind,],wts[ind])
      ws[j]<- sum(wts[ind]) # ws is the within slice sum of weights
      M <- M + ws[j]*IminusC %*% IminusC
      }
    M <- M/sum(ws)
    return (list (M=M,slice.info=slices))
}


#####################################################################
#     pHd, pHdy and pHdres
#####################################################################

dr.M.phdy <- function(...) {dr.M.phd(...)}
dr.M.mphd <- function(...) stop("Multivariate pHd not implemented!")

dr.M.phdres <- function(...) {dr.M.phd(...)}
dr.M.mphdres <- function(...) stop("Multivariate pHd not implemented!")

dr.M.phd <-function(object,...) {
    wts <- dr.wts(object)
    z <- dr.z(object)
    y <- dr.y(object)
    M<- (t(apply(z,2,"*",wts*y)) %*% z) / sum(wts)
    return(list(M=M))
}

dr.M.mphd <- function(...) stop("Multivariate pHd not implemented!")


dr.y.phdy <- function(object) {y <- object$y ; y - mean(y)}
dr.y.phdres <- function(object) { 
   y <- object$y
   sw <- sqrt(object$weights)
   qr.resid(object$qr,sw*(y-mean(y)))}
dr.y.phd <- function(object) {dr.y.phdres(object)}

dr.test.phd<-function(object,nd) {
# Modified by Jorge de la Vega, February, 2001
#compute the phd asymptotic test statitics under restrictions for the
#first nd directions, based on response = OLS residuals
# order the absolute eigenvalues
    e<-sort(abs(object$evalues))
    p<-length(object$evalues)
# get the response
    resi<-dr.y(object)
    varres<-2*var(resi)
    n<-object$cases
    st<-df<-pv<-0
    nt<-min(p,nd)
    for (i in 0:nt-1)
# the statistic is partial sums of squared absolute eigenvalues
      {st[i+1]<-n*(p-i)*mean(e[seq(1,p-i)]^2)/varres
# compute the degrees of freedom
       df[i+1]<-(p-i)*(p-i+1)/2
# use asymptotic chi-square distrbution for p-values.  Requires normality.
       pv[i+1]<-1-pchisq(st[i+1],df[i+1])
      }
# compute the test for complete independence
    indep <- dr.indep.test.phdres(object,st[1])
# compute tests that do not require normal theory (linear combination of 
# chi-squares):
    lc <- dr.test2.phdres(object,st)
# report results
    z<-data.frame(cbind(st,df,pv,c(indep[[2]],rep(NA,length(st)-1)),lc[,2]))
    rr<-paste(0:(nt-1),"D vs >= ",1:nt,"D",sep="")
    cc<-c("Stat","df","Normal theory","Indep. test","General theory")
    dimnames(z)<-list(rr,cc)
    z
}

dr.test.phdres<-function(object,nd){dr.test.phd(object,nd)}

dr.test.phdy<-function(object,nd) {
#compute the phd test statitics for the first nd directions
#based on response = y.  According to Li (1992), this requires normal
#predictors.
# order the absolute eigenvalues
    e<-sort(abs(object$evalues))
    p<-length(object$evalues)
# get the response
    resi<-dr.y(object)
    varres<-2*var(resi)
    n<-object$cases
    st<-df<-pv<-0
    nt<-min(p,nd)
    for (i in 0:nt-1)
# the statistic is partial sums of squared absolute eigenvalues
      {st[i+1]<-n*(p-i)*mean(e[seq(1,p-i)]^2)/varres
# compute the degrees of freedom
       df[i+1]<-(p-i)*(p-i+1)/2
# use asymptotic chi-square distrbution for p-values.  Requires normality.
       pv[i+1]<-1-pchisq(st[i+1],df[i+1])
      }
# report results
    z<-data.frame(cbind(st,df,pv))
    rr<-paste(0:(nt-1),"D vs >= ",1:nt,"D",sep="")
    cc<-c("Stat","df","p-value")
    dimnames(z)<-list(rr,cc)
    z}


#####################################################################
# phdq, Reference:  Li (1992, JASA)
# Function to build quadratic form in the full quadratic fit.
# Corrected phdq by Jorge de la Vega 7/10/01
#####################################################################

dr.M.phdq<-function(object,...){
 pars <- fullquad.fit(object)$coef
 k<-length(pars)
 p<-(-3+sqrt(9+8*(k-1)))/2 #always k=1+2p+p(p-1)/2
 mymatrix <- diag(pars[round((p+2):(2*p+1))]) # doesn't work in S without 'round'
 pars <- pars[-(1:(2*p+1))]
 for (i in 1:(p - 1)) {
      mymatrix[i,(i+1):p] <- pars[1:(p - i)]/2
      mymatrix[(i+1):p,i] <- pars[1:(p - i)]/2
      pars <- pars[-(1:(p - i))]
 }
 return(list(M=mymatrix))
 }

fullquad.fit <-function(object) {
 x <- dr.z(object)
 y <- object$y
 w <- dr.wts(object)
 z <- cbind(x,x^2)
 p <- NCOL(x)
 for (j in 1:(p-1)){
   for (k in (j+1):p) { z <- cbind(z,matrix(x[,j]*x[,k],ncol=1))}}
 lm(y~z,weights=w)
 }

 
dr.y.phdq<-function(object){residuals(fullquad.fit(object),type="pearson")}
dr.test.phdq<-function(object,nd){dr.test.phd(object,nd)}
dr.M.mphdq <- function(...) stop("Multivariate pHd not implemented!")

#####################################################################
#     Helper methods
#####################################################################


#####################################################################
#     Recover the direction vectors
#####################################################################

dr.directions <- function(object, which, norm, x) {UseMethod("dr.direction")}
dr.direction  <- function(object, which, norm, x) {UseMethod("dr.direction")}

dr.direction.default <- 
  function(object, which=1:object$numdir,norm=FALSE,x=dr.x(object)) {
    ans <- (apply(x,2,function(x){x-mean(x)}) %*% object$evectors)[,which]
    ans <- if (norm) apply(ans,2,function(x){x/(sqrt(sum(x^2)))}) else ans
    if (length(which) > 1) dimnames(ans) <-
                       list ( attr(x,"dimnames")[[1]],paste("Dir",which,sep=""))
    ans
  }
  
#####################################################################
#     Plotting methods
#####################################################################

plot.dr <- function
      (x,which=1:x$numdir,mark.by.y=FALSE,plot.method=pairs,...) {
 d <- dr.direction(x,which)
 if (mark.by.y == FALSE) {
    plot.method(cbind(dr.y(x),d),labels=c(dr.y.name(x),colnames(d)),...)
    }
 else {plot.method(d,labels=colnames(d),col=markby(dr.y(x)),...)}
 }

#############This does not work in Splus
# point marking in S is not easy because the arg col is expected to be an
# integer rather than a list
#########################################################################
#plot.dr <- function(object,which=1:object$numdir,mark.by.y=FALSE,
#                panel=points.col,colors,...) {
# d <- dr.direction(object,which)
# if (mark.by.y == FALSE) {
#    pairs(cbind(dr.y(object),d),labels=c(dr.y.name(object),colnames(d)),
#                 ,panel,colors=colors,...)
#    }
# else {pairs(d,labels=colnames(d),col=markby(dr.y(object)),...)}
# }

#points.col <- function(x,y,colors,...){
#  if (version$language == "R") {points(x,y,col=colors,...)}
#  else if (length(colors) == 1) points(x,y,col=colors,...)
#   else {
#    unique.colors <- unique(colors)
#    for (j in 1:length(unique.colors)){
#     sel<- colors == unique.colors[j]
#     points(x[sel],y[sel],col=unique.colors[j],...)}}}

dr.coplot <- function(object,which=1:object$numdir,mark.by.y=FALSE,...) {
 d <- data.frame(dr.direction(object,which))
 if (mark.by.y == FALSE){
      d$yvar <- object$y
      coplot(yvar~Dir1|Dir2,data=d,...)}
    else
     {coplot(Dir1~Dir2|Dir3,data=d,col=markby(object$y),...)}
}

givens.rotation <- function(theta,p=2,which=c(1,2)){
 m <- matrix(rep(0,p^2),nrow=p)
 m[which,which]<-
        matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),nrow=2)
 m}

#
# The following function rotplot gives static views of a 3D rotating
# plot.  A call might be 
#   rotplot(dr.directions(m1,1:2),dr.y(m1),number=16)
# Does not seem very useful...

rotplot <- function(x,y,number=9,theta=seq(0,pi/2,length=number),...){
 z<-NULL
 for (j in 1:number) { y
   z<- rbind(z, cbind(y, x %*% givens.rotation(theta[j]), theta[j]))}
 z <- data.frame(z)
 formula <- if (version$language == "R") y~V2|V4 else z.1~z.2|z.4
 coargs <- co.intervals(z[,4],number=number,overlap=0)
 coplot(formula,data=z,given.values=coargs,
          xlab=c("Linear combination","Angle (radians)"),
          ylab=deparse(substitute(y)),...)}
  
dr.persp<-function(object,which=1:2,h=c(.1,.1),...){
 require(sm) # uses the sm library
 if (length(which) == 2){
  d1<-dr.direction(object,which,norm=TRUE)
  y<-dr.y(object)
  sm.regression(d1,y,h=h,
                     xlab=dimnames(d1)[[2]][1],
                     ylab=dimnames(d1)[[2]][2],
                     zlab=names(object$model[1]))
  }
  else
  print("This method requires specifying two directions")
 }

markby <-
function(z,use="color",values=NULL,color.fn=rainbow,na.action="na.use") {
 u <- unique(z)
 lu <- length(u)
 ans <- 0
 vals <- if (use == "color") 
      {if (!is.null(values) && length(values) == lu)
                values else color.fn(lu)}
   else
      {if (!is.null(values) && length(values) == lu)
                   values else 1:lu}
 for (j in 1:lu)
      if (is.na(u[j])) 
         {ans[which(is.na(z))] <- 
         if(na.action=="na.use") vals[j] else NA} else
         {ans[z==u[j]] <- vals[j]}
 ans}
   

###################################################################
#
#  basic print method for dimension reduction
#
###################################################################
"print.dr" <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("Eigenvectors:\n")
    evectors<-x$evectors
    print.default(evectors)
    cat("Eigenvalues:\n")
    print.default(x$evalues)
    cat("\n")
    invisible(x)
}

###################################################################
#  basic summary method for dimension reduction
###################################################################
"summary.dr" <- function (object, ...)
{   z <- object
    ans <- z[c("call")]
    nd <- min(z$numdir,length(which(abs(z$evalues)>1.e-15)))
    ans$evectors <- z$evectors[,1:nd]
    ans$method <- z$method
    ans$nslices <- z$slice.info$nslices
    ans$sizes <- z$slice.info$slice.sizes
    ans$weights <- dr.wts(z)
    sw <- sqrt(ans$weights)
    y <- z$y #dr.y(z)
    ans$n <- z$cases #NROW(z$model)
    ols.fit <- qr.fitted(object$qr,sw*(y-mean(y)))
    angles <- cosangle(dr.direction(object),ols.fit)
    angles <- if (is.matrix(angles)) angles[,1:nd] else angles[1:nd]
    if (is.matrix(angles)) dimnames(angles)[[1]] <- z$y.name
    angle.names <- if (!is.matrix(angles)) "R^2(OLS|dr)" else
                        paste("R^2(",dimnames(angles)[[1]],"|dr)",sep="")
    ans$evalues <-rbind (z$evalues[1:nd],angles)
    dimnames(ans$evalues)<-
     list(c("Eigenvalues",angle.names),
          paste("Dir", 1:NCOL(ans$evalues), sep=""))
    ans$test <- dr.test(object,nd)
    class(ans) <- "summary.dr"
    ans
}

###################################################################
#
# basic print.summary method for dimension reduction
#
###################################################################
"print.summary.dr" <-
function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")#S: ' ' instead of '\n'
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    cat("Method:\n")#S: ' ' instead of '\n'
    if(is.null(x$nslices)){
       cat(paste(x$method, ", n = ", x$n,sep=""))
       if(!is.null(x$weights)) cat(", using weights.\n") else cat(".\n")}
       else {
         cat(paste(x$method," with ",x$nslices, " slices, n = ",
                   x$n,sep=""))
         if(!is.null(x$weights)) cat(", using weights.\n") else cat(".\n")
         cat("\nSlice Sizes:\n")#S: ' ' instead of '\n'
         cat(x$sizes,"\n")}
    cat("\nEigenvectors:\n")
    print(x$evectors,digits=digits)
    cat("\n")
    print(x$evalues,digits=digits)
    if (length(x$omitted) > 0){
      cat("\nModel matrix is not full rank.  Deleted columns:\n")
      cat(x$omitted,"\n")}
    if (!is.null(x$test)){
      cat("\nAsymp. Chi-square tests for dimension:\n")
      print(as.matrix(x$test),digits=digits)}
    invisible(x)
}



###################################################################3
##
##  Translation of methods in Arc for testing with pHd to R
##  Original lisp functions were mostly written by R. D. Cook
##  Translation to R by S. Weisberg, February, 2001
##
###################################################################3
# this function is a translation from Arc.  It computes the matrices W and
# eW described in Sec. 12.3.1 of Cook (1998), Regression Graphics.
# There are separate versions of this function for R and for Splus because
cov.ew.matrix <- function(object,scaled=FALSE) {
  mat.normalize <- function(a){apply(a,2,function(x){x/(sqrt(sum(x^2)))})}
  n <- dim(dr.x(object))[1]
  TEMPwts <- object$weights
  sTEMPwts <- sqrt(TEMPwts)
  v <- sqrt(n)* mat.normalize(
         apply(scale(dr.x(object),center=TRUE,scale=FALSE),2,"*",sTEMPwts) %*% 
               object$evectors)
  y <- dr.y(object) # get the response
  y <- if (scaled) y-mean(sTEMPwts*y) else 1 # a multiplier in the matrix
  p <- dim(v)[2]
  ew0 <- NULL
  for (i in 1:p){
   for (j in i:p){
    ew0 <- cbind(ew0, if (i ==j) y*(v[,j]^2-1) else y*sqrt(2)*v[,i]*v[,j])}}
  wmean <- function (x,w)  { sum(x * w) / sum (w) }
  tmp <- apply(ew0,2,function(x,w,wmean){sqrt(w)*(x-wmean(x,w))},TEMPwts,wmean)
  ans<-(1/sum(TEMPwts)) * t(tmp) %*% tmp
  ans} 

#translation of :general-pvalues method for phd in Arc
dr.test2.phdres <- function(object,stats){
  covew <- cov.ew.matrix(object,scaled=TRUE)
  C <- .5/var(dr.y(object))
  p <- length(stats)
  pval <- NULL
  d2 <-dim(dr.x(object))[2]
  start <- -d2
  end <- dim(covew)[2]
  for (i in 1:p) {
   start <- start + d2-i+2
   evals <- 
     eigen(as.numeric(C)*covew[start:end,start:end],only.values=TRUE)$values
   pval<-c(pval,wood.pvalue(evals,stats[i]))}
# report results
    z<-data.frame(cbind(stats,pval))
    rr<-paste(0:(p-1),"D vs >= ",1:p,"D",sep="")
    dimnames(z)<-list(rr,c("Stat","p-value"))
    z}
   
dr.indep.test.phdres <- function(object,stat) {
  eval <- eigen(cov.ew.matrix(object,scaled=FALSE),only.values=TRUE)
  pval<-wood.pvalue(.5*eval$values,stat)
# report results
    z<-data.frame(cbind(stat,pval))
    dimnames(z)<-list(c("Test of independence"),c("Stat","p-value"))
    z}


wood.pvalue <- function (coef, f, tol=0.0, print=FALSE){
#Returns an approximation to P(coef'X > f) for X=(X1,...,Xk)', a vector of iid
#one df chi-squared rvs.  coef is a list of positive coefficients. tol is used
#to check for near-zero conditions.
#See Wood (1989), Communications in Statistics, Simulation 1439-1456.
#Translated from Arc function wood-pvalue.
#  error checking
  if (min(coef) < 0) stop("negative eigenvalues")
  if (length(coef) == 1)
     {pval <- 1-pchisq(f/coef,1)} else
     {k1 <-     sum(coef)
      k2 <- 2 * sum(coef^2)
      k3 <- 8 * sum(coef^3)
      t1 <- 4*k1*k2^2 + k3*(k2-k1^2)
      t2 <- k1*k3 - 2*k2^2
    if ((t2 <= tol) && (tol < t2) ){
        a1 <- k1^2/k2
    b  <- k1/k2
    pval <- 1 - pgamma(b*f,a1)
        if (print) 
      print(paste("Satterthwaite-Welsh Approximation =", pval))}
      else if( (t1 <= tol) && (tol < t2)){
        a1 <-2 + (k1/k2)^2
    b  <- (k1*(k1^2+k2))/k2
    pval <- if (f < tol) 1.0 else 1 - pgamma(b/f,a1)
        if (print) print(paste("Inverse gamma approximation =",pval))}
      else if ((t1 > tol) && (t2 > tol)) {
        a1 <- (2*k1*(k1*k3 + k2*k1^2 -k2^2))/t1
     b <- t1/t2
    a2 <- 3 + 2*k2*(k2+k1^2)/t2
    pval <- 1-pbeta(f/(f+b),a1,a2)
        if (print) print(paste(
          "Three parameter F(Pearson Type VI) approximation =", pval))}
      else {
        pval <- -1
        if (print) print("Wood's Approximation failed")}}
   pval}


#########################################################################
#
# permutation tests for dimenison reduction
#
#########################################################################

dr.permutation.test <- function(object,npermute=50,numdir=object$numdir,
                                permute.weights=TRUE) {
   call <- object$call
   call[[1]] <- as.name("dr.compute")
   call$formula <- call$data <- call$subset <- call$na.action <- NULL
   x <- dr.directions(object) # predictors with a 'better' basis
   call$y <- object$y
   weights <- object$weights
# nd is the number of dimensions to test
   nd <- min(numdir,length(which(abs(object$evalues)>1.e-8))-1)
   nt <- nd + 1
# observed value of the test statistics = obstest
   obstest<-dr.permutation.test.statistic(object,nt)
# count and val keep track of the results and are initialized to zero
   count<-rep(0,nt)
   val<-rep(0,nt)
# main loop
   for (j in 1:npermute){   #repeat npermute times
    perm<-sample(1:object$cases)
    call$weights<- if (permute.weights) weights[perm] else weights    
# inner loop
    for (col in 0:nd){
# z gives a permutation of x.  For a test of dim = col versus dim >= col+1,
# all columns of x are permuted EXCEPT for the first col columns.
        call$x<-if (col==0) x[perm,] else   cbind(x[,(1:col)],x[perm,-(1:col)])
        iperm <- eval(call)
        val[col+1]<- dr.permutation.test.statistic(iperm,col+1)[col+1]
        }  # end of inner loop
# add to counter if the permuted test exceeds the obs. value
     count[val>obstest]<-count[val>obstest]+1
   }
# summarize
   pval<-(count)/(npermute+1)
   ans1 <- data.frame(cbind(obstest,pval))
   dimnames(ans1) <-list(paste(0:(nt-1),"D vs >= ",1:nt,"D",sep=""),
                         c("Stat","p-value"))
   ans<-list(summary=ans1,npermute=npermute)
   class(ans) <- "dr.permutation.test"
   ans
   }

"print.dr.permutation.test" <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
   cat("\nPermutation tests\nNumber of permutations:\n")
   print.default(x$npermute)
   cat("\nTest results:\n")
   print(x$summary,digits=digits) 
   invisible(x)
}

"summary.dr.permutation.test" <- function(...)
              {print.dr.permutation.test(...)}


#########################################################################
#
# dr.permutation.test.statistic method
#
#########################################################################

dr.permutation.test.statistic <- function(object,nd)
  {UseMethod("dr.permutation.test.statistic")}

dr.permutation.test.statistic.default <- function(object,nd){
   object$cases*rev(cumsum(rev(object$evalues)))[1:nd]}

dr.permutation.test.statistic.phdy <- function(object,nd){
       dr.permutation.test.statistic.phd(object,nd)}
dr.permutation.test.statistic.phdres <- function(object,nd){
       dr.permutation.test.statistic.phd(object,nd)}
dr.permutation.test.statistic.phd <- function(object,nd){
   (.5*object$cases*rev(cumsum(rev(object$evalues^2)))/var(dr.y(object)))[1:nd]}

#####################################################################
#
#     dr.slices returns non-overlapping slices based on y
#     y is either a list of n numbers or an n by p matrix
#     nslices is either the total number of slices, or else a
#     list of the number of slices in each dimension
#
#####################################################################

dr.slices <- function(y,nslices) {
  p <- if (is.matrix(y)) dim(y)[2] else 1
  h <- if (length(nslices) == p) nslices else rep(ceiling(nslices^(1/p)),p)
  a <- dr.slice.1d( if(is.matrix(y)) y[,1] else y, h[1])
  if (p > 1){
    for (col in 2:p) {
       ns <- 0
       for (j in unique(a$slice.indicator)) {
         b <- dr.slice.1d(y[a$slice.indicator==j,col],h[col])
         a$slice.indicator[a$slice.indicator==j] <- 
                a$slice.indicator[a$slice.indicator==j] + 10^(col-1)*b$slice.indicator
         ns <- ns + b$nslices}
       a$nslices <- ns }
#recode unique values to 1:nslices and fix up slice sizes
    v <- unique(a$slice.indicator)
    L <- NULL
    for (i in 1:length(v)) {
       sel <- a$slice.indicator==v[i]
       a$slice.indicator[sel] <- i
       L <- c(L,length(a$slice.indicator[sel]))}
    a$slice.sizes <- L }
  a}

# modified 11/8/03 for a categorical predictor, but this will still fail in
# some circumstances if length(z) > h if the first value of y has too many
# observations.
dr.slice.1d <- function(y,h) {
 z<-unique(y)
 if (length(z) >= h) dr.slice2(y,h) else dr.slice1(y,length(z),sort(z))}
# dr.slice2(y,h)}

dr.slice1 <- function(y,h,u){
  z <- sizes <- 0
  for (j in 1:length(u)) {
      temp <- which(y==u[j])
      z[temp] <- j
      sizes[j] <- length(temp) } 
  list(slice.indicator=z, nslices=length(u), slice.sizes=sizes)
  }

dr.slice2<-function(y,h)
{
  or <- order(y)
  n <- length(y)
  m<-floor(n/h)
  r<-n-m*h
  start<-sp<-ans<-0
  j<-1
  while((start+m)<n)
    { if (r==0)
        start<-start
      else 
        {start<-start+1
         r<-r-1
        }
       while (y[or][start+m]==y[or][start+m+1])
          start<-start+1
       sp[j]<-start+m
       start<-sp[j]
       j<-j+1
     }
# next line added 6/17/02 to assure that the last slice has at least 2 obs.
  if (sp[j-1] == n-1) j <- j-1
  sp[j]<-n
  ans[or[1:sp[1]]] <- 1
  for (k in 2:j){ans[ or[(sp[k-1]+1):sp[k] ] ] <- k}
  list(slice.indicator=ans, nslices=j, slice.sizes=c(sp[1],diff(sp)))
}

#####################################################################
#
#     Misc. Auxillary functions: cosine of the angle between two 
#     vectors and between a vector and a subspace.
#
#####################################################################

#
# angle between a vector vec and a subspace span(mat)
#
cosangle <- function(mat,vecs){
 ans <-NULL
 if (!is.matrix(vecs)) ans<-cosangle1(mat,vecs) else {
   for (i in 1:dim(vecs)[2]){ans <-rbind(ans,cosangle1(mat,vecs[,i]))}
   dimnames(ans) <- list(colnames(vecs),NULL) }
 ans}

cosangle1 <- function(mat,vec) {
  ans <- 
   cumsum((t(qr.Q(qr(mat))) %*% scale(vec)/sqrt(length(vec)-1))^2) 
# R returns a row vector, but Splus returns a column vector.  The next line
# fixes this difference
  if (version$language == "R") ans else t(ans)
}


#####################################################################
# R Functions for reweighting for elliptical symmetry
# modified from reweight.lsp for Arc
# Sanford Weisberg, sandy@stat.umn.edu
# March, 2001, rev June 2004

# Here is an outline of the function:
#   1.  Estimates of the mean m and covariance matrix S are obtained.  The
#       function cov.rob in the MASS package is used for this purpose.  
#   2.  The matrix X is replaced by Z = (X - 1t(m))S^{-1/2}.  If the columns
#       of X were from a multivariate normal N(m,S), then the rows of Z are
#       multivariate normal N(0, I).
#   3.  Randomly sample from N(0, sigma*I), where sigma is a tuning
#       parameter.  Find the row of Z closest to the random data, and increase
#       its weight by 1
#   4.  Return the weights divided by nsamples and multiplied by n.
#
#     dr.weights
#
#####################################################################
dr.weights <-
function (formula, data = list(), subset, na.action=na.fail, 
   covmethod="mve",sigma=1,nsamples=NULL,...)
{
# Create design matrix from the formula and call dr.estimate.weights
    mf1 <- match.call(expand.dots=FALSE)
    mf1$... <- NULL # ignore ...
    mf1$covmethod  <- mf1$nsamples <- NULL
    mf1[[1]] <- as.name("model.frame")
    mf <- eval(mf1, sys.frame(sys.parent()))
    mt <- attr(mf,"terms")
    x <- model.matrix(mt, mf)
    int <- match("(Intercept)", dimnames(x)[[2]], nomatch=0)
    if (int > 0) x <- x[, -int, drop=FALSE] # drop the intercept from X
    z <- robust.center.scale(x,method=covmethod,...)
    n <- dim(z)[1]   # number of obs
    p <- dim(z)[2]   # number of predictors
    ns <- if (is.null(nsamples)) 10*n else nsamples
    dist <- wts <- rep(0,n)  # initialize distances and weights
    for (i in 1:ns) {
       point <- rnorm(p) * sigma      # Random sample from a normal N(0, sigma^2I)
       dist <- apply(z,1,function(x,point){sum((point-x)^2)},point) 
#Dist to each point 
       sel <- dist == min(dist)               # Find closest point(s)
       wts[sel]<-wts[sel]+1/length(wts[sel])} # Increase weights for those points
    w <- n*wts/ns
    if (missing(subset)) return(w)
    if (is.null(subset)) return(w) else {
# find dimension of mf without subset specified
    mf1$subset <- NULL
    w1 <- rep(NA,length=dim(eval(mf1))[1])
    w1[subset] <- w
    return(w1)}}


robust.center.scale<- function (x,method, ...) { 
# This function takes an n by p matrix x, finds a robust center m and scale
# S, and then returns (x-1t(m)) %*% S^(-1/2).  
# This method uses the cov.rob method in the MASS package
# Additional args to the function are passed to cov.rob.
#
# make sure the library is loaded
 if(version$language == "R")
   {if(!require(MASS)){ 
      stop("MASS package not in library, but is required for this function")}} else
   {if(!exists("covRob")) library(MASS)}
 ans <- cov.rob(x,method=method, ...) 
 m <- if (method == "classical") apply(x,2,"median") else ans$center
 s<-svd(ans$cov)
 sweep(x,2,m) %*% s$u %*% diag(1/sqrt(s$d))}
 
#####################################################################
##
##  Add functions to Splus that are built-in to R
##
#####################################################################

if (version$language != "R") {

"is.empty.model" <- function (x)
{
    tt <- terms(x)
    (length(attr(tt, "factors")) == 0) & (attr(tt, "intercept")==0)
}
"NROW" <-
function(x) if(is.array(x)||is.data.frame(x)) nrow(x) else length(x)
"NCOL" <-
function(x) if(is.array(x)||is.data.frame(x)) ncol(x) else as.integer(1)

"colnames" <-
function(x, do.NULL = TRUE, prefix = "col")
{
    dn <- dimnames(x)
    if(!is.null(dn[[2]]))
    dn[[2]]
    else {
    if(do.NULL) NULL else paste(prefix, seq(length=NCOL(x)), sep="")
    }
}

"getOption" <- function(x) options(x)[[1]]
# end of special functions
}

# END OF FILE
