#####################################################################
#
#     Dimension reduction methods for R and Splus
#     Written in July, 2000 by Sandy Weisberg
#     copyright 2001, Sanford Weisberg
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
# needed for the various engines (I have not tested Splus 3.0)
# whichengine provides a switch for the various engines.
# The code reads the variable thisengine, to get whichengine.
# The default is whichengine <- "R"
#
#####################################################################

# choices for whichengine are:	
# "R" for the R system
# "s6" for Splus 6  (works as of August 2001)
# "s2k" for splus 2000 
# "s5" for Splus 5 

whichengine <- "R"

#####################################################################
#
#     Splus class compatibility, not for R
#
#####################################################################

if (whichengine == "s5") {
  setOldClass(c("sir","dr"))
  setOldClass(c("save","dr"))
  setOldClass(c("phd","dr"))
  setOldClass(c("phdy","dr"))
  setOldClass(c("phdres","dr"))
  setOldClass(c("phdq","dr"))
  setOldClass(c("msir","sir","mdr"))
  setOldClass(c("msave","save","mdr"))
  setOldClass(c("mphdy","mdr"))
  setOldClass(c("mphdq","mdr"))
  setOldClass(c("mphdres","mdr"))}

#####################################################################
#
#     dr is the primary function
#
#####################################################################
dr <-
function (formula, data = list(), subset, weights=NULL, na.action=na.omit,
      method = "sir", contrasts = NULL, offset = NULL, ...)
{
#this first section is copied from the lm function
    mf <- match.call(expand.dots=FALSE) #puts args in a standard form
    mf$... <- NULL   # drop ...
    mf$method <- NULL # do not pass method to the model.frame
    mf$contrasts <- NULL # do not pass the contrasts to the model.frame
    mf[[1]] <- as.name("model.frame")
    mf$drop.unused.levels <- TRUE  #remove lin. dep. columns
    mf <- eval(mf, sys.frame(sys.parent()))
    mt <- attr(mf,"terms")
    xvars <- as.character(attr(mt, "variables"))[-1]
    if((yvar <- attr(mt, "response")) > 0) xvars <- xvars[-yvar]
    xlev <- if(length(xvars) > 0) {
                       xlev <- lapply(mf[xvars], levels)
                       xlev[!sapply(xlev, is.null)] }
    y <- model.extract(mf, "response")
    ynames <- colnames(y)
    offset <- mf$"(offset)"
    if(!is.null(offset) && length(offset) != NROW(y))
     stop(paste("Length of the offset", length(offset),
         ", must equal", NROW(y), " the number of cases"))
    if (is.empty.model(mt)) stop(paste("No model specified!"))
#Set the class name
    classname<- if (is.matrix(y))
                    c(paste("m",method,sep=""),method) else method
    genclassname<-"dr"
# define z to be an object of class classname, so model matrix will work
# correctly
    z <- list()
    switch(whichengine,
       R=class(z) <- c(classname,genclassname),
       s2k=class(z) <- c(classname,genclassname),
       s6=class(z) <- c(classname,genclassname),
       s5=oldClass(z) <- classname)
#set up the model matrix, x
    x <- model.matrix(z, mt, mf, contrasts)
#check lengths
    if (NROW(y) != nrow(x))
     stop("The response and predictors have differing number of observations")
#work with the weights
    w <- mf$"(weights)"
#if any weights are equal to zero, set to NA, issue a warning, and continue
    if(!is.null(w)){
        w<- length(w)*w/sum(w)  # scale weights to length n
        wsel<- w <= 0
        if (any(wsel,na.rm=TRUE)){
            w[wsel] <- NA # set zero weights to missing
            warning("Zero weight cases set to NA")}}
#scale weights to add to n
    wts <- if(is.null(w)) rep(1,NROW(y)) else NROW(y) * w/sum(w) 
#set the offset to be zero if null
    off <- if(is.null(offset)) rep(0,NROW(y)) else offset
#dr.z centers, rotates, and get the rank.  z$z is full column rank
    z <- dr.z(x,wts,center=TRUE,rotate=TRUE,decomp="qr") 
#cols gives the columns that give a full-rank parameterization
    cols <- z$QR$pivot[1:z$QR$rank]
#fitval are the WLS fitted values.  
#The fitval are correct even with zero case weights
    ols.coef<-qr.coef(z$QR,sqrt(wts)*(y-off))
    fitval <-off+sum(wts*y-off)/sum(wts)+dr.z(x,wts,rotate=FALSE)$z %*% ols.coef
#update the object and then call the fit method
    z <- list(formula=formula,contrasts= attr(x, "contrasts"),
              xlevels=xlev, call = match.call (),
          ols.coef=ols.coef, ols.fit=fitval,
          weights=wts,cols.used=cols, offset=offset,
          terms=mt, method=method, response.name=ynames,
          model=mf, cases= NROW(y))
#set the class again, as it has been lost
    switch(whichengine,
       R=class(z) <- c(classname,genclassname),
       s2k=class(z) <- c(classname,genclassname),
       s6=class(z) <- c(classname,genclassname),
       s5=oldClass(z) <- classname)
    z1 <- dr.fit(object=z,...)
#assign the results to names
    z <- c(z, z1)
#since z has been reassigned, assign its class again
    switch(whichengine,
       R=class(z) <- c(classname,genclassname),
       s2k=class(z) <- c(classname,genclassname),
       s6=class(z) <- c(classname,genclassname),
       s5=oldClass(z) <- classname)
#return the object
    z }

#####################################################################
#
#     dr.weights
#
#####################################################################
dr.weights <-
function (formula, data = list(), subset, weights=NULL, na.action=na.omit,
      method = "sir", contrasts = NULL, offset = NULL, ...)
{
#this first section is copied from the lm function
    mf <- match.call(expand.dots=FALSE) #puts args in a standard form
    mf$... <- NULL   # drop ...
    mf$method <- NULL # do not pass method to the model.frame
    mf$contrasts <- NULL # do not pass the contrasts to the model.frame
    mf[[1]] <- as.name("model.frame")
    mf$drop.unused.levels <- TRUE  #remove lin. dep. columns
    mf <- eval(mf, sys.frame(sys.parent()))
    mt <- attr(mf,"terms")
    xvars <- as.character(attr(mt, "variables"))[-1]
    if((yvar <- attr(mt, "response")) > 0) xvars <- xvars[-yvar]
    xlev <- if(length(xvars) > 0) {
                       xlev <- lapply(mf[xvars], levels)
                       xlev[!sapply(xlev, is.null)] }
    y <- model.extract(mf, "response")
    ynames <- colnames(y)
    offset <- mf$"(offset)"
    if(!is.null(offset) && length(offset) != NROW(y))
     stop(paste("Length of the offset", length(offset),
         ", must equal", NROW(y), " the number of cases"))
    if (is.empty.model(mt)) stop(paste("No model specified!"))
#Set the class name
    classname<- if (is.matrix(y))
                    c(paste("m",method,sep=""),method) else method
    genclassname<-"dr"
# define z to be an object of class classname, so model matrix will work
# correctly
    z <- list()
    switch(whichengine,
       R=class(z) <- c(classname,genclassname),
       s2k=class(z) <- c(classname,genclassname),
       s6=class(z) <- c(classname,genclassname),
       s5=oldClass(z) <- classname)
#set up the model matrix, x
    x <- model.matrix(z, mt, mf, contrasts)
#check lengths
    if (NROW(y) != nrow(x))
     stop("The response and predictors have differing number of observations")
#estimate weights
     w<-dr.estimate.weights(x,...) # compute the estimated weights
#if any weights are equal to zero, set to NA, issue a warning, and continue
    if(!is.null(w)){
        w<- length(w)*w/sum(w)  # scale weights to length n
        wsel<- w <= 0
        if (any(wsel,na.rm=TRUE)){
            w[wsel] <- NA # set zero weights to missing
            warning("Zero weight cases set to NA")}}
#if weights.only=T, return only the value of w and exit; otherwise, continue
    {return(w)}}

model.matrix.dr <- function(object,...) {
    mat <- model.matrix(...)
    int <- match("(Intercept)", dimnames(mat)[[2]], nomatch=0)
    if (int > 0) mat <- mat[, -int, drop=FALSE]
    mat}       


#####################################################################
#
#     Fitting function
#
#####################################################################
dr.fit <-function(object,numdir=4,tol=1.e-7,...){
    decomp <- "svd"
# Use singular value decomposition for all calculations
    x <- dr.x(object)
    z <- dr.z(x,object$weights,center=TRUE,rotate=TRUE,decomp="svd")
    cols <- object$cols.used
# compute M, with a different function for each method
    yvar <- dr.fit.y(object)
    M <- dr.fit.M(object=object,z=z$z,y=yvar,w=object$weights,...)
# compute singular values and vectors of M
    D <- svd(M$M,nrow(M$M),0)
    evalues <- if (nrow(M$M)==ncol(M$M)) D$d else (D$d)^2
    or <- rev(order(abs(evalues)))  # order absolute eigenvalues
    evalues<-evalues[or]
    raw.evectors <- D$u[,or]
    evectors<-switch(decomp,
      svd=apply( (z$v %*% diag(1/z$d) %*% D$u)[,or],2,
                          function(x){x/(sqrt(sum(x^2)))}),
      qr=apply( solve(qr.R(z$QR)[cols,cols],D$u)[,or],2,
                          function(x){x/(sqrt(sum(x^2)))}))
# assign names to eigen arrays
    dimnames(evectors)<-
         list(colnames(x)[cols], paste("Dir", 1:NCOL(evectors), sep=""))
    aa<-c( list(evectors=evectors,evalues=evalues, 
                numdir=min(numdir,dim(evectors)[2]),
                raw.evectors=raw.evectors,decomp=decomp), M)
    return(aa)
}

#####################################################################
#
#     dr.fit.M functions, one for each fitting method
#     dr.fit.y functions, one for each fitting method
#
#####################################################################

dr.fit.M <- switch(whichengine,
            s5=function(object, ...){UseMethod("dr.fit.M")},
            s2k=function(object, ...){UseMethod("dr.fit.M")},
            s6=function(object, ...){UseMethod("dr.fit.M")},
            R= function(object, ...){UseMethod("dr.fit.M",object)})
dr.fit.y <- switch(whichengine,
            s5=function(object){UseMethod("dr.fit.y")},
            s2k=function(object){UseMethod("dr.fit.y")},
            s6=function(object){UseMethod("dr.fit.y")},
            R= function(object){UseMethod("dr.fit.y",object)})

# default is used for sir and save
dr.fit.y.default<-function(object) { 
  dr.y(object) - if (is.null(object$offset)) 0 else object$offset}

#####################################################################
#     Sliced Inverse Regression
#####################################################################

dr.fit.M.sir <-function(object,z,y,w=NULL,nslices=NULL,
                        slice.info=NULL,...) {
# get slice information
    h <- if (!is.null(nslices)) nslices else max(8, NCOL(z)+3)
    slices<- if(is.null(slice.info)) dr.slices(y,h) else slice.info
# initialize slice means matrix
    zmeans <- matrix(0,slices$nslices,NCOL(z))
    slice.weight <- rep(0,slices$nslices)  # NOT rep(0,NCOL(z))
# make sure weights add to n
    wts <- if(is.null(w)) rep(1,NROW(z)) else NROW(z) * w /sum(w)
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

dr.fit.M.msir <-function(...) {dr.fit.M.sir(...)}


#####################################################################
#     Sliced Average Variance Estimation
#####################################################################

dr.fit.M.save <-function(object,z,y,w=NULL,nslices=NULL,
                        slice.info=NULL,...) {
# get slice information
    h <- if (!is.null(nslices)) nslices else max(8, ncol(z)+3)
    slices<- if(is.null(slice.info)) dr.slices(y,h) else slice.info
# initialize M
    M <- matrix(0,NCOL(z),NCOL(z))
    ws <-rep(0,slices$nslices)
# Compute weighted within-slice covariance matrices, skipping any slice with
# total weight smaller than 1
    wts <- if(is.null(w)) rep(1,NROW(z)) else NROW(z) * w /sum(w)
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

dr.fit.y.phdy <- function(object){y <- dr.y(object) ; y - mean(y)}
dr.fit.M.phdy <- function(...) {dr.fit.M.phd(...)}
dr.fit.M.mphd <- function(...) stop("Multivariate pHd not implemented!")

dr.fit.y.phdres <- function(...) {dr.fit.y.phd(...)}
dr.fit.M.phdres <- function(...) {dr.fit.M.phd(...)}
dr.fit.M.mphdres <- function(...) stop("Multivariate pHd not implemented!")

dr.fit.y.phd <- function(object){
  y<-dr.y(object) ; ols.fit <- object$ols.fit; y-ols.fit-mean(y-ols.fit)}
dr.fit.M.phd <-function(object,z,y,w=NULL,...) {
# get the weights, and be sure they add to n
    wts <- if(is.null(w)) rep(1,NROW(z)) else {NROW(z) * w / sum(w)}
# compute M 
    M<- (t(apply(z,2,"*",wts*y)) %*% z) / sum(wts)
    return(list(M=M))
}
dr.fit.M.mphd <- function(...) stop("Multivariate pHd not implemented!")

#####################################################################
# Function to build quadratic form in the full quadratic fit.
# Corrected phdq by Jorge de la Vega 7/10/01
#####################################################################
# Reference:  Li (1992, JASA)

fullquad.fit <-function(x,y,w,...) {
 z <- cbind(1,x,x^2)
 p <- NCOL(x)
 for (j in 1:(p-1)){
   for (k in (j+1):p) { z <- cbind(z,matrix(x[,j]*x[,k],ncol=1))}}
 if(whichengine == "R")
  {if (is.null(w)) lm.fit(z,y,...) else lm.wfit(z,y,w,method="qr",...)} else
  {if (is.null(w)) lm.fit(z,y) else lm.wfit(z,y,w,method="qr")}}

make.Mhat.phdq<-function (pars)
{   k<-length(pars)
    p<-(-3+sqrt(9+8*(k-1)))/2 #always k=1+2p+p(p-1)/2
    mymatrix <- diag(pars[(p+2):(2*p+1)])
    pars <- pars[-(1:(2*p+1))]
    for (i in 1:(p - 1)) {
      mymatrix[i,(i+1):p] <- pars[1:(p - i)]/2
      mymatrix[(i+1):p,i] <- pars[1:(p - i)]/2
      pars <- pars[-(1:(p - i))]
    }
    mymatrix
}

fullquad.residuals <-function(...){fullquad.fit(...)$residuals}

dr.fit.M.phdq<-function(object,z,y,w=NULL,...){
 x <- dr.x(object)
 y <- dr.y(object)
 M<-make.Mhat.phdq(fullquad.fit(x,y,offset=object$offset,
                                 w=object$weights)$coef)
 return(list(M=M))
 }

dr.fit.y.phdq<-function(object){
 x <- dr.x(object)
 y <- dr.y(object)
 fullquad.residuals(x,y,offset=object$offset,w=object$weights)}

dr.test.phdq<-function(object,nd){dr.test.phd(object,nd)}
dr.fit.M.mqphd <- function(...) stop("Multivariate pHd not implemented!")

#####################################################################
#
#     Accessor functions for dimension reduction
#
#####################################################################

# recover and x and y data

dr.x <- function(object)
  {model.matrix(object,object$terms,object$model,object$contrasts)[,object$cols.used]}
dr.y <- function(object){ model.extract(object$model, "response")}

dr.y.name <- function(object){
  which <- attr(attr(object$model, "terms"),"response")
  sel <- switch(whichengine,R=which+1,s5=which,s2k=which,s6=which)
  as.character(attr(attr(object$model, "terms"), "variables")[sel])
  }

dr.x.omitted <- function(z) 
   {dimnames(model.matrix(z,z$terms,z$model,z$contrasts))[[2]][-z$cols.used]}


dr.z <- function(x,weights=NULL,center=TRUE,rotate=TRUE,decomp="svd"){
   if (is.null(class(x))) {z<-x; wts<-weights; n<-dim(x)[1]} else 
      if (class(x)== "matrix") {z<-x; wts<-weights; n<-dim(x)[1]} else 
      if (class(x)== "model.matrix") {z<-x; wts<-weights; n<-dim(x)[1]} else 
         {z <- dr.x(x); wts<-x$weights; n<-NROW(z)}
   wts <- if (is.null(wts)) rep(1,dim(z)[1]) else n*wts/sum(wts)
# always scale the weights to add to n
   swts <- sqrt(wts)
   ssumwts <- sqrt(sum(wts))
   wcenter <- function(x,wts) {x - sum(x*wts)/sum(wts)}
   z <- if (center) apply(z,2,wcenter,wts) else x
   if (rotate){
# Singular value decomposition
    if (decomp == "svd"){
      SVD <- svd(z * swts,nu=0) # temporary variable
      InvSqrtSigma <- SVD$v %*% diag(1/SVD$d) * ssumwts
# The next statement sets Z = W^{1/2}(X - 1t(m))Sigma^{-1/2}
      z <- apply(z,2,"*",swts) %*% InvSqrtSigma  # centered AND rotated.
      list(z=z,v=SVD$v,d=SVD$d,decomp="svd")} else {
# QR factorization
      QR <- qr(z * swts)
      rank <-QR$rank
      cols <-QR$pivot[1:rank]  
      z <- qr.Q(QR)[,cols]*ssumwts
      list(z=z,QR=QR, decomp="qr")}}  else
    {list(z=z,decomp=NULL)}} 


# recover the direction vectors

dr.directions <- function(object, ...) {UseMethod("dr.direction")}
dr.direction <- function(object, ...) {UseMethod("dr.direction")}

dr.direction.default <- 
  function(object, which=1:object$numdir,norm=FALSE,x=dr.x(object)) {
    ans <- (apply(x,2,function(x){x-mean(x)}) %*% object$evectors)[,which]
    ans <- if (norm) apply(ans,2,function(x){x/(sqrt(sum(x^2)))}) else ans
    if (length(which) > 1) dimnames(ans) <-
                       list ( attr(x,"dimnames")[[1]],paste("Dir",which,sep=""))
    ans
  }

  
#####################################################################
#
#     Plotting methods
#
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
#  if (whichengine == "R") {points(x,y,col=colors,...)}
#  else if (length(colors) == 1) points(x,y,col=colors,...)
#   else {
#    unique.colors <- unique(colors)
#    for (j in 1:length(unique.colors)){
#     sel<- colors == unique.colors[j]
#     points(x[sel],y[sel],col=unique.colors[j],...)}}}

dr.coplot <- function(object,which=1:object$numdir,mark.by.y=FALSE,...) {
 d <- data.frame(dr.direction(object,which))
 if (mark.by.y == FALSE){
      d$yvar <- dr.y(object)
      coplot(yvar~Dir1|Dir2,data=d,...)}
    else
     {coplot(Dir1~Dir2|Dir3,data=d,col=markby(dr.y(object)),...)}
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
 formula <- if (whichengine == "R") y~V2|V4 else z.1~z.2|z.4
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
#
#  basic summary method for dimension reduction
#
###################################################################
"summary.dr" <- function (object, ...)
{
#   z <- switch(whichengine,R=,.Alias(object),s5=,object,s2k=,object,s6=,object)
    z <- object
    ans <- z[c("call", "terms")]
    nd <- min(z$numdir,length(which(abs(z$evalues)>1.e-15)))
    ans$evectors <- z$evectors[,1:nd]
    ans$method <- z$method
    ans$nslices <- z$slice.info$nslices
    ans$sizes <- z$slice.info$slice.sizes
    ans$n <- NROW(z$model)
    angles <- cosangle(dr.direction(object),object$ols.fit)
    angles <- if (is.matrix(angles)) angles[,1:nd] else angles[1:nd]
    if (is.matrix(angles)) dimnames(angles)[[1]] <- z$response.name
    angle.names <- if (!is.matrix(angles)) "R^2(OLS|dr)" else
                        paste("R^2(",dimnames(angles)[[1]],"|dr)",sep="")
    ans$omitted <- dr.x.omitted(z)
    ans$evalues <-rbind (z$evalues[1:nd],angles)
    ans$weights <- z$weights 
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
    cat("Terms:\n")#S: ' ' instead of '\n'
    if (whichengine == "R"){
      cat(paste(deparse(x$terms), sep="\n", collapse = "\n"), "\n\n", sep="")}
    else if (whichengine == "s5"){
      print(attr(x$terms,"formula"))
      cat("\n")}
    else if (whichengine == "s2k"){
      print(attr(x$terms,"formula"))
      cat("\n")}
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
      print(as.matrix(x$test),digits)}
    invisible(x)
}

###################################################################
#
# Asymptotic test methods for dimenison reduction.  Separate 
# functions for each computing method.
#
###################################################################
switch(whichengine,
  R= dr.test <- function(object, ...){ UseMethod("dr.test",object)},
  s5=dr.test <- function(object, ...){ UseMethod("dr.test")},
  s2k=dr.test <- function(object, ...){ UseMethod("dr.test")},
  s6=dr.test <- function(object, ...){ UseMethod("dr.test")}
)

dr.test.default <-function(object, nd) {NULL}

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

#########################################################################
#
# Asymptotic test for phdres (no test for phdy
# Modified by Jorge de la Vega, February, 2001
#
#########################################################################

dr.test.phd<-function(object,nd) {
#compute the phd asymptotic test statitics under restrictions for the
#first nd directions, based on response = OLS residuals
# order the absolute eigenvalues
    e<-sort(abs(object$evalues))
    p<-length(object$evalues)
# get the response
    resi<-dr.fit.y(object)
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
# chi-squares:
    lc <- dr.test2.phdres(object,st)
# report results
    z<-data.frame(cbind(st,df,pv,c(indep[[2]],rep(NA,length(st)-1)),lc[,2]))
    rr<-paste(0:(nt-1),"D vs >= ",1:nt,"D",sep="")
    cc<-c("Stat","df","Normal theory","Indep. test","General theory")
    dimnames(z)<-list(rr,cc)
#   z<-data.frame(cbind(c(indep[[1]],st),c(NA,df),c(indep[[2]],pv)))
#   rr<-c("Test of Independence",paste(0:(nt-1),"D vs >= ",1:nt,"D",sep=""))
#   dimnames(z)<-list(rr,c("Stat","df","p-value"))
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
    resi<-dr.fit.y(object)
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
# I do not know how to make a local function get an argument from the
# function without passing an argument.
cov.ew.matrix <- function(object,scaled=FALSE) {
  n <- dim(dr.x(object))[1]
  TEMPwts <- if (is.null(object$weights)) rep(1,n) else 
                  {n*object$weights/sum(object$weights)}
  sTEMPwts <- sqrt(TEMPwts)
# The products Z %*% theta where theta are the eigenvectors in the Z scale
# are the same as the products X %*% evectors in the original scale, after
# centering X and then scaling the result so columns have length sqrt(wts).
# This avoids computing two singular value decompositions
  v <- sqrt(n)* mat.normalize(
         apply(dr.z(object,rotate=FALSE)$z,2,"*",sTEMPwts) %*% object$evectors)
  y <- dr.fit.y(object) # get the response
  y <- if (scaled) y-mean(y) else 1 # a multiplier in the matrix
  p <- dim(v)[2]
  ew0 <- NULL
  for (i in 1:p){
   for (j in i:p){
    ew0 <- cbind(ew0, if (i ==j) y*(v[,j]^2-1) else y*sqrt(2)*v[,i]*v[,j])}}
  if(whichengine != "R") {TEMPwts <<- TEMPwts}
  if(whichengine != "R") {sTEMPwts <<- sTEMPwts}
  if(whichengine != "R")
        {wmean <<- function (x) { sum(x * TEMPwts) / sum (TEMPwts) }} else
        {wmean <- function (x) { sum(x * TEMPwts) / sum (TEMPwts) }}
  tmp <- apply(ew0,2,function(x){sTEMPwts*(x-wmean(x))})
  ans<-(1/sum(TEMPwts)) * t(tmp) %*% tmp
  if(whichengine != "R"){ rm(TEMPwts);rm(sTEMPwts);rm(wmean)}
  ans} 

#translation of :general-pvalues method for phd in Arc
dr.test2.phdres <- function(object,stats){
  covew <- cov.ew.matrix(object,scaled=TRUE)
  C <- .5/var(dr.fit.y(object))
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
# nd is the number of dimensions to test for
   nd <- min(numdir,length(which(abs(object$evalues)>1.e-8))-1)
   nt <- nd + 1
# y is generally the response variable, but it is residuals for phdres
   y <- dr.fit.y(object)
   n <- NROW(y)
# observed value of the test statistics = obstest
   obstest<-dr.permutation.test.statistic(object,object$evalues,nt,
                          object$cases,y)
# make sure weights are a list and not just NULL
   wts <- if(is.null(object$weights)){rep(1,n)} else {object$weights}
# z2 is a matrix like X consisting of the principal directions.  z2 and
# x span the same subspace and are therefore equivalent.
   z2<- dr.direction(object,which=1:(NCOL(dr.x(object))))
# count and val keep track of the results and are initialized to zero
   count<-rep(0,nt)
   val<-rep(0,nt)
# main loop
   for (j in 1:npermute){                         #repeat npermute times
    perm<-permute(n)                #perm is a permutation of subscripts
# The weights are permuted if permute.weights is TRUE; otherwise, they
# are not permuted.
    wp<- if (permute.weights) wts[perm] else wts    
# inner loop
    for (col in 0:nd){
# z gives a permutation of z2.  For a test of dim = col versus dim >= col+1,
# all columns of z2 are permuted EXCEPT for the first col columns.
        z<-if (col==0) z2[perm,] else     
              cbind(z2[,(1:col)],z2[perm,-(1:col)]) 
# normalize (center, scale and rotate) the permuted data
        z<-dr.z(z,wp)$z      # uses singular value decomposition 
# the eigenvalues/singular values of M determine the test statistic
        M <- dr.fit.M(object,z,y,wp,method=object$method,
                        slice.info=object$slice.info)$M 
# get the eigenvalues
        evalues <- if(NROW(M)==NCOL(M)) rev(sort(abs(svd(M,0,0)$d)))
                      else rev(sort(abs(svd(M,0,0)$d^2))) #get eigenvalues
        val[col+1]<-
            dr.permutation.test.statistic(object,evalues,col+1,n,y)[col+1]
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

dr.permutation.test.statistic <- 
  switch(whichengine,
  R= function(object,...)
                     {UseMethod("dr.permutation.test.statistic",object)},
  s5=function(object,...)
                     {UseMethod("dr.permutation.test.statistic")},
  s2k=function(object,...)
                     {UseMethod("dr.permutation.test.statistic")},
  s6=function(object,...)
                     {UseMethod("dr.permutation.test.statistic")})

dr.permutation.test.statistic.default <- function(object,values,nd,n,...){
   n*rev(cumsum(rev(values)))[1:nd]}

dr.permutation.test.statistic.phdy <- function(object,...){
       dr.permutation.test.statistic.phd(object,...)}
dr.permutation.test.statistic.phdres <- function(object,...){
       dr.permutation.test.statistic.phd(object,...)}
dr.permutation.test.statistic.phd <- function(object,values,nd,n,y){
   (.5*n*rev(cumsum(rev(values^2)))/var(y))[1:nd]}

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

dr.slice.1d <- function(y,h) {
# z<-unique(y)
# if (length(z) >= h) dr.slice2(y,h) else dr.slice1(y,length(z),sort(z))}
  dr.slice2(y,h)}

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
#     Auxillary functions
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
# For reasons I do not understand, qr.qty does not seem to work.
# I have used a work-around that actually requires forming Q explictly
# cumsum(qr.qty(qr(mat),scale(vec)/sqrt(length(vec)-1))^2) 
  ans <- 
   cumsum((t(qr.Q(qr(mat))) %*% scale(vec)/sqrt(length(vec)-1))^2) 
# R returns a row vector, but Splus returns a column vector.  The next line
# fixes this difference
  if (whichengine == "R") ans else t(ans)
}

#
# random permutation
#

permute <- function(n){order(runif(n))}

mat.normalize <- function(a){apply(a,2,function(x){x/(sqrt(sum(x^2)))})}

#####################################################################
##
##  Add functions to Splus that are built-in to R
##
#####################################################################

if (whichengine != "R") {

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
# end
}
#####################################################################
##
##  End of functions in R added to Splus
##
#####################################################################

# R Functions for reweighting for elliptical symmetry
# modified from reweight.lsp for Arc
# Sanford Weisberg, sandy@stat.umn.edu
# March, 2001

# Here is an outline of the function:
#   1.  Estimates of the mean m and covariance matrix S are obtained.  The
#       function cov.rob in the lqs package is used for this purpose in R.  
#       For Splus, the function covRob in the robust package is used All
#       the args of this function in ... are passed to cov.rob or covRob.
#   2.  The matrix X is replaced by Z = (X - 1t(m))S^{-1/2}.  If the columns
#       of X were from a multivariate normal N(m,S), then the rows of Z are
#       multivariate normal N(0, I).
#   3.  Randomly sample from N(0, sigma*I), where sigma is a tuning
#       parameter.  Find the row of Z closest to the random data, and increase
#       its weight by 1
#   4.  Return the weights divided by nsamples and multiplied by n.

dr.estimate.weights <- function(object,sigma=1,
                          covmethod="mve",nsamples=NULL,...) {
 x <- if(is.null(class(object))) object else {
          if (class(object) == "data.frame") as.matrix(object) else
	  if (class(object) == "model.matrix") as.matrix(object) else
          dr.x(object)}
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
 n*wts/ns}                                 # Return weights, normalized to n

robust.center.scale<- function (x,method, ...) {
# This function takes an n by p matrix x, finds a robust center m and scale
# S, and then returns (x-1t(m)) %*% S^(-1/2).  
# In R, this method uses the cov.rob method in the lqs package
# In Splus, this method uses the covRob method in the robust package
# Additional args to the function are passed to cov.rob.
#
# make sure the library is loaded
 if(whichengine == "R")
   {if(!require(lqs)){ 
      stop("lqs package not in library, but is required for this function")}} else
   {if(!exists("covRob")) library(robust)}
 ans <- if(whichengine == "R") {cov.rob(x,method=method, ...) } else
                               {covRob(x,...)}
 m <- if (method == "classical") apply(x,2,"median") else ans$center
 s<-svd(ans$cov)
 sweep(x,2,m) %*% s$u %*% diag(1/sqrt(s$d))}

# END OF FILE
