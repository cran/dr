#####################################################################
#
#     Inverse regression methods for R and Splus
#     Written in July, 2000 by Sandy Weisberg
#     sandy@stat.umn.edu
#     ir.permutation test revised and corrected Feb. 15, 2001
#
#####################################################################
#
# Since R and Splus are not identical, different code is sometimes
# needed for the various engines (I have not tested Splus 3.0)
# whichengine provides a switch for the various engines.
# The code reads the variable thisengine, to get whichengine.
# Possible values are "R" and "s5".  If thisengine is not set, ask
# the user to set it.
#
#####################################################################

whichengine <- "s5"
#if (!exists("thisengine"))
# {
#  print(cat("Which S engine are you using?\n"))
#  thisengine<-switch(menu(c("R","Splus 5 (Unix)", "Splus 2000"))+1, 
#                            "R", "R", "s5", "s2k") 
# } else thisengine

#####################################################################
#
#     Splus class compatibility, not for R
#
#####################################################################

if (whichengine == "s5") {
  setOldClass(c("sir","ir"))
  setOldClass(c("save","ir"))
  setOldClass(c("phd","ir"))
  setOldClass(c("phdy","ir"))
  setOldClass(c("phdres","ir"))
  setOldClass(c("msir","mir"))
  setOldClass(c("msave","mir"))
  setOldClass(c("mphd","mir"))
  setOldClass(c("mphdy","mir"))
  setOldClass(c("mphdres","mir"))}

#####################################################################
#
#     ir is the primary function
#
#####################################################################
ir <-
function (formula, data = list(), subset, weights, na.action=na.omit,
      method = "sir", contrasts = NULL, offset = NULL, ...)
{
#this first section is copied from the lm function
    mf <- match.call(expand.dots=F) #puts args in a standard form
    mf$... <- NULL   # drop ...
    mf$method <- NULL # don't pass method to the model.frame
    mf$contrasts <- NULL # don't pass the contrasts to the model.frame
    mf[[1]] <- as.name("model.frame")
    mf$drop.unused.levels <- TRUE  #remove lin. dep. columns
    mf <- eval(mf, sys.frame(sys.parent()))
    mt <- attr(mf,"terms")
    xvars <- as.character(attr(mt, "variables"))[-1]
    if((yvar <- attr(mt, "response")) > 0) xvars <- xvars[-yvar]
    xlev <- if(length(xvars) > 0) {
                       xlev <- lapply(mf[xvars], levels)
                       xlev[!sapply(xlev, is.null)] }
    y <- model.response(mf, "numeric")
    w <- model.weights(mf)
    offset <- model.offset(mf)
    if(!is.null(offset) && length(offset) != NROW(y))
     stop(paste("Length of the offset", length(offset),
         ", must equal", NROW(y), " the number of cases"))
    if (is.empty.model(mt)) stop(paste("No model specified!"))
#Set the class name
    classname<-paste(if(is.matrix(y)) "m",method,sep="")
    genclassname<-paste(if(is.matrix(y)) "m","ir",sep="")
    z <- list(formula=formula)
    switch(whichengine,
       R=,class(z) <- c(classname,genclassname),
       s2k=,class(z) <- c(classname,genclassname),
       s5=,oldClass(z) <- classname)
#set up the model matrix
    x <- model.matrix(z,mt, mf, contrasts) 
    if (NROW(y) != nrow(x))
     stop("The response and predictors have differing number of observations")
#call the fit method
    z <- ir.fit(object=z,x=x,y=y,w=w,offset=offset,...)
#assign the results to names
    z <- c(z, list(offset=offset, contrasts = attr(x, "contrasts"),
                   xlevels = xlev, call = match.call(),
                   terms = mt, method = method,
                   model = mf, cases = NROW(y)))
#since z has been reassigned, assign its class again
    switch(whichengine,
       R=,class(z) <- c(classname,genclassname),
       s2k=,class(z) <- c(classname,genclassname),
       s5=,oldClass(z) <- classname)
#return the object
    z }


#####################################################################
#
#     a model.matrix method for ir
#
#####################################################################
model.matrix.ir<-function(x,...) {
    mat <- model.matrix(...)
    int <- match("(Intercept)", dimnames(mat)[[2]], nomatch=0)
    if (int > 0) mat <- mat[, -int, drop=F]
    mat}

#####################################################################
#
#     Fitting function
#
#####################################################################
ir.fit <-function(object,x,y,w,offset,numdir=3,
                    tol=1.e-7,decomp="qr",...) {
    wts <- if(is.null(w)) rep(1,NROW(y)) else w
    off <- if(is.null(offset)) rep(0,NROW(y)) else offset
# fit via (weighted) least squares to get fitted values, determine
# the rank of the model matrix, and get the qr factorization.
    z <- scale(x,center=T,scale=F) # z is centered
    if (whichengine == "R") {
             lsfit <- lm.wfit(z,y,wts,off,tol=tol)
             rank <-lsfit$qr$rank 
             cols <-lsfit$qr$pivot[1:rank]  }
    else if (whichengine == "s5" || whichengine == "s2k") {
             swts <- sqrt(wts)
             QR <- qr(z * swts)
             rank <-QR$rank
             cols <-QR$pivot[1:rank]  }
    if (decomp == "svd"){ 
       z <- z[,cols] # z has redundant cols removed
       D <- svd(sweep(z,1,sqrt(wts),"*"),nu=0) # temporary variable
       InvSqrtSigma <- D$v %*% diag(1/D$d) * sqrt(sum(wts))
       z <- z %*% InvSqrtSigma} # z is rotated.
    else { #  Using QR
       SqrtSigma <- qr.R(switch(whichengine,R=,lsfit$qr,s5=,QR,s2k=,QR))[cols,cols]*sqrt(sum(wts))
       z <- qr.Q(switch(whichengine,R=,lsfit$qr,s5=,QR,s2k=,QR))[,cols]*sqrt(sum(wts))
       }
# compute M, with a different function for each method
    fitval <- switch(whichengine,s5=,qr.fitted(QR,swts*(y-off)),
                                ,s2k=,qr.fitted(QR,swts*(y-off)),
                                 R=,lsfit$fitted.values) + mean(y)
    yvar <- switch(whichengine,
                   R=,ir.fit.y(object=object,y=y,offset=off,ols.fit=fitval),
                   s5=,ir.fit.y(object,y,off,fitval),
                   s2k=,ir.fit.y(object,y,off,fitval))
    M <- ir.fit.M(object=object,z=z,y=yvar,w=wts,...)
# compute singular values and vectors of M
    D <- svd(M$M,nrow(M$M),0)
    evalues <- if (nrow(M$M)==ncol(M$M)) D$d else (D$d)^2
    or <- rev(order(abs(evalues)))  # order absolute eigenvalues
    evalues<-evalues[or]
    raw.evectors <- D$u[,or]
    if (decomp == "svd") {
         evectors<-apply( (InvSqrtSigma %*% D$u)[,or],2,
                          function(x){x/(sqrt(sum(x^2)))})}
       else {
         evectors<-apply( solve(SqrtSigma,D$u)[,or],2,
                          function(x){x/(sqrt(sum(x^2)))})}
# assign names to eigen arrays
    dimnames(evectors)<-
         list(colnames(x)[cols], paste("Dir", 1:NCOL(evectors), sep=""))
    aa<-c( list(evectors=evectors,evalues=evalues, numdir=numdir,
                ols.coef=switch(whichengine,R=,lsfit$coef[-1],
                                 s5=,qr.coef(QR,swts*(y-off)),
                                 s2k=,qr.coef(QR,swts*(y-off))),
                offset=offset,weights=w,raw.evectors=raw.evectors,
                cols.used=cols,decomp=decomp, ols.fit=fitval), M)
    return(aa)
}

#####################################################################
#
#     ir.fit.M functions, one for each fitting method
#     ir.fit.y functions, one for each fitting method
#
#####################################################################

ir.fit.M <- switch(whichengine,
            s5=,function(object, ...){UseMethod("ir.fit.M")},
            s2k=,function(object, ...){UseMethod("ir.fit.M")},
            R= ,function(object, ...){UseMethod("ir.fit.M",object)})
ir.fit.y <- switch(whichengine,
            s5=,function(object, ...){UseMethod("ir.fit.y")},
            s2k=,function(object, ...){UseMethod("ir.fit.y")},
            R= ,function(object, ...){UseMethod("ir.fit.y",object)})

# default is used for sir and save
ir.fit.y.default<-function(object,y,offset,ols.fit) { 
  y - if (is.null(offset)) 0 else offset}

#####################################################################
#     Sliced Inverse Regression
#####################################################################

ir.fit.M.sir <-function(object,z,y,w=NULL,nslices=NULL,
                        slice.info=NULL,...) {
# get slice information
    h <- if (!is.null(nslices)) nslices else max(8, NCOL(z)+3)
    slices<- if(is.null(slice.info)) ir.slices(y,h) else slice.info
# get slice means
    xmeans <- matrix(0,slices$nslices,ncol(z))
    for (j in 1:slices$nslices){
      xmeans[j,]<- apply(z[slices$slice.indicator==j,],2,mean)}
# get M matrix for sir
    M <- t(xmeans) %*% diag(slices$slice.sizes) %*% xmeans / NROW(y)
    return (list (M=M,slice.info=slices))
}

#####################################################################
#     Sliced Average Variance Estimation
#####################################################################

ir.fit.M.save <-function(object,z,y,w=NULL,nslices=NULL,
                        slice.info=NULL,...) {
# get slice information
    h <- if (!is.null(nslices)) nslices else max(8, ncol(z)+3)
    slices<- if(is.null(slice.info)) ir.slices(y,h) else slice.info
# get M = sum (ns-1)(I-C)^2/sum(ns-1)
    M <- matrix(0,NCOL(z),NCOL(z))
    for (j in 1:slices$nslices){
      ns <- slices$slice.sizes[j]
      IminusC <- diag(rep(1,NCOL(z))) - var(z[slices$slice.indicator==j,])
      M <- M + (ns-1)*IminusC %*% IminusC
      }
    M <- M/(NROW(z)-slices$nslices)
    return (list (M=M,slice.info=slices))
}

#####################################################################
#     pHd, pHdy and pHdres
#####################################################################

ir.fit.y.phdy <- function(object,y,offset,ols.fit){y - mean(y)}
ir.fit.M.phdy <- function(...) {ir.fit.M.phd(...)}

ir.fit.y.phdres <- function(...) {ir.fit.y.phd(...)}
ir.fit.M.phdres <- function(...) {ir.fit.M.phd(...)}

ir.fit.y.phd <- function(object,y=ir.y(object),offset=ojbect$offset,
         ols.fit=object$ols.fit) {y-ols.fit}
ir.fit.M.phd <-function(object,z,y,w=NULL,...) {
    wts <- if(is.null(w)) rep(1,NROW(z)) else w
# compute M
    M<- (t(z*(wts*y)) %*% z) / sum(wts)
    return(list(M=M))
}

#####################################################################
#
#     Accessor functions for inverse regressions
#
#####################################################################

# recover and x and y data

ir.x <- function(z) {model.matrix(z,z$terms,z$model,z$contrasts)[,z$cols.used]}
ir.y <- function(object){ model.response(object$model, "numeric")}

ir.y.name <- function(object){
  which <- attr(attr(object$model, "terms"),"response")
  sel <- switch(whichengine,R=,which+1,s5=,which,s2k=,which)
  as.character(attr(attr(object$model, "terms"), "variables")[sel])
  }

ir.x.omitted <- function(z) 
   {dimnames(model.matrix(z,z$terms,z$model,z$contrasts))[[2]][-z$cols.used]}

# recover the direction vectorss

ir.directions <- function(object, ...) {UseMethod("ir.direction")}
ir.direction <- function(object, ...) {UseMethod("ir.direction")}

ir.direction.default <- 
  function(object, which=1:object$numdir,norm=F,x=ir.x(object)) {
    ans <- (scale(x,center=T,scale=F) %*% object$evectors)[,which]
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


plot.ir <- function(object,which=1:object$numdir,mark.by.y=F,...) {
 d <- ir.direction(object,which)
 if (mark.by.y == F) {
    pairs(cbind(ir.y(object),d),labels=c(ir.y.name(object),colnames(d)),...)
    }
 else {pairs(d,labels=colnames(d),col=markby(ir.y(object)),...)}
 }

#############This does not work in Splus
# point marking in S is not easy because the arg col is expected to be an
# integer rather than a list
#########################################################################
#plot.ir <- function(object,which=1:object$numdir,mark.by.y=F,
#                panel=points.col,colors,...) {
# d <- ir.direction(object,which)
# if (mark.by.y == F) {
#    pairs(cbind(ir.y(object),d),labels=c(ir.y.name(object),colnames(d)),
#                 ,panel,colors=colors,...)
#    }
# else {pairs(d,labels=colnames(d),col=markby(ir.y(object)),...)}
# }

points.col <- function(x,y,colors,...){
  if (whichengine == "R") {points(x,y,col=colors,...)}
  else if (length(colors) == 1) points(x,y,col=colors,...)
   else {
    unique.colors <- unique(colors)
    for (j in 1:length(unique.colors)){
     sel<- colors == unique.colors[j]
     points(x[sel],y[sel],col=unique.colors[j],...)}}}

ir.coplot <- function(object,which=1:object$numdir,mark.by.y=F,...) {
 d <- data.frame(ir.direction(object,which))
 if (mark.by.y == F){
      d$yvar <- ir.y(object)
      coplot(yvar~Dir1|Dir2,data=d,...)}
    else
     {coplot(Dir1~Dir2|Dir3,data=d,col=markby(ir.y(object)),...)}
}

givens.rotation <- function(theta,p=2,which=c(1,2)){
 m <- matrix(rep(0,p^2),c(p,p))
 m[which,which]<-
        matrix(c(cos(theta),-sin(theta),sin(theta),cos(theta)),c(2,2))
 m}

#
# The following function rotplot gives static views of a 3D rotating
# plot.  A call might be 
#   rotplot(ir.directions(m1,1:2),ir.y(m1),number=16)
# Doesn't seem very useful...
rotplot <- function(x,y,number=9,theta=seq(0,pi/2,length=number),...){
 z<-NULL
 for (j in 1:number) { 
   z<- rbind(z, cbind(y, x %*% givens.rotation(theta[j]), theta[j]))}
   dimnames(z)[[2]]<-
           c( deparse(substitute(y)), "Linear Combination", "LC2", "Angle")
   coplot(z[,1]~z[,2]|z[,4],number=number,overlap=0,...)}
  
"panel.smooth2" <-
function (x, y, col = par("col"), pch = par("pch"),
    span = 2/3, iter = 0, ...) 
{
    points(x, y, pch = pch, col = col)
    lin.col <- unique(col)
    for (i in 1:length(unique(col))){
        ok<- if(length(col)==1)
              is.finite(x) & is.finite(y) else
              is.finite(x) & is.finite(y) & col==lin.col[i]
        if (any(ok)) {
            or <- order(x[ok])
            a<-data.frame(cbind(x[ok][or],y[ok][or]))
            l<-loess(a[,2]~a[,1],...)
            l2<-loess(l$residuals^2~a[,1],data=a,...)
            lines(l$x,l$fitted,col=lin.col[i],...)
            lines(l$x,l$fitted+sqrt(max(0,l2$fitted)),col=lin.col[i],...)
            lines(l$x,l$fitted-sqrt(max(0,l2$fitted)),col=lin.col[i],...)
    }}
}


ir.persp<-function(object,which=1:2,h=c(.1,.1),...){
 if (length(which) == 2){
  d1<-ir.direction(object,which,norm=T)
  y<-ir.y(object)
  sm.regression(d1,y,h=h,
                     xlab=dimnames(d1)[[2]][1],
                     ylab=dimnames(d1)[[2]][2],
                     zlab=names(object$model[1]))
  }
  else
  print("This method requires specifying two directions")
 }

markby <- function(z,use="color",values=NULL,color.fn=rainbow) {
 u <- unique(z)
 lu <- length(u)
 ans <- 0
 vals <- if (use == "color") 
      {if (!is.null(values) && length(values) == lu)
                values else color.fn(lu)}
   else
      {if (!is.null(values) && length(values) == lu)
                   values else 1:lu}
 for (j in 1:lu){ans[z==u[j]] <- vals[j]}
 ans}
   

###################################################################
#
#  basic print method for inverse regression
#
###################################################################
"print.ir" <-
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
#  basic summary method for inverse regression
#
###################################################################
"summary.ir" <-
function (object)
{
    z <- switch(whichengine,R=,.Alias(object),s5=,object,s2k=,object)
    ans <- z[c("call", "terms")]
    nd <- min(z$numdir,length(which(abs(z$evalues)>1.e-8)))
    ans$evectors <- z$evectors[,1:nd]
    ans$method <- z$method
    ans$nslices <- z$slice.info$nslices
    ans$sizes <- z$slice.info$slice.sizes
    ans$n <- NROW(z$model)
    ans$omitted <- ir.x.omitted(z)
    ans$evalues <-rbind (z$evalues[1:nd],
                         cosangle(ir.direction(object),object$ols.fit)[1:nd])
    dimnames(ans$evalues)<-
     list(c("Eigenvalues","R^2(OLS|ir)"),
          paste("Dir", 1:NCOL(ans$evalues), sep=""))
    ans$test <- ir.test(object,nd)
    class(ans) <- "summary.ir"
    ans
}

###################################################################
#
# basic print.summary method for inverse regression
#
###################################################################
"print.summary.ir" <-
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
    if(is.null(x$nslices))
       cat(paste(x$method, ", n = ", x$n, ".\n",sep=""))
       else {
         cat(paste(x$method," with ",x$nslices, " slices, n = ",
                   x$n,".\n",sep=""))
         cat("\nSlice Sizes:\n")#S: ' ' instead of '\n'
         cat(x$sizes,"\n")}
    cat("\nEigenvectors:\n")
    print(round(x$evectors,digits))
    cat("\n")
    print(round(x$evalues,digits))
    if (length(x$omitted) > 0){
      cat("\nModel matrix is not full rank.  Deleted columns:\n")
      cat(x$omitted,"\n")}
    if (!is.null(x$test)){
      cat("\nAsymp. Chi-square tests for dimension:\n")
      print(round(as.matrix(x$test),digits))}
    invisible(x)
}

###################################################################
#
# Asymptotic test methods for inverse regression.  Separate 
# functions for each computing method.
#
###################################################################
switch(whichengine,
  R=, ir.test <- function(object, ...){ UseMethod("ir.test",object)},
  s5=,ir.test <- function(object, ...){ UseMethod("ir.test")},
  s2k=,ir.test <- function(object, ...){ UseMethod("ir.test")}
)

ir.test.default <-function(object, ...) {NULL}

ir.test.sir<-function(object,nd) {
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
# Written by Jorge de la Vega, February, 2001
#
#########################################################################

ir.test.phd<-function(object,nd) {
#compute the phd asymptotic test statitics under restrictions for the
#first nd directions, based on response = OLS residuals
# order the absolute eigenvalues
    e<-sort(abs(object$evalues))
    p<-length(object$evalues)
# get the response
    resi<-ir.fit.y(object,ir.y(object),object$offset,object$ols.fit)
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
    indep <- ir.indep.test.phdres(object,st[1])
# compute tests that do not require normal theory (linear combination of 
# chi-squares:
    lc <- ir.test2.phdres(object,st)
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


###################################################################3
##
##  Translation of methods in Arc for testing with pHd to R
##  Original lisp functions were mostly written by R. D. Cook
##  Translation to R by S. Weisberg, February, 2001
##
###################################################################3
# this function is a translation from Arc.  It computes the matrices W and
# eW described in Sec. 12.3.1 of Cook (1998), Regression Graphics.
cov.ew.matrix <- function(object,scaled=F){
  n <- dim(ir.x(object))[1]
  wts <- if (is.null(object$weights)) rep(1,n) else object$weights
# The products Z %*% theta where theta are the eigenvectors in the Z scale
# are the same as the products X %*% evectors in the original scale, after
# centering X and then scaling the result so columns have length sqrt(wts).
# This avoids computing two singular value decompositions
# v <- ir.normalize.z(ir.x(object),wts) %*%
#           object$raw.evectors # eigenvectors in z scale 
  v <- sqrt(sum(wts))* mat.normalize(
         sweep(scale(ir.x(object),center=T,scale=F),1,sqrt(wts),"*") %*%
	 object$evectors)
  y <- ir.fit.y(object) # get the response
  y <- if (scaled) y-mean(y) else 1 # a multiplier in the matrix
  p <- dim(v)[2]
  ew0 <- NULL
  for (i in 1:p){
   for (j in i:p){
    ew0 <- cbind(ew0, if (i ==j) y*(v[,j]^2-1) else y*sqrt(2)*v[,i]*v[,j])}}
  var(apply(ew0,2,"*",sqrt(wts)))
  }

#translation of :general-pvalues method for phd in Arc
ir.test2.phdres <- function(object,stats){
  covew <- cov.ew.matrix(object,scaled=T)
  C <- .5/var(ir.fit.y(object))
  p <- length(stats)
  pval <- NULL
  d2 <-dim(ir.x(object))[2]
  start <- -d2
  end <- dim(covew)[2]
  for (i in 1:p) {
   start <- start + d2-i+2
   evals <- eigen(C*covew[start:end,start:end],only.values=TRUE)$values
   pval<-c(pval,wood.pvalue(evals,stats[i]))}
# report results
    z<-data.frame(cbind(stats,pval))
    rr<-paste(0:(p-1),"D vs >= ",1:p,"D",sep="")
    dimnames(z)<-list(rr,c("Stat","p-value"))
    z}
   
ir.indep.test.phdres <- function(object,stat) {
  eval <- eigen(cov.ew.matrix(object,scaled=F),only.values=T)
  pval<-wood.pvalue(.5*eval$values,stat)
# report results
    z<-data.frame(cbind(stat,pval))
    dimnames(z)<-list(c("Test of independence"),c("Stat","p-value"))
    z}


wood.pvalue <- function (coef, f, tol=0.0, print=F){
#Returns an approximation to P(coef'X > f) for X=(X1,...,Xk)', a vector of iid
#one df chi-squared rv's.  coef is a list of positive coefficients. tol is used
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
# permutation tests for inverse regression
#
#########################################################################

ir.permutation.test <- function(object,npermute=50,numdir=object$numdir,
                                permute.weights=TRUE) {
# nd is the number of dimensions to test for
   nd <- min(numdir,length(which(abs(object$evalues)>1.e-8))-1)
   nt <- nd + 1
# y is generally the response variable, but it is residuals for phdres
   y <- ir.fit.y(object,y=ir.y(object),
                        if(is.null(object$offset)) 0 else object$offset,
                        object$ols.fit)
   n <- NROW(y)
# observed value of the test statistics = obstest
   obstest<-ir.permutation.test.statistic(object,object$evalues,nt,
                          object$cases,y)
# make sure weights are a list and not just NULL
   wts <- if(is.null(object$weights)){rep(1,n)} else {object$weights}
# z2 is a matrix like X consisting of the principal directions.  z2 and
# x span the same subspace and are therefore equivalent.
   z2<- ir.direction(object,which=1:(NCOL(ir.x(object))))
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
	z<-ir.normalize.z(z,wp)      # uses singular value decomposition 
# the eigenvalues/singular values of M determine the test statistic
        M <- ir.fit.M(object,z,y,wp,method=object$method,
                        slice.info=object$slice.info)$M 
# get the eigenvalues
        evalues <- if(NROW(M)==NCOL(M)) rev(sort(abs(svd(M,0,0)$d)))
                      else rev(sort(abs(svd(M,0,0)$d^2))) #get eigenvalues
        val[col+1]<-
            ir.permutation.test.statistic(object,evalues,col+1,n,y)[col+1]
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
   class(ans) <- "ir.permutation.test"
   ans
   }

"print.ir.permutation.test" <-
function(ans, digits = max(3, getOption("digits") - 3), ...)
{
   cat("\nPermutation tests\nNumber of permutations:\n")
   print.default(ans$npermute)
   cat("\nTest results:\n")
   print(ans$summary,digits=digits) 
   invisible(ans)
}

"summary.ir.permutation.test" <- function(...)
              {print.ir.permutation.test(...)}


#########################################################################
#
# ir.permutation.test.statistic method
#
#########################################################################

ir.permutation.test.statistic <- 
  switch(whichengine,
  R=, function(object,...)
                     {UseMethod("ir.permutation.test.statistic",object)},
  s5=,function(object,...)
                     {UseMethod("ir.permutation.test.statistic")},
  s2k=,function(object,...)
                     {UseMethod("ir.permutation.test.statistic")})

ir.permutation.test.statistic.default <- function(object,values,nd,n,...){
   n*rev(cumsum(rev(values)))[1:nd]}

ir.permutation.test.statistic.phdy <- function(object,...){
       ir.permutation.test.statistic.phd(object,...)}
ir.permutation.test.statistic.phdres <- function(object,...){
       ir.permutation.test.statistic.phd(object,...)}
ir.permutation.test.statistic.phd <- function(object,values,nd,n,y){
   (.5*n*rev(cumsum(rev(values^2)))/var(y))[1:nd]}

#####################################################################
#
#     ir.slices returns non-overlapping slices based on y
#
#####################################################################

ir.slices <- function(y,h) {
  z<-unique(y)
  if (length(z) >= h) ir.slice2(y,h) else ir.slice1(y,length(z),sort(z))}

ir.slice1 <- function(y,h,u){
  z <- sizes <- 0
  for (j in 1:length(u)) {
      temp <- which(y==u[j])
      z[temp] <- j
      sizes[j] <- length(temp) } 
  list(slice.indicator=z, nslices=length(u), slice.sizes=sizes)
  }

ir.slice2<-function(y,h)
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
  sp[j]<-n
  ans[or[1:sp[1]]] <- 1
  for (k in 2:j){ans[ or[(sp[k-1]+1):sp[k] ] ] <- k}
  list(slice.indicator=ans, nslices=j, slice.sizes=c(sp[1],diff(sp)))
}

ir.normalize.z <- function(x,wts,center=TRUE){
      z<- if (center) scale(x,center=T,scale=FALSE) else x
      SVD <- svd(sweep(z,1,sqrt(wts),"*"),nu=0) # temporary variable
      InvSqrtSigma <- SVD$v %*% diag(1/SVD$d) * sqrt(sum(wts))
      z %*% InvSqrtSigma  # centered AND rotated.
      }



#####################################################################
#
#     Auxillary functions
#
#####################################################################

#
# angle between a vector vec and a subspace span(mat)
#

cosangle <- function(mat,vec) {
# For reasons I don't understand, qr.qty does not seem to work.
# I have used a work-around that actually requires forming Q explictly
# cumsum(qr.qty(qr(mat),scale(vec)/sqrt(length(vec)-1))^2) 
  cumsum((t(qr.Q(qr(mat))) %*% scale(vec)/sqrt(length(vec)-1))^2) 
}

#
# random permutation
#

permute <- function(n){order(runif(n))}

#####################################################################
##
##  Add functions to Splus that are built-in to R
##
#####################################################################

if (whichengine != "R") {

"model.response" <- function (data, type = "any")
{
 model.extract(data,"response")
}
"model.weights" <- function(x) x$"(weights)"
"model.offset" <- function(x) x$"(offset)"
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
