
#############################################################################
# partial sir --- see Chiaromonte, Cook and Li (2002).  The implemetation
# here follows Shao, Y. (2007) Topics on dimension reduction, unpublished PhD
# Dissertation, School of Statistics, University of Minnesota.
#############################################################################
 
dr.fit.psir <-function(object,numdir=4,...){ 
    M <- dr.M(object,...)  
    D <- eigen(M$M)
    or <- rev(order(abs(D$values)))
    evalues <- D$values[or]
    raw.evectors <- D$vectors[,or] 
    "%^%"<-function(A,n) {
            if (dim(A)[2]==1) {A^n}
            else {eg<-eigen(A)
                    (eg$vectors) %*% diag(abs(eg$values)^n) %*% t(eg$vectors)}
                    } 
    evectors <- M$Sigma.pool %^% (-1/2) %*% raw.evectors
    evectors <- if (is.matrix(evectors)) evectors else matrix(evectors,ncol=1)
    evectors <- apply(evectors,2,function(x) x/sqrt(sum(x^2)))
    names <- colnames(dr.x(object))[1:object$qr$rank]
    dimnames(evectors)<-
         list(names, paste("Dir", 1:NCOL(evectors), sep=""))
    aa<-c( object, list(evectors=evectors,evalues=evalues, 
                numdir=min(numdir,dim(evectors)[2],dim(M$M)[1]),
                raw.evectors=raw.evectors), M)
    class(aa) <- class(object)
    return(aa)
}
  
dr.M.psir <- function(object,nslices=2,slice.function=dr.slices,pool=FALSE,...){
    y <- dr.y(object)
    z <- dr.x(object) # centered/rotated z computed in this function
    w <- dr.wts(object)
    n <- length(w)
    p <- dim(z)[2]
    if(is.null(object$group)) object$group <- rep(1,dim(z)[1])
    group.names <- unique(as.factor(object$group))
    nw <- table(object$group)
    if (any(nw < p) ) stop("At least one group has too few cases")
# compute Sigma.pool
    Sigma.pool <- matrix(0,p,p)
    Sigma <- array(0,c(p,p,length(group.names)))
    wt.means <- matrix(0,length(group.names),p)
    for (j in 1:length(group.names)){
      sel <- object$group == group.names[j]
      wt.means[j,] <- apply(z[sel,],2,function(x) sum(x*w[sel])/sum(w[sel]))
      z[sel,] <- t(apply(z[sel,],1,function(x) x-wt.means[j,]))# centers
      Sigma[,,j] <- t(z[sel,]) %*% apply(z[sel,],2,"*",w[sel])/ sum(w[sel])
      Sigma.pool <- Sigma.pool + sum(w[sel]) *Sigma[,,j]/ n
      }
    "%^%"<-function(A,n) {
            if (dim(A)[2]==1) {A^n}
            else {eg<-eigen(A)
                    (eg$vectors) %*% diag(abs(eg$values)^n) %*% t(eg$vectors)}
                    }
# Center/scale z.  Chiaromonte et al. used Sigma.pooled; Shao uses separate
# Sigma for each group.  Compute slice info and slice means.
    slices <- slice.function(rep(1,length(w)),1)
    slices$nslices <- NULL
    zmeans <- NULL
    numslices <- 0
    for (j in 1:length(group.names)){ 
      sel <- object$group==group.names[j]
      z[sel,] <- z[sel,] %*% if (pool==TRUE)
                 Sigma.pool %^% (-1/2) else Sigma[,,j] %^% (-1/2)
      s <- slice.function(y[sel],nslices)
      slices$slice.indicator[sel] <- s$slice.indicator+numslices
      slices$nslices <- c(slices$nslices,s$nslices)
      slices$slice.sizes <- if( numslices == 0) s$slice.sizes else
            c(slices$slice.sizes,s$slice.sizes)
      numslices <- sum(slices$nslices)
      for (h in 1:s$nslices){
       slice.sel <- s$slice.indicator==h
       zmeans <- rbind(zmeans,
                        sqrt(sum(w[slice.sel])/sum(w)) * 
                         apply(z[slice.sel,],2,
                          function(x) sum(x*w[slice.sel])/sum(w[slice.sel])))
       }}
      return(list(M = t(zmeans) %*% zmeans, Sigma.pool=Sigma.pool,
        slice.info=slices))
       }
 
"summary.psir" <- function (object, ...)
{   z <- object
    ans <- z[c("call")]
    nd <- min(z$numdir,length(which(abs(z$evalues)>1.e-15)))
    ans$evectors <- z$evectors[,1:nd]
    ans$method <- z$method
    ans$nslices <- sum(z$slice.info$nslices)
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
    class(ans) <- "summary.psir"
    ans
}


"print.summary.psir" <-
function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")#S: ' ' instead of '\n'
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    cat("Method:\n")#S: ' ' instead of '\n'
    cat("partial sir with",x$nslices, "slices, n =",x$n) 
    cat("\nSlice Sizes:\n")#S: ' ' instead of '\n'
    cat(x$sizes,"\n")
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


# The function dr.temp.psir() is the core function which does the actual fitting and testing.
# dr.fit.psir() and dr.coordinate.test.sir() are just wrappers.
# The function takes a data set as input.  The outputs are:
# 1. directions: the estimated partial central subspace.
# 2. test(): the function to test marginal dimensional hypothesis.
# 3. coordinate.test(): the function to test marginal coordinate hypothesis.
dr.temp.psir=function(dat) {
    Y=dat$Y
    X=dat$X
    n=dim(X)[1]
    p=dim(X)[2]

    group.names <- unique(as.factor(dat$G))
    nG=length(group.names)
    G=numeric(n)
    for (j in 1:nG) G[dat$G==group.names[j]]=j

    # If no slices info is given, then we assume each group has two slices.
    slice=rep(2,nG)
    if (!is.null(dat$slice)) slice=rep(dat$slice,nG)
    zmeans=array(0,c(p,sum(slice)))
    SigmaPool=matrix(0,p,p)

    "%^%"<-function(A,n) {
            if (dim(A)[2]==1) {
                    A^n
            }
            else {
                    eg<-eigen(A)
                    (eg$vectors) %*% diag(abs(eg$values)^n) %*% t(eg$vectors)
            }}

    # This is a simplied SIR function for each group k.
    # Function is private.
    # The returning function stat() computes the marginal coorindate test statistic for the group,
    # as well the weights of the chi-squared distribution of the test statistic.
    temp.sir=function(k) {
        h=slice[k]
        ind=(G==k)
        Z=X[ind,]
        y=Y[ind]
        slices=dr.slices(y,h)
        n=length(y)
        Sigma=cov(Z)
        z=(Z-matrix(1,n,1)%*%apply(Z,2,mean)) %*% (Sigma %^% (-1/2))  
        zmeans=array(0,c(p,h))
        for (j in 1:h) {
            ind=(slices$slice.indicator==j)
            if (sum(ind)<1) mu=rep(0,p)
            else mu=apply(z[ind,],2,mean)
            zmeans[,j]=sqrt(sum(ind)/n)*mu
        }

        # The function stat() computes the marginal coorindate test statistic for the group,
        stat=function(gamma) {
            r <- p - dim(gamma)[2]
            gamma <- Sigma %^% (1/2) %*% gamma
            H <- as.matrix(qr.Q(qr(gamma), complete = TRUE)[, (p - r + 1):p])
            st=sum((t(H)%*%zmeans)^2)*n
            wts <- rep(1,h-1)
            wts[1:min(p,h-1)]=1-svd(zmeans)$d[1:min(p,h-1)]^2
            return(list(st=st,wts=wts))
        }

        return(list(zmeans=zmeans,stat=stat))
    }

    # The following code computes the zmeans matrix.
    for (k in 1:nG) {
        ind=(G==k)
        Z=X[ind,]
        if (sum(ind)<1) 
            zmeans[,(sum(slice[1:k])-slice[k]+1):sum(slice[1:k])]=rep(0,p)
        else {
            Sigma=cov(Z)
            SigmaPool=SigmaPool+Sigma
            temp=temp.sir(k)$zmeans
            zmeans[,(sum(slice[1:k])-slice[k]+1):sum(slice[1:k])]=temp
        }
    }

    # The following code estimates the partial central subspace.
    D=svd(zmeans,p)
    directions=SigmaPool%^%(-1/2)%*%D$u

    # The function test() does the sequential marginal dimension testing from 0 to d.
    # Function is public.
    test=function(d) {
        h=sum(slice)
        d=min(d,h-nG)
        st=df=pv=0
        for (i in 0:(d-1)) {
            st[i+1] <- sum(D$d[(i+1):min(p,h-nG)]^2)*n
            df[i+1] <- (h-i-nG)*(p-i) 
            pv[i+1] <- 1 - pchisq(st[i+1], df[i+1])
        }
        return(data.frame(Test=st,df=df,Pvalue=pv))
    }

    # The function coordinate.test() tests the marginal coordinate hypothesis.
    # The test statistic is constructed by summing up the test statistics in each group.
    # The weights are contructed by combining all the weights in each group.
    # Function is public.
    coordinate.test=function(H) {
        r=p-dim(H)[2]
        st=0
        wts=0
        for (i in 1:nG) { 
            tp=temp.sir(i)$stat(H)
            st=st+tp$st
            wts=cbind(wts,tp$wts)
        }
        wts=rep(wts,r)
        testr=dr.pvalue(wts[wts>1e-5],st,a=chi2approx)
        df=testr$df.adj
        pv=testr$pval.adj
        return(data.frame(cbind(Test=st,df=df,Pvalue=pv)))
    }

    return(list(directions=data.frame(directions),test=test,coordinate.test=coordinate.test))
}

dr.coordinate.test.psir<-function(object,hypothesis,d=NULL,
    chi2approx=object$chi2approx,...){
    gamma <- if (class(hypothesis) == "formula")
        coord.hyp.basis(object, hypothesis)
        else as.matrix(hypothesis)
    p<-length(object$evalues)
    n<-object$cases
    z <- dr.z(object)
    ev<-sort(object$evalues,decreasing=TRUE)
    d <- if(is.null(d)) length(ev) else d
    slices<-object$slice.info
    h<-slices$nslices
    M<-object$M
    r<-p-dim(gamma)[2]
    H<- (dr.R(object)) %*% gamma  
    H <- qr.Q(qr(H),complete=TRUE) # a p times p matrix
    QH<- H[,1:(p-r)] # first p-r columns
    H <- H[,(p-r+1):p] # last r columns
    st<-n*sum(ev[1:d])-n*sum(sort(eigen(t(QH)%*%M%*%QH)$values,
               decreasing=TRUE)[1:min(d,p-r)])
    wts <- 1-ev[1:min(d,h-1)]
# Bentler-Xie modified statistic, assuming conditions C1, C2, C3
# each eigenvalue occurs r times.  
    testr <- dr.pvalue(rep(wts,r),st,chi2approx=chi2approx) 
# general test
    epsilon<-array(0,c(n,h))
    zmeans<-array(0,c(p,h)) 
    for (i in 1:h) {
        sel<-(slices$slice.indicator==i)
        f_k<-sum(sel)/n
        zmeans[,i]<-apply(z[sel,],2,mean)
        epsilon[,i]<-(sel-f_k-z%*%zmeans[,i]*f_k)/sqrt(f_k)
        }
    HZ<- z%*%H
    Phi<-svd(zmeans,nv=h)$v[,1:d]
    epsilonHZ<-array(0,c(n,r*d))
    for (j in 1:n) epsilonHZ[j,]<-t(Phi)%*%t(t(epsilon[j,]))%*%t(HZ[j,])
    wts <- eigen(((n-1)/n)*cov(epsilonHZ))$values
    testg<-dr.pvalue(wts[wts>0],st,chi2approx=chi2approx)
    testg<-
    z <- data.frame(cbind(st, testr$pval.adj, testg$pval.adj))
    dimnames(z) <- list("Test", c("Statistic", "p.val(Res)", 
        "p.val(Gen)"))
    z
}
