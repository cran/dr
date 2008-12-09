
# The function dr.temp.psave() is the core function which does the actual fitting and testing.
# dr.fit.psave() and dr.coordinate.test.save() are just wrappers.
# The function takes a data set as input.  The outputs are:
# 1. directions: the estimated partial central subspace.
# 2. test(): the function to test marginal dimensional hypothesis.
# 3. coordinate.test(): the function to test marginal coordinate hypothesis.
dr.temp.psave=function(dat) {
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

    M=matrix(0,p,p)
    A=array(0,c(sum(slice),p,p))
    SigmaPool=matrix(0,p,p)

    "%^%"<-function(A,n) {
            if (dim(A)[2]==1) {
                    A^n
            }
            else {
                    eg<-eigen(A)
                    (eg$vectors) %*% diag(abs(eg$values)^n) %*% t(eg$vectors)
            }
    }

    # This is a simplied SAVE function for each group.
    # But here the predictors are already normalized.
    # This function is private.
    temp.save=function(z,y,h) {
        slices=dr.slices(y,h)
        n=length(y)
        M=matrix(0,p,p)
        A=array(0,c(h,p,p))
        for (j in 1:h) {
            ind=(slices$slice.indicator==j)
            if (sum(ind)<3) IminusC=diag(rep(1,p))
            else IminusC=diag(rep(1,p))-cov(z[ind,]) 
            M=M+sum(ind)/n* IminusC%*%IminusC
            A[j,,]=sqrt(sum(ind)/n)*IminusC
        }
        return(list(M=M,A=A))
    }

    # The following code computes the M matrix and the A matrices.
    for (k in 1:nG) {
        ind=(G==k)
        Z=X[ind,]
        if (sum(ind)<3) 
            A[(sum(slice[1:k])-slice[k]+1):sum(slice[1:k]),,]=diag(rep(0,p))
        else {
            Sigma=cov(Z)
            SigmaPool=SigmaPool+Sigma*sum(ind)/n
            Z=(X[ind,]-matrix(1,sum(ind),1)%*%apply(X[ind,],2,mean)) %*% (Sigma %^% (-1/2))
            temp=temp.save(Z,Y[ind],slice[k])
            M=M+temp$M*sum(ind)/n
            A[(sum(slice[1:k])-slice[k]+1):sum(slice[1:k]),,]=temp$A*sqrt(sum(ind)/n)
        }
    }

    # The following code computes the eigenvectors of M, 
    # and transform the eigenvectors into X scale.
    D <- eigen(M)
    or <- rev(order(abs(D$values)))
    evectors <- D$vectors[, or]
    directions=SigmaPool%^%(-1/2)%*%evectors

    # The function temp.test() tests the marginal coordinate hypothesis.
    # The argument H is a matrix transormed from the actual hypothesis.
    # Function is private.
    temp.test=function(H) {
        r <- dim(H)[2]
        st <- 0
        h <- sum(slice)
        for (j in 1:h) {
            st <- st + sum((t(H) %*% A[j, , ] %*% H)^2) * n/2
        }
        df <- (h - nG) * r * (r + 1)/2
        pv <- 1 - pchisq(st, df)
        return(data.frame(Test=st,df=df,Pvalue=pv))
    }

    # The function coordinate.test() tests the marginal coordinate hypothesis.
    # The argument gamma is the actual hypothesis.
    # Function is public.
    coordinate.test=function(gamma) {
        r <- p - dim(gamma)[2]
        gamma <- SigmaPool %^% (1/2) %*% gamma
        H <- as.matrix(qr.Q(qr(gamma), complete = TRUE)[, (p - r + 1):p])
        return(temp.test(H))
    }

    # The function test1() tests the marginal dimension hypothesis that dimension=d.
    # Function is private.
    test1=function(d) {
        return(temp.test(as.matrix(evectors[,(d+1):p])))
    }

    # The function test() does the sequential marginal dimension testing from 0 to d.
    # Function is public.
    test=function(d) {
        tp=data.frame()
        for (i in 1:d) 
            tp=rbind(tp,test1(i-1))
        rr<-paste(0:(d-1),"D vs >= ",1:d,"D",sep="")
        dimnames(tp)<-list(rr,c("Stat","df","p-values"))
        return(tp)
    }

    return(list(directions=data.frame(directions),test=test,coordinate.test=coordinate.test))
}

dr.fit.psave <-function(object,numdir=4,nslices=NULL,...){ 
    dat=list()
    dat$X=dr.x(object)
    dat$Y=dr.y(object)
    dat$G=object$group
    if (is.null(nslices)) dat$slice=nslices
    else dat$slice=2
    object$temp=dr.temp.psave(dat)
    object$result=object$temp$directions[,1:numdir]
    object$test=object$temp$test(numdir)
    class(object)=c("dr","psave")
    return(object)
}

dr.coordinate.test.psave=function(object,hypothesis,d=NULL,...) {
    gamma <- if (class(hypothesis) == "formula")
        coord.hyp.basis(object, hypothesis)
        else as.matrix(hypothesis)
    object$temp$coordinate.test(gamma)
}

"print.psave"=function(object,...) {
    print(object$call)
    print(object$result)
    cat("Asymp. Chi-square tests for dimension:\n")
    print(object$test)
}

"summary.psave"=function(object,...) {
    ans=object["call"]
    ans$result=object$result
    ans$test=object$test
    class(ans)="summary.psave"
    ans
}
