drop1.dr <-
function (object, scope=NULL, update=TRUE, test = "general",...) 
{
    keep <- if(is.null(scope)) NULL else
         attr(terms(update(object$terms,scope)),"term.labels")
    all <- attr(object$terms,"term.labels")
    candidates <- setdiff(all,keep)
    if(length(candidates) == 0) stop("Error---nothing to drop")
    ans <- NULL
    for (label in candidates){
     ans <- rbind(ans,
      dr.coordinate.test(object,as.formula(paste("~.-",label,sep="")),...))
    }
    row.names(ans) <- paste("-",candidates)
    ncols <- ncol(ans)
    or <- order(-ans[,if(test=="general") ncols else (ncols-1)])
    form <- formula(object) 
    attributes(form) <- NULL
    fout <- deparse(form,width.cutoff=50)
    for (f in fout) cat("\n",f)
    cat("\n")
    print(ans[or,]) 
    if (is.null(object$stop)) { object$stop <- 0 }
    stopp <- if(ans[or[1],if(test=="general") ncols else (ncols-1)] <
                  object$stop) TRUE else FALSE 
    if (stopp == TRUE){ 
     cat("\nStopping Criterion Met\n")
    object$stop <- TRUE; object } else
     if(update==TRUE) {
      update(object,as.formula(paste("~.",row.names(ans)[or][1],sep="")))} 
       else invisible(ans)
    }
  
dr.step <-
function(object,scope=NULL,d=NULL,numdir=object$numdir,stop=0,...) {
   if(is.null(object$stop)){ object$stop <- stop }
   if(object$stop == TRUE) object else {
    numdir <- object$numdir
    keep <- if(is.null(scope)) NULL else
         attr(terms(update(object$terms,scope)),"term.labels")
    all <- attr(object$terms,"term.labels")
    if (length(keep) >= length(all) | length(all) <= numdir) object else {
    if(dim(object$x)[2] <= numdir) object else {
      obj1 <- drop1(object,scope=scope,d=d,...)
      dr.step(obj1,scope=scope,d=d,stop=stop,...)}}}}
