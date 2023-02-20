

read.tsv=function(t,verbose=FALSE,stringsAsFactors=FALSE,...) {

  readit=function(file,...)
    if (requireNamespace("data.table",quietly = TRUE) && summary(file(file))$class!="gzfile") {
      as.data.frame(data.table::fread(file = file,stringsAsFactors=stringsAsFactors,check.names=FALSE,data.table = FALSE,...))
    } else {
      read.delim(file = file,stringsAsFactors=stringsAsFactors,check.names=FALSE,...)
    }

  if (suppressWarnings(requireNamespace("RCurl",quietly = TRUE)) && RCurl::url.exists(t)) {
    file <- tempfile()
    if (verbose) cat("Downloading file...\n")
    download.file(url, file, quiet=!verbose)
    if (verbose) cat("Reading file...\n")
    t=readit(file,...)
    if (verbose) cat("Deleting temporary file...\n")
    unlink(file)
  } else {
    if (verbose) cat("Reading file...\n")
    t=readit(t,...)
  }

  if (stringsAsFactors==TRUE) t=as.data.frame(lapply(t,function(c) if (is.character(c)) factor(c,levels=unique(c)) else c),check.names=FALSE)
  t
}

confint.nls.lm=function (object, parm, level = 0.95, ...)
{
  cc <- coef(object)
  if (missing(parm))
    parm <- seq_along(cc)
  levs <- c((1 - level)/2, 0.5 + level/2)
  dfval <- (length(object$fvec) - length(object$par))
  tdist <- qt(levs[2], dfval)
  m1 <- outer(sqrt(diag(vcov(object))), c(-1, 1)) * tdist
  m2 <- sweep(m1, cc, MARGIN = 1, "+")
  colnames(m2) <- paste(100 * levs, "%", sep = "")
  m2[parm, ]
}

equal = function(a,b) length(a)==length(b) && all(a==b)

#' Defer calling a function
#'
#' This generates a function with one mandatory parameter (and additional optional parameters)
#' that, when called, (i) also receives the parameters given when calling \code{Defer}, and (ii)
#' after calling it each element of the \code{add} list is appended by \code{+}. When no optional parameters
#' are given, the result is cached.
#'
#' @param FUN the function to be deferred
#' @param ... additional parameters to be used when the deferred function is called
#' @param add list containing additional elements to be added \code{+} to the result of the deferred function
#' @param cache use caching mechanism
#'
#' @return a function that can be called
#' @export
#'
#' @details The following expressions are very similar: \code{f <- function(d) Heavy.function(d)} and \code{f <- Defer(Heavy.function)}. In
#' both cases, you get a function \code{f} that you can call for some \code{d}, which in turn calls \code{Heavy.function}. The only
#' difference is that in the second case, the result is cached: \code{Heavy.function} is called only once when first calling \code{f},
#' if \code{f} is called a second time, the previous result is returned. This makes sense if the parameter \code{d} is constant (like a grandR object)
#' and if \code{Heavy.function} is deterministic.
#'
#' @details If additional parameters are provided to \code{f}, caching is disabled. If any of these additional parameters has the same name as the parameters
#' given to \code{Defer()}, the parameters given to \code{Defer()} are overwritten. Be careful if \code{Heavy.function} is not deterministic (see examples).
#'
#' @details Use case scenario: You want to produce a heatmap from a grandR object to be used as \code{plot.static} in the shiny web interface.
#' \code{\link{PlotHeatmap}} takes some time, and the resulting object is pretty large in memory. Saving the heatmap object to disk is very
#' inefficient (the Rdata file will be huge, especially with many heatmaps). Deferring the call without caching also is bad, because whenever
#' the user clicks onto the heatmap, it is regenerated.
#'
#' @examples
#' Heavy.function <- function(data) rnorm(5,mean=data)
#' f1=Defer(Heavy.function)
#' f2=function(d) Heavy.function(d)
#' f2(4)
#' f2(4) # these are not equal, as rnorm is called twice
#' f1(4)
#' f1(4) # these are equal, as the result of rnorm is cached
#'
#' @concept helper
Defer=function(FUN,...,add=NULL, cache=TRUE) {
  param=list(...)
  value=NULL
  function(data,...) {
    pp=list(...)
    if (length(pp)>0) {
      re=do.call(FUN,c(list(data),utils::modifyList(param,pp)))
      if (!is.null(add)) for (e in if (methods::is(add,"gg")) list(add) else add) re=re+e
      return(re)
    }

    if (is.null(value)) {
      value<<-do.call(FUN,c(list(data),param))
      if (!is.null(add)) for (e in if (methods::is(add,"gg")) list(add) else add) value<<-value+e
    }
    value
  }
}


#' Convert a structure into a vector
#'
#' The structure is supposed to be a list. Flattening is done by extracting the given fields (\code{return.fields})
#' and applying the additional function (\code{return.extra}). This is mainly to be used within \code{sapply} and similar.
#'
#' @param d the data structure
#' @param return.fields which fields should be extracted directly (may be NULL)
#' @param return.extra apply a function returning a flat list or vector (may be NULL)
#'
#'
#' @return the data flattened into a vector
#' @export
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' sars <- Normalize(sars)
#' fit <- FitKineticsGeneLeastSquares(sars,"SRSF6")$Mock
#' print(fit)
#' kinetics2vector(fit)
#'
#' @concept helper
structure2vector=function(d,return.fields=NULL,return.extra=NULL) {
  r=list()
  if (!is.null(return.fields)) r=c(r,d[return.fields])
  if (!is.null(return.extra)) r=c(r,return.extra(d))
  unlist(r)
}
#' @describeIn structure2vector Convert the output of the FitKinetics methods into a vector
#' @param condition if the original grandR object had \code{\link{Condition}} set, which condition to extract (NULL otherwise)
#' @export
kinetics2vector=function(d,condition=NULL,return.fields=c("Synthesis","Half-life"),return.extra=NULL) structure2vector(if (is.null(condition)) d else d[[condition]],return.fields=return.fields,return.extra=return.extra)

my.precision=function (x)
{
  x <- unique(x)
  x <- x[is.finite(x)]
  if (length(x) <= 1) {
    return(1)
  }
  smallest_diff <- min(diff(sort(x)))
  if (smallest_diff < sqrt(.Machine$double.eps)) {
    1
  }
  else {
    pmin(10^(floor(log10(smallest_diff))), 1)
  }
}

sd.from.hessian=function(H) {
  m=-H
  keep=rep(T,nrow(m))
  while(sum(keep>0) && rcond(m[keep,keep,drop=FALSE])<.Machine$double.eps) {
    s=apply(abs(m),1,sum)
    keep=keep&s>min(s[keep])
  }
  d=rep(0,nrow(m))
  d[keep]=sqrt(diag(solve(m[keep,keep])))
  d
}


cnt=function(m) {
  m=as.matrix(m)
  mode(m) <- "integer"
  m
}

#' Estimate dispersion parameters for a count matrix using DESeq2
#'
#' @param ss the count matrix
#'
#' @return a vector of dispersion parameters (to be used as size=1/dispersion for Xnbinom functions)
#' @export
#'
#' @concept helper
estimate.dispersion=function(ss) {
  checkPackages("DESeq2")

  dds=DESeq2::DESeqDataSetFromMatrix(countData = cnt(ss),colData=data.frame(rep(1,ncol(ss))),design = ~1)
  dds=DESeq2::estimateSizeFactors(dds)
  dds=DESeq2::estimateDispersions(dds,quiet=TRUE)
  disp=DESeq2::dispersions(dds)
  disp[is.na(disp)]=0.1 # ss=0 0 0 ...
  disp
}

GetField=function(name,field,sep=".") sapply(strsplit(as.character(name),sep,fixed=TRUE),function(v) v[field])


#' Parallel (s/l)apply
#'
#' Depending on whether \link{SetParallel} has been called, execute in parallel or not.
#'
#' @param ... forwarded to lapply or parallel::mclapply
#' @param seed Seed for the random number generator
#'
#' @details If the code uses random number specify the seed to make it deterministic
#'
#' @return a vector (psapply) or list (plapply)
#' @export
#'
#' @concept helper
psapply=function(...,seed=NULL) {simplify2array(plapply(...,seed=seed))}
#' @rdname psapply
#' @export
plapply=function(...,seed=NULL) {
  if (!IsParallel()) return(opt$lapply(...))

  if (!is.null(seed)) {
    rng=RNGkind()[1]
    RNGkind("L'Ecuyer-CMRG") # make it deterministic!
    set.seed(seed)
  }

  re=opt$lapply(...)

  if (!is.null(seed)) {
   RNGkind(rng)
  }

  re
}


opt <- new.env()
opt$lapply=function(...) lapply(...)
opt$sapply=function(...) simplify2array(opt$lapply(...))
opt$nworkers=0
opt$singlemsg=list()

singleMessage=function(msg) {
  if (is.null(opt$singlemsg[[msg]])) {
    cat(sprintf("%s\n (This message is only shown once per session)\n",msg))
    opt$singlemsg[[msg]] = 1
  }
}


#' Set up parallel execution
#'
#' Set the number of cores for parallel execution.
#'
#' @param cores number of cores
#'
#' @details Whenever \link{psapply} or \link{plapply} are used, they are executed in parallel.
#'
#' @return No return value, called for side effects
#' @export
#'
#' @concept helper
SetParallel=function(cores=max(1,parallel::detectCores()-2)) {
  opt$nworkers=cores
  if (cores>1) {
    if (.Platform$OS.type!="unix") {
      warning("Parallelism is not supported under windows. Will use a single thread!")
      opt$nworkers=0
    } else {
      opt$lapply<-function(...) parallel::mclapply(...,mc.cores=cores)
    }
  } else {
    opt$lapply<-function(...) lapply(...)
  }
}

#' Checks for parallel execution
#'
#' @return whether or not parallelism is activated
#' @export
#' @concept helper
IsParallel=function() {opt$nworkers>1}


checkPackages=function(pp,error=TRUE,warn=TRUE) {
  missing = c();
  for (p in pp) {
    if (!suppressWarnings(requireNamespace(p,quietly = TRUE))) {
      missing = c(missing,p)
    }
  }
  if (length(missing)>0) {
    msg = sprintf("For this function, you need additional packages. Missing package%s: %s. Please install them e.g. using \n BiocManager::install(c('%s'))\n",
                  if (length(missing)>1) "s" else "",paste(missing,collapse=","),paste(missing,collapse="','"))
    if (error) stop(msg)
    if (warn) warning(msg)
    return(FALSE)
  }
  return(TRUE)
}

