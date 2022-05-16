

read.tsv=function(t,verbose=FALSE,stringsAsFactors=TRUE,...) {

  if (has.package("RCurl") && RCurl::url.exists(t)) {
    file <- tempfile()
    if (verbose) cat("Downloading file...\n")
    download.file(url, file, quiet=!verbose)
    if (verbose) cat("Reading file...\n")
    t=read.delim(file,stringsAsFactors=FALSE,check.names=FALSE,...)
    if (verbose) cat("Deleting temporary file...\n")
    unlink(file)
  } else {
    if (verbose) cat("Reading file...\n")
    t=read.delim(t,stringsAsFactors=FALSE,check.names=FALSE,...)
  }

  if (stringsAsFactors==TRUE) t=as.data.frame(lapply(t,function(c) if (is.character(c)) factor(c,levels=unique(c)) else c),check.names=FALSE)
  t
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
structure2vector=function(d,return.fields=NULL,return.extra=NULL) {
  r=list()
  if (!is.null(return.fields)) r=c(r,d[return.fields])
  if (!is.null(return.extra)) r=c(r,return.extra(d))
  unlist(r)
}
#' @describeIn structure2vector Convert the output of the FitKinetics methods into a vector
#' @param condition if the original grandR object had \code{\link{Condition}} set, which condition to extract (NULL otherwise)
#' @export
kinetics2vector=function(d,condition=NULL,return.fields=c("Synthesis","Half-life","rmse"),return.extra=NULL) structure2vector(if (is.null(condition)) d else d[[condition]],return.fields=return.fields,return.extra=return.extra)

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
estimate.dispersion=function(ss) {
  dds=DESeq2::DESeqDataSetFromMatrix(countData = cnt(ss),colData=data.frame(rep(1,ncol(ss))),design = ~1)
  dds=DESeq2::estimateSizeFactors(dds)
  dds=DESeq2::estimateDispersions(dds,quiet=TRUE)
  disp=DESeq2::dispersions(dds)
  disp[is.na(disp)]=0.1 # ss=0 0 0 ...
  disp
}

has.package=function(name) name %in% rownames(installed.packages())

GetField=function(name,field,sep=".") sapply(strsplit(as.character(name),sep,fixed=TRUE),function(v) v[field])


psapply=function(...,seed=NULL) simplify2array(plapply(...,seed=seed))
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
opt$lapply=function(...) lapply2(...)
opt$sapply=function(...) simplify2array(opt$lapply(...))
opt$nworkers=0
SetParallel=function(cores=max(1,parallel::detectCores()-2)) {
  opt$nworkers=cores
  if (cores>1) {
    if (.Platform$OS.type!="unix") stop("Parallelism is not supported under windows!")
    opt$lapply<-function(...) parallel::mclapply(...,mc.cores=cores)
  } else {
    opt$lapply<-function(...) lapply2(...)
  }
}
IsParallel=function() opt$nworkers>1


#' Wrapper around lapply to track progress
#'
#' @param X         a vector (atomic or list) or an expressions vector. Other
#'                  objects (including classed objects) will be coerced by
#'                  as.list
#' @param FUN       the function to be applied to
#' @param ...       optional arguments to FUN
#' @param progress track progress?
#' @param style    style of progress bar (see txtProgressBar)
#'
#' @examples
#' x <- lapply2(1:1000, function(i, y) Sys.sleep(0.01))
#' x <- lapply2(1:3, function(i, y) Sys.sleep(1))
#'
#' dat <- lapply(1:10, function(x) rnorm(100))
#' func <- function(x, arg1) mean(x)/arg1
#' lapply2(dat, func, arg1=10)
lapply2 <- function(..., progress=FALSE,style=3)
{
  if (!progress) return(lapply(...))

  l=list(...)
  X=if ("X" %in% names(l)) l[["X"]] else l[[1]]

  env=environment()
  s=0
  progress=txtProgressBar(min = 0, max = length(X), style = style)

  FUN2=function(...){
    v=get("s", envir = env)
    assign("s", v+1 ,envir=env)
    setTxtProgressBar(get("progress", envir=env), v+1)
    FUN(...)
  }
  re=lapply(X, FUN2, ...)
  close(progress)
  re
}

#' Wrapper around mclapply to track progress
#'
#' Based on http://stackoverflow.com/questions/10984556
#'
#' @param X         a vector (atomic or list) or an expressions vector. Other
#'                  objects (including classed objects) will be coerced by
#'                  as.list
#' @param FUN       the function to be applied to
#' @param ...       optional arguments to FUN
#' @param mc.preschedule see mclapply
#' @param mc.set.seed see mclapply
#' @param mc.silent see mclapply
#' @param mc.cores see mclapply
#' @param mc.cleanup see mclapply
#' @param mc.allow.recursive see mclapply
#' @param progress track progress?
#' @param style    style of progress bar (see txtProgressBar)
#'
#' @examples
#' x <- mclapply2(1:1000, function(i, y) Sys.sleep(0.01))
#' x <- mclapply2(1:3, function(i, y) Sys.sleep(1), mc.cores=1)
#'
#' dat <- lapply(1:10, function(x) rnorm(100))
#' func <- function(x, arg1) mean(x)/arg1
#' mclapply2(dat, func, arg1=10, mc.cores=2)
mclapply2 <- function(X, FUN, ...,
                      mc.preschedule = TRUE, mc.set.seed = TRUE,
                      mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
                      mc.cleanup = TRUE, mc.allow.recursive = TRUE,
                      progress=FALSE, style=3)
{
  if (!is.vector(X) || is.object(X)) X <- as.list(X)

  if (progress && Sys.getenv("RSTUDIO") == "1") {
    if(getOption("RStudio.progress.warning",TRUE)) {
      warning("Showing progress bars for parallelized  work does not work under rstudio! This message will only be shown once in this session!")
      options("RStudio.progress.warning"=FALSE)
    }
    progress=FALSE
  }


  if (progress) {
    f <- fifo(tempfile(), open="w+b", blocking=T)
    p <- parallel:::mcfork()
    pb <- txtProgressBar(0, length(X), style=style)
    setTxtProgressBar(pb, 0)
    progress <- 0
    if (inherits(p, "masterProcess")) {
      while (progress < length(X)) {
        readBin(f, "double")
        progress <- progress + 1
        setTxtProgressBar(pb, progress)
      }
      cat("\n")
      parallel:::mcexit()
    }
  }
  tryCatch({
    result <- parallel::mclapply(X, ..., function(...) {
      res <- FUN(...)
      if (progress) writeBin(1, f)
      res
    },
    mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
    mc.silent = mc.silent, mc.cores = mc.cores,
    mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
    )

  }, finally = {
    if (progress) close(f)
  })
  result
}
