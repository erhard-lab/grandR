

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
    pmin(10^(ceiling(log10(smallest_diff))), 1)
  }
}


cnt=function(m) {
  m=as.matrix(m)
  mode(m) <- "integer"
  m
}

get.varargs.names=function(...) sapply(as.list(substitute(list(...)))[-1L],as.character)

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
opt$lapply=function(...) lapply(...)
opt$sapply=function(...) simplify2array(opt$lapply(...))
SetParallel=function(cores=max(1,parallel::detectCores()-2)) {
  if (cores>1) {
    if (.Platform$OS.type!="unix") stop("Parallelism is not supported under windows!")
    opt$lapply<-function(...) parallel::mclapply(...,mc.cores=cores)
  } else {
    opt$lapply<-function(...) lapply(...)
  }
}
IsParallel=function() grepl("mclapply",capture.output(print(opt$lapply))[1])
