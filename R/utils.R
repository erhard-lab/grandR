

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

has.package=function(name) name %in% rownames(installed.packages())

GetField=function(name,field,sep=".") sapply(strsplit(as.character(name),sep,fixed=TRUE),function(v) v[field])

opt <- new.env()
opt$lapply=function(...) lapply(...)
opt$sapply=function(...) simplify2array(opt$lapply(...))
SetParallel=function(cores=max(1,parallel::detectCores()-2)) {
  if (cores>1) {
    opt$lapply<-function(...) parallel::mclapply(...,mc.cores=cores)
  } else {
    opt$lapply<-function(...) lapply(...)
  }
}

