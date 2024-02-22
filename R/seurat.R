

#' Create Seurat object from a grandR object
#'
#' @param data a grandR object
#' @param modalities vector defining modalities to include in the Seurat object (see details)
#' @param hls half-lives for computing previous RNA, only required for "prev" modality (see details)
#' @param time labeling time, only required for "prev" modality (see details)
#' @param mode how to integrate modalities into seurat object (see details)
#'
#' @details Modalities must be a named character vector. The only allowed elements are "total" (total counts),
#' "new" (new counts), "old" (old counts), "prev" (estimated previous time point counts). The names of the elements are further
#' used depending on mode.
#'
#' @details To compute the previous time point counts, a vector of half lives and the labeling time is required. The half-lives
#' must be given in the correct order (same as in the grandR object).
#'
#' @details The mode parameter defines how the defined modalities are represented in the Seurat object. "assay" means
#' that for each modality, the Seurat object will contain an assay (named according to the corresponding name in modalities).
#' "cells" means that cells will be copied for each modality and cell names are prefixed by the corresponding name in modalities
#'  (i.e., if the grandR object has 1000 cells named c1,...,c1000, and modalities=c(RNA="total",newRNA="new"), the Seurat object
#'  will have 2000 cells named RNA.c1,...,RNA.c1000,newRNA.c1,...,newRNA.c1000). "genes" means that genes fill be copied for each
#'  modality and gene names are prefixed by the corresponding name in modalities. "list" means that instead of a single Seurat object,
#'  a list of Seurat objects is returned.
#'
#' @return a Seurat object
#'
#' @exportS3Method  Seurat::as.Seurat grandR
#' @export as.Seurat.grandR
#' @export
#'
#' @concept load
as.Seurat.grandR=function(data,modalities=c(RNA="total",newRNA="new"),hls=NULL,time=NULL,mode=c("assay","cells","genes","list")) {

  checkPackages(c("Seurat"))

  if (length(modalities)==0) stop("No modality specified!")

  mats=list(total=GetMatrix(data,mode.slot='count'))
  rows=Matrix::rowSums(mats$total)>0
  cols=Matrix::colSums(mats$total)>0
  mats$total=mats$total[rows,cols]
  mats$ntr=GetMatrix(data,mode.slot='ntr')[rows,cols]
  if (any(c("old","new","prev") %in% modalities)) {
    mats$new=round(mats$total*mats$ntr)
    if (any(c("old","prev") %in% modalities)) mats$old=mats$total-mats$new
    if ("prev" %in% modalities) {
      stopifnot(!is.null(hls));
      mats$prev=mats$old*exp(time*log(2)/pmin(pmax(hls,0.25),24))
      rownames(mats$prev)=rownames(mats$total)
      colnames(mats$prev)=colnames(mats$total)
    }
  }
#  if (any(c("old.lower","new.upper") %in% modalities)) {
#    if(is.null(d$data$upper)) stop("You need to load Grand3 data including CIs! set read.CI=TRUE when calling ReadGrand3!")
#    mats$new.upper=round(mats$total*d$data$upper[rows,cols])
#    if ("old.lower" %in% modalities) mats$old.lower=mats$total-mats$new.upper
#  }
#  if (any(c("new.lower","old.upper") %in% modalities)) {
#    if(is.null(d$data$lower)) stop("You need to load Grand3 data including CIs! set read.CI=TRUE when calling ReadGrand3!")
#    mats$new.lower=round(mats$total*d$data$lower[rows,cols])
#    if ("old.upper" %in% modalities) mats$old.upper=mats$total-mats$new.lower
#  }
  if (!all(modalities %in% names(mats))) stop("Modalities unknown! Can be any of total,new,old,ntr,prev!")

  mats=mats[modalities]
  if (!is.null(names(modalities))) names(mats)=names(modalities)

  re=switch(mode[1],assay={
    s=Seurat::CreateSeuratObject(counts=mats[[1]],assay=names(mats)[1],project=basename(data$prefix))
    if (length(mats)>1) for (i in 2:length(mats)) s[[names(mats)[i]]]=Seurat::CreateAssayObject(mats[[i]])
    s
  },cells={
    for (i in 1:length(mats)) colnames(mats[[i]])=paste(colnames(mats[[i]]),names(mats)[i],sep=".")
    mat=do.call("cbind",mats)
    mode=do.call("c",lapply(1:length(mats),function(i) rep(names(mats)[i],ncol(mats[[i]]))))
    s=Seurat::CreateSeuratObject(counts=mat,project=basename(data$prefix))
    s[["mode"]]=factor(mode,levels=unique(mode))
    s
  },genes={
    for (i in 1:length(mats)) rownames(mats[[i]])=paste(rownames(mats[[i]]),names(mats)[i],sep=".")
    mat=do.call("rbind",mats)
    Seurat::CreateSeuratObject(counts=mat,project=basename(data$prefix))
  },list={
    re=lapply(1:length(mats),function(i) Seurat::CreateSeuratObject(counts=mats[[i]],project=names(mats)[i]))
    names(re)=names(mats)
    re
  })
  append.meta=function(s) {
    for (i in seq_len(ncol(data$coldata))) s[[names(data$coldata)[i]]]=rep(data$coldata[cols,i],nrow(s[[]])/nrow(data$coldata[cols,]))
    s
  }
  re = if ("Seurat" %in% class(re)) append.meta(re) else lapply(re,append.meta)
  invisible(re)
}

as.Seurat.legacy.grandR=function(d,old=TRUE,new=TRUE,ntr=FALSE,prev=FALSE,hls=NULL,mode=c("assay","cells","genes","list")) {

  total.mat=d$data$count
  rows=rowSums(total.mat)>0
  cols=colSums(total.mat)>0
  total.mat=total.mat[rows,cols]
  ntr.mat=d$data$ntr[rows,cols]
  if (old||new||prev) {
    new.mat=round(total.mat*ntr.mat)
    old.mat=total.mat-new.mat
    if (prev) {
      stopifnot(!is.null(hls));
      prev.mat=old.mat*exp(2*log(2)/pmin(pmax(hls,0.25),24))
    }
  }
  re=switch(mode[1],assay={
    s=Seurat::CreateSeuratObject(counts=total.mat,project=basename(d$prefix))
    if (ntr) s[["NTR"]]=Seurat::CreateAssayObject(ntr.mat)
    if (new) s[["newRNA"]]=Seurat::CreateAssayObject(new.mat)
    if (old) s[["oldRNA"]]=Seurat::CreateAssayObject(old.mat)
    if (prev) s[["prevRNA"]]=Seurat::CreateAssayObject(prev.mat)
    s
  },cells={
    colnames(total.mat)=paste(colnames(total.mat),"total",sep=".")
    mode=rep(c("total"),ncol(total.mat))
    if (ntr) {
      colnames(ntr.mat)=paste(colnames(ntr.mat),"new",sep=".")
      total.mat=cbind(total.mat,ntr.mat)
      mode=c(mode,rep(c("ntr"),ncol(ntr.mat)))
    }
    if (new) {
      colnames(new.mat)=paste(colnames(new.mat),"new",sep=".")
      total.mat=cbind(total.mat,new.mat)
      mode=c(mode,rep(c("new"),ncol(new.mat)))
    }
    if (old) {
      colnames(old.mat)=paste(colnames(old.mat),"old",sep=".")
      total.mat=cbind(total.mat,old.mat)
      mode=c(mode,rep(c("old"),ncol(old.mat)))
    }
    if (prev) {
      colnames(prev.mat)=paste(colnames(prev.mat),"prev",sep=".")
      total.mat=cbind(total.mat,prev.mat)
      mode=c(mode,rep(c("prev"),ncol(prev.mat)))
    }
    s=Seurat::CreateSeuratObject(counts=total.mat,project=basename(d$prefix))
    s[["mode"]]=factor(mode,levels=unique(mode))
    s
  },genes={
    if (ntr) {
      rownames(ntr.mat)=paste(rownames(ntr.mat),"ntr",sep=".")
      total.mat=rbind(total.mat,ntr.mat)
    }
    if (new) {
      rownames(new.mat)=paste(rownames(new.mat),"new",sep=".")
      total.mat=rbind(total.mat,new.mat)
    }
    if (old) {
      rownames(old.mat)=paste(rownames(old.mat),"old",sep=".")
      total.mat=rbind(total.mat,old.mat)
    }
    if (prev) {
      rownames(prev.mat)=paste(rownames(prev.mat),"prev",sep=".")
      total.mat=rbind(total.mat,prev.mat)
    }
    Seurat::CreateSeuratObject(counts=total.mat,project=basename(d$prefix))
  },list={
    s=list(RNA=Seurat::CreateSeuratObject(counts=total.mat,project="RNA"))
    if (ntr) s[["NTR"]]=Seurat::CreateSeuratObject(ntr.mat,project="NTR")
    if (new) s[["newRNA"]]=Seurat::CreateSeuratObject(new.mat,project="newRNA")
    if (old) s[["oldRNA"]]=Seurat::CreateSeuratObject(old.mat,project="oldRNA")
    if (prev) s[["prevRNA"]]=Seurat::CreateSeuratObject(prev.mat,project="prevRNA")
    s
  })
  append.meta=function(s) {
    for (i in seq_len(ncol(d$coldata))) s[[names(d$coldata)[i]]]=rep(d$coldata[cols,i],nrow(s[[]])/nrow(d$coldata[cols,]))
    s
  }
  re = if ("Seurat" %in% class(re)) append.meta(re) else lapply(re,append.meta)
  invisible(re)
}
#' Create Pseudobulk Table from a Seurat object
#'
#' @param data a Seurat object
#' @param name.column name of the metadata column containing the sample/cell names. Default "Name".
#' @param pseudobulk.column name of the metadata column containing the Pseudobulk names. Default "Condition".
#'
#' @details This function returns a table which can be used as input for GRAND3
#'
#' @return a table with two columns "Cell" and "Pseudobulk"
#'
#' @concept data
CreatePseudobulkTable <- function(data,name.column="Name",pseudobulk.column="Condition") {
  table = data[[]][,c(pseudobulk.column, name.column)]
  rownames(table) = seq(nrow(table))
  colnames(table) = c("Pseudobulk","Cell")
  return(table)
}

#' Create Convolution Table from a Seurat object
#'
#' @param data a Seurat object
#' @param n.neighbors the number of neighbors to be convoluted
#'
#' @details This function returns a table which can be used as input for GRAND3. Note that a data set contatining multiple time points should be split before convolution.
#'
#' @return a table with two columns "Cell" and "Pseudobulk"
#'
#' @concept data
CreateConvolutionTable<- function(data,n.neighbors=20) {
  checkPackages(c("Seurat"))

  data <- Seurat::FindNeighbors(data, dims = 1:10, k.param = n.neighbors)
  knn <- data@graphs$RNA_nn
  knn <- data.frame(knn)
  tab <- matrix(nrow=0, ncol=2)
  colnames(tab) <- c("Cell", "Pseudobulk")

  for(x in (1:length(knn[1,]))){
    for(y in colnames(knn[x,(knn[x,] == 1)])){
      tab <- rbind(tab, c(y, rownames(knn[x,])))
    }
  }
  return(tab)
}
