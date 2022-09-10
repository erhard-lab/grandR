

as.Seurat.grandR=function(d,modalities=c(RNA="total",newRNA="new"),hls=NULL,mode=c("assay","cells","genes","list")) {

  if (length(modalities)==0) stop("No modality specified!")

  mats=list(total=GetSparseMatrix(d,mode.slot='count'))
  rows=Matrix::rowSums(mats$total)>0
  cols=Matrix::colSums(mats$total)>0
  mats$total=mats$total[rows,cols]
  mats$ntr=GetSparseMatrix(d,mode.slot='ntr')[rows,cols]
  if (any(c("old","new","prev") %in% modalities)) {
    mats$new=round(mats$total*mats$ntr)
    if (any(c("old","prev") %in% modalities)) mats$old=mats$total-mats$new
    if ("prev" %in% modalities) {
      stopifnot(!is.null(hls));
      mats$prev=mats$old*exp(2*log(2)/pmin(pmax(hls,0.25),24))
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
    s=Seurat::CreateSeuratObject(counts=mats[[1]],assay=names(mats)[1],project=basename(d$prefix))
    if (length(mats)>1) for (i in 2:length(mats)) s[[names(mats)[i]]]=Seurat::CreateAssayObject(mats[[i]])
    s
  },cells={
    for (i in 1:length(mats)) colnames(mats[[i]])=paste(colnames(mats[[i]]),names(mats)[i],sep=".")
    mat=do.call("cbind",mats)
    mode=do.call("c",lapply(1:length(mats),function(i) rep(names(mats)[i],ncol(mats[[i]]))))
    s=Seurat::CreateSeuratObject(counts=mat,project=basename(d$prefix))
    s[["mode"]]=factor(mode,levels=unique(mode))
    s
  },genes={
    for (i in 1:length(mats)) rownames(mats[[i]])=paste(rownames(mats[[i]]),names(mats)[i],sep=".")
    mat=do.call("rbind",mats)
    Seurat::CreateSeuratObject(counts=mat,project=basename(d$prefix))
  },list={
    re=lapply(1:length(mats),function(i) Seurat::CreateSeuratObject(counts=mats[[i]],project=names(mats)[i]))
    names(re)=names(mats)
    re
  })
  append.meta=function(s) {
    for (i in seq_len(ncol(d$coldata))) s[[names(d$coldata)[i]]]=rep(d$coldata[cols,i],nrow(s[[]])/nrow(d$coldata[cols,]))
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
