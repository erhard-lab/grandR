
library(Matrix)

IsSparse=function(prefix) !file.exists(paste0(prefix,".targets/data.tsv.gz"))

ReadGrand3=function(prefix,...) if (IsSparse(prefix)) ReadGrand3_sparse(prefix,...) else ReadGrand3_dense(prefix,...)

ReadGrand3_sparse=function(prefix, verbose=FALSE, design=c(Design$Library,Design$Sample,Design$Barcode), label="4sU",estimator="Binom", read.CI=FALSE) {
  
  cols=readLines(paste0(prefix, ".targets/barcodes.tsv.gz"))
  conds=strsplit(cols,".",fixed=TRUE)[[1]]
  if (length(conds)!=length(design)) stop(paste0("Design parameter is incompatible with input data: ",paste(conds,collapse=".")))
  
  
  if (verbose) cat("Reading count matrix...\n")
  count=as(readMM(paste0(prefix, ".targets/matrix.mtx.gz")),Class = "dgCMatrix")
  gene.info=setNames(read.delim(paste0(prefix, ".targets/features.tsv.gz"),header = FALSE,stringsAsFactors = FALSE),c("Gene","Symbol","Mode","Category","Length"))
    
  if (anyDuplicated(gene.info$Gene)) {
    dupp=table(gene.info$Gene)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Gene=make.unique(gene.info$Gene)
  }
  if (anyDuplicated(gene.info$Symbol)) {
    dupp=table(gene.info$Symbol)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene symbols (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Symbol=make.unique(gene.info$Symbol)
  }

  if (any(make.names(gene.info$Symbol)!=gene.info$Symbol)) {
  	ill=gene.info$Symbol[which(make.names(gene.info$Symbol)!=gene.info$Symbol)]
    warning(sprintf("Illegal identifiers among the gene symbols (e.g. %s), making legal!",paste(head(ill),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Symbol=make.names(gene.info$Symbol)
  }

  colnames(count)=cols
  rownames(count)=gene.info$Symbol
  re=list()
  re$count=count
  
  if (verbose) cat("Reading NTRs...\n")
  ntr=as(readMM(sprintf("%s.targets/%s.%s.ntr.mtx.gz",prefix,label,estimator)),Class = "dgCMatrix")
  colnames(ntr)=cols
  rownames(ntr)=gene.info$Symbol
  re$ntr=ntr
  
  if (read.CI && file.exists(sprintf("%s.targets/%s.%s.lower.mtx.gz",prefix,label,estimator)) && file.exists(sprintf("%s.targets/%s.%s.upper.mtx.gz",prefix,label,estimator))) {
    if (verbose) cat("Reading CIs...\n")
    lower=as(readMM(sprintf("%s.targets/%s.%s.lower.mtx.gz",prefix,label,estimator)),Class = "dgCMatrix")
    colnames(lower)=cols
    rownames(lower)=gene.info$Symbol
    re$lower=lower
    
    upper=as(readMM(sprintf("%s.targets/%s.%s.upper.mtx.gz",prefix,label,estimator)),Class = "dgCMatrix")
    colnames(upper)=cols
    rownames(upper)=gene.info$Symbol
    re$upper=upper
  }
  gene.info$Type=gsub(".*\\(","",gsub(")","",gene.info$Category,fixed=TRUE))
  gene.info$Type=factor(gene.info$Type,levels=unique(gene.info$Type))
  
  
  # use make coldata instead. add no4sU column!
  coldata=data.frame(Name=cols)
  spl=strsplit(as.character(coldata$Name),".",fixed=TRUE)
  for (i in 1:length(design)) coldata=cbind(coldata,factor(sapply(spl,function(v) v[i]),levels=unique(sapply(spl,function(v) v[i]))))
  names(coldata)[-1]=design
  rownames(coldata)=coldata$Name
  
  invisible(grandR(prefix,gene.info,re,coldata))
}
 
ReadGrand3_dense=function(prefix, verbose=FALSE, design=c(design$Condition,Design$Replicate)) {
  
}

as.Seurat.grandR=function(d,modalities=c(RNA="total",newRNA="new"),hls=NULL,mode=c("assay","cells","genes","list")) {
	
  if (length(modalities)==0) stop("No modality specified!")

  mats=list(total=d$data$count)
  rows=rowSums(mats$total)>0
  cols=colSums(mats$total)>0
  mats$total=mats$total[rows,cols]
  mats$ntr=d$data$ntr[rows,cols]
  if (any(c("old","new","prev") %in% modalities)) {
    mats$new=round(mats$total*mats$ntr)
    if (any(c("old","prev") %in% modalities)) mats$old=mats$total-mats$new
    if ("prev" %in% modalities) {
      stopifnot(!is.null(hls));
      prev.mat=old.mat*exp(2*log(2)/pmin(pmax(hls,0.25),24))
    }
  }
  if (any(c("old.lower","new.upper") %in% modalities)) {
    if(is.null(d$data$upper)) stop("You need to load Grand3 data including CIs! set read.CI=TRUE when calling ReadGrand3!")
    mats$new.upper=round(mats$total*d$data$upper[rows,cols])
    if ("old.lower" %in% modalities) mats$old.lower=mats$total-mats$new.upper  
  }
  if (any(c("new.lower","old.upper") %in% modalities)) {
    if(is.null(d$data$lower)) stop("You need to load Grand3 data including CIs! set read.CI=TRUE when calling ReadGrand3!")
    mats$new.lower=round(mats$total*d$data$lower[rows,cols])
    if ("old.upper" %in% modalities) mats$old.upper=mats$total-mats$new.lower
  }
  if (!all(modalities %in% names(mats))) stop("Modalities unknown! Can be any of total,new,old,ntr,prev,new.lower,new.upper,old.lower,old.upper!")

  mats=mats[modalities]

  re=switch(mode[1],assay={
    s=CreateSeuratObject(counts=mats[[1]],assay=names(modalities)[1],project=basename(d$prefix))
    if (length(mats)>1) for (i in 2:length(mats)) s[[names(modalities)[i]]]=CreateAssayObject(mats[[i]])
    s
  },cells={
    for (i in 1:length(mats)) colnames(mats[[i]])=paste(colnames(mats[[i]]),names(mats)[i],sep=".")
    mat=do.call("cbind",mats)
    mode=do.call("c",lapply(1:length(mats),function(i) rep(names(mats)[i],ncol(mats[[i]]))))
    s=CreateSeuratObject(counts=mat,project=basename(d$prefix))
    s[["mode"]]=factor(mode,levels=unique(mode))
    s
  },genes={
    for (i in 1:length(mats)) rownames(mats[[i]])=paste(rownames(mats[[i]]),names(mats)[i],sep=".")
    mat=do.call("rbind",mats)
    CreateSeuratObject(counts=mat,project=basename(d$prefix))
  },list={
    re=lapply(1:length(mats),function(i) CreateSeuratObject(counts=mats[[i]],project=names(mats)[i]))
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
    s=CreateSeuratObject(counts=total.mat,project=basename(d$prefix))
    if (ntr) s[["NTR"]]=CreateAssayObject(ntr.mat)
    if (new) s[["newRNA"]]=CreateAssayObject(new.mat)
    if (old) s[["oldRNA"]]=CreateAssayObject(old.mat)
    if (prev) s[["prevRNA"]]=CreateAssayObject(prev.mat)
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
    s=CreateSeuratObject(counts=total.mat,project=basename(d$prefix))
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
    CreateSeuratObject(counts=total.mat,project=basename(d$prefix))
  },list={
    s=list(RNA=CreateSeuratObject(counts=total.mat,project="RNA"))
    if (ntr) s[["NTR"]]=CreateSeuratObject(ntr.mat,project="NTR")
    if (new) s[["newRNA"]]=CreateSeuratObject(new.mat,project="newRNA")
    if (old) s[["oldRNA"]]=CreateSeuratObject(old.mat,project="oldRNA")
    if (prev) s[["prevRNA"]]=CreateSeuratObject(prev.mat,project="prevRNA")
    s
  })
  append.meta=function(s) {
    for (i in seq_len(ncol(d$coldata))) s[[names(d$coldata)[i]]]=rep(d$coldata[cols,i],nrow(s[[]])/nrow(d$coldata[cols,]))
    s
  }
  re = if ("Seurat" %in% class(re)) append.meta(re) else lapply(re,append.meta)
  invisible(re)
}


    

ReadNewTotal=function(genes, cells, new.matrix, total.matrix, detection.rate=1,verbose=FALSE) {
  
  cols=read.csv(cells,check.names = FALSE,stringsAsFactors = FALSE)
  gene.info=setNames(read.csv(genes,check.names = FALSE,stringsAsFactors = FALSE),c("Gene","Biotype","Symbol"))
  
  if (verbose) cat("Reading total count matrix...\n")
  count=as(readMM(total.matrix),Class = "dgCMatrix")
  
  if (anyDuplicated(gene.info$Gene)) {
    dupp=table(gene.info$Gene)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Gene=make.unique(gene.info$Gene)
  }
  if (anyDuplicated(gene.info$Symbol)) {
    dupp=table(gene.info$Symbol)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene symbols (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Symbol=make.unique(gene.info$Symbol)
  }
  
  if (any(make.names(gene.info$Symbol)!=gene.info$Symbol)) {
    ill=gene.info$Symbol[which(make.names(gene.info$Symbol)!=gene.info$Symbol)]
    warning(sprintf("Illegal identifiers among the gene symbols (e.g. %s), making legal!",paste(head(ill),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Symbol=make.names(gene.info$Symbol)
  }
  
  colnames(count)=cols[,1]
  rownames(count)=gene.info$Symbol
  
  if (verbose) cat("Reading new count matrix...\n")
  new=as(readMM(new.matrix),Class = "dgCMatrix")
  
  if (verbose) cat("Computing NTRs...\n")
  new=new/detection.rate
  ntr=new/count
  ntr@x[ntr@x>1]=1
  ntr@x[is.nan(ntr@x)]=0
  ntr=as(ntr,Class = "dgCMatrix")
  
  colnames(ntr)=colnames(count)
  rownames(ntr)=rownames(count)
  
  re=list()
  re$count=count
  re$ntr=ntr
  
  invisible(grandR("",gene.info,re,cols))
}
