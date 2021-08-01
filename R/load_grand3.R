
library(Matrix)

IsSparse=function(prefix) !file.exists(paste0(prefix,".targets/data.tsv.gz"))

ReadGrand3=function(prefix,...) if (IsSparse(prefix)) ReadGrand3_sparse(prefix,...) else ReadGrand3_dense(prefix,...)

ReadGrand3_sparse=function(prefix, verbose=FALSE, design=c(Design$Library,Design$Sample,Design$Barcode), label="4sU",estimator="Binom") {
  
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

  if (verbose) cat("Reading NTRs...\n")
  ntr=as(readMM(sprintf("%s.targets/%s.%s.ntr.mtx.gz",prefix,label,estimator)),Class = "dgCMatrix")
  colnames(ntr)=cols
  rownames(ntr)=gene.info$Symbol
  
  gene.info$Type=gsub(".*\\(","",gsub(")","",gene.info$Category,fixed=TRUE))
  gene.info$Type=factor(gene.info$Type,levels=unique(gene.info$Type))
  
  re=list()
  re$count=count
  re$ntr=ntr
  
  coldata=data.frame(Name=cols)
  spl=strsplit(as.character(coldata$Name),".",fixed=TRUE)
  for (i in 1:length(design)) coldata=cbind(coldata,factor(sapply(spl,function(v) v[i]),levels=unique(sapply(spl,function(v) v[i]))))
  names(coldata)[-1]=design
  rownames(coldata)=coldata$Name
  
  invisible(grandR(prefix,gene.info,re,coldata))
}
 
ReadGrand3_dense=function(prefix, verbose=FALSE, design=c(design$Condition,Design$Replicate)) {
  
}

as.Seurat.grandR=function(d,old=TRUE,new=TRUE,ntr=FALSE,prev=FALSE,hls=NULL,mode=c("mode","cells","genes")) {

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
  s=switch(mode[1],mode={
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
  })
  for (i in seq_len(ncol(d$coldata))) s[[names(d$coldata)[i]]]=rep(d$coldata[cols,i],nrow(s[[]])/nrow(d$coldata[cols,]))
  invisible(s)
}


    
