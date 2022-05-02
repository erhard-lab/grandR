
comp.hl=function(p,time=1) ifelse(p==0,Inf,log(2)/(-1/time*log(1-p)))

comp.tpm=function(cmat,lengths,subset=NULL) {
  zerolen=lengths==0
  lengths[zerolen]=1
  rpk=cmat/(lengths/1000)
  rpk[zerolen,]=NA
  scale=colSums(if(!is.null(subset)) rpk[subset,] else rpk,na.rm=T)/1E6
  re=t(t(rpk)/scale)
  re
}

comp.fpkm=function(cmat,lengths,subset=NULL) {
  scale=colSums(if(!is.null(subset)) cmat[subset,] else cmat,na.rm=T)/1E6
  rpm=t(t(cmat)/scale)

  zerolen=lengths==0
  lengths[zerolen]=1
  re=rpm/(lengths/1000)
  re[zerolen,]=NA
  re
}
comp.rpm=function(cmat,subset=NULL) {
  scale=colSums(if(!is.null(subset)) cmat[subset,] else cmat,na.rm=T)/1E6
  re=t(t(cmat)/scale)
  re
}


ComputeNtrCI=function(data,CI.size=0.95,name.lower="lower",name.upper="upper") {
  data=ComputeNtrPosteriorLower(data=data,CI.size=CI.size,name=name.lower)
  data=ComputeNtrPosteriorUpper(data=data,CI.size=CI.size,name=name.upper)
  data
}
ComputeNtrPosteriorLower=function(data,CI.size=0.95,name="lower") ComputeNtrPosteriorQuantile(data=data,quantile=(1-CI.size)/2,name=name)
ComputeNtrPosteriorUpper=function(data,CI.size=0.95,name="upper") ComputeNtrPosteriorQuantile(data=data,quantile=1-(1-CI.size)/2,name=name)

ComputeNtrPosteriorQuantile=function(data,quantile,name) {
  a=as.matrix(GetTable(data,type="alpha",name.by = "Gene"))
  b=as.matrix(GetTable(data,type="beta",name.by = "Gene"))
  v=qbeta(quantile,a,b)
  AddSlot(data,name,v)
}

ComputeSteadyStateHalfLives=function(data,time=Design$dur.4sU,name, max.HL=48, CI.size=0.95, compute.CI=FALSE, as.analysis=FALSE) {
  if (is.character(time) && length(time)==1) time=data$coldata[[time]]

  ntrs=as.matrix(GetTable(data,type="ntr",name.by = "Gene"))
  if (length(time)==1) time=rep(time,ncol(ntrs))
  stopifnot(ncol(ntrs)==length(time))

  if (compute.CI) {
    as.analysis=TRUE
    data=ComputeNtrCI(data,CI.size = CI.size)
    lower=as.matrix(GetTable(data,type="lower",name.by = "Gene"))
    upper=as.matrix(GetTable(data,type="upper",name.by = "Gene"))
    hls = do.call(cbind,lapply(1:length(time),function(i) cbind(
                                              pmin(comp.hl(p = upper[,i],time = time[i]),max.HL),
                                              pmin(comp.hl(p = ntrs[,i],time = time[i]),max.HL),
                                              pmin(comp.hl(p = lower[,i],time = time[i]),max.HL)
                                              )))
    colnames(hls)=paste0(rep(c("lower.","tMAP.","upper."),ncol(ntrs)),rep(colnames(ntrs),each=3))
  }
  else {
    hls = sapply(1:length(time),function(i) comp.hl(p = ntrs[,i],time = time[i]))
    colnames(hls)=colnames(ntrs)
  }

  if (as.analysis) {
    AddAnalysis(data,name = name,table = as.data.frame(hls))
  } else {
    AddSlot(data,name,hls)
  }
}


ComputeAbsolute=function(data,dilution=2E6,volume=200,name="absolute") {
  fd=data.frame(gene_short_name=rownames(data$data$tpm))
  rownames(fd)=rownames(data$data$tpm)
  ercc=data$data$tpm[data$gene.info$Type=='ERCC',]

  cds <- newCellDataSet(data$data$tpm,
                        featureData = new ("AnnotatedDataFrame",fd),
                        lowerDetectionLimit = 0.1,
                        expressionFamily = tobit(Lower = 0.1))
  rpc_matrix <- monocle::relative2abs(cds, method = "num_genes", ERCC_controls=ercc, ERCC_annotation=monocle::spike_df,dilution=dilution, volume=volume)
  rpc_matrix[data$data$tpm==0]=0

  colnames(rpc_matrix)=colnames(data$data$tpm)
  rownames(rpc_matrix)=rownames(data$data$tpm)
  rpc_matrix[,grepl("test",colnames(rpc_matrix))]=NA
  AddSlot(data,name,rpc_matrix)
}



Normalize=function(data,sizeFactors=NULL,genes=NULL,name="norm",slot="count",return.sf=FALSE,set.to.default=TRUE) {
  stopifnot(is.grandR(data))
  if (is.null(genes)) genes=1:nrow(data)
  mat=as.matrix(GetTable(data,type=slot,genes=genes,ntr.na = FALSE,name.by = "Gene"))
  if (is.null(sizeFactors)) sizeFactors=DESeq2::estimateSizeFactorsForMatrix(mat)
  if (return.sf) return(sizeFactors)

  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  data=AddSlot(data,name,matrix = t(t(mat)/sizeFactors),set.to.default=set.to.default)
  data
}
NormalizeFPKM=function(data,genes=Genes(data),tlen=data$gene.info$Length,name="fpkm",slot="count",set.to.default=TRUE) {
  genes=ToIndex(data,genes)
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  data=AddSlot(data,name,comp.fpkm(mat,tlen,subset = genes),set.to.default=set.to.default)
  data
}
NormalizeRPM=function(data,genes=Genes(data),name="rpm",slot="count",set.to.default=TRUE) {
  genes=ToIndex(data,genes)
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  data=AddSlot(data,name,comp.rpm(mat,subset = genes),set.to.default=set.to.default)
  data
}
NormalizeTPM=function(data,genes=Genes(data),tlen=data$gene.info$Length,name="tpm",slot="count",set.to.default=TRUE) {
  genes=ToIndex(data,genes)
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  data=AddSlot(data,name,comp.tpm(mat,tlen,subset = genes),set.to.default=set.to.default)
  data
}

NormalizeBaseline=function(data,baseline=FindReferences(data,reference=Condition==levels(Condition)[1],as.list=TRUE),name="baseline",slot=DefaultSlot(data),set.to.default=FALSE,LFC.fun=PsiLFC,...) {
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  mat=sapply(names(baseline),function(n) LFC.fun(mat[,n],rowMeans(mat[,baseline[[n]]]),...))
  data=AddSlot(data,name,mat,set.to.default=set.to.default)
  data
}


Scale=function(data,name="scaled",slot=DefaultSlot(data),set.to.default=FALSE,group=NULL,...) {
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  if (is.null(group)) {
    mat=t(scale(t(mat)))
  } else {
    gr=Coldata(data)[[group]]
    for (c in unique(gr)) mat[,gr==c]=t(scale(t(mat[,gr==c])))
  }
  data=AddSlot(data,name,mat,set.to.default=set.to.default)
  data
}



# minsamp: GetSummarizeMatrix(d,subset=NULL,average=FALSE), if all replicates of a sample exceed minval
FilterGenes=function(data,type='count',minval=100,mincol=ncol(data)/2,min.cond=NULL,use=NULL,keep=NULL,return.genes=FALSE) {
  if (!is.null(use) & !is.null(keep)) stop("Do not specify both use and keep!")

  if (is.null(use)) {
    t=GetTable(data,type=type)
    if (!is.null(min.cond)) {
      m=GetSummarizeMatrix(data,columns=NULL,average=FALSE)
      use=rowSums(sapply(1:ncol(m),function(i) apply(t[,m[,i]>0]>=minval,1,all)))>=min.cond
    } else {
      use=apply(t,1,function(v) sum(v>=minval,na.rm=TRUE)>=mincol)
    }
    if (!is.null(keep)) use = use | rownames(t) %in% rownames(t[keep,])
  } else {
    use=ToIndex(data,use)
  }

  if (return.genes) return(data$gene.info$Gene[use])
  return(data.apply(data,function(t) t[use,],fun.gene.info = function(t) t[use,]))
}


ComputeExpressionPercentage=function(data,name,genes,mode.slot=DefaultSlot(data),multiply.by.100=TRUE) {
  gof=colSums(GetTable(data,type=mode.slot,ntr.na = FALSE,genes = genes))
  total=colSums(GetTable(data,type=mode.slot,ntr.na = FALSE))
  percentage=gof/total
  if (multiply.by.100) percentage=percentage*100
  Coldata(data,name)=percentage
  data
}
