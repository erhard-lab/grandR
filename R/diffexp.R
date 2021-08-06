

CreateDESeq2<-function(data,Total=TRUE,New=FALSE,Old=FALSE, conditions=colnames(data$data$count), formula=~Sample*Type) {
	cbindcheck=function(a,b) if (is.null(a)) b else cbind(a,b)
	coldata=NULL
	cnt=NULL
	if (Total) {
		coldata=rbind(coldata,cbind(data$coldata[conditions,],data.frame(Type="Total")))
		cnt=cbindcheck(cnt,GetData(data,"count",table = TRUE,conditions=conditions))
	}
	if (New) {
		coldata=rbind(coldata,cbind(data$coldata[conditions,],data.frame(Type="New")))
		cnt=cbindcheck(cnt,GetData(data,"new.count",table = TRUE,conditions=conditions))
	}
	if (Old) {
		coldata=rbind(coldata,cbind(data$coldata[conditions,],data.frame(Type="Old")))
		cnt=cbindcheck(cnt,GetData(data,"old.count",table = TRUE,conditions=conditions))
	}
	allzero=apply(cnt,2,function(v) all(v==0))
	cnt=cnt[,!allzero]
	coldata=droplevels(coldata[!allzero,])
	for (i in 1:ncol(coldata)) {
		if (is.factor(coldata[,i])) levels(coldata[,i]) = gsub("+","_P_",gsub("-","_M_",levels(coldata[,i]),fixed=TRUE),fixed=TRUE)
	}
	DESeq(DESeqDataSetFromMatrix(countData = cnt ,colData = coldata,design = formula))
}

CreateDESeq2NewOld<-function(data,...) CreateDESeq2(data,Total=FALSE,New=TRUE,Old=TRUE,...)
CreateDESeq2New<-function(data,...) CreateDESeq2(data,Total=FALSE,New=TRUE,Old=FALSE,...)
CreateDESeq2Old<-function(data,...) CreateDESeq2(data,Total=FALSE,New=FALSE,Old=TRUE,...)



cnt=function(m) {
	mode(m) <- "integer"
	m
}


TestGenesLRT=function(data,target=~Sample,background=~1,name="lrt",verbose=FALSE,subset=!data$coldata$no4sU) {

	
	colData=droplevels(data$coldata[subset,])
	countData=data$data$count[,subset]
	ntrData=data$data$ntr[,subset]
	
	if (verbose) cat("Testing total...\n")
	dds.tot <- DESeq(DESeqDataSetFromMatrix(countData = cnt(countData),colData=colData,design = target),test="LRT", reduced=background)
	if (verbose) cat("Testing new...\n")
	dds.new <- DESeq(DESeqDataSetFromMatrix(countData = cnt(countData*ntrData),colData=colData,design = target),test="LRT", reduced=background)
	if (verbose) cat("Testing old...\n")
	dds.old <- DESeq(DESeqDataSetFromMatrix(countData = cnt(countData*(1-ntrData)),colData=colData,design = target),test="LRT", reduced=background)

	res.tot <- DESeq2::results(dds.tot)
	res.new <- DESeq2::results(dds.new)
	res.old <- DESeq2::results(dds.old)

	
	df=data.frame(
		total.p=res.tot$pvalue,
		total.q=p.adjust(res.tot$pvalue,method="BH"),
		new.p=res.new$pvalue,
		new.q=p.adjust(res.new$pvalue,method="BH"),
		old.p=res.old$pvalue,
		old.q=p.adjust(res.old$pvalue,method="BH")
	)
	names(df)=paste0(name,".",names(df))

	if (is.null(data$data$diffexp)) data$data$diffexp = df else data$data$diffexp=cbind(data$data$diffexp[,!(names(data$data$diffexp) %in% names(df))],df)

	invisible(data)
}

GetFdrTab=function(data) {
	d=data$data$diffexp[,grepl(".q$",names(data$data$diffexp))]
	ord=order(apply(d,1,function(v) sum(log(v+1E-8),na.rm=T)),decreasing=FALSE)
	cbind(data$gene.info[,!(names(data$gene.info) %in% c("Length"))],d)[ord,]
}

GetSummarizeMatrix <- function (x, ...) {
  UseMethod("GetSummarizeMatrix", x)
}
GetSummarizeMatrix.grandR=function(data,...) GetSummarizeMatrix.default(data$coldata,...)
GetSummarizeMatrix.default=function(coldata,column="Sample",subset=!coldata$no4sU,average=TRUE) {
	re=NULL
	for (v in unique(coldata[,column])) re=cbind(re,ifelse(coldata[,column]==v,1,0))
	rownames(re)=rownames(coldata)
	colnames(re)=unique(coldata[,column])
	re[-which(subset),]=0
	re=re[,colSums(re)>0]
	if (average) re=t(t(re)/colSums(re))
	re
}

GetContrasts <- function (x, ...) {
  UseMethod("GetContrasts", x)
}
GetContrasts.grandR=function(data,subset=!data$coldata$no4sU,...) GetContrasts.default(coldata=data$coldata,subset=subset,...)
GetContrasts.default=function(names=NULL,design=NULL,coldata=MakeColData(names,design),contrast,covariate=NULL,subset=NULL) {
  if (length(contrast)==1) contrast=c(contrast,levels(coldata[,contrast])[1:2])
  if (length(contrast)!=3 || !contrast[1]%in%names(coldata) || !all(contrast[2:3] %in% coldata[,contrast[1]])) stop("Illegal contrasts (either a name from design, or a vector of length 3 (name from design vector and two levels)")
  contr=function(use=TRUE) {
    re=rep(0,nrow(coldata))
    re[coldata[,contrast[1]]==contrast[2] & use]=1
    re[coldata[,contrast[1]]==contrast[3] & use]=-1
    if (!is.null(subset)) {
      save=re
      re=rep(0,nrow(coldata))
      re[subset]=save[subset]
    }
    re
  }
  if (is.null(covariate))  return(contr())
  covvec=interaction(coldata[covariate],drop=FALSE,sep=".")
  as.data.frame(setNames(lapply(levels(covvec),function(bin) contr(use=covvec==bin)),levels(covvec)),check.names=FALSE)
}

LFC=function(data,type="total.count",contrasts,name.suffix=paste0(".",type),LFC.fun=PsiLFC,...) {
  mat=as.matrix(GetData(data,type=type,table=TRUE,keep.ntr.na=FALSE))
  re=as.data.frame(setNames(lapply(contrasts,function(ctr) {
    LFC.fun(rowSums(mat[,ctr==1]),rowSums(mat[,ctr==-1]),...)
  }),colnames(contrasts)),check.names=FALSE)
  names(re)=paste0(names(re),name.suffix)
  re
}

