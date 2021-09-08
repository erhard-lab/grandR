

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

DropDiffExp=function(data) {
  data$diffexp=NULL
  invisible(data)
}
AddDiffExp=function(data,name,mode,table) {
  if (is.null(data$diffexp)) data$diffexp=list()
  if (is.null(data$diffexp[[name]])) data$diffexp[[name]]=list()
  if (is.null(data$diffexp[[name]][[mode]])) {
    data$diffexp[[name]][[mode]]=table
  } else {
    for (n in names(table)) data$diffexp[[name]][[mode]][[n]]=table[[n]]
  }
  invisible(data)
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

	data=AddDiffExp(data,name,"Total",
	                data.frame(
	                  M=res.tot$baseMean,
	                  S=res.tot$stat,
	                  P=res.tot$pvalue,
	                  Q=p.adjust(res.tot$pvalue,method="BH"))
	)
	data=AddDiffExp(data,name,"New",
	                data.frame(
	                  M=res.new$baseMean,
	                  S=res.new$stat,
	                  P=res.new$pvalue,
	                  Q=p.adjust(res.new$pvalue,method="BH"))
	)
	data=AddDiffExp(data,name,"Old",
	                data.frame(
	                  M=res.old$baseMean,
	                  S=res.old$stat,
	                  P=res.old$pvalue,
	                  Q=p.adjust(res.old$pvalue,method="BH"))
	)
	
	invisible(data)
}

LFC=function(data,contrasts,LFC.fun=PsiLFC,...) {
  comp=function(type,name) {
    mat=as.matrix(GetData(data,type=type,table=TRUE,keep.ntr.na=FALSE))
    l=setNames(lapply(contrasts,function(ctr) {
      LFC.fun(rowSums(mat[,ctr==1,drop=FALSE]),rowSums(mat[,ctr==-1,drop=FALSE]),...)
    }),colnames(contrasts))
    for (n in names(l)) data=AddDiffExp(data,n,name,data.frame(LFC=l[[n]]))
    data
  }
  data=comp("total.count","Total")
  data=comp("new.count","New")
  data=comp("old.count","Old")
  invisible(data)
}


GetDiffExpTable=function(data,cols=NA) {
  if (!is.null(data$data$diffexp) && is.null(data$diffexp)) data$diffexp=data$data$diffexp
  re=data$gene.info[,!(names(data$gene.info) %in% c("Length")),drop=FALSE]
  ord=rep(0,nrow(re))
  for (name in names(data$diffexp)) {
    for (mode in names(data$diffexp[[name]])) {
      t=data$diffexp[[name]][[mode]]
      if (!is.na(cols)) t=t[,intersect(names(t),cols)]
      names(t)=paste(name,mode,names(t),sep=".")
      re=cbind(re,t)
      if (!is.null(t$q)) ord=re+t$q
    }
  }
	re[order(ord),]
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
GetContrasts.default=function(names=NULL,design=NULL,coldata=MakeColData(names,design),contrast,covariate=NULL,name="LFC",subset=NULL) {
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
  if (is.null(covariate))  return(setNames(data.frame(contr()),name))
  covvec=interaction(coldata[covariate],drop=FALSE,sep=".")
  re=as.data.frame(setNames(lapply(levels(covvec),function(bin) contr(use=covvec==bin)),levels(covvec)),check.names=FALSE)
  re=re[,!apply(re==0,2,all),drop=FALSE]
  re
}

VulcanoPlot=function(data,name=names(data$diffexp)[1],mode="Total",aest=aes(),p.cutoff=0.05,lfc.cutoff=1,label.numbers=TRUE) {
  df=data$diffexp[[name]][[mode]]
  aes=modifyList(aes(LFC,-log10(Q),color=density2d(LFC,-log10(Q))),aest)
  g=ggplot(df,mapping=aes)+
    geom_point(size=0.1)+
    scale_color_viridis_c(guide=FALSE)+
    xlab(bquote(log[2]~FC))+
    ylab(bquote("-"~log[10]~FDR))+
    geom_hline(yintercept=-log10(p.cutoff),linetype=2)+
    geom_vline(xintercept=c(-lfc.cutoff,lfc.cutoff),linetype=2)+
    ggtitle(paste0(name," (",mode,")"))
  if (label.numbers) {
    n=table(cut(df$LFC,breaks=c(-Inf,-lfc.cutoff,lfc.cutoff,Inf)),df$Q>p.cutoff)
    g=g+annotate("label",x=c(-Inf,0,Inf,-Inf,0,Inf),y=c(Inf,Inf,Inf,-Inf,-Inf,-Inf),label=paste0("n=",as.numeric(n)),hjust=c(-0.1,0.5,1.1,-0.1,0.5,1.1),vjust=c(1.1,1.1,1.1,-0.1,-0.1,-0.1))
  }
  g
}



MAPlot=function(data,name=names(data$diffexp)[1],mode="Total",aest=aes(),p.cutoff=0.05,lfc.cutoff=1) {
  df=data$diffexp[[name]][[mode]]
  aes=modifyList(aes(M+1,LFC,color=ifelse(Q<p.cutoff,"Sig.","NS")),aest)
  g=ggplot(df,mapping=aes)+
    geom_point(size=0.5)+
    scale_x_log10()+
    scale_color_manual(values=c(Sig.="black",NS="grey30"),guide=FALSE)+
    ylab(bquote(log[2]~FC))+
    xlab("Total expression")+
    geom_hline(yintercept=c(-lfc.cutoff,lfc.cutoff),linetype=2)+
    ggtitle(paste0(name," (",mode,")"))
  g
}
