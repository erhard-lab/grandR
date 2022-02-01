

cnt=function(m) {
	mode(m) <- "integer"
	m
}

DropDiffExp=function(data,pattern=NULL) {
  if (is.null(pattern)) {
    data$diffexp=NULL
  } else {
    data$diffexp=data$diffexp[!grepl(pattern,names(data$diffexp))]
  }
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

TestGenesLRT=function(data,target=~Condition,background=~1,name="lrt",verbose=FALSE,columns=!data$coldata$no4sU,total=TRUE,new=TRUE,old=TRUE) {
	colData=droplevels(Coldata(data)[columns,])
	countData=data$data$count[,columns]
	ntrData=data$data$ntr[,columns]

	colData=as.data.frame(lapply(colData,function(c) {
	  if (is.factor(c)) levels(c)=make.names(levels(c),unique = TRUE)
	  c
	}))

	if (total) {
  	if (verbose) cat("Testing total...\n")
  	dds.tot <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(countData = cnt(countData),colData=colData,design = target),test="LRT", reduced=background)
  	res.tot <- DESeq2::results(dds.tot)
  	data=AddDiffExp(data,name,"Total",
  	                data.frame(
  	                  M=res.tot$baseMean,
  	                  S=res.tot$stat,
  	                  P=res.tot$pvalue,
  	                  Q=p.adjust(res.tot$pvalue,method="BH"))
  	)
	}

	if(new) {
  	if (verbose) cat("Testing new...\n")
  	dds.new <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(countData = cnt(countData*ntrData),colData=colData,design = target),test="LRT", reduced=background)
  	res.new <- DESeq2::results(dds.new)
  	data=AddDiffExp(data,name,"New",
  	                data.frame(
  	                  M=res.new$baseMean,
  	                  S=res.new$stat,
  	                  P=res.new$pvalue,
  	                  Q=p.adjust(res.new$pvalue,method="BH"))
  	)
	}


	if (old) {
  	if (verbose) cat("Testing old...\n")
  	dds.old <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(countData = cnt(countData*(1-ntrData)),colData=colData,design = target),test="LRT", reduced=background)
  	res.old <- DESeq2::results(dds.old)
  	data=AddDiffExp(data,name,"Old",
  	                data.frame(
  	                  M=res.old$baseMean,
  	                  S=res.old$stat,
  	                  P=res.old$pvalue,
  	                  Q=p.adjust(res.old$pvalue,method="BH"))
  	)
	}

	invisible(data)
}

TestPairwise=function(data,contrasts,total=TRUE,new=TRUE,old=TRUE) {
  comp=function(type,name) {
    mat=as.matrix(GetTable(data,type=type,ntr.na=FALSE))

    l=setNames(lapply(contrasts,function(ctr) {
      B=mat[,ctr==1,drop=FALSE]
      A=mat[,ctr==-1,drop=FALSE] # reverse compared to LFC! (does not matter in the end, only the q val is extracted)
      coldata=data.frame(comparison=c(rep("A",ncol(A)),rep("B",ncol(B))))
      dds <- DESeqDataSetFromMatrix(countData = cbind(A,B),
                                    colData = coldata,
                                    design= ~ comparison-1)
      results(DESeq(dds))
    }),colnames(contrasts))
    for (n in names(l)) data=AddDiffExp(data,n,name,data.frame(
      M=l[[n]]$baseMean,
      S=l[[n]]$stat,
      P=l[[n]]$pvalue,
      Q=l[[n]]$padj
    ))
    data
  }
  if (total) data=comp("total.count","Total")
  if (new) data=comp("new.count","New")
  if (old) data=comp("old.count","Old")
  invisible(data)
}

LFC=function(data,contrasts,LFC.fun=PsiLFC,slot="count",total=TRUE,new=TRUE,old=TRUE,compute.test=FALSE,...) {
  comp=function(type,name) {
    mat=as.matrix(GetTable(data,type=type,ntr.na=FALSE))
    l=setNames(lapply(contrasts,function(ctr) {
      LFC.fun(rowSums(mat[,ctr==1,drop=FALSE]),rowSums(mat[,ctr==-1,drop=FALSE]),...)
    }),colnames(contrasts))
    for (n in names(l)) data=AddDiffExp(data,n,name,data.frame(LFC=l[[n]]))
    data
  }
  if (total) data=comp(paste0("total.",slot),"Total")
  if (new) data=comp(paste0("new.",slot),"New")
  if (old) data=comp(paste0("old.",slot),"Old")
  if (compute.test) data = TestPairwise(data=data,contrasts=contrasts,total=total,new=new,old=old)
  invisible(data)
}


GetDiffExpTable=function(data,gene.info=TRUE,rownames="Symbol",names=NULL,modes=NULL,cols=NULL,sort=FALSE) {
  if (!is.null(data$data$diffexp) && is.null(data$diffexp)) {
    data$diffexp=data$data$diffexp
    if (!is.null(cols)) cols=c(cols,tolower(cols))
  }
  re=data$gene.info
  if (!is.null(rownames)) {
    rownames(re)=if (rownames %in% names(data$gene.info)) data$gene.info[[rownames]] else data$gene.info[,1]
  }
  ord=rep(0,nrow(re))
  sintersect=function(a,b) if (is.null(b)) a else intersect(a,b)
  for (name in sintersect(names(data$diffexp),names)) {
    for (mode in sintersect(names(data$diffexp[[name]]),modes)) {
      t=data$diffexp[[name]][[mode]]
      if (!is.null(cols)) t=t[,intersect(names(t),cols),drop=FALSE]
      if (!is.null(t$Q)) ord=ord+log(t$Q)
      names(t)=paste(name,mode,names(t),sep=".")
      re=cbind(re,t)
    }
  }
  if (sort) re=re[order(ord),]
  if (is.logical(gene.info) && !gene.info) re=re[,(ncol(data$gene.info)+1):ncol(re)]
  if (is.character(gene.info)) re=re[,-which(!names(data$gene.info) %in% gene.info)]
  re
}

GetSummarizeMatrix <- function (x, ...) {
  UseMethod("GetSummarizeMatrix", x)
}
GetSummarizeMatrix.grandR=function(data,...) GetSummarizeMatrix.default(data$coldata,...)
GetSummarizeMatrix.default=function(coldata,group=Design$Condition,columns=!coldata$no4sU,average=TRUE) {
	re=NULL
	for (v in unique(coldata[,group])) re=cbind(re,ifelse(coldata[,group]==v,1,0))
	rownames(re)=rownames(coldata)
	colnames(re)=unique(coldata[,group])
	if (!is.null(columns)) {
  	save=re[columns,]
  	re[,]=0
  	re[columns,]=save
	}
	re=re[,colSums(re)>0]
	if (average) re=t(t(re)/colSums(re))
	re
}

GetContrasts <- function (x, ...) {
  UseMethod("GetContrasts", x)
}
GetContrasts.grandR=function(data,columns=NULL,...) GetContrasts.default(coldata=data$coldata,columns=columns,...)
GetContrasts.default=function(names=NULL,design=NULL,coldata=MakeColdata(names,design),contrast,name=NULL,group=NULL,columns=NULL) {
  # either level 1 against level 2 of some coldata column (3 entries in contrast)
  # or each level against one specific level (2 entries)
  # or each all pairwise comparison (1 entry)
  # in both cases, groups would just subset by the group levels first, and then do it for each subset separately
  if (!(length(contrast) %in% 1:3) || !contrast[1]%in%names(coldata) || (length(contrast)>1 && !all(contrast[2:length(contrast)] %in% coldata[,contrast[1]]))) stop("Illegal contrasts (either a name from design (all pairwise comparisons), a name and a reference level (all comparisons vs. the reference), or a name and two levels (exactly this comparison))")

  make.col=function(contrast,use=TRUE) {
    re=rep(0,nrow(coldata))
    re[coldata[,contrast[1]]==contrast[2] & use]=1
    re[coldata[,contrast[1]]==contrast[3] & use]=-1
    if (!is.null(columns)) {
      save=re
      re=rep(0,nrow(coldata))
      re[columns]=save[columns]
    }
    matrix(re)
  }

  contr=if (length(contrast)==3) {
    function(use=TRUE) make.col(contrast,use)
  } else if (length(contrast)==2) {
    function(use=TRUE) {
      if (is.null(columns)) columns=TRUE
      ll=if (is.factor(coldata[,contrast[1]])) levels(droplevels(coldata[columns,contrast[1]])) else unique(coldata[columns,contrast[1]])
      ll=setdiff(ll,contrast[2])
      re=sapply(ll,function(l) make.col(c(contrast[1],l,contrast[2]),use))
      #re=matrix(0,nrow=nrow(coldata),ncol=length(ll))
      #re[coldata[,contrast[1]]==contrast[2] & use,]=-1
      #for (c in 1:length(ll)) re[coldata[,contrast[1]]==ll[c] & use,c]=1
      #if (!is.null(columns)) {
      #  save=re
      #  re=matrix(0,nrow=nrow(coldata),ncol=length(ll))
      #  re[columns,]=save[columns,]
      #}
      colnames(re)=ll
      re
    }
  } else {
    function(use=TRUE) {
      if (is.null(columns)) columns=TRUE
      ll=if (is.factor(coldata[,contrast[1]])) levels(droplevels(coldata[columns,contrast[1]])) else unique(coldata[columns,contrast[1]])
      re=combn(ll,2,FUN=function(v) make.col(c(contrast,v),use)[,1])
      colnames(re)=combn(ll,2,FUN=paste,collapse=" vs ")
      re
    }
  }
  re=if (is.null(group)) {
    as.data.frame(contr())
  } else {
    covvec=interaction(coldata[group],drop=FALSE,sep=".")
    names=levels(covvec)
    re=lapply(levels(covvec),function(bin) contr(use=covvec==bin))
    re=if (!is.null(colnames(re))) setNames(re,paste0(names,names(re),sep=".")) else setNames(re,names)
    as.data.frame(re,check.names=FALSE)
  }
  if (ncol(re)==1 && names(re)=="V1") {
    names(re)=if (is.null(name)) contrast[1] else name;
    name=NULL;
  }
  if (!is.null(name)) {
    names(re)= if (length(name)==ncol(re)) name else paste(name,names(re),sep=".")
  }
  re=re[,!apply(re==0,2,all),drop=FALSE]
  rownames(re)=rownames(coldata)

  remove=apply(re>=0,2,all) | apply(re<=0,2,all)
  if (sum(remove)>0) {
    warning(sprintf("Removed columns without matching experiment: %s",paste(colnames(re)[remove],collapse = ",")))
    re=re[,!remove]
  }
  re
}

VulcanoPlot=function(data,name=names(data$diffexp)[1],mode="Total",aest=aes(),p.cutoff=0.05,lfc.cutoff=1,label.numbers=TRUE) {
  df=data$diffexp[[name]][[mode]]
  aes=modifyList(aes(LFC,-log10(Q),color=density2d(LFC,-log10(Q))),aest)
  g=ggplot(df,mapping=aes)+
    geom_point(size=0.1)+
    scale_color_viridis_c(guide='none')+
    xlab(bquote(log[2]~FC))+
    ylab(bquote("-"~log[10]~FDR))+
    geom_hline(yintercept=-log10(p.cutoff),linetype=2)+
    geom_vline(xintercept=c(-lfc.cutoff,lfc.cutoff),linetype=2)+
    ggtitle(paste0(name," (",mode,")"))
  if (label.numbers) {
    n=table(cut(df$LFC,breaks=c(-Inf,-lfc.cutoff,lfc.cutoff,Inf)),factor(df$Q>p.cutoff,levels=c("FALSE","TRUE")))
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





# Normalization of NTRs such that: median logFC new RNA vs. new RNA is 0, there is no correlation of this logFC vs the NTR
NormalizeEffectiveLabeling=function(data,reference=colnames(data),slot="norm",verbose=FALSE) {
  w=rowMeans(GetTable(data,type=slot,columns=reference),na.rm=TRUE)
  ntr=rowMeans(GetTable(data,type="ntr",columns=reference),na.rm=TRUE)
  w=w*ntr
  use=!is.na(w) && w>0
  w=w[use]

  for (s in colnames(data)[!data$coldata$no4sU])  {
    if (verbose) cat(sprintf("Fitting model for %s...\n",s))
    d.ntr=data$data$ntr[use,s]
    d=GetTable(data,type=slot,columns=s)[use,1]
    df=data.frame(ntr=ntr,lfc=log2(d.ntr*d/w))
    df=df[!is.na(df$lfc) & !is.infinite(df$lfc),]
    fit=quantreg::lprq(df$ntr,df$lfc,h = 0.05)
    fn=splinefun(fit$xx,fit$fv)
    data$data$ntr[,s]=pmax(pmin(1,data$data$ntr[,s]/2^fn(ntr)),0)
  }

  data
}
