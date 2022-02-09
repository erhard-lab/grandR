

cnt=function(m) {
	mode(m) <- "integer"
	m
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

PairwiseRegulation=function(data,name,contrasts,slot=DefaultSlot(data),time=Design$dur.4sU,steady.state.columns,N=10000,seed=NULL) {

  if(!is.null(seed)) set.seed(seed)

  alpha=as.matrix(GetTable(data,type="alpha"))
  beta=as.matrix(GetTable(data,type="beta"))

  ss=as.matrix(GetTable(data,type="count",columns=colnames(data)[steady.state.columns]))
  ss.norm=as.matrix(GetTable(data,type=slot,columns=colnames(data)[steady.state.columns]))

  ApplyContrasts(data,name=name,contrasts=contrasts,mode.slot=slot,analysis="Regulation",FUN=function(mat,A,B) {
    count.A=mat[,A,drop=FALSE]
    count.B=mat[,B,drop=FALSE]

    alpha.A=alpha[,A,drop=FALSE]
    alpha.B=alpha[,B,drop=FALSE]
    beta.A=beta[,A,drop=FALSE]
    beta.B=beta[,B,drop=FALSE]

    t=unique(c(Coldata(data)[[time]][A],Coldata(data)[[time]][B]))

    if (ncol(ss)==1) {
      #cannot estimate dispersion, set them to (very conservative) 1
      disp=rep(1,ncol(count.A))
    } else {
      # estimate dispersions and size factors ( must be done even if normalized values were already computed, no way to tell whether any slot is valid data for DESeq2)
      dds=DESeq2::DESeqDataSetFromMatrix(countData = cnt(ss),colData=data.frame(rep(1,ncol(ss))),design = ~1)
      dds=DESeq2::estimateSizeFactors(dds)
      dds=DESeq2::estimateDispersions(dds,quiet=TRUE)
      disp=DESeq2::dispersions(dds)
      disp[is.na(disp)]=1 # ss=0 0 0 ...
    }

    re=plapply(1:nrow(count.A),function(i) {
      # for each column: sample from the posterior of s and d:
    #for (i in 1:nrow(count.A)) {
      f0=rnbinom(N,size=1/disp[i],mu=mean(ss.norm[i,]))
      sample.A=lapply(1:ncol(alpha.A), function(c){
        if (count.A[i,c]==0) return(cbind(s=rep(1,N),d=rep(0,N)))
        ntr.A=rbeta(N,alpha.A[i,c],beta.A[i,c])
        F.A=count.A[i,c]*(1-ntr.A)/f0
        d=-1/t*log(F.A)
        s=-1/t*count.A[i,c]*ntr.A * ifelse(F.A>=1,-1,ifelse(is.infinite(F.A),0,log(F.A)/(1-F.A)))
        cbind(s,d)
      })
      sample.B=lapply(1:ncol(alpha.B), function(c){
        if (count.B[i,c]==0) return(cbind(s=rep(1,N),d=rep(0,N)))
        ntr.B=rbeta(N,alpha.B[i,c],beta.B[i,c])
        F.B=count.B[i,c]*(1-ntr.B)/f0
        d=-1/t*log(F.B)
        s=-1/t*count.B[i,c]*ntr.B * ifelse(F.B>=1,-1,ifelse(is.infinite(F.B),0,log(F.B)/(1-F.B)))
        cbind(s,d)
      })
      s.A=do.call("c",lapply(sample.A,function(m) m[,'s']))
      d.A=do.call("c",lapply(sample.A,function(m) m[,'d']))
      s.B=do.call("c",lapply(sample.B,function(m) m[,'s']))
      d.B=do.call("c",lapply(sample.B,function(m) m[,'d']))

      mean.and.sig=function(a,b) {
        use=a>0&b>0
        all=log2(a[use]/b[use])
        LFC=mean(c(all,rep(0,sum(!use))))
        c(LFC,max(0,sum(sign(all)==sign(LFC))/length(a)*2-1))
      }

      setNames(c(mean.and.sig(s.A,s.B),mean.and.sig(d.B,d.A)),c("LFC.s","Prob.s","LFC.HL","Prob.HL"))  # that's right d.B/d.A is HL.A/HL.B
    },seed=seed)

    as.data.frame(t(simplify2array(re)))
  })


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


ApplyContrasts=function(data,name,contrasts,mode.slot="count",FUN,analysis=attr(FUN,"analysis"),...) {
  mat=as.matrix(GetTable(data,type=mode.slot,ntr.na=FALSE))
  mode.slot=get.mode.slot(data,mode.slot)
  for (n in names(contrasts)) {
    re.df=FUN(mat,contrasts[[n]]==1,contrasts[[n]]==-1,...)
    data=AddAnalysis(data,description=MakeAnalysis(name = paste0(name,".",n),analysis = analysis,mode = mode.slot$mode,slot=mode.slot$slot,columns = colnames(mat)[contrasts[[n]]==1|contrasts[[n]]==-1]),table = re.df)
  }
  data
}

LFC=function(data,name,contrasts,LFC.fun=PsiLFC,mode.slot="count",...) {
  ApplyContrasts(data,name=name,contrasts=contrasts,mode.slot=mode.slot,analysis="LFC",FUN=function(mat,A,B) {
    lfcs=LFC.fun(rowSums(mat[,A,drop=FALSE]),rowSums(mat[,B,drop=FALSE]),...)
    if (is.data.frame(lfcs)) lfcs else data.frame(LFC=lfcs)
  })
}

#' Create a summarize matrix
#'
#' If this matrix is multiplied with a count table (e.g. obtained by \code{\link{GetTable}}),
#' either the average (average=TRUE) or the sum (average=FALSE) of all columns (samples or cells)
#' belonging to the same \code{\link{Condition}} is computed.
#'
#' @param data A grandR object
#' @param no4sU Use no4sU columns (TRUE) or not (FALSE)
#' @param columns logical vector of which columns (samples or cells) to use (or NULL: use all)
#' @param average matrix to compute the average (TRUE) or the sum (FALSE)
#' @param v a named vector (the names indicate the sample names, the value the conditions to be summarized)
#' @param subset logical vector of which elements of the vector v to use (or NULL: use all)
#'
#' @return A matrix to be multiplied with a count table
#'
#' @details The method for grandR object simply calls the general method
#'
#' @seealso \link{GetTable}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#'
#' GetSummarizeMatrix(sars)
#' head(as.matrix(GetTable(sars)) %*% GetSummarizeMatrix(sars))   # average by matrix multiplication
#' head(GetTable(sars,summarize = TRUE))                          # shortcut, does the same
#'
#' GetSummarizeMatrix(c(A=1,B=1,C=2,D=3))
#'
#' @export
#'
GetSummarizeMatrix <- function (x, ...) {
  UseMethod("GetSummarizeMatrix", x)
}
#' @rdname GetSummarizeMatrix
#' @export
GetSummarizeMatrix.grandR=function(data,no4sU=FALSE,columns=NULL,average=TRUE) {
  if (is.null(Condition(data))) stop("Does not have conditions!")
  columns=if (is.null(columns)) no4sU | !Coldata(data)$no4sU else columns&(no4sU | !Coldata(data)$no4sU)
  GetSummarizeMatrix.default(setNames(Condition(data),colnames(data)),subset=columns,average=average)
}
#' @rdname GetSummarizeMatrix
#' @export
GetSummarizeMatrix.default=function(v,subset=NULL,average=TRUE) {
	re=NULL
	for (e in unique(v)) re=cbind(re,ifelse(v==e,1,0))
	rownames(re)=names(v)
	colnames(re)=unique(v)
	if (!is.null(subset)) {
  	save=re[subset,]
  	re[,]=0
  	re[subset,]=save
	}
	re=re[,colSums(re)>0]
	if (average) re=t(t(re)/colSums(re))
	re
}

#' Create a contrast matrix
#'
#' Each column of a contrast matrix represents a pairwise comparison of all samples or cells of
#' a grandR object (or a column annotation table). Elements being 1 are contrasted vs. elements being -1
#' (and all 0 are irrelevant for this comparison).
#'
#' @param data A grandR object
#' @param coldata A column annotation table
#' @param contrast A vector describing what should be contrasted
#' @param no4sU Use no4sU columns (TRUE) or not (FALSE)
#' @param columns logical vector of which columns (samples or cells) to use (or NULL: use all)
#' @param group Split the samples or cells according to this column of the column annotation table (and adapt the of the output table)
#' @param name.format Format string for generating the column from the contrast vector (see details)
#'
#' @return A contrast matrix to be used in \code{\link{ApplyContrasts}}, \code{\link{LFC}}, \code{\link{TestPairwise}}
#'
#' @details To compare one specific factor level \emph{A} against another level \emph{B} in
#' a particular column \emph{COL} of the column annotation table, specify contrast=c("COL","A","B")
#'
#' @details To compare all levels against a specific level \emph{A} in
#' a particular column \emph{COL} of the column annotation table, specify contrast=c("COL","A")
#'
#' @details To perform all pairwise comparisons of all levels from
#' a particular column \emph{COL} of the column annotation table, specify contrast=c("COL")
#'
#' @details If the column \emph{COL} only has two levels, all three are equivalent.
#'
#' @details In all cases, if groups is not NULL, the columns annotation table is first split and contrasts are applied within all samples or cells
#' with the same \emph{group} factor level.
#'
#' @details The format string specifies the column name in the generated contrast matrix (which is used as the \emph{Analysis} name when calling
#' \code{\link{ApplyContrasts}}, \code{\link{LFC}}, \code{\link{TestPairwise}}, etc.). The keywords \emph{$COL}, \emph{$A} and \emph{$B} are substituted
#' by the respective elements of the contrast vector.
#'
#' @details The method for grandR objects simply calls the general method
#'
#' @seealso \code{\link{ApplyContrasts}}, \code{\link{LFC}}, \code{\link{TestPairwise}}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition","Time",Design$Replicate))
#'
#' GetContrasts(sars,contrast="Condition")                    # Compare all Mock vs. all SARS
#' GetContrasts(sars,contrast=c("Condition","SARS","Mock"))   # This direction of the comparison is more reasonable
#' GetContrasts(sars,contrast=c("Condition","SARS","Mock"),group="Time")   # Compare SARS vs Mock per time point
#' GetContrasts(sars,contrast=c("Time","no4sU"), group="Condition",name.format="$COL ($A vs $B)")   # Compare each sample against the respective no4sU sample
#'
#' @export
#'
GetContrasts <- function (x, ...) {
  UseMethod("GetContrasts", x)
}
GetContrasts.grandR=function(data,contrast="Condition",no4sU=FALSE,columns=NULL,group=NULL,name.format="$A vs $B") {
  columns=if (is.null(columns)) no4sU | !Coldata(data)$no4sU else columns&(no4sU | !Coldata(data)$no4sU)
  GetContrasts.default(coldata=Coldata(data),contrast=contrast,group=group,columns=columns,name.format=name.format)
}
GetContrasts.default=function(coldata,contrast,columns=NULL,group=NULL,name.format="$A vs $B") {
  if (!(length(contrast) %in% 1:3) || !contrast[1]%in%names(coldata) || (length(contrast)>1 && !all(contrast[2:length(contrast)] %in% coldata[,contrast[1]]))) stop("Illegal contrasts (either a name from design (all pairwise comparisons), a name and a reference level (all comparisons vs. the reference), or a name and two levels (exactly this comparison))")

  make.name=function(contrast) gsub("$COL",contrast[1], gsub("$A",contrast[2], gsub("$B",contrast[3], name.format,fixed=TRUE),fixed=TRUE),fixed=TRUE)
  make.col=function(contrast,use=TRUE) {
    re=rep(0,nrow(coldata))
    re[coldata[,contrast[1]]==contrast[2] & use]=1
    re[coldata[,contrast[1]]==contrast[3] & use]=-1
    if (!is.null(columns)) {
      save=re
      re=rep(0,nrow(coldata))
      re[columns]=save[columns]
    }
    setNames(data.frame(re,check.names=FALSE),make.name(contrast))
  }

  contr=if (length(contrast)==3) {
    function(use=TRUE) make.col(contrast,use)
  } else if (length(contrast)==2) {
    function(use=TRUE) {
      if (is.null(columns)) columns=TRUE
      ll=if (is.factor(coldata[,contrast[1]])) levels(droplevels(coldata[columns,contrast[1]])) else unique(coldata[columns,contrast[1]])
      ll=setdiff(ll,contrast[2])
      as.data.frame(lapply(ll,function(l) make.col(c(contrast[1],l,contrast[2]),use)),check.names=FALSE)
    }
  } else {
    function(use=TRUE) {
      if (is.null(columns)) columns=TRUE
      ll=if (is.factor(coldata[,contrast[1]])) levels(droplevels(coldata[columns,contrast[1]])) else unique(coldata[columns,contrast[1]])
      re=combn(ll,2,FUN=function(v) make.col(c(contrast,v),use)[,1])
      colnames(re)=combn(ll,2,FUN=function(v) names(make.col(c(contrast,v),use))[1])
      as.data.frame(re,check.names=FALSE)
    }
  }
  re=if (is.null(group)) {
    contr()
  } else {
    covvec=interaction(coldata[group],drop=FALSE,sep=".")
    names=levels(covvec)
    as.data.frame(lapply(levels(covvec),function(bin) {
      r=contr(use=covvec==bin)
      setNames(r,paste0(names(r),".",bin))
    }),check.names=FALSE)
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
