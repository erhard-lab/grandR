
Kinetics=function(a0,HL,steady) KineticsRates(a0,steady*log(2)/HL,log(2)/HL)
KineticsRates=function(a0,s,d) function(t) (a0-s/d)*exp(-t*d)+s/d


density2d=function(x, y, facet=NULL, n=100) {
    if (is.null(facet)) {
        use=is.finite(x+y)
        d=MASS::kde2d(x[use], y[use], n=n)
        r=rep(NA,length(x))
        r[use]=d$z[cbind(findInterval(x[use], d$x),findInterval(y[use], d$y))]
	r=r/max(r,na.rm=T)
        return(r)
    }
    re<-rep(NA,length(x))
    for (f in unique(facet)) re[f==facet]=density2d(x[f==facet],y[f==facet],n=n)
    re
}



PlotPCA=function(data,type="count",ntop=500,aest=aes(color=Sample),x=1,y=2,subset=NULL) {
  
  my.cbind=function(...) {
    l=list(...)
    l=lapply(1:length(l),function(i) as.matrix(setNames(l[[i]],paste0(type[i],".",names(l[[i]])))))
    do.call("cbind",l)
  }
  
  mat=do.call("my.cbind",lapply(type,GetTable,data=data))
	#mat=cnt(switch(type[1],Total=data$data$count,New=data$data$count*data$data$ntr,Old=data$data$count*(1-data$data$ntr)))
	cd=do.call("rbind",lapply(1:length(type), function(i) cbind(data$coldata,data.frame(type=type[i]))))

	if (!is.null(subset)) {
		mat=mat[,subset]
		cd=cd[subset,]
	}
	rm.na=!apply(is.na(mat),2,sum)==nrow(mat)
	mat=mat[,rm.na]
	cd=cd[rm.na,]
	vsd <- DESeq2::vst(mat)

  	rv <- matrixStats::rowVars(vsd)
  	select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  	pca <- prcomp(t(vsd[select,]))
	percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
	d <- as.data.frame(pca$x)
	names(d)=paste0("PC",1:dim(d)[2])	
	d=cbind(d, cd)

	ggplot(d,modifyList(aes_string(paste0("PC",x),paste0("PC",y)),aest))+ geom_point(size=3)+xlab(paste0("PC",x,": ",round(percentVar[x] * 100),"% variance"))+ylab(paste0("PC",y,": ",round(percentVar[y] * 100),"% variance"))+coord_fixed()
}

Transform.Z=function(mat) {re=t(scale(t(mat))); attr(re,"label")="z score"; re}
Transform.VST=function(mat) {re=vst(mat); attr(re,"label")="VST"; re}
Transform.logFC=function(LFC.fun=PsiLFC,lfc.reference=NULL, contrasts=NULL,...) {
  function(m) {
    if (!is.null(contrasts)) {
      re=as.matrix(as.data.frame(setNames(lapply(contrasts,function(ctr) {
        LFC.fun(rowSums(m[,ctr==1,drop=FALSE]),rowSums(m[,ctr==-1,drop=FALSE]),...)
      }),colnames(contrasts))))
    }
    else {
      if (is.null(lfc.reference)) lfc.reference=1:ncol(m)
      ref=rowMeans(m[,lfc.reference,drop=F])
      re=apply(m,2,function(v) LFC.fun(v,ref))
    }
    attr(re,"label")="log FC"
    re
  }
}

PlotHeatmap=function(data,
                     genes=Genes(data),
                     type=c("Total","New","Old"),
                     slot=DefaultSlot(data),
                     summarize=NULL,
                     transform=Transform.Z,
                     subset=NULL,
                     cluster.genes=NULL,
                     label.genes=length(genes)<=50,
                     breaks=NULL, # either actual breaks, or null (50 and 99 quantile), or the quantiles (number matching to colors); if values centere around 0, then the center color is always 0
                     colors=RColorBrewer::brewer.pal(5,"RdBu"),
                     verbose=FALSE,
                     sample.names=colnames(data),
                     title=NULL,return.matrix=FALSE,...) {
  
  if (length(genes)==1 && genes %in% names(data$diffexp) && "Q" %in% names(data$diffexp[[genes]][[type[1]]])) {
    n=genes
    genes=data$diffexp[[genes]][[type[1]]]$Q<0.05
    if (verbose) cat(sprintf("Selected %d genes significant in %s\n",sum(genes),n))
    if (is.null(cluster.genes)) cluster.genes=TRUE
  }
  
  if (is.null(cluster.genes)) cluster.genes=length(genes)==length(Genes(data)) && all(genes==Genes(data))
  
  if (tolower(type[1])=="new") slot=paste0("new.",slot)
  if (tolower(type[1])=="old") slot=paste0("old.",slot)
  if (tolower(type[1])!="total" && is.null(subset)) subset=!data$coldata$no4sU
  
  mat=as.matrix(GetTable(data,type=slot,gene = genes,subset=subset,summarize = summarize,keep.ntr.na = FALSE))
  mat=transform(mat)

  name=attr(mat,"label")
  if (quantile(mat,0.25)<0) {
    nq=(length(colors)-1)/2
    if (is.null(breaks) || (length(breaks)%%2==1 && length(breaks)==nq)) {
      quant=if (is.null(breaks)) c(seq(0,1,length.out=nq+1)[c(-1,-nq-1)]*100,95) else breaks
      upper=quantile(mat[mat>0],quant/100)
      lower=quantile(-mat[mat<0],quant/100)
      breaks=c(-rev(pmax(upper,lower)),0,pmax(upper,lower))
    }
  } else {
    nq=length(colors)
    if (is.null(breaks) || length(breaks)==nq) {
      quant=if (is.null(breaks)) c(5,seq(0,1,length.out=nq)[c(-1,-nq)]*100,95) else breaks
      breaks=quantile(mat,quant/100)
    }
  }
  colnames(mat)=if (is.null(subset)) sample.names else sample.names[subset]
  col=circlize::colorRamp2(breaks = breaks,colors=colors)
  hm=ComplexHeatmap::Heatmap(mat,name=name,
                          cluster_rows = cluster.genes,
                          cluster_columns = FALSE,
                          show_row_names = label.genes,
                          col = col,
                          column_title=if (is.null(title)) "" else title,
                          row_title=sprintf("n=%d",nrow(mat)),
                          ...)  
  
  if (return.matrix) return(list(Matrix=mat,Heatmap=hm))
  hm
}

PlotTestOverlap=function(data,name="lrt",alpha=0.05,type=c("venn","euler")) {
  mat=GetDiffExpTable(f,gene.info=FALSE,name=name,cols='Q')
	df=setNames(as.data.frame(mat<alpha & !is.na(mat)),gsub(".Q$","",names(mat)))
	pl=switch(type[1],euler=eulerr::euler(df),venn=eulerr::venn(df))
	plot(pl,main=name)
}


PlotScatter=function(...)  {
	UseMethod('PlotScatter',list(...)[[1]])  
}

PlotScatter.default=function(...) {
  l=list(...)
  df=data.frame(l[[1]],l[[2]])
  names(df)=names(l)[1:2]
  l=c(list(df=df),l[-1:-2])
  do.call("PlotScatter.data.frame",l)
}

PlotScatter.grandR=function(data,cx,cy,type=DefaultSlot(data),type.x=type,type.y=type,...) {
  a=GetData(data,type.x,conditions=cx,table=T)[,1]
	b=GetData(data,type.y,conditions=cy,table=T)[,1]
	df=data.frame(a,b)
	names(df)=c(cx,cy)
	PlotScatter(df,...)
}

PlotScatter.data.frame=function(df,xcol=1,ycol=2,log=FALSE,log.x=log,log.y=log) {
  df=df[,c(xcol,ycol)]
	dfnames=names(df)
	names(df)=c("A","B")
	xmean=mean(df[,1])
	ymean=mean(df[,2])
	xlim=boxplot.stats(df[,1])$stats[c(1,5)] #(quantile(df[,1]-xmean,pnorm(c(-2,2)))*1.5)+xmean
	ylim=boxplot.stats(df[,2])$stats[c(1,5)] #(quantile(df[,2]-ymean,pnorm(c(-2,2)))*1.5)+ymean
	clip=function(v,lim) ifelse(v<lim[1],-Inf,ifelse(v>lim[2],Inf,v))
	df$A=clip(df$A,xlim)
	df$B=clip(df$B,ylim)
	
	trans.x=if(log.x) function(x) log(x) else function(x) x
	trans.y=if(log.y) function(y) log(y) else function(y) y
	
	g=ggplot(df,aes(A,B,color=density2d(trans.x(A), trans.y(B), n = 100)))+
			geom_point()+
			scale_color_viridis_c(name = "Density",guide=FALSE)+
			xlab(dfnames[1])+ylab(dfnames[2])+
			coord_cartesian(ylim=ylim,xlim=xlim)
	if (log.x) g=g+scale_x_log10()
	if (log.y) g=g+scale_y_log10()
	g
}

#PlotToxicityTest=function(data,w4sU,no4sU,ylim=c(-1,1),LFC.fun=PsiLFC,hl.quantile=0.8) {
#	w=GetData(data,"count",conditions=w4sU,table=T)[,1]
#	n=if (is.numeric(no4sU)) no4sU[data$gene.info$Gene] else GetData(data,"count",conditions=no4sU,table=T)[,1]
#	ntr=GetData(data,"ntr",conditions=w4sU,table=T)[,1]
#	use=!is.na(w+n+ntr)
#	w=w[use]
#	n=n[use]
#	ntr=ntr[use]
#
#	phl=comp.hl(ntr)	
#	df=data.frame(lfc=LFC.fun(w,n),PHL=phl)[ntr<1,]
#	df=df[df$PHL<quantile(df$PHL[is.finite(df$PHL)],hl.quantile),]
#	ggplot(df,aes(PHL,lfc,color=density2d(PHL, lfc, n = 100)))+
#			scale_color_viridis_c(name = "Density",guide=FALSE)+
#			geom_point(alpha=1)+
#			geom_hline(yintercept=0)+
##			geom_smooth(method="loess")+
#			xlab("RNA half-life")+ylab("log FC 4sU/no4sU")+
#			scale_x_continuous(breaks=c())+
#			coord_cartesian(ylim=ylim)
#}



PlotExpressionTest=function(data,w4sU,no4sU,ylim=c(-1,1),LFC.fun=PsiLFC,hl.quantile=0.8) {
	w=GetData(data,"count",conditions=w4sU,table=T)[,1]
	n=if (is.numeric(no4sU)) no4sU[data$gene.info$Gene] else GetData(data,"count",conditions=no4sU,table=T)[,1]
	use=!is.na(w+n)
	w=w[use]
	n=n[use]

	df=data.frame(lfc=LFC.fun(w,n),M=(log10(w+1)+log10(n+1))/2)
	ggplot(df,aes(M,lfc,color=density2d(M, lfc, n = 100)))+
			scale_color_viridis_c(name = "Density",guide=FALSE)+
			geom_point(alpha=1)+
			geom_hline(yintercept=0)+
#			geom_smooth(method="loess")+
			xlab("Mean expression")+ylab("log FC 4sU/no4sU")+
			coord_cartesian(ylim=ylim)
}

PlotTypeDistribution=function(data,type="tpm",relative=FALSE) {
	df=GetData(d,type=type,table=T)
	df=sapply(levels(data$gene.info$Type),function(type) colSums(df[ data$gene.info$Type==type,]))
	df=df[,colSums(df)>0]
	if (relative) {
		df=df/rowSums(df)*100
		type=sprintf("%s [%%]",type)
	}
	df=melt(df,varnames=c("Condition","Type"))
	ggplot(df,aes(Condition,value,fill=Type))+geom_bar(stat="Identity")+scale_fill_viridis_d()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab(type)+xlab(NULL)
}

PlotGeneOldVsNew=function(data,gene,slot=DefaultSlot(data),show.CI=FALSE,aest=aes(color=Condition,shape=Replicate)) {
  new=paste0("new.",slot)
  old=paste0("old.",slot)
  df=GetData(data,gene=gene,type=c(old,new),melt=F,coldata=T,keep.ntr.na = FALSE)
  g=ggplot(df,modifyList(aes_string(old,new),aest))+
    geom_point(size=2)+
    scale_x_log10("Old RNA")+
    scale_y_log10("New RNA")
  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df=cbind(df,GetData(data,gene=gene,type=c("lower","upper",slot),melt=F,coldata=T,keep.ntr.na = FALSE)[,c("lower","upper",slot)])
    g=g+geom_errorbar(data=df,mapping=aes_string(ymin=paste0("lower*",slot),ymax=paste0("upper*",slot)))
    g=g+geom_errorbarh(data=df,mapping=aes_string(xmin=paste0("(1-upper)*",slot),xmax=paste0("(1-lower)*",slot)))
  }
  g
}

PlotGeneTotalVsNtr=function(data,gene,slot=DefaultSlot(data),show.CI=FALSE,aest=aes(color=Condition,shape=Replicate)) {
  df=GetData(data,gene=gene,type=c("ntr",slot),melt=F,coldata=T,keep.ntr.na = FALSE)
  g=ggplot(df,modifyList(aes_string(slot,"ntr"),aest))+
    geom_point(size=2)+
    scale_x_log10("Total RNA")+
    scale_y_continuous("NTR")
  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df=cbind(df,GetData(data,gene=gene,type=c("lower","upper"),melt=F,coldata=T,keep.ntr.na = FALSE)[,c("lower","upper")])
    g=g+geom_errorbar(data=df,mapping=aes(ymin=lower,ymax=upper))
  }
  g
}

PlotGeneGroupsPoints=function(data,gene,group="Condition",slot=DefaultSlot(data),type="total",show.CI=FALSE,aest=aes(color=Condition,shape=Replicate)) {
  df=GetData(data,gene=gene,type=c(slot,"ntr"),melt=F,coldata=T,keep.ntr.na = FALSE)
  df$value=switch(type[1],total=df[[slot]],new=df[[slot]]*df[["ntr"]],old=df[[slot]]*(1-df[["ntr"]]),stop(paste0(type," unknown!")))
  g=ggplot(df,modifyList(aes_string(group,"value"),aest))+
    geom_point(size=2,position=if(show.CI) position_dodge(width=0.4) else "identity")+
    xlab(NULL)+
    scale_y_log10(paste0(toupper(substr(type,1,1)),substr(type,2,nchar(type))," RNA (",slot,")"))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df=cbind(df,GetData(data,gene=gene,type=c("lower","upper"),melt=F,coldata=T,keep.ntr.na = FALSE)[,c("lower","upper")])
    g=switch(type[1],
             total=g,
             new=g+geom_errorbar(data=df,mapping=aes_string(ymin=paste0("lower*",slot),ymax=paste0("upper*",slot)),width=0,position=position_dodge(width=0.4)),
             old=g+geom_errorbar(data=df,mapping=aes_string(ymin=paste0("(1-upper)*",slot),ymax=paste0("(1-lower)*",slot)),width=0,position=position_dodge(width=0.4)),
             stop(paste0(type," unknown!"))
    )
  }
  g
}

Plot=function(fun=NULL,...,gg=NULL) {
  function(data,gene) {
    if (is.null(fun)) return(NULL)
    re=fun(data=data,gene=gene,...)
    if (!is.null(gg)) re=re+gg
    re
  }
}


# Caching will work as long as it is the same data object that is used!
DPlot=function(FUN,...,height=7,width=7,add=NULL) {
  param=list(...)
  value=NULL
  function(data,...) {
    pp=list(...)
    if (length(pp)>0) {
      re=do.call(FUN,c(list(data),param,pp))
      if (!is.null(add)) re=re+add
      return(re)
    }
    
    if (is.null(value)) {
      value<<-do.call(FUN,c(list(data),param))
      if (!is.null(add))
        value<<-value+add
    }
    value
  }
}

