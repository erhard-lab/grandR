
density2d=function(x, y, facet=NULL, n=100, margin='n') {
  bandwidth.nrd.ex=function (x)
  {
    r <- range(x)
    h <- (r[2L] - r[1L])/1.34
    4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
  }

  if (is.null(facet)) {
        use=is.finite(x+y)
        bw.x=MASS::bandwidth.nrd(x[use])
        if (bw.x==0) bw.x=bandwidth.nrd.ex(x[use])
        bw.y=MASS::bandwidth.nrd(y[use])
        if (bw.y==0) bw.y=bandwidth.nrd.ex(y[use])
        d=MASS::kde2d(x[use], y[use], h=c(bw.x,bw.y), n=n)
        if (margin=='x') d$z=d$z/apply(d$z,1,max)
        else if (margin=='y') d$z=t(t(d$z)/apply(d$z,2,max))
        else d$z=d$z/max(d$z)
        r=rep(NA,length(x))
        r[use]=d$z[cbind(findInterval(x[use], d$x),findInterval(y[use], d$y))]
	r=r/max(r,na.rm=T)
        return(r)
    }
    re<-rep(NA,length(x))
    for (f in unique(facet)) re[f==facet]=density2d(x[f==facet],y[f==facet],n=n)
    re
}



PlotPCA=function(data,type="count",ntop=500,aest=aes(color=Condition),x=1,y=2,columns=NULL) {

  my.cbind=function(...) {
    l=list(...)
    l=lapply(1:length(l),function(i) as.matrix(setNames(l[[i]],paste0(type[i],".",names(l[[i]])))))
    do.call("cbind",l)
  }

  mat=do.call("my.cbind",lapply(type,GetTable,data=data))
	#mat=cnt(switch(type[1],Total=data$data$count,New=data$data$count*data$data$ntr,Old=data$data$count*(1-data$data$ntr)))
	cd=do.call("rbind",lapply(1:length(type), function(i) cbind(data$coldata,data.frame(type=type[i]))))

	if (!is.null(columns)) {
		mat=mat[,columns]
		cd=cd[columns,]
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

Transform=function(name,label=NULL) function(mat) {
  re=get(paste0("Transform.",name))(mat)
  if (!is.null(label)) attr(re,"label")=label
  re
}
Transform.noop=function(mat) {re=mat; attr(re,"label")="Expression"; re}
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
                     columns=NULL,
                     cluster.genes=NULL,
                     label.genes=length(genes)<=50,
                     breaks=NULL, # either actual breaks, or null (50 and 99 quantile), or the quantiles (number matching to colors); if values centere around 0, then the center color is always 0
                     colors=RColorBrewer::brewer.pal(5,"RdBu"),
                     verbose=FALSE,
                     col.names=colnames(data),
                     title=NULL,return.matrix=FALSE,...) {

  if (length(genes)==1 && genes %in% Analyses(data) && "Q" %in% names(GetAnalysisTable(data,analyses=genes,regex=FALSE))) {
    n=genes
    genes=GetAnalysisTable(data,analyses=genes,columns="Q")$Q<0.05
    if (verbose) cat(sprintf("Selected %d genes significant in %s\n",sum(genes),n))
    if (is.null(cluster.genes)) cluster.genes=TRUE
  }

  if (is.null(cluster.genes)) cluster.genes=is.logical(genes) || (length(genes)==length(Genes(data)) && all(genes==Genes(data)))

  if (tolower(type[1])=="new") slot=paste0("new.",slot)
  if (tolower(type[1])=="old") slot=paste0("old.",slot)
  if (tolower(type[1])!="total" && is.null(columns)) columns=!data$coldata$no4sU

  mat=as.matrix(GetTable(data,type=slot,genes = genes,columns=columns,summarize = summarize,ntr.na = FALSE))
  mat=transform(mat)

  name=attr(mat,"label")
  if (quantile(mat,0.25)<0) {
    nq=(length(colors)-1)/2
    if (is.null(breaks) || (length(colors)%%2==1 && length(breaks)==nq)) {
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
  colnames(mat)=if (is.null(columns)) col.names else col.names[columns]
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

PlotTestOverlap=function(data,names=NULL,alpha=0.05,type=c("venn","euler")) {
  mat=GetAnalysisTable(data,gene.info=FALSE,genes=names,columns='^Q$')
	df=setNames(as.data.frame(mat<alpha & !is.na(mat)),gsub(".Q$","",names(mat)))
	pl=switch(type[1],euler=eulerr::euler(df),venn=eulerr::venn(df))
	plot(pl,main=name)
}


PlotScatter=function(df,xcol=1,ycol=2,x=NULL,y=NULL,log=FALSE,log.x=log,log.y=log,color=NA,remove.outlier=1.5,size=0.3,xlim=NULL,ylim=NULL, highlight=NULL, label=NULL, columns=NULL) {
  df=as.data.frame(df)
  if (!is.data.frame(df)) stop("df must be a data frame (or at least coercable into a data frame)")
  adaptInf=function(df,rx,ry) {
    # workaround to also "brush" infinity points at the border of the plane
    if (log.x) rx=log10(rx)
    rx=c(rx[1]-0.04*(rx[2]-rx[1]),rx[2]+0.04*(rx[2]-rx[1]))
    if (log.x) rx=10^rx
    df$A[is.infinite(df$A) & df$A>0]=rx[2]
    df$A[is.infinite(df$A) & df$A<0]=rx[1]

    if (log.y) ry=log10(ry)
    ry=c(ry[1]-0.04*(ry[2]-ry[1]),ry[2]+0.04*(ry[2]-ry[1]))
    if (log.y) ry=10^ry
    df$B[is.infinite(df$B) & df$B>0]=ry[2]
    df$B[is.infinite(df$B) & df$B<0]=ry[1]

    df
  }

  x=substitute(x)
  y=substitute(y)

  dfnames=c(xcol,ycol)
  rn=rownames(df)

  A=if (is.null(x)) df[,xcol] else {
    eval(x,df,parent.frame())
  }
  B=if (is.null(y)) df[,ycol] else {
    eval(y,df,parent.frame())
  }

  dfnames[1]=if (is.null(x)) names(df[,xcol,drop=F]) else deparse(x)
  dfnames[2]=if (is.null(y)) names(df[,ycol,drop=F]) else deparse(y)

  df$A=A
  df$B=B
  #df=data.frame(A=A,B=B)
  rownames(df)=rn
  if (!is.null(columns)) df=df[columns,]

  df$A.trans=if(log.x) log10(df$A) else df$A
  df$B.trans=if(log.y) log10(df$B) else df$B

  set.coord=remove.outlier!=FALSE || !is.null(xlim) || !is.null(ylim)
  if (set.coord) {
    if (is.null(xlim) && remove.outlier) {
      xlim=boxplot.stats(df$A.trans[is.finite(df$A.trans)],coef=remove.outlier)$stats[c(1,5)] #(quantile(df[,1]-xmean,pnorm(c(-2,2)))*1.5)+xmean
      xlim=c(df$A[which(df$A.trans==xlim[1])[1]],df$A[which(df$A.trans==xlim[2])[1]])
    }
    if (is.null(ylim) && remove.outlier) {
      ylim=boxplot.stats(df$B.trans[is.finite(df$B.trans)],coef=remove.outlier)$stats[c(1,5)] #(quantile(df[,2]-ymean,pnorm(c(-2,2)))*1.5)+ymean
      ylim=c(df$B[which(df$B.trans==ylim[1])[1]],df$B[which(df$B.trans==ylim[2])[1]])
    }
    if (is.null(xlim)) xlim=range(df$A[!is.infinite(df$A)])
    if (is.null(ylim)) ylim=range(df$B[!is.infinite(df$B)])
    clip=function(v,ch,lim,minus) ifelse(ch<lim[1],minus,ifelse(ch>lim[2],Inf,v))
    df$A.trans=clip(df$A.trans,df$A,xlim,-Inf)
    df$B.trans=clip(df$B.trans,df$B,ylim,-Inf)
    df$A=clip(df$A,df$A,xlim,if (log.x) 0 else -Inf)
    df$B=clip(df$B,df$B,ylim,if (log.y) 0 else -Inf)
  }
  else {
    xlim=range(df$A[!is.infinite(df$A)])
    ylim=range(df$B[!is.infinite(df$B)])
  }
  if (is.na(color)) {
    if (nrow(df)>1000) {
      df$color=density2d(df$A.trans, df$B.trans, n = 100)
      colorscale=scale_color_viridis_c(name = "Density",guide="none")
    } else {
      df$color=rep(1,nrow(df))
      colorscale=scale_color_manual(values="black",guide="none")
    }
  } else {
    df$color=color
    cmap=setNames(unique(color),unique(color))
    colorscale=scale_color_manual(name = NULL,values=cmap)
  }
  g=ggplot(df,aes(A,B,color=color))+
    geom_point(size=size)+
    colorscale+
    xlab(dfnames[1])+ylab(dfnames[2])
  if (!is.null(highlight)) {
    if (is.list(highlight)){
      for (col in names(highlight)) {
        g=g+geom_point(data=df[highlight[[col]],],color=col,size=size*2)
      }
    } else {
      g=g+geom_point(data=df[highlight,],color='red',size=size*2)
    }
  }
  if (!is.null(label)) {
    df2=df
    df2$label=""
    df2[label,"label"]=rownames(df2)[label]
    g=g+ggrepel::geom_label_repel(data=df2,mapping=aes(label=label),show.legend = FALSE)
  }
  if (set.coord) g=g+coord_cartesian(ylim=ylim,xlim=xlim)
  if (log.x) g=g+scale_x_log10()
  if (log.y) g=g+scale_y_log10()

  attr(g, 'df') <- adaptInf(df,xlim,ylim)
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
#			scale_color_viridis_c(name = "Density",guide='none')+
#			geom_point(alpha=1)+
#			geom_hline(yintercept=0)+
##			geom_smooth(method="loess")+
#			xlab("RNA half-life")+ylab("log FC 4sU/no4sU")+
#			scale_x_continuous(breaks=c())+
#			coord_cartesian(ylim=ylim)
#}



PlotExpressionTest=function(data,w4sU,no4sU,ylim=c(-1,1),LFC.fun=PsiLFC,hl.quantile=0.8) {
	w=GetTable(data,type="count",columns=w4sU)[,1]
	n=if (is.numeric(no4sU)) no4sU[data$gene.info$Gene] else GetTable(data,type="count",columns=no4sU)[,1]
	use=!is.na(w+n)
	w=w[use]
	n=n[use]

	df=data.frame(lfc=LFC.fun(w,n),M=(log10(w+1)+log10(n+1))/2)
	ggplot(df,aes(M,lfc,color=density2d(M, lfc, n = 100)))+
			scale_color_viridis_c(name = "Density",guide='none')+
			geom_point(alpha=1)+
			geom_hline(yintercept=0)+
#			geom_smooth(method="loess")+
			xlab("Mean expression")+ylab("log FC 4sU/no4sU")+
			coord_cartesian(ylim=ylim)
}

PlotAnalyses=function(data,plot.fun,analyses=Analyses(data),add=NULL,...) {
  lapply(analyses,function(analysis) {
    re=plot.fun(data,analysis=analysis,...)
    if (!is.null(add)) for (e in if (is.list(add)) add else list(add)) re=re+e
    re
  })
}

VulcanoPlot=function(data,analysis=Analyses(data)[1],p.cutoff=0.05,lfc.cutoff=1,label.numbers=TRUE,...) {
  df=GetAnalysisTable(data,analyses=analysis,regex=FALSE,columns=c("LFC|Q"),gene.info = FALSE)
  names(df)=gsub(".*.Q","Q",gsub(".*.LFC","LFC",names(df)))
  g=PlotScatter(df,x=LFC,y=-log10(Q),remove.outlier = FALSE,...)+
    xlab(bquote(log[2]~FC))+
    ylab(bquote("-"~log[10]~FDR))+
    geom_hline(yintercept=-log10(p.cutoff),linetype=2)+
    geom_vline(xintercept=c(-lfc.cutoff,lfc.cutoff),linetype=2)+
    ggtitle(analysis)


  if (label.numbers) {
    n=table(cut(df$LFC,breaks=c(-Inf,-lfc.cutoff,lfc.cutoff,Inf)),factor(df$Q>p.cutoff,levels=c("FALSE","TRUE")))
    g=g+annotate("label",x=c(-Inf,0,Inf,-Inf,0,Inf),y=c(Inf,Inf,Inf,-Inf,-Inf,-Inf),label=paste0("n=",as.numeric(n)),hjust=c(-0.1,0.5,1.1,-0.1,0.5,1.1),vjust=c(1.1,1.1,1.1,-0.1,-0.1,-0.1))
  }
  g
}



MAPlot=function(data,analysis=Analyses(data)[1],aest=aes(),p.cutoff=0.05,lfc.cutoff=1,label.numbers=TRUE,highlight=NULL,label=NULL,repel=1) {
  df=GetAnalysisTable(data,analyses=analysis,regex=FALSE,columns=c("M|LFC|Q"),gene.info = FALSE)
  if (is.numeric(analysis)) analysis=Analyses(data)[analysis]
  names(df)=gsub(".*.Q","Q",gsub(".*.LFC","LFC",gsub(".*.M","M",names(df))))
  aes=modifyList(aes(M+1,LFC,color=ifelse(Q<p.cutoff,"Sig.","NS")),aest)
  g=ggplot(df,mapping=aes)+
    geom_point(size=1)+
    scale_x_log10()+
    scale_color_manual(values=c(Sig.="black",NS="grey50"),guide='none')+
    ylab(bquote(log[2]~FC))+
    xlab("Total expression")+
    geom_hline(yintercept=c(-lfc.cutoff,lfc.cutoff),linetype=2)+
    ggtitle(analysis)

  if (!is.null(highlight)) {
    if (is.list(highlight)){
      for (col in names(highlight)) {
        g=g+geom_point(data=df[ToIndex(data,highlight[[col]]),],color=col,size=1.5)
      }
    } else {
      g=g+geom_point(data=df[ToIndex(data,highlight),],color='red',size=1.5)
    }
  }
  if (!is.null(label)) {
    if (label=="auto") label=abs(df$LFC)>lfc.cutoff & df$Q<p.cutoff & !is.na(df$LFC) & !is.na(df$Q)
    df2=df
    df2$label=""
    df2[label,"label"]=rownames(df2)[label]
    g=g+ggrepel::geom_label_repel(data=df2,mapping=aes(label=label),show.legend = FALSE,force=repel)
  }
  if (label.numbers) {
    n=c(sum(df$LFC>lfc.cutoff & df$Q<p.cutoff,na.rm = TRUE),sum(df$LFC< -lfc.cutoff & df$Q<p.cutoff,na.rm = TRUE))
    g=g+annotate("label",x=c(Inf,Inf),y=c(Inf,-Inf),label=paste0("n=",n),hjust=c(1.1,1.1),vjust=c(1.1,-0.1))
  }
  g
}

PlotTypeDistribution=function(data,mode.slot=DefaultSlot(data),relative=FALSE) {
	df=GetTable(data,type=mode.slot)
	df=sapply(levels(data$gene.info$Type),function(type) colSums(df[ data$gene.info$Type==type,]))
	df=df[,colSums(df)>0]
	if (relative) {
		df=df/rowSums(df)*100
		type=sprintf("%s [%%]",type)
	}
	df=melt(df,varnames=c("Condition","Type"))
	ggplot(df,aes(Condition,value,fill=Type))+geom_bar(stat="Identity")+scale_fill_viridis_d()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab(mode.slot)+xlab(NULL)
}

PlotGeneOldVsNew=function(data,gene,slot=DefaultSlot(data),show.CI=FALSE,aest=aes(color=Condition,shape=Replicate)) {
  new=paste0("new.",slot)
  old=paste0("old.",slot)
  df=GetData(data,genes=gene,mode.slot=c(old,new),melt=F,coldata=T,ntr.na = FALSE)
  g=ggplot(df,modifyList(aes_string(old,new),aest))+
    geom_point(size=2)+
    scale_x_log10("Old RNA")+
    scale_y_log10("New RNA")
  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df=cbind(df,GetData(data,genes=gene,mode.slot=c("lower","upper",slot),melt=F,coldata=T,ntr.na = FALSE)[,c("lower","upper",slot)])
    g=g+geom_errorbar(data=df,mapping=aes_string(ymin=paste0("lower*",slot),ymax=paste0("upper*",slot)))
    g=g+geom_errorbarh(data=df,mapping=aes_string(xmin=paste0("(1-upper)*",slot),xmax=paste0("(1-lower)*",slot)))
  }
  g
}

PlotGeneTotalVsNtr=function(data,gene,slot=DefaultSlot(data),show.CI=FALSE,aest=aes(color=Condition,shape=Replicate)) {
  df=GetData(data,genes=gene,mode.slot=c("ntr",slot),melt=F,coldata=T,ntr.na = FALSE)
  g=ggplot(df,modifyList(aes_string(slot,"ntr"),aest))+
    geom_point(size=2)+
    scale_x_log10("Total RNA")+
    scale_y_continuous("NTR")
  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df=cbind(df,GetData(data,genes=gene,mode.slot=c("lower","upper"),melt=F,coldata=T,ntr.na = FALSE)[,c("lower","upper")])
    g=g+geom_errorbar(data=df,mapping=aes(ymin=lower,ymax=upper))
  }
  g
}

PlotGeneGroupsPoints=function(data,gene,group="Condition",slot=DefaultSlot(data),type="total",log=TRUE,show.CI=FALSE,aest=aes(color=Condition,shape=Replicate)) {
  df=GetData(data,genes=gene,mode.slot=c(slot,"ntr"),melt=F,coldata=T,ntr.na = FALSE)
  df$value=switch(type[1],total=df[[slot]],new=df[[slot]]*df[["ntr"]],old=df[[slot]]*(1-df[["ntr"]]),stop(paste0(type," unknown!")))
  g=ggplot(df,modifyList(aes_string(group,"value"),aest))+
    geom_point(size=2,position=if(show.CI) position_dodge(width=0.4) else "identity")+
    xlab(NULL)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  if (log) {
    g=g+scale_y_log10(paste0(toupper(substr(type,1,1)),substr(type,2,nchar(type))," RNA (",slot,")"))
  } else {
    g=g+scale_y_continuous(paste0(toupper(substr(type,1,1)),substr(type,2,nchar(type))," RNA (",slot,")"))
  }

  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df=cbind(df,GetData(data,genes=gene,mode.slot=c("lower","upper"),melt=F,coldata=T,ntr.na = FALSE)[,c("lower","upper")])
    g=switch(type[1],
             total=g,
             new=g+geom_errorbar(data=df,mapping=aes_string(ymin=paste0("lower*",slot),ymax=paste0("upper*",slot)),width=0,position=position_dodge(width=0.4)),
             old=g+geom_errorbar(data=df,mapping=aes_string(ymin=paste0("(1-upper)*",slot),ymax=paste0("(1-lower)*",slot)),width=0,position=position_dodge(width=0.4)),
             stop(paste0(type," unknown!"))
    )
  }
  g
}

PlotGeneGroupsBars=function(data,gene,slot=DefaultSlot(data),show.CI=FALSE) {
  df=GetData(data,genes=gene,mode.slot=paste0(c("new.","old."),slot),melt=T,coldata=T,ntr.na = FALSE)
  g=ggplot(df,aes(Name,Value,fill=Type))+
    geom_bar(stat="identity",position=position_stack())+
    scale_fill_manual(NULL,values = c('red','gray'))+
    xlab(NULL)+ylab(slot)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df2=GetData(data,genes=gene,mode.slot=c("lower","upper",slot),melt=F,coldata=T,ntr.na = FALSE)
    g=g+geom_errorbar(data=df2,mapping=aes_string(y=slot,fill=NULL,ymin=paste0("(1-upper)*",slot),ymax=paste0("(1-lower)*",slot)),width=0)
  }
  g
}

PlotGeneTimeCourse=function(data,gene,group="Condition",time=Design$dur.4sU,slot=DefaultSlot(data),aest=aes(color=Condition,shape=Replicate),average.lines=TRUE,log.y=TRUE, show.CI=FALSE) {
  df=GetData(data,genes=gene,mode.slot=slot,melt=F,coldata=T,ntr.na = FALSE)
  aes=modifyList(aes_string(time,"Value",group=group),aest)
  g=ggplot(df,mapping=aes)+
    geom_point(size=2)+
    xlab(NULL)+
    ylab(paste0(" RNA (",slot,")"))
  if (log.y) g=g+scale_y_log10()
  if (average.lines) {
    # compute average line:
    print(df)
    ddf=as.data.frame(lapply(aes,function(col) rlang::eval_tidy(col,data=df)))
    print(ddf)
    ddf=ddply(ddf,.(x,colour,group),function(s) c(Value=mean(s$y,na.rm=TRUE)))
    ddf[[group]]=ddf$group
    print(ddf)
    g=g+geom_line(data=ddf,mapping=aes(x,Value,color=colour,group=interaction(colour,group)),inherit.aes=F)
  }
  if (show.CI) {
    if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
    df2=GetData(data,genes=gene,mode.slot=c("lower","upper",slot),melt=F,coldata=T,ntr.na = FALSE)
    if (slot=="ntr") df2$ntr=1
    g=g+geom_errorbar(data=df2,mapping=aes_string(y=slot,fill=NULL,ymin=paste0("lower*",slot),ymax=paste0("upper*",slot)),width=0)
  }
  g
}

Plot=function(fun=NULL,...,gg=NULL) {
  function(data,gene) {
    if (is.null(fun)) return(NULL)
    re=fun(data=data,genes=gene,...)
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
      if (!is.null(add)) for (e in if (class(add)=="list") add else list(add)) re=re+e
      return(re)
    }

    if (is.null(value)) {
      value<<-do.call(FUN,c(list(data),param))
      if (!is.null(add)) for (e in if (class(add)=="list") add else list(add)) value<<-value+e
    }
    value
  }
}

