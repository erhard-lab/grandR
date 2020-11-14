
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



PlotPCA=function(data,type=c("Total","New","Old"),ntop=500,aest=aes(color=Sample),x=1,y=2) {
	mat=cnt(switch(type[1],Total=data$data$count,New=data$data$count*data$data$ntr,Old=data$data$count*(1-data$data$ntr)))
	vsd <- vst(mat)

  	rv <- rowVars(vsd)
  	select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  	pca <- prcomp(t(vsd[select,]))
	percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
	d <- as.data.frame(pca$x)
	names(d)=paste0("PC",1:dim(d)[2])	
	d=cbind(d, data$coldata)

	ggplot(d,modifyList(aes_string(paste0("PC",x),paste0("PC",y)),aest))+ geom_point(size=3)+xlab(paste0("PC",x,": ",round(percentVar[x] * 100),"% variance"))+ylab(paste0("PC",y,": ",round(percentVar[y] * 100),"% variance"))+coord_fixed()
}

PlotTestOverlap=function(data,name="lrt",alpha=0.01,type=c("venn","euler")) {
	mat=data$data$diffexp[,grepl(paste0("^",name,".*p$"),names(data$data$diffexp))]
	df=setNames(as.data.frame(mat<alpha & !is.na(mat)),c("Total","New","Old"))
	pl=switch(type[1],euler=euler(df),venn=venn(df))
	plot(pl,main=name)
}


PlotScatter=function(x,...)  {
	UseMethod('PlotScatter',x)  
}



PlotScatter.grandR=function(data,cx,cy,type="tpm",type.x=type,type.y=type,...) {
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
	xlim=(quantile(df[,1]-xmean,pnorm(c(-2,2)))*1.5)+xmean
	ylim=(quantile(df[,2]-ymean,pnorm(c(-2,2)))*1.5)+ymean

	g=ggplot(df,aes(A,B,color=densit2d(log(1+A), log(1+B), n = 100)))+
			geom_point()+
			scale_color_viridis_c(name = "Density",guide=FALSE)+
			xlab(dfnames[1])+ylab(dfnames[2])+
			coord_cartesian(ylim=ylim,xlim=xlim)
	if (log.x) g=g+scale_x_log10()
	if (log.y) g=g+scale_y_log10()
	g
}

PlotToxicityTest=function(data,w4sU,no4sU,ylim=c(-1,1),LFC.fun=PsiLFC,hl.quantile=0.8) {
	w=GetData(data,"count",conditions=w4sU,table=T)[,1]
	n=if (is.numeric(no4sU)) no4sU[data$gene.info$Gene] else GetData(data,"count",conditions=no4sU,table=T)[,1]
	ntr=GetData(data,"ntr",conditions=w4sU,table=T)[,1]
	use=!is.na(w+n+ntr)
	w=w[use]
	n=n[use]
	ntr=ntr[use]

	phl=comp.hl(ntr,1)	
	df=data.frame(lfc=LFC.fun(w,n),PHL=phl)[ntr<1,]
	df=df[df$PHL<quantile(df$PHL[is.finite(df$PHL)],hl.quantile),]
	ggplot(df,aes(PHL,lfc,color=density2d(PHL, lfc, n = 100)))+
			scale_color_viridis_c(name = "Density",guide=FALSE)+
			geom_point(alpha=1)+
			geom_hline(yintercept=0)+
#			geom_smooth(method="loess")+
			xlab("RNA half-life")+ylab("log FC 4sU/no4sU")+
			scale_x_continuous(breaks=c())+
			coord_cartesian(ylim=ylim)
}

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


