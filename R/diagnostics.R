

PlotMismatchPositionForCondition=function(x,...)  {
	UseMethod('PlotMismatchPositionForCondition',x)
}

PlotMismatchPositionForCondition.grandR=function(data,...) {
	t=read.tsv(paste0(data$prefix,".mismatch.position.tsv"))
	c=read.tsv(paste0(data$prefix,".clip.tsv"),header=F)
	c=setNames(c[,2],c[,1])
	PlotMismatchPositionForCondition(tab=t,clip=c,...)
}



PlotMismatchPositionForCondition.data.frame=function(tab,condition,clip=NA,Sense=NA,Category=NA,Corrected=NA) {
	if (!is.na(Sense)) tab=tab[tab$Sense==Sense,]
	if (!is.na(Category)) tab=tab[tab$Category==Category,]
	if (!is.na(Corrected)) tab=tab[tab$Corrected==Corrected,]
	tab$`Mismatch frequency`=tab[,condition]
	max1=max(c(0,tab$Position[tab$First.read==1]))
	max2=max(c(0,tab$Position[tab$First.read==0]))
	offset=round(max1*1.2)
	tab$x=ifelse(tab$First.read==1,tab$Position,max2-tab$Position+offset)
	tab$Corrected=factor(c("Uncorrected","Retained")[tab$Corrected+1],levels=c("Retained","Uncorrected"))
	tab$Sense=factor(c("Antisense","Sense")[tab$Sense+1],levels=c("Sense","Antisense"))

	sep=tab[tab$First.read==1&tab$Position==1,]
	sep$`Mismatch frequency`=NA
	sep$x=max(tab$Position[tab$First.read==1])+1
	tab=droplevels(rbind(tab,sep))
	tab$SenseCat=paste(tab$Sense,tab$Category,sep=",")

	lab1=labeling::extended(1,max1,5)
	lab2=labeling::extended(1,max2,5)

	ymax=quantile(tab$`Mismatch frequency`,0.99,na.rm=T)*1.2
	g=ggplot(tab,aes(x,`Mismatch frequency`,color=SenseCat,linetype=factor(Corrected),group=interaction(SenseCat,Corrected)))+cowplot::theme_cowplot()
	if (!is.na(clip["Used 5p1"]) || !is.na(clip["Used 3p1"])) {
		xmin=if (!is.na(clip["Used 5p1"])) clip["Used 5p1"]+1 else 1
		xmax=if (!is.na(clip["Used 3p1"])) max1-clip["Used 3p1"] else max1
		g=g+geom_rect(xmin=xmin,xmax=xmax,ymin=0,ymax=Inf,color="gray",fill="gray")
	}
	if (max2>0 && (!is.na(clip["Used 5p2"]) || !is.na(clip["Used 3p2"]))) {
		xmin=if (!is.na(clip["Used 3p2"])) clip["Used 3p2"]+offset else offset
		xmax=if (!is.na(clip["Used 5p2"])) max2+offset-clip["Used 5p2"]-1 else max2+offset
		g=g+geom_rect(xmin=xmin,xmax=xmax,ymin=0,ymax=Inf,color="gray",fill="gray")
	}
	if (!is.na(clip["Inferred 5p1"])) g=g+geom_vline(xintercept=clip["Inferred 5p1"]+1)
	if (!is.na(clip["Inferred 3p1"])) g=g+geom_vline(xintercept=max1-clip["Inferred 3p1"])
	if (max2>0 && !is.na(clip["Inferred 5p2"])) g=g+geom_vline(xintercept=max2-clip["Inferred 5p2"]+offset-1)
	if (max2>0 && !is.na(clip["Inferred 3p2"])) g=g+geom_vline(xintercept=offset+clip["Inferred 3p2"])
	g=g+geom_line()+facet_grid(Read~Genomic)+scale_linetype_discrete(NULL,guide=if(length(unique(tab$Corrected))<=1) F else guide_legend())+scale_color_brewer(NULL,palette="Dark2")+scale_x_continuous(breaks=c(lab1,max2-lab2+offset),labels=c(lab1,lab2))+coord_cartesian(ylim=c(0,ymax))+xlab(NULL)+ggtitle(condition)
	g
}


PlotMismatchPositionForType=function(x,...)  {
	UseMethod('PlotMismatchPositionForType',x)
}

PlotMismatchPositionForType.grandR=function(data,...) {
	t=read.tsv(paste0(data$prefix,".mismatch.position.tsv"))
	PlotMismatchPositionForType(t,...)
}

PlotMismatchPositionForType.data.frame=function(tab,genomic,read,clip=c(Inferred1=NA,Inferred2=NA,Used1=NA,Used2=NA),Sense=NA,Category=NA,Corrected=NA) {
	tab=tab[tab$Genomic==genomic & tab$Read==read,]
	if (!is.na(Sense)) tab=tab[tab$Sense==Sense,]
	if (!is.na(Category)) tab=tab[tab$Category==Category,]
	if (!is.na(Corrected)) tab=tab[tab$Corrected==Corrected,]
	max1=max(c(0,tab$Position[tab$First.read==1]))
	max2=max(c(0,tab$Position[tab$First.read==0]))
	offset=round(max1*1.2)
	tab$Corrected=factor(c("Uncorrected","Retained")[tab$Corrected+1],levels=c("Retained","Uncorrected"))
	tab$Sense=factor(c("Antisense","Sense")[tab$Sense+1],levels=c("Sense","Antisense"))

	tab=reshape2::melt(tab,value.name="Mismatch frequency",variable.name="Sample",id.vars=names(tab)[1:7])
	tab$SenseCat=paste(tab$Sense,tab$Category,sep=",")
	tab$x=ifelse(tab$First.read==1,tab$Position,max2-tab$Position+offset)
	sep=tab[tab$First.read==1&tab$Position==1,]
	sep$`Mismatch frequency`=NA
	sep$x=max(tab$Position[tab$First.read==1])+1
	tab=droplevels(rbind(tab,sep))


	lab1=labeling::extended(1,max1,5)
	lab2=labeling::extended(1,max2,5)

	ymax=max(plyr::ddply(tab,plyr::.(Sample),function(s) quantile(s$`Mismatch frequency`,0.95,na.rm=TRUE)*1.2)[,2])
	g=ggplot(tab,aes(x,`Mismatch frequency`,color=Sample,linetype=factor(Corrected)))+cowplot::theme_cowplot()
	if (!is.na(clip["Used 5p1"]) || !is.na(clip["Used 3p1"])) {
		xmin=if (!is.na(clip["Used 5p1"])) clip["Used 5p1"]+1 else 1
		xmax=if (!is.na(clip["Used 3p1"])) max1-clip["Used 3p1"] else max1
		g=g+geom_rect(xmin=xmin,xmax=xmax,ymin=0,ymax=Inf,color="gray",fill="gray")
	}
	if (max2>0 && (!is.na(clip["Used 5p2"]) || !is.na(clip["Used 3p2"]))) {
		xmin=if (!is.na(clip["Used 3p2"])) clip["Used 3p2"]+offset else offset
		xmax=if (!is.na(clip["Used 5p2"])) max2+offset-clip["Used 5p2"]-1 else max2+offset
		g=g+geom_rect(xmin=xmin,xmax=xmax,ymin=0,ymax=Inf,color="gray",fill="gray")
	}
	if (!is.na(clip["Inferred 5p1"])) g=g+geom_vline(xintercept=clip["Inferred 5p1"]+1)
	if (!is.na(clip["Inferred 3p1"])) g=g+geom_vline(xintercept=max1-clip["Inferred 3p1"])
	if (max2>0 && !is.na(clip["Inferred 5p2"])) g=g+geom_vline(xintercept=max2-clip["Inferred 5p2"]+offset-1)
	if (max2>0 && !is.na(clip["Inferred 3p2"])) g=g+geom_vline(xintercept=offset+clip["Inferred 3p2"])
	g=g+geom_line()+facet_grid(Category~Sense)+scale_linetype_discrete(NULL,guide=if(length(unique(tab$Corrected))<=1) F else guide_legend())+scale_x_continuous(breaks=c(lab1,max2-lab2+offset),labels=c(lab1,lab2))+coord_cartesian(ylim=c(0,ymax))+xlab(NULL)+ggtitle(paste0(genomic,">",read))
	g
}


PlotMismatchFreq=function(x,...)  {
	UseMethod('PlotMismatchFreq',x)
}

PlotMismatchFreq.grandR=function(data,...) {
	t=read.tsv(paste0(data$prefix,".mismatch.freq.tsv"))
	subr=read.tsv(paste0(data$prefix,".subread.tsv"))
	subr$Semantic=factor(as.character(subr$Semantic),levels=subr$Semantic)
	t=merge(tab,subr,by="Subread")

	PlotMismatchFreq(t,...)
}

PlotMismatchFreq.data.frame=function(tab,category,ncond.boxplot=120) {
	tab=tab[tab$Category==category,]
	ncond=nrow(tab)/nrow(unique(data.frame(tab$Genomic,tab$Read,tab$Semantic)))

	if (ncond>=ncond.boxplot) {
		max=max(plyr::ddply(tab,plyr::.(Category,Condition,Semantic,Genomic,Read),function(s) c(max=quantile(s$Frequency,0.99)))$max)
		if (length(unique(tab$Condition))<ncond/100) {
			ggplot(tab,aes(paste0(Genomic,"->",Read),Frequency,fill=Condition))+cowplot::theme_cowplot()+
		    geom_hline(yintercept=0,linetype=2)+geom_boxplot(width=0.4,outlier.size = 0.1)+facet_grid(Semantic~.,scales="free_y")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab(NULL)+ggtitle(category)+scale_fill_brewer(NULL,palette="Set2")+coord_cartesian(ylim=c(0,max))
		} else {
			ggplot(tab,aes(paste0(Genomic,"->",Read),Frequency))+cowplot::theme_cowplot()+
		    geom_hline(yintercept=0,linetype=2)+geom_boxplot(width=0.8)+facet_grid(Semantic~.,scales="free_y")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab(NULL)+ggtitle(category)+coord_cartesian(ylim=c(0,max))
		}
	} else {
		ggplot(tab,aes(paste0(Genomic,"->",Read),Frequency,color=Condition))+cowplot::theme_cowplot()+
	    geom_hline(yintercept=0,linetype=2)+geom_point(position=position_dodge(width=0.7))+facet_grid(Semantic~.,scales="free_y")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab(NULL)+ggtitle(category)
	}

}




PlotModelNtr=function(x,...)  {
	UseMethod('PlotModelNtr',x)
}
PlotModelNtr.grandR=function(data,label="4sU",...) {
	t=read.tsv(paste0(data$prefix,".model.parameters.tsv"))
	PlotModelNtr(t,label,...)
}
PlotModelNtr.data.frame=function(tab,label="4sU",estimator="Separate",model=c("Binom","TB-Binom")) {
	tab=tab[tab$Label==label & tab$Estimator==estimator,]

	param=paste(model[1],"ntr")
	lower=paste("Lower",model[1],"ntr")
	upper=paste("Upper",model[1],"ntr")
	qq=function(s) paste0("`",s,"`")

	if (lower %in% names(tab)) {
		ggplot(tab,aes_string("Condition",qq(param),color="Subread",ymin=qq(lower),ymax=qq(upper)))+
	    cowplot::theme_cowplot()+
	    geom_errorbar(width=0.1,position=position_dodge(w=0.2))+
	    geom_point(position=position_dodge(w=0.2))+
	    coord_cartesian(ylim=c(0,max(tab[,upper])))+
	    ylab("Global NTR")+xlab(NULL)+
	    scale_color_brewer(NULL,palette="Dark2")+
	    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	}else{
		ggplot(tab,aes_string("Condition",qq(param),color="Subread"))+
	    cowplot::theme_cowplot()+
	    geom_point(position=position_dodge(w=0.2))+
	    coord_cartesian(ylim=c(0,max(tab[,param])))+
	    ylab("Global NTR")+xlab(NULL)+
	    scale_color_brewer(NULL,palette="Dark2")+
	    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	}
}


PlotModelConv=function(x,...)  {
	UseMethod('PlotModelConv',x)
}
PlotModelConv.grandR=function(data,label="4sU",...) {
	t=read.tsv(paste0(data$prefix,".model.parameters.tsv"))
	PlotModelConv(t,label,...)
}
PlotModelConv.data.frame=function(tab,label="4sU",estimator="Separate",model=c("Binom","TB-Binom")) {
	tab=tab[tab$Label==label & tab$Estimator==estimator,]

	if (model[1]=="TB-Binom") {
		param=paste(model[1],"p.mconv")
		lower=paste("Lower",model[1],"p.mconv")
		upper=paste("Upper",model[1],"p.mconv")
		qq=function(s) paste0("`",s,"`")

		tab2=tab
		tab2$`TB-Binom p.mconv`=etbeta(tab2$`TB-Binom p.err`,tab2$`TB-Binom p.mconv`,exp(tab$`TB-Binom shape`),exp(-tab$`TB-Binom shape`))
		if (lower %in% names(tab)) {
			tab2$`Lower TB-Binom p.mconv`=etbeta(tab2$`Lower TB-Binom p.err`,tab2$`Lower TB-Binom p.mconv`,exp(tab2$`Lower TB-Binom shape`),exp(-tab2$`Lower TB-Binom shape`))
			tab2$`Upper TB-Binom p.mconv`=etbeta(tab2$`Upper TB-Binom p.err`,tab2$`Upper TB-Binom p.mconv`,exp(tab2$`Upper TB-Binom shape`),exp(-tab2$`Upper TB-Binom shape`))
			tab=rbind(cbind(Type="Max",tab),cbind(Type="Mean",tab2))
			ggplot(tab,aes_string("Condition",qq(param),color="Subread",shape="Type",ymin=qq(lower),ymax=qq(upper)))+cowplot::theme_cowplot()+
			  geom_errorbar(width=0.1,position=position_dodge(w=0.2))+geom_point(position=position_dodge(w=0.2))+coord_cartesian(ylim=c(0,max(tab[,upper])))+ylab("p.conv")+xlab(NULL)+scale_color_brewer(NULL,palette="Dark2")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_shape_manual(NULL,values=c(Max=1,Mean=19))
		} else {
			tab=rbind(cbind(Type="Max",tab),cbind(Type="Mean",tab2))
			ggplot(tab,aes_string("Condition",qq(param),color="Subread",shape="Type"))+cowplot::theme_cowplot()+
			  geom_point(position=position_dodge(w=0.2))+coord_cartesian(ylim=c(0,max(tab[,param])))+ylab("p.conv")+xlab(NULL)+scale_color_brewer(NULL,palette="Dark2")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+scale_shape_manual(NULL,values=c(Max=1,Mean=19))
		}
	} else {
		param=paste(model[1],"p.conv")
		lower=paste("Lower",model[1],"p.conv")
		upper=paste("Upper",model[1],"p.conv")
		qq=function(s) paste0("`",s,"`")

		if (lower %in% names(tab)) {
			ggplot(tab,aes_string("Condition",qq(param),color="Subread",ymin=qq(lower),ymax=qq(upper)))+cowplot::theme_cowplot()+
		    geom_errorbar(width=0.1,position=position_dodge(w=0.2))+geom_point(position=position_dodge(w=0.2))+coord_cartesian(ylim=c(0,max(tab[,upper])))+ylab("p.conv")+xlab(NULL)+scale_color_brewer(NULL,palette="Dark2")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
		} else {
			ggplot(tab,aes_string("Condition",qq(param),color="Subread"))+cowplot::theme_cowplot()+
		    geom_point(position=position_dodge(w=0.2))+coord_cartesian(ylim=c(0,max(tab[,param])))+ylab("p.conv")+xlab(NULL)+scale_color_brewer(NULL,palette="Dark2")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
		}
	}
}

PlotModelErr=function(x,...)  {
	UseMethod('PlotModelErr',x)
}
PlotModelErr.grandR=function(data,label="4sU",...) {
	t=read.tsv(paste0(data$prefix,".model.parameters.tsv"))
	PlotModelErr(t,label,...)
}
PlotModelErr.data.frame=function(tab,label="4sU",estimator="Separate",model=c("Binom","TB-Binom")) {
	tab=tab[tab$Label==label & tab$Estimator==estimator,]

	param=paste(model[1],"p.err")
	lower=paste("Lower",model[1],"p.err")
	upper=paste("Upper",model[1],"p.err")
	qq=function(s) paste0("`",s,"`")

	if (lower %in% names(tab)) {
		ggplot(tab,aes_string("Condition",qq(param),color="Subread",ymin=qq(lower),ymax=qq(upper)))+cowplot::theme_cowplot()+
	    geom_errorbar(width=0.1,position=position_dodge(w=0.2))+geom_point(position=position_dodge(w=0.2))+coord_cartesian(ylim=c(0,max(tab[,upper])))+ylab("p.err")+xlab(NULL)+scale_color_brewer(NULL,palette="Dark2")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	}else{
		ggplot(tab,aes_string("Condition",qq(param),color="Subread"))+cowplot::theme_cowplot()+
	    geom_point(position=position_dodge(w=0.2))+coord_cartesian(ylim=c(0,max(tab[,param])))+ylab("p.err")+xlab(NULL)+scale_color_brewer(NULL,palette="Dark2")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	}
}


PlotModelShape=function(x,...)  {
	UseMethod('PlotModelShape',x)
}
PlotModelShape.grandR=function(data,label="4sU",...) {
	t=read.tsv(paste0(data$prefix,".model.parameters.tsv"))
	PlotModelShape(t,label,...)
}
PlotModelShape.data.frame=function(tab,label="4sU",estimator="Separate") {
	tab=tab[tab$Label==label & tab$Estimator==estimator,]
	if ("Lower TB-Binom shape" %in% names(tab)) {
		ggplot(tab,aes(Condition,`TB-Binom shape`,color=Subread,ymin=`Lower TB-Binom shape`,ymax=`Upper TB-Binom shape`))+cowplot::theme_cowplot()+
	    geom_errorbar(width=0.1,position=position_dodge(w=0.2))+geom_point(position=position_dodge(w=0.2))+ylab("shape")+xlab(NULL)+scale_color_brewer(NULL,palette="Dark2")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	}else{
		ggplot(tab,aes(Condition,`TB-Binom shape`,color=Subread))+cowplot::theme_cowplot()+
	    geom_point(position=position_dodge(w=0.2))+ylab("shape")+xlab(NULL)+scale_color_brewer(NULL,palette="Dark2")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	}
}


PlotModelLabelTimeCourse=function(x,...)  {
	UseMethod('PlotModelLabelTimeCourse',x)
}
PlotModelLabelTimeCourse.grandR=function(data,label="4sU",...) {
	t=read.tsv(paste0(data$prefix,".model.parameters.tsv"))
	e=read.tsv(paste0(data$prefix,".experimentalDesign.tsv"))
	PlotModelLabelTimeCourse(t,e,label,...)
}
PlotModelLabelTimeCourse.data.frame=function(tab,e,label="4sU",estimator="Separate") {
	tab=tab[tab$Label==label & tab$Estimator==estimator,]
	e$Condition=ifelse(is.na(e$Sample),as.character(e$Library),paste(e$Library,e$Sample,sep="."))
	e=unique(cbind(Condition=e$Condition,e[,grep(label,names(e))]))
	names(e)=gsub(paste0(label," "),"",names(e))
	e=e[!is.na(e$concentration),]

	x=seq(0,1,length.out=100)
	make.df=function(row) data.frame(x=x,`Labeling efficiency`=qtbeta(x,tab$`TB-Binom p.err`[row],tab$`TB-Binom p.mconv`[row],exp(tab$`TB-Binom shape`[row]),exp(-tab$`TB-Binom shape`[row])),Condition=tab$Condition[row],Subread=tab$Subread[row],check.names=FALSE)
	df=do.call("rbind",lapply(1:dim(tab)[1],make.df))

	if (length(unique(df$Condition))!=length(e$Condition)) stop("Conditions in experimental design and parameter file are inconsistent!")
	df=merge(df,e)
	df$concentration=factor(df$concentration,levels=sort(unique(df$concentration)))
	df$Time=df$x*df$duration+(1-df$x)*df$chase
	ggplot(df,aes(Time,`Labeling efficiency`,color=Condition,linetype=Subread))+cowplot::theme_cowplot()+
	  geom_line()+scale_color_viridis_d(NULL,)+scale_linetype_discrete(NULL)
}


PlotModelCompareErrPrior=function(x,...)  {
	UseMethod('PlotModelCompareErrPrior',x)
}
PlotModelCompareErrPrior.grandR=function(data,label="4sU",...) {
	t=read.tsv(paste0(data$prefix,".model.parameters.tsv"))
	PlotModelCompareErrPrior(t,e,label,...)
}
PlotModelCompareErrPrior.data.frame=function(tab,label="4sU",estimator="Separate",model=c("Binom","TB-Binom")) {
	tab=tab[tab$Label==label & tab$Estimator==estimator,]

	param=paste(model[1],"p.err")
	qq=function(s) paste0("`",s,"`")

	ggplot(tab,aes_string("(`Lower prior p.err`+`Upper prior p.err`)/2",qq(param),color="Subread",xmin="`Lower prior p.err`",xmax="`Upper prior p.err`"))+cowplot::theme_cowplot()+
	  geom_point()+geom_errorbarh()+geom_abline()+scale_color_brewer(NULL,palette="Dark2")+xlab("Prior p.err")
}

PlotModelCompareNtr=function(x,...)  {
	UseMethod('PlotModelCompareNtr',x)
}
PlotModelCompareNtr.grandR=function(data,label="4sU",...) {
	t=read.tsv(paste0(data$prefix,".model.parameters.tsv"))
	PlotModelCompareNtr(t,e,label,...)
}
PlotModelCompareNtr.data.frame=function(tab,label="4sU",estimator="Separate") {
	tab=tab[tab$Label==label & tab$Estimator==estimator,]
	ggplot(tab,aes(`Binom ntr`,`TB-Binom ntr`,color=Subread))+cowplot::theme_cowplot()+
	  geom_point()+geom_abline()+scale_color_brewer(NULL,palette="Dark2")
}

PlotModelCompareErr=function(x,...)  {
	UseMethod('PlotModelCompareErr',x)
}
PlotModelCompareErr.grandR=function(data,label="4sU",...) {
	t=read.tsv(paste0(data$prefix,".model.parameters.tsv"))
	PlotModelCompareErr(t,e,label,...)
}
PlotModelCompareErr.data.frame=function(tab,label="4sU",estimator="Separate") {
	tab=tab[tab$Label==label & tab$Estimator==estimator,]
	ggplot(tab,aes(`Binom p.err`,`TB-Binom p.err`,color=Subread))+cowplot::theme_cowplot()+
	  geom_point()+geom_abline()+scale_color_brewer(NULL,palette="Dark2")
}


PlotModelCompareConv=function(x,...)  {
	UseMethod('PlotModelCompareConv',x)
}
PlotModelCompareConv.grandR=function(data,label="4sU",...) {
	t=read.tsv(paste0(data$prefix,".model.parameters.tsv"))
	PlotModelCompareConv(t,e,label,...)
}
PlotModelCompareConv.data.frame=function(tab,label="4sU",estimator="Separate") {
	tab=tab[tab$Label==label & tab$Estimator==estimator,]
	ggplot(tab,aes(`Binom p.conv`,etbeta(`TB-Binom p.err`,`TB-Binom p.mconv`,exp(`TB-Binom shape`),exp(-`TB-Binom shape`)),color=Subread))+cowplot::theme_cowplot()+
	  geom_point()+geom_abline()+scale_color_brewer(NULL,palette="Dark2")+ylab("Mean TB-Binom p.conv")
}


PlotModelCompareLL=function(x,...)  {
	UseMethod('PlotModelCompareLL',x)
}
PlotModelCompareLL.grandR=function(data,label="4sU",...) {
	t=read.tsv(paste0(data$prefix,".model.parameters.tsv"))
	PlotModelCompareLL(t,e,label,...)
}
PlotModelCompareLL.data.frame=function(tab,label="4sU",estimator="Separate") {
	tab=tab[tab$Label==label & tab$Estimator==estimator,]
	ggplot(tab,aes(Condition,`Binom log likelihood`-`TB-Binom log likelihood`,color=Subread))+cowplot::theme_cowplot()+
	  geom_point(position=position_dodge(w=0.2))+coord_cartesian(ylim=c(min(tab$`Binom log likelihood`-tab$`TB-Binom log likelihood`),0))+ylab("log likelihood ratio Binom/TB-Binom")+xlab(NULL)+scale_color_brewer(NULL,palette="Dark2")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_hline(yintercept=0)+geom_hline(yintercept=-qchisq(0.95,df=1),linetype=2)
}



PlotProfileLikelihood=function(x,...)  {
	UseMethod('PlotProfileLikelihood',x)
}

PlotProfileLikelihood.grandR=function(data,condition,...) {
	t=read.tsv(paste0(data$prefix,".model.profile.tsv"))
	PlotProfileLikelihood(t,condition,...)
}

PlotProfileLikelihood.data.frame=function(tab,condition,subread,label="4sU") {
	tab=tab[tab$Condition==condition & tab$Subread==subread & tab$Label==label,]

	plot1=function(x,y,ylab=y,xlab=NULL) {
		t=tab[tab$Parameter==x,]
		g=ggplot(t,aes_string(x,y))+cowplot::theme_cowplot()+
		  ylab(ylab)+xlab(xlab)+theme(axis.text.x = element_text(angle = 40, vjust = 1, hjust=1))
		if (y=="deltaLL") g=g+geom_hline(yintercept=0,color="grey",linetype=2)
		g+geom_line()
	}

	g=plot_grid(
		plot1("ntr","deltaLL",ylab="log likelihood ratio",xlab="ntr"),
		plot1("p.err","deltaLL",ylab="log likelihood ratio",xlab="p.err"),
		plot1("p.mconv","deltaLL",ylab="log likelihood ratio",xlab="p.mconv"),
		plot1("shape","deltaLL",ylab="log likelihood ratio",xlab="shape"),

		ggdraw(),
		plot1("p.err","ntr"),
		plot1("p.mconv","ntr"),
		plot1("shape","ntr"),

		plot1("ntr","p.err"),
		ggdraw(),
		plot1("p.mconv","p.err"),
		plot1("shape","p.err"),

		plot1("ntr","p.mconv"),
		plot1("p.err","p.mconv"),
		ggdraw(),
		plot1("shape","p.mconv"),

		plot1("ntr","shape"),
		plot1("p.err","shape"),
		plot1("p.mconv","shape"),
		ggdraw(),

		align="hv",axis="b",ncol=4)

	title=ggdraw()+draw_label(sprintf("%s (%s) - %s",condition,subread,label), size=15)
	plot_grid(title, g, ncol=1, rel_heights=c(0.1, 1))
}


