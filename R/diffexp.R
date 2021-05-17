

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

	res.tot <- results(dds.tot)
	res.new <- results(dds.new)
	res.old <- results(dds.old)

	
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
