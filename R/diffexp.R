

cnt=function(m) {
	mode(m) <- "integer"
	m
}


TestGenesLRT=function(data,target=~Combined,background=~1,name="lrt",verbose=FALSE,subset=!data$coldata$no4sU) {

	
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
