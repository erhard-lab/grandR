
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)


comp.hl=function(p,time) ifelse(p==0,Inf,log(2)/(-1/time*log(1-p)))

comp.tpm=function(cmat,lengths,subset=NULL) {
	rpk=cmat/(lengths/1000)
	scale=colSums(if(!is.null(subset)) rpk[subset,] else rpk)/1E6
	t(t(rpk)/scale)
}
tomat=function(m,names,cnames){
	m=as.matrix(m)
	m[is.na(m)]=0
	colnames(m)=gsub(" .*","",cnames)
	rownames(m)=names
	m
}

GeneType=list(
	mito=function(data) grepl("^MT-",data$Symbol,ignore.case=TRUE),
	ERCC=function(data) grepl("ERCC-[0-9]{5}",data$Gene),
	Cellular=function(data) grepl("^ENS.*G\\d+$",data$Gene),
	Unknown=function(data) rep(T,dim(data)[1])
)

Design=list(conc.4sU="concentration.4sU",dur.4sU="duration.4sU","Replicate"="Replicate")

# TODO: param should be folder to read numis, rates, etc.
ReadGRAND=function(prefix, verbose=FALSE, classify.genes=GeneType,design=c("Condition",Design$Replicate),Unknown=NA,rename.samples=NULL) {


checknames=function(a,b){
	if (!all(colnames(a)==colnames(b))) stop("Column names do not match!")
	if (!all(rownames(a)==rownames(b))) stop("Row names do not match!")
}


	file=if (file.exists(prefix)) prefix else paste0(prefix,".tsv")
	prefix=gsub(".tsv$","",file)
	
	if (verbose) cat("Checking file...\n")
	con <- file(file,"r")
	header <- strsplit(readLines(con,n=1),"\t")[[1]]
	close(con)

	if (header[1]!="Gene" || header[2]!="Symbol" || !grepl("Readcount$",header[3])) stop("File is not a GRAND-SLAM output file!")
	conds=strsplit(gsub(" Readcount","",header[3]),".",fixed=TRUE)[[1]]
	if (length(conds)!=length(design)) stop(paste0("Design parameter is incompatible with input data: ",paste(conds,collapse=".")))


	if (verbose) cat("Reading files...\n")
	data=read.delim(file,stringsAsFactors=FALSE,check.names=FALSE)
	if (anyDuplicated(data$Gene)) {
		warning("Duplicate gene names present, making unique!",call. = FALSE,immediate. = TRUE)
		data$Gene=make.unique(data$Gene)
	}
	if (anyDuplicated(data$Symbol)) {
		warning("Duplicate gene symbols present, making unique!",call. = FALSE,immediate. = TRUE)
		data$Symbol=make.unique(data$Symbol)
	}
	if (!is.null(rename.samples)) {
		if (verbose) cat("Renaming samples...\n")
		for (from in names(rename.samples)) {
			names(data)=gsub(from,rename.samples[from],names(data))
		}
	}
	if (verbose) cat("Processing...\n")

	ntr.mean = grepl("Mean",names(data))
	ntr = grepl("MAP",names(data))
	alpha = grepl("alpha",names(data))
	beta = grepl("beta",names(data))
	count = grepl("Readcount",names(data))

	if (!is.na(Unknown)) names(classify.genes)[names(classify.genes)=="Unknown"]=Unknown
	gene.info = data.frame(Gene=as.character(data$Gene),Symbol=as.character(data$Symbol),Length=data$Length,stringsAsFactors=FALSE)
	gene.info$Type=NA
	for (i in length(classify.genes):1) gene.info$Type[classify.genes[[i]](gene.info)]=names(classify.genes)[i]
	gene.info$Type=factor(gene.info$Type,levels=names(classify.genes))

	re=list()
	re$count=tomat(data[,count],data$Gene,names(data)[count])
	re$tpm=comp.tpm(re$count,gene.info$Length)
	checknames(re$count,re$tpm)
	re$ntr.mean=tomat(data[,alpha]/(data[,alpha]+data[,beta]),data$Gene,names(data)[ntr.mean])
	re$ntr=tomat(data[,ntr],data$Gene,names(data)[ntr])
	re$alpha=tomat(data[,alpha],data$Gene,names(data)[alpha])
	re$beta=tomat(data[,beta],data$Gene,names(data)[beta])
	
	no4sU.cols=!(colnames(re$count) %in% colnames(re$ntr))
	if (sum(no4sU.cols)>0) {
		correctmat=function(m) {
			r=matrix(NA,ncol=dim(re$count)[2],nrow=dim(re$count)[1])
			rownames(r)=rownames(re$count)
			colnames(r)=colnames(re$count)
			r[,colnames(r)[!no4sU.cols]]=m[,colnames(r)[!no4sU.cols]]
			r
		}
		re$ntr.mean=correctmat(re$ntr.mean)
		re$ntr=correctmat(re$ntr)
		re$alpha=correctmat(re$alpha)
		re$beta=correctmat(re$beta)
	}

	checknames(re$count,re$ntr)
	checknames(re$count,re$ntr.mean)
	checknames(re$count,re$alpha)
	checknames(re$count,re$beta)

	coldata=data.frame(Name=colnames(re$count))
	spl=strsplit(as.character(coldata$Name),".",fixed=TRUE)
	for (i in 1:length(design)) coldata=cbind(coldata,factor(sapply(spl,function(v) v[i]),levels=unique(sapply(spl,function(v) v[i]))))
	names(coldata)[-1]=design
	rownames(coldata)=coldata$Name
	coldata$Sample=interaction(coldata[ !(names(coldata) %in% c("Name","Replicate"))],drop=TRUE)
	coldata$no4sU=no4sU.cols

	invisible(grandR(prefix,gene.info,re,coldata))
}

grandR=function(prefix,gene.info,data,coldata) {
	info=list()
	info$prefix=prefix
	info$gene.info=gene.info
	info$data=data
	info$coldata=coldata
	class(info)="grandR"
	invisible(info)
}
dim.grandR=function(d) c(NumRows(d),NumCond(d))
is.grandR <- function(x) inherits(x, "grandR")

data.apply=function(data,fun,gene.info=TRUE,...) {
	re=list()
	for (l1 in names(data$data)) {
		re[[l1]]=fun(data$data[[l1]],...)
	}
	ngene.info=if (gene.info) fun(data$gene.info,...) else data$gene.info
	invisible(grandR(data$prefix,ngene.info,re,data$coldata))
}


GetAbsolute=function(data,dilution=2E6,volume=200) {
	fd=data.frame(gene_short_name=rownames(data$data$tpm))
	rownames(fd)=rownames(data$data$tpm)
	ercc=data$data$tpm[data$gene.info$Type=='ERCC',]

	cds <- newCellDataSet(data$data$tpm,
			featureData = new ("AnnotatedDataFrame",fd),
		        lowerDetectionLimit = 0.1,
		        expressionFamily = tobit(Lower = 0.1))
	rpc_matrix <- monocle::relative2abs(cds, method = "num_genes", ERCC_controls=ercc, ERCC_annotation=monocle::spike_df,dilution=dilution, volume=volume)
	rpc_matrix[data$data$tpm==0]=0

	colnames(rpc_matrix)=colnames(data$data$tpm)
	rownames(rpc_matrix)=rownames(data$data$tpm)
	rpc_matrix[,grepl("test",colnames(rpc_matrix))]=NA

	data$data$absolute=rpc_matrix
	invisible(data)
}

ToIndex=function(data,gene) {
	if (is.numeric(gene)) return(gene)
	m.gene=merge(cbind(Index=1:dim(data$gene.info)[1],data$gene.info),data.frame(Gene=gene))
	m.symbol=merge(cbind(Index=1:dim(data$gene.info)[1],data$gene.info),data.frame(Symbol=gene))

	if (dim(m.gene)[1]>dim(m.symbol)[1]) m.gene$Index else m.symbol$Index
}

NumRows=function(data) dim(data$gene.info)[1]
NumCond=function(data) dim(data$coldata)[1]


GetData=function(data,type,conditions=colnames(data$data$count),gene=1:NumRows(data),melt=FALSE,coldata=TRUE,table=FALSE,keep.ntr.na=TRUE) {
	if (table) { 
		melt=FALSE
		coldata=FALSE
		if (length(type)!=1) stop("Only one type for table output!")
	}
	og=gene
	gene=ToIndex(data,gene)
	og=data$gene.info$Symbol[gene]
	uno=function(type) {
		tno="t"
		spl=strsplit(type,".",fixed=TRUE)[[1]]
		if (length(spl)>1) {tno=spl[1]; type=spl[2];}
		mf = switch(tolower(substr(tno,1,1)),t=1,n=data$data$ntr[gene,conditions],o=1-data$data$ntr[gene,conditions],stop(paste0(type," unknown!"))) 
		if (!keep.ntr.na) {
			mf[is.na(mf)]=if(tolower(substr(tno,1,1))=="n") 0 else 1
		}
		conv=if (type=="count") function(m) {mode(m) <- "integer";m} else function(m) m

		if (!(type %in% names(data$data))) stop(paste0(type," unknown!"))
		if (length(gene)==1) data.frame(conv(data$data[[tolower(type)]][gene,conditions]*mf)) else as.data.frame(conv(t(data$data[[tolower(type)]][gene,conditions]*mf)))
	}
	re=as.data.frame(lapply(type,uno))
	if(length(type)==1 && length(gene)==1) names(re)="Value" else if (length(type)==1) names(re)=og else if (length(gene)==1) names(re)=type else names(re)=paste0(rep(og,length(type))," ",rep(type,each=length(og)))
	if (coldata) re = cbind(data$coldata[conditions,],re)
	if (melt && (length(gene)>1 || length(type)>1)) {
		re = melt(re,id.vars=if(coldata) names(data$coldata) else c(),value.name="Value")
		if (length(type)==1) names(re)[dim(re)[2]-1]="Gene" else if (length(gene)==1) names(re)[dim(re)[2]-1]="Type" else {
			re=cbind(re[,c(1:(dim(re)[2]-2))],setNames(as.data.frame(t(as.data.frame(strsplit(as.character(re$variable)," ")))),c("Gene","Type")),Value=re$Value)
		}
	}
	rownames(re)=NULL
	if (table) {
		names(re)=og
		re=as.data.frame(t(re))
		names(re)=conditions
	}
	re
}

Normalize=function(data,sizeFactors=NULL,name="norm") {
	if (is.null(sizeFactors)) sizeFactors=DESeq2::estimateSizeFactorsForMatrix(data$data$count)
	data$data[[name]] = t(t(data$data$count)/sizeFactors)
	data
}

FilterGenes=function(data,type='tpm',minval=1,mincond=NumCond(data)/2,use=NULL,keep=NULL) {
	if (!is.null(use) & !is.null(keep)) stop("Do not specify both use and keep!")

	if (is.null(use)) {
		t=GetData(data,type=type,table=TRUE)
		use=apply(t,1,function(v) sum(v>=minval,na.rm=TRUE)>=mincond)
		if (!is.null(keep)) use = use | rownames(t) %in% rownames(t[keep,])
	}
	return(data.apply(data,function(t) t[use,]))
}

SwapSamples=function(data,s1,s2) {
	i1=if(is.numeric(s1)) s1 else which(rownames(data$coldata)==s1)
	i2=if(is.numeric(s2)) s2 else which(rownames(data$coldata)==s2)
	return(data.apply(data,function(t) {
		tmp=t[,i1]
		t[,i1]=t[,i2]
		t[,i2]=tmp
		t
	},gene.info=FALSE))
}

TestFilterGenes=function(data,type='tpm',minval='1',mincond=NumCond(data)/2,use=NULL) {
	if (is.null(use)) {
		t=GetData(data,type=type,table=TRUE)
		use=apply(t,1,function(v) sum(v>=minval,na.rm=TRUE)>=mincond)
	}
	return(sum(use))
}






