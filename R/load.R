
library(ggplot2)
library(cowplot)
library(reshape2)
library(plyr)
library(lfc)
library(rclipboard)

comp.hl=function(p,time=1) ifelse(p==0,Inf,log(2)/(-1/time*log(1-p)))

comp.tpm=function(cmat,lengths,subset=NULL) {
	zerolen=lengths==0
	lengths[zerolen]=1
	rpk=cmat/(lengths/1000)
	rpk[zerolen,]=NA
	scale=colSums(if(!is.null(subset)) rpk[subset,] else rpk,na.rm=T)/1E6
	re=t(t(rpk)/scale)
	re
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


semantics.noop=function(s,name) setNames(data.frame(s),name)
semantics.time=function(s,name) {
  time=rep(NA,length(s))
  
  no4sU=c("no4sU","-")
  time[s %in% no4sU]=0
  
  h=grepl("[0-9.]+h(p.)?",s)
  time[h]=as.numeric(gsub("([0-9.]+)h(p.)?","\\1",s[h]))
  
  min=grepl("[0-9.]+min",s)
  time[min]=as.numeric(substr(s[min],1,nchar(s[min])-3))/60
  
  if (any(is.na(time))) stop(paste0("Time semantics cannot be used for this: ",paste(s[is.na(time)],collapse=",")))
  setNames(data.frame(time),name)
}


Design=list(
  has.4sU="has.4sU",
  conc.4sU="concentration.4sU",
  dur.4sU="duration.4sU",
  Replicate="Replicate",
  Condition="Condition",
  hpi="hpi",
  hps="hps",
  Library="Library",
  Sample="Sample",
  Barcode="Barcode",
  Origin="Origin"
  )
# Sample is without replicate! Condition might have multiple samples, but not the other way around; Condition is used all over the place (e.g. kinetic fitting) to identify
Design.Semantics=list(
  hpi=semantics.time,
  hps=semantics.time,
  dur.4sU=semantics.time
)
GetField=function(name,field,sep=".") sapply(strsplit(as.character(name),sep,fixed=TRUE),function(v) v[field])

ConvFields=function(v) {
  if (is.data.frame(v) || is.matrix(v)) {
    re=as.data.frame(lapply(as.data.frame(v),ConvFields))
    colnames(re)=colnames(v)
    rownames(re)=rownames(v)
    return(re)
  }
  v=as.character(v)
  if (sum(is.na(suppressWarnings(as.logical(v))))==0) {
    as.logical(v)
  } else if (sum(is.na(suppressWarnings(as.integer(v))))==0) {
    as.integer(v)
  } else if (sum(is.na(suppressWarnings(as.double(v))))==0) {
    as.double(v)
  } else {
    factor(v,levels=unique(v))
  }
}


MakeColdata=function(names,design,semantics=Design.Semantics,rownames=TRUE) {
  coldata=data.frame(Name=names,check.names=FALSE,stringsAsFactors = TRUE)
  spl=strsplit(as.character(coldata$Name),".",fixed=TRUE)
  if (any(sapply(spl, length)!=length(design))) stop(paste0("Design parameter is incompatible with input data (e.g., ",paste(coldata$Name[which(sapply(spl, length)!=length(design))[1]]),")"))
  
  for (i in 1:length(design)) if (!is.na(design[i])) coldata=cbind(coldata,ConvFields(sapply(spl,function(v) v[i])))
  names(coldata)[-1]=design[!is.na(design)]
  if (rownames) rownames(coldata)=coldata$Name
  coldata$Sample=interaction(coldata[ !(names(coldata) %in% c("Name","Replicate"))],drop=TRUE)
  
  for (sname in names(semantics)) {
    if (sname %in% design) {
      s=as.character(coldata[[sname]])
      df=semantics[[sname]](s,sname)
      coldata=cbind(coldata[!names(coldata) %in% names(df)],df)
    }
  }
  
  
  coldata
}


ReadGRAND=function(prefix, verbose=FALSE, classify.genes=GeneType,design=c(Design$Condition,Design$Replicate),Unknown=NA,rename.samples=NULL,read.percent.conv=FALSE) {


checknames=function(a,b){
	if (!all(colnames(a)==colnames(b))) stop("Column names do not match!")
	if (!all(rownames(a)==rownames(b))) stop("Row names do not match!")
}


	file=if (file.exists(prefix)) prefix else paste0(prefix,".tsv")
	if (!file.exists(file) && file.exists(paste0(file,".gz"))) file = paste0(file,".gz")
	prefix=gsub(".tsv(.gz)?$","",file)
	
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
	  dupp=table(data$Gene)
	  dupp=names(dupp)[which(dupp>1)]
	  warning(sprintf("Duplicate gene names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
		data$Gene=make.unique(data$Gene)
	}
	if (anyDuplicated(data$Symbol)) {
	  dupp=table(data$Symbol)
	  dupp=names(dupp)[which(dupp>1)]
	  warning(sprintf("Duplicate gene symbols (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
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
	conv = grepl("Conversions",names(data))
	cove = grepl("Coverage",names(data)) & !grepl("Double-Hit Coverage",names(data))
	
	if (!is.na(Unknown)) names(classify.genes)[names(classify.genes)=="Unknown"]=Unknown
	gene.info = data.frame(Gene=as.character(data$Gene),Symbol=as.character(data$Symbol),Length=data$Length,stringsAsFactors=FALSE)
	gene.info$Type=NA
	for (i in length(classify.genes):1) gene.info$Type[classify.genes[[i]](gene.info)]=names(classify.genes)[i]
	gene.info$Type=factor(gene.info$Type,levels=names(classify.genes))

	re=list()
	re$count=tomat(data[,count],data$Gene,names(data)[count])
#	re$tpm=comp.tpm(re$count,gene.info$Length)
#	checknames(re$count,re$tpm)
	#re$ntr.mean=tomat(data[,alpha]/(data[,alpha]+data[,beta]),data$Gene,names(data)[ntr.mean])
	re$ntr=tomat(data[,ntr],data$Gene,names(data)[ntr])
	re$alpha=tomat(data[,alpha],data$Gene,names(data)[alpha])
	re$beta=tomat(data[,beta],data$Gene,names(data)[beta])
	
	if (read.percent.conv) {
	  re$percent_conv=tomat(data[,conv]/data[,cove],data$Gene,names(data)[conv])
	}
	
	no4sU.cols=!(colnames(re$count) %in% colnames(re$ntr))
	# remove no4sU conditions from ntr relevant matrices
	if (sum(no4sU.cols)>0) {
		correctmat=function(m) {
			r=matrix(NA,ncol=dim(re$count)[2],nrow=dim(re$count)[1])
			rownames(r)=rownames(re$count)
			colnames(r)=colnames(re$count)
			r[,colnames(r)[!no4sU.cols]]=m[,colnames(r)[!no4sU.cols]]
			r
		}
		#re$ntr.mean=correctmat(re$ntr.mean)
		re$ntr=correctmat(re$ntr)
		re$alpha=correctmat(re$alpha)
		re$beta=correctmat(re$beta)
	}

	checknames(re$count,re$ntr)
	checknames(re$count,re$ntr.mean)
	checknames(re$count,re$alpha)
	checknames(re$count,re$beta)
	
	#coldata=data.frame(Name=colnames(re$count))
	#spl=strsplit(as.character(coldata$Name),".",fixed=TRUE)
	#for (i in 1:length(design)) coldata=cbind(coldata,factor(sapply(spl,function(v) v[i]),levels=unique(sapply(spl,function(v) v[i]))))
	#names(coldata)[-1]=design
	#rownames(coldata)=coldata$Name
	#coldata$Sample=interaction(coldata[ !(names(coldata) %in% c("Name","Replicate"))],drop=TRUE)
	coldata=MakeColdata(colnames(re$count),design)
	coldata$no4sU=no4sU.cols

	re=grandR(prefix=prefix,gene.info=gene.info,data=re,coldata=coldata,metadata=list())
	DefaultSlot(re)="count"
	invisible(re)
}


grandR=function(prefix=parent$prefix,gene.info=parent$gene.info,data=parent$gene.info,coldata=parent$coldata,metadata=parent$metadata,parent=NULL) {
	info=list()
	info$prefix=prefix
	info$gene.info=gene.info
	info$data=data
	info$coldata=coldata
	info$metadata=metadata
	class(info)="grandR"
	invisible(info)
}

VersionString=function(data) {
  "grandR v0.1.0"
}

Title=function(data) {
  x=strsplit(data$prefix,"/")[[1]]
  x[length(x)]
}

`DefaultSlot<-` <- function(data, value) {
  data$metadata$default.slot=value
  data
}
DefaultSlot <- function(data) {
  data$metadata$default.slot
}

Slots=function(data) {
  names(data$data)
}

Genes=function(data, use.symbols=TRUE) data$gene.info[[if (use.symbols) "Symbol" else "Gene"]]

AddSlot=function(data,name,matrix) {
  if (!all(colnames(matrix)==colnames(data$data$count))) stop("Column names do not match!")
  if (!all(rownames(matrix)==rownames(data$data$count))) stop("Row names do not match!")
  data$data[[name]]=matrix
  data
}

dim.grandR=function(data) c(dim(data$gene.info)[1],dim(data$coldata)[1])
is.grandR <- function(x) inherits(x, "grandR")
dimnames.grandR=function(data) dimnames(data$data$count)

data.apply=function(data,fun,fun.gene.info=NULL,fun.coldata=NULL,...) {
	re=list()
	for (l1 in names(data$data)) {
		re[[l1]]=fun(data$data[[l1]],...)
	}
	ngene.info=if (!is.null(fun.gene.info)) fun.gene.info(data$gene.info,...) else data$gene.info
	ncoldata=if (!is.null(fun.coldata)) fun.coldata(data$coldata,...) else data$coldata
	invisible(grandR(parent=data,gene.info=ngene.info,data=re,coldata=ncoldata))
}

reorder.grandR=function(data,columns) {
  r=subset.grandR(data,columns)
  r$coldata=ConvFields(r$coldata)
  r
}
subset.grandR=function(data,columns) {
  keep=rownames(data$coldata)[columns]
  data.apply(data,function(m) m[,intersect(keep,colnames(m))],fun.coldata = function(t) droplevels(t[columns,]))
}

split.grandR=function(data,column.name=Design$Condition) {
  col=as.factor(data$coldata[[column.name]])
  setNames(lapply(levels(col),function(c) {re=subset(data,col==c); re$coldata[[Design$Origin]]=c; re }),levels(col))
}

`Condition<-` <- function(data, value) {
  data$coldata$Condition <- interaction(data$coldata[value])
  data
}

RenameSamples=function(data,map=NULL,fun=NULL) {
  if (!is.null(fun)) {
    map=setNames(sapply(colnames(data),fun),colnames(data))
  }
  names=rownames(data$coldata)
  names[names %in% names(map)]=unlist(map[names[names %in% names(map)]])
  rownames(data$coldata)=names
  data$coldata$Name=names
  data.apply(data,function(m) {colnames(m)=names; m})
}

merge.grandR=function(...,list=NULL,column.name=Design$Origin) {
  list=c(list(...),list)
  if (length(list)==1) return(list[[1]])
  
  re=list[[1]]
  if (!is.null(names(list))) re$coldata[[column.name]]=names(list)[1]
  for (i in 2:length(list)) {
    add=list[[i]]
    if (!is.null(names(list))) add$coldata[[column.name]]=names(list)[i]
    if (any(colnames(add) %in% colnames(re))) stop("Sample names must be unique!")
    if (any(rownames(add)!=rownames(re))) stop("Data sets must have the same genes!")
    if (any(colnames(add$coldata)!=colnames(re$coldata))) stop("Data sets must have the coldata columns!")
    if (!all(names(add$data) %in% names(re$data))) stop("Data sets must have the same data tables!")
    re$coldata=rbind(re$coldata,add$coldata)
    
    for (n in names(re$data)) re$data[[n]]=cbind(re$data[[n]],add$data[[n]])
    
  }
  re  
}

ComputeNtrCI=function(data,CI.size=0.95,name.lower="lower",name.upper="upper") {
  data=ComputeNtrPosteriorLower(data=data,CI.size=CI.size,name=name.lower)
  data=ComputeNtrPosteriorUpper(data=data,CI.size=CI.size,name=name.upper)
  data
}
ComputeNtrPosteriorLower=function(data,CI.size=0.95,name="lower") ComputeNtrPosteriorQuantile(data=data,quantile=(1-CI.size)/2,name=name)
ComputeNtrPosteriorUpper=function(data,CI.size=0.95,name="upper") ComputeNtrPosteriorQuantile(data=data,quantile=1-(1-CI.size)/2,name=name)

ComputeNtrPosteriorQuantile=function(data,quantile,name) {
  a=as.matrix(GetTable(data,type="alpha",name.by = "Gene"))
  b=as.matrix(GetTable(data,type="beta",name.by = "Gene"))
  v=qbeta(quantile,a,b)
  AddSlot(data,name,v)
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
  if (is.logical(gene) && length(gene)==nrow(data)) return(which(gene))
  if (all(gene %in% data$gene.info$Gene)) return(setNames(1:nrow(data),data$gene.info$Gene)[gene])
  if (all(gene %in% data$gene.info$Symbol)) return(setNames(1:nrow(data),data$gene.info$Symbol)[gene])
	stop("Could not find all genes!")
}



GetTable=function(data,type,subset=NULL,gene=Genes(data),keep.ntr.na=TRUE,summarize=NULL,prefix=NULL,name.by="Symbol") {
  r=GetData(data,type,subset=subset,gene,keep.ntr.na = keep.ntr.na,coldata=FALSE, table=TRUE,name.by = name.by)
  if (!is.null(summarize)) {
    if (is.logical(summarize) && length(summarize)==1) summarize=GetSummarizeMatrix(data)
    r=as.data.frame(as.matrix(r) %*% summarize)
  }
  if (!is.null(prefix)) colnames(r)=paste0(prefix,colnames(r))
  r
}
GetData=function(data,type=DefaultSlot(data),subset=NULL,gene=Genes(data),melt=FALSE,coldata=TRUE,table=FALSE,keep.ntr.na=TRUE,name.by="Symbol") {
	if (table) { 
		melt=FALSE
		coldata=FALSE
		if (length(type)!=1) stop("Only one type for table output!")
	}
  if (is.null(subset)) subset=colnames(data)
	og=gene
	gene=ToIndex(data,gene)
	og=data$gene.info[[name.by]][gene]
	uno=function(type) {
		tno="t"
		spl=strsplit(type,".",fixed=TRUE)[[1]]
		if (length(spl)>1) {tno=spl[1]; type=spl[2];}
		mf = switch(tolower(substr(tno,1,1)),t=1,n=data$data$ntr[gene,subset],o=1-data$data$ntr[gene,subset],stop(paste0(type," unknown!"))) 
		if (!keep.ntr.na) {
			mf[is.na(mf)]=if(tolower(substr(tno,1,1))=="n") 0 else 1
		}
		conv=if (type=="count") function(m) {mode(m) <- "integer";m} else if (type=="ntr" && !keep.ntr.na) function(m) {m[is.na(m)]=0; m} else function(m) m

		if (!(type %in% names(data$data))) stop(paste0(type," unknown!"))
		if (length(gene)==1) data.frame(conv(data$data[[type]][gene,subset]*mf)) else as.data.frame(conv(t(data$data[[type]][gene,subset]*mf)))
	}
	re=as.data.frame(lapply(type,uno))
	if(length(type)==1 && length(gene)==1) names(re)="Value" else if (length(type)==1) names(re)=og else if (length(gene)==1) names(re)=type else names(re)=paste0(rep(og,length(type))," ",rep(type,each=length(og)))
	if (coldata) re = cbind(data$coldata[subset,],re)
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
		names(re)=data$coldata[subset,"Name"]
	}
	re
}

Normalize=function(data,sizeFactors=NULL,genes=NULL,name="norm",slot="count",return.sf=FALSE,set.to.default=TRUE) {
  if (is.null(genes)) genes=1:nrow(data)
  mat=as.matrix(GetTable(data,slot,gene=genes,keep.ntr.na = FALSE,name.by = "Gene"))
	if (is.null(sizeFactors)) sizeFactors=DESeq2::estimateSizeFactorsForMatrix(mat)
  if (return.sf) return(sizeFactors)
  data=AddSlot(data,name,t(t(mat)/sizeFactors))
	if (set.to.default) DefaultSlot(data)=name
	data
}
NormalizeTPM=function(data,tlen=data$gene.info$Length,name="tpm",slot="count",set.to.default=TRUE) {
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,slot,keep.ntr.na = FALSE,name.by = "Gene"))
  data=AddSlot(data,name,comp.tpm(mat,tlen))
  if (set.to.default) DefaultSlot(data)=name
  data
}

NormalizeBaseline=function(data,baseline=FindReferences(data,reference=Condition==levels(Condition)[1]),name="baseline",slot=DefaultSlot(data),set.to.default=TRUE,LFC.fun=PsiLFC,...) {
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,slot,keep.ntr.na = FALSE,name.by = "Gene"))
  mat=sapply(names(baseline),function(n) LFC.fun(mat[,n],rowMeans(mat[,baseline[[n]]]),...))
  data=AddSlot(data,name,mat)
  if (set.to.default) DefaultSlot(data)=name
  data
}


FindReferences=function(data,reference, covariate=NULL) {
  stopifnot(is.grandR(data))
  df=data$coldata
  df$group=if(is.null(covariate)) 1 else interaction(df[covariate],drop=FALSE,sep=".")
  e=substitute(reference)
  map=dlply(df,.(group),function(s) as.character(s$Name[eval(e,s,parent.frame())]))
  pairs=setNames(lapply(df$group,function(g) map[[g]]),df$Name)
  pairs
}


# minsamp: GetSummarizeMatrix(d,subset=NULL,average=FALSE), if all replicates of a sample exceed minval
FilterGenes=function(data,type='count',minval=100,mincond=ncol(data)/2,minsamp=NULL,use=NULL,keep=NULL,return.genes=FALSE) {
	if (!is.null(use) & !is.null(keep)) stop("Do not specify both use and keep!")

	if (is.null(use)) {
		t=GetData(data,type=type,table=TRUE)
		if (!is.null(minsamp)) {
		  m=GetSummarizeMatrix(d,subset=NULL,average=FALSE)
		  use=rowSums(sapply(1:ncol(m),function(i) apply(t[,m[,i]>0]>=minval,1,all)))>=minsamp
		} else {
      use=apply(t,1,function(v) sum(v>=minval,na.rm=TRUE)>=mincond)
	  }
		if (!is.null(keep)) use = use | rownames(t) %in% rownames(t[keep,])
	} else {
	  use=ToIndex(data,use)
	}
  
  if (return.genes) return(data$gene.info$Gene[use])
	return(data.apply(data,function(t) t[use,],fun.gene.info = function(t) t[use,]))
}

SwapSamples=function(data,s1,s2) {
	i1=if(is.numeric(s1)) s1 else which(rownames(data$coldata)==s1)
	i2=if(is.numeric(s2)) s2 else which(rownames(data$coldata)==s2)
	return(data.apply(data,function(t) {
		tmp=t[,i1]
		t[,i1]=t[,i2]
		t[,i2]=tmp
		t
	}))
}



