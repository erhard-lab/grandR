




#' Build the type column for the gene info table.
#'
#' A list of functions. These are called one after the other by \code{\link{ReadGRAND}} or \code{\link{ReadGRAND3}} to build the \emph{Type} column in the \code{\link{GeneInfo}} table. The name of the first function returning TRUE for a row in the table is used as its \emph{type}.
#'
#' @details This is by default given to the \emph{classify.genes} parameter of \code{\link{ReadGRAND}} and \code{\link{ReadGRAND3}}.
#' It will assign the type \emph{mito} if the gene symbol starts with \emph{mt-}, \emph{ERCC} if it starts with \emph{ERCC} and \emph{Cellular} if the gene name starts with \emph{ENS}.
#'
#' @seealso \link{ReadGRAND}, \link{ReadGRAND3}
#' @examples
#' classi <- c(GeneType,
#'            `SARS-CoV-2`=function(gene.info) gene.info$Symbol %in%
#'                     c('ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF10','ORF1ab','S')
#'            )
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate),
#'                   classify.genes=classi,
#'                   verbose=TRUE)
#' table(GeneInfo(sars)$Type)
#'
#' @export
GeneType=list(
	mito=function(gene.info) grepl("^MT-",gene.info$Symbol,ignore.case=TRUE),
	ERCC=function(gene.info) grepl("ERCC-[0-9]{5}",gene.info$Gene),
	Cellular=function(gene.info) grepl("^ENS.*G\\d+$",gene.info$Gene)
)



#' A list of predefined names for design vectors
#'
#' These predefined names mainly are implemented here to harmonize analyses.
#' It is good practise to use these names if sensible.
#'
#' @export
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


#' Add additional columns to the \code{\link{Coldata}} table.
#'
#' Design.Semantics is a list of functions that is supposed to be used as \code{semantics} parameter when calling \code{\link{MakeColdata}}.
#' For each design vector element matching a name of this list the corresponding function is called by \link{MakeColdata} to add additional columns.
#' By default, \code{duration.4sU} is mapped to \code{\link{Semantics.time}}.
#'
#' @name Design.Semantics
#'
#' @param s the original column in the \code{Coldata} table column
#' @param name the name of the column in the \code{Coldata} table column
#'
#' @details
#' \itemize{
#'   \item{\strong{Semantics.noop:} Just add the column as is}
#'   \item{\strong{Semantics.time:} Contents such as 3h or 30min are converted into a numerical value (in hours), and no4sU is converted into 0.}
#'   }
#'
#' @format
#' @seealso \link{MakeColdata}
#' @examples
#'
#' Semantics.time(c("5h","30min","no4sU"),"Test")
#'
#'
#' sema <- list(duration.4sU=
#'      function(s,name) {
#'         r<-Semantics.time(s,name)
#'         cbind(r,data.frame(hpi=paste0(r$duration.4sU+3,"h")))
#'      })
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=function(names) MakeColdata(names,c("Cell",Design$dur.4sU,Design$Replicate),semantics=sema),
#'                   classify.genes=classi,
#'                   verbose=TRUE)
#' Coldata(sars)
#'
NULL

#' @rdname Design.Semantics
#' @export
Semantics.noop=function(s,name) setNames(data.frame(s),name)
#' @rdname Design.Semantics
#' @export
Semantics.time=function(s,name) {
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

#' @rdname Design.Semantics
#' @export
Design.Semantics=list(
  duration.4sU=Semantics.time
)


#' Extract an annotation table from a formatted names vector
#'
#' If columns (i.e. sample or cell) follow a specific naming pattern, this can be used to conveniently set up an annotation table.
#'
#' @param names Formatted names vector (see details)
#' @param design Titles for the columns of the annotation table
#' @param semantics Additional semantics to apply to given annotations (see details)
#' @param rownames Add rownames to the annotation table
#' @param keep.originals To not discard the original values for all annotations where semantics were applied
#'
#' @return A data frame representing the annotation table
#'
#' @details The names have to contain dots (.) to separate the fields for the column annotation table.
#' E.g. the name \emph{Mock.4h.A} will be split into the fields \emph{Mock}, \emph{4h} and  \emph{A}.
#' For such names, a design vector of length 3 has to be given, that describes the meaning of each field.
#' A reasonable design vector for the example would be \code{c("Treatment","Time","Replicate")}.
#' Some names are predefined in the list \link{Design}.
#'
#' @details The names given in the design vector might even have additional semantics:
#' E.g. for the name \emph{duration.4sU} the values are interpreted (e.g. 4h is converted into the number 4, or 30min into 0.5, or no4sU into 0).
#'
#' @details Semantics can be user-defined via the \emph{semantics} list:
#' For each name in the design vector matching to a name in this list, the corresponding function in the list is run.
#' Functions must accept 2 parameters, the first is the original column in the annotation table, the second the original name.
#' The function must return a data.frame with the number of rows matching to the annotation table.
#' In most cases it is easier to manipulate the returned data frame instead of changing the semantics.
#' However, the build-in semantics provide a convenient way to reduce this kind of manipulation in most cases.
#'
#' @export
#'
#' @seealso \link{ReadGRAND},\link{Design.Semantics},\link{Coldata}
#'
#' @examples
#' coldata <- MakeColdata(c("Mock.0h.A","Mock.0h.B","Mock.2h.A","Mock.2h.B"), design=c("Cell",Design$dur.4sU,Design$Replicate))
#'
MakeColdata=function(names,design,semantics=Design.Semantics,rownames=TRUE,keep.originals=TRUE) {
  coldata=data.frame(Name=names,check.names=FALSE,stringsAsFactors = TRUE)
  spl=strsplit(as.character(coldata$Name),".",fixed=TRUE)
  if (any(sapply(spl, length)!=length(design))) stop(paste0("Design parameter is incompatible with input data (e.g., ",paste(coldata$Name[which(sapply(spl, length)!=length(design))[1]]),")"))

  for (i in 1:length(design)) if (!is.na(design[i])) coldata=cbind(coldata,type.convert(sapply(spl,function(v) v[i]),as.is=FALSE))
  names(coldata)[-1]=design[!is.na(design)]
  if (rownames) rownames(coldata)=coldata$Name

  for (sname in names(semantics)) {
    if (sname %in% design) {
      s=coldata[[sname]]
      df=semantics[[sname]](as.character(s),sname)
      if (keep.originals) df=cbind(df,setNames(data.frame(s),paste0(sname,".original")))
      coldata=cbind(coldata[!names(coldata) %in% names(df)],df)
    }
  }

  coldata
}



#' Read the output of GRAND-SLAM 2.0 into a grandR object.
#'
#' Metabolic labeling - nucleotide conversion RNA-seq data (such as generated by SLAM-seq,TimeLapse-seq or TUC-seq)
#' must be carefully analyzed to remove bias due to incomplete labeling. GRAND-SLAM is a software package that
#' employs a binomial mixture modeling approach to obtain precise estimates of the new-to-total RNA ratio (NTR) per gene and sample (or cell).
#' This function directly reads the output of GRAND-SLAM 2.0 into a grandR object.
#'
#' @param prefix Can either be the prefix used to call GRAND-SLAM with, or the main output file ($prefix.tsv.gz);
#' if the RCurl package is installed, this can also be a URL
#' @param design Either a design vector (see details), or a data.frame providing metadata for all columns (samples/cells),
#' or a function that is called with the condition name vector and is supposed to return this data.frame.
#' @param classify.genes A list of functions that is used to add the \emph{type} column to the gene annotation table
#' @param Unknown If no function from \emph{classify.genes} to a row in the gene annotation table, this is used as the \emph{type}
#' @param read.percent.conv Should the percentage of convertions also be read?
#' @param verbose Verbose status outputs
#'
#' @return A grandR object containing the read counts, NTRs, information on the NTR posterior distribution (alpha,beta)
#' and potentially additional information of all genes detected by GRAND-SLAM
#'
#' @details If columns (samples/cells) are named systematically in a particular way, the design vector provides
#' a powerful and easy way to create the column annotations.
#'
#' @details The column names have to contain dots (.) to separate the fields for the column annotation table.
#' E.g. the name \emph{Mock.4h.A} will be split into the fields \emph{Mock}, \emph{4h} and  \emph{A}.
#' For such names, a design vector of length 3 has to be given, that describes the meaning of each field.
#' A reasonable design vector for the example would be \code{c("Treatment","Time","Replicate")}.
#' Some names are predefined in the list \link{Design}.
#'
#' @details The names given in the design vector might even have additional semantics:
#' E.g. for the name \emph{duration.4sU} the values are interpreted (e.g. 4h is converted into the number 4,
#' or 30min into 0.5, or no4sU into 0). Semantics can be user-defined by calling \code{\link{MakeColdata}}
#' and using the return value as the design parameter, or a function that calls MakeColdata.
#' In most cases it is easier to manipulate the \code{\link{Coldata}} table after loading data instead of using this mechanism;
#' the build-in semantics simply provide a convenient way to reduce this kind of manipulation in most cases.
#'
#' @seealso \link{GeneType},\link{MakeColdata},\link{Design.Semantics}
#'
#' @examples
#' sars <- ReadGRAND("https://zenodo.org/record/5834034/files/sars.tsv.gz", design=c("Cell",Design$dur.4sU,Design$Replicate), verbose=TRUE)
#'
#' @export
#'
ReadGRAND=function(prefix, design=c(Design$Condition,Design$Replicate),classify.genes=GeneType,Unknown=NA,read.percent.conv=FALSE, verbose=FALSE) {

  tomat=function(m,names,cnames){
    m=as.matrix(m)
    m[is.na(m)]=0
    colnames(m)=gsub(" .*","",cnames)
    rownames(m)=names
    m
  }

  checknames=function(a,b){
  	if (!all(colnames(a)==colnames(b))) stop("Column names do not match!")
  	if (!all(rownames(a)==rownames(b))) stop("Row names do not match!")
  }
  do.callback=function() {}

  url=NULL
  if (has.package("RCurl")) {
    if (RCurl::url.exists(prefix)) {
      url=prefix
    } else if (RCurl::url.exists(paste0(prefix,".tsv"))) {
      url=paste0(prefix,".tsv")
    } else if (RCurl::url.exists(paste0(prefix,".tsv.gz"))) {
      url=paste0(prefix,".tsv.gz")
    } else {
      url=NULL
    }
    if (!is.null(url)) {
      file <- tempfile()
      if (verbose) cat(sprintf("Downloading file (destination: %s) ...\n",file))
      download.file(url, file, quiet=!verbose)
      prefix=gsub(".tsv(.gz)?$","",url)
      do.callback=function() {
        if (verbose) cat("Deleting temporary file...\n")
        unlink(file)
      }
    }
  }

  if (is.null(url)) {
  	file=if (file.exists(prefix)) prefix else paste0(prefix,".tsv")
  	if (!file.exists(file) && file.exists(paste0(file,".gz"))) file = paste0(file,".gz")
  	prefix=gsub(".tsv(.gz)?$","",file)
  }

  if (!file.exists(file)) stop("File not found; If you want to access non-local files directly, please install the RCurl package!")


	if (verbose) cat("Checking file...\n")
	con <- file(file,"r")
	header <- strsplit(readLines(con,n=1),"\t")[[1]]
	close(con)

	if (header[1]!="Gene" || header[2]!="Symbol" || !grepl("Readcount$",header[3])) stop("File is not a GRAND-SLAM output file!")
	conds=gsub(" Readcount","",header[grepl(" Readcount",header)])
	terms=strsplit(conds[1],".",fixed=TRUE)[[1]]

	if (is.data.frame(design)) {
	  if (length(conds)!=nrow(design)) stop(paste0("Design parameter (table) is incompatible with input data: ",paste(terms,collapse=".")))
	  coldata=design
	} else if (is.function(design)) {
	  coldata=design(conds)
	  if (length(conds)!=nrow(coldata)) stop(paste0("Design parameter (function) is incompatible with input data: ",paste(terms,collapse=".")))
	} else {
  	if (length(terms)!=length(design)) stop(paste0("Design parameter is incompatible with input data: ",paste(terms,collapse=".")))
  	coldata=MakeColdata(conds,design)
	}


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
	if (verbose) cat("Processing...\n")

	ntr.mean = grepl("Mean",names(data))
	ntr = grepl("MAP",names(data))
	alpha = grepl("alpha",names(data))
	beta = grepl("beta",names(data))
	count = grepl("Readcount",names(data))
	conv = grepl("Conversions",names(data))
	cove = grepl("Coverage",names(data)) & !grepl("Double-Hit Coverage",names(data))


	classify.genes=c(classify.genes,Unknown=function(gene.info) rep(T,dim(gene.info)[1]))
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

	coldata$no4sU=no4sU.cols

	do.callback()

	# insert name no rep and set this as default condition, if there is no condition field!
	re=grandR(prefix=prefix,gene.info=gene.info,slots=re,coldata=coldata,metadata=list(Description="GRAND-SLAM 2.0 data"))
	DefaultSlot(re)="count"
	re
}





# TODO: do not use read.delim or .csv, nor file.exists (to allow using urls!)
isSparse=function(prefix) !file.exists(paste0(prefix,".targets/data.tsv.gz"))

ReadGRAND3=function(prefix,...) if (isSparse(prefix)) ReadGRAND3_sparse(prefix,...) else ReadGRAND3_dense(prefix,...)

ReadGRAND3_sparse=function(prefix, verbose=FALSE, design=c(Design$Library,Design$Sample,Design$Barcode), label="4sU",estimator="Binom", read.CI=FALSE) {

  cols=readLines(paste0(prefix, ".targets/barcodes.tsv.gz"))
  conds=strsplit(cols,".",fixed=TRUE)[[1]]
  if (length(conds)!=length(design)) stop(paste0("Design parameter is incompatible with input data: ",paste(conds,collapse=".")))


  if (verbose) cat("Reading count matrix...\n")
  count=as(readMM(paste0(prefix, ".targets/matrix.mtx.gz")),Class = "dgCMatrix")
  gene.info=read.delim(paste0(prefix, ".targets/features.tsv.gz"),header = FALSE,stringsAsFactors = FALSE)
  if (ncol(gene.info)==4) gene.info$Length=1
  gene.info=setNames(gene.info,c("Gene","Symbol","Mode","Category","Length"))

  if (anyDuplicated(gene.info$Gene)) {
    dupp=table(gene.info$Gene)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Gene=make.unique(gene.info$Gene)
  }
  if (anyDuplicated(gene.info$Symbol)) {
    dupp=table(gene.info$Symbol)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene symbols (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Symbol=make.unique(gene.info$Symbol)
  }

  if (any(make.names(gene.info$Symbol)!=gene.info$Symbol)) {
    ill=gene.info$Symbol[which(make.names(gene.info$Symbol)!=gene.info$Symbol)]
    warning(sprintf("Illegal identifiers among the gene symbols (e.g. %s), making legal!",paste(head(ill),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Symbol=make.names(gene.info$Symbol)
  }

  colnames(count)=cols
  rownames(count)=gene.info$Symbol
  re=list()
  re$count=count

  if (verbose) cat("Reading NTRs...\n")
  ntr=as(readMM(sprintf("%s.targets/%s.%s.ntr.mtx.gz",prefix,label,estimator)),Class = "dgCMatrix")
  colnames(ntr)=cols
  rownames(ntr)=gene.info$Symbol
  re$ntr=ntr

  if (read.CI && file.exists(sprintf("%s.targets/%s.%s.lower.mtx.gz",prefix,label,estimator)) && file.exists(sprintf("%s.targets/%s.%s.upper.mtx.gz",prefix,label,estimator))) {
    if (verbose) cat("Reading CIs...\n")
    lower=as(readMM(sprintf("%s.targets/%s.%s.lower.mtx.gz",prefix,label,estimator)),Class = "dgCMatrix")
    colnames(lower)=cols
    rownames(lower)=gene.info$Symbol
    re$lower=lower

    upper=as(readMM(sprintf("%s.targets/%s.%s.upper.mtx.gz",prefix,label,estimator)),Class = "dgCMatrix")
    colnames(upper)=cols
    rownames(upper)=gene.info$Symbol
    re$upper=upper
  }
  gene.info$Type=gsub(".*\\(","",gsub(")","",gene.info$Category,fixed=TRUE))
  gene.info$Type=factor(gene.info$Type,levels=unique(gene.info$Type))


  # use make coldata instead. add no4sU column!
  coldata=data.frame(Name=cols)
  spl=strsplit(as.character(coldata$Name),".",fixed=TRUE)
  for (i in 1:length(design)) coldata=cbind(coldata,factor(sapply(spl,function(v) v[i]),levels=unique(sapply(spl,function(v) v[i]))))
  names(coldata)[-1]=design
  rownames(coldata)=coldata$Name

  grandR(prefix,gene.info,slots=re,coldata)
}

ReadGRAND3_dense=function(prefix, verbose=FALSE, design=c(design$Condition,Design$Replicate)) {

}


ReadNewTotal=function(genes, cells, new.matrix, total.matrix, detection.rate=1,verbose=FALSE) {

  cols=read.csv(cells,check.names = FALSE,stringsAsFactors = FALSE)
  gene.info=setNames(read.csv(genes,check.names = FALSE,stringsAsFactors = FALSE),c("Gene","Biotype","Symbol"))

  if (verbose) cat("Reading total count matrix...\n")
  count=as(readMM(total.matrix),Class = "dgCMatrix")

  if (anyDuplicated(gene.info$Gene)) {
    dupp=table(gene.info$Gene)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Gene=make.unique(gene.info$Gene)
  }
  if (anyDuplicated(gene.info$Symbol)) {
    dupp=table(gene.info$Symbol)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene symbols (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Symbol=make.unique(gene.info$Symbol)
  }

  if (any(make.names(gene.info$Symbol)!=gene.info$Symbol)) {
    ill=gene.info$Symbol[which(make.names(gene.info$Symbol)!=gene.info$Symbol)]
    warning(sprintf("Illegal identifiers among the gene symbols (e.g. %s), making legal!",paste(head(ill),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Symbol=make.names(gene.info$Symbol)
  }

  colnames(count)=cols[,1]
  rownames(count)=gene.info$Symbol

  if (verbose) cat("Reading new count matrix...\n")
  new=as(readMM(new.matrix),Class = "dgCMatrix")

  if (verbose) cat("Computing NTRs...\n")
  new=new/detection.rate
  ntr=new/count
  ntr@x[ntr@x>1]=1
  ntr@x[is.nan(ntr@x)]=0
  ntr=as(ntr,Class = "dgCMatrix")

  colnames(ntr)=colnames(count)
  rownames(ntr)=rownames(count)

  re=list()
  re$count=count
  re$ntr=ntr

  grandR("",gene.info,slots=re,cols)
}


