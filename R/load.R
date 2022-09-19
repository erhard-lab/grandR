




#' Build the type column for the gene info table.
#'
#' Returns a function to be used as \code{classify.genes} parameter for \code{\link{ReadGRAND}}.
#'
#' @param use.default if TRUE, use the default type inference (priority after the user defined ones); see details
#' @param drop.levels if TRUE, drop unused types from the factor that is generated
#' @param name.unknown the type to be used for all genes where no type was identified
#' @param ... additional functions to define types (see details)
#'
#' @details This function returns a function. Usually, you do not use it yourself but \code{ClassifyGenes} is usually as \code{classify.genes} parameter
#' for  \code{\link{ReadGRAND}} to build the \emph{Type} column in the \code{\link{GeneInfo}} table. See the example
#' to see how to use it directly.
#'
#' @details Each ... parameter must be a function that receives the gene info table and must return a logical vector, indicating for each row
#' in the gene info table, whether it matches to a specific type. The name of the parameter is used as the type name.
#'
#' @details If a gene matches to multiple type, the first function returning TRUE for a row in the table is used.
#'
#' @details By default, this function will recognize mitochondrial genes (MT prefix of the gene symbol), ERCC spike-ins,
#' and Ensembl gene identifiers (which it will call "cellular"). These three are the last functions to be checked (in case a user defined type via ...) also
#' matches to, e.g., an Ensembl gene).
#'
#' @return a function that takes the original \link{GeneInfo} table and adds the Type column
#'
#' @seealso \link{ReadGRAND}
#' @examples
#'
#' viral.genes <- c('ORF3a','E','M','ORF6','ORF7a','ORF7b','ORF8','N','ORF10','ORF1ab','S')
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate),
#'                   classify.genes=ClassifyGenes(`SARS-CoV-2`=
#'                              function(gene.info) gene.info$Symbol %in% viral.genes),
#'                   verbose=TRUE)
#' table(GeneInfo(sars)$Type)
#'
#' fun<-ClassifyGenes(viral=function(gene.info) gene.info$Symbol %in% viral.genes)
#' table(fun(GeneInfo(sars)))
#'
#' @export
#'
#' @concept load
ClassifyGenes=function(...,use.default=TRUE, drop.levels=TRUE, name.unknown="Unknown") {
  ll=list(...)
  ll=c(ll,list(
  	mito=function(gene.info) grepl("^MT-",gene.info$Symbol,ignore.case=TRUE),
  	ERCC=function(gene.info) grepl("ERCC-[0-9]{5}",gene.info$Gene),
  	Cellular=function(gene.info) grepl("^ENS.*G\\d+$",gene.info$Gene)
  ))
  ll=c(ll,setNames(list(function(gene.info) rep(T,dim(gene.info)[1])),name.unknown))

  function(gene.info) {
    gene.info$Type=NA
    for (i in length(ll):1) gene.info$Type[ll[[i]](gene.info)]=names(ll)[i]
    re=factor(gene.info$Type,levels=names(ll))
    if (drop.levels) droplevels(re) else re
  }
}



#' A list of predefined names for design vectors
#'
#' These predefined names mainly are implemented here to harmonize analyses.
#' It is good practise to use these names if sensible.
#'
#' @export
#' @concept load
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


#' Build the design semantics list
#'
#' This is used to add additional columns to the \code{\link{Coldata}} table by giving additional semantics to existing columns.
#'
#' @param ... named parameter list of functions (see details)
#'
#' @details DesignSemantics returns a list of functions that is supposed to be used as \code{semantics} parameter when calling \code{\link{MakeColdata}}.
#' For each design vector element matching a name of this list the corresponding function is called by \link{MakeColdata} to add additional columns.
#'
#' @details Each function takes two parameters, the first being the original column in the \code{Coldata} table column, the second being its name.
#'
#' @details Semantics.time is such a predefined function: Contents such as 3h or 30min are converted into a numerical value (in hours), and no4sU is converted into 0.
#'
#' @details By default, this is used for the names duration.4sU and Experimental.time
#'
#' @return a named list; the names should correspond to column names in the \link{Coldata} table,
#' and the values are functions to add semantics to this table
#'
#' @seealso \link{MakeColdata}
#' @examples
#'
#' Semantics.time(c("5h","30min","no4sU"),"Test")
#'
#'
#' myfun <- function(s,name) {
#'         r<-Semantics.time(s,name)
#'         cbind(r,data.frame(hpi=paste0(r$duration.4sU+3,"h")))
#' }
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=function(names)
#'                     MakeColdata(names,c("Cell",Design$dur.4sU,Design$Replicate),
#'                   semantics=DesignSemantics(duration.4sU=myfun)),
#'                   verbose=TRUE)
#' Coldata(sars)
#'
#' @export
#' @concept load
DesignSemantics=function(...) {
  ll=list(...)
  if (!"duration.4sU" %in% names(ll)) ll=c(ll,list(duration.4sU=Semantics.time))
  if (!"Experimental.time" %in% names(ll)) ll=c(ll,list(Experimental.time=Semantics.time))
  ll
}


#' Semantics for time columns
#'
#' Defines additional semantics for columns representing temporal dimensions
#'
#' @param s original column
#' @param name the column name
#'
#' @return a data frame with a single numeric column, where <x>h from s is replaced by x, <x>min is replaced by
#' x/60, and no4sU is replaced by 0
#'
#' @export
#'
#' @concept load
Semantics.time=function(s,name) {
  time=rep(NA,length(s))

  no4sU=c("nos4U","no4sU","-")
  time[s %in% no4sU]=0
  s=gsub("_",".",s)

  h=grepl("[0-9.]+h(p.)?",s)
  time[h]=as.numeric(gsub("([0-9.]+)h(p.)?","\\1",s[h]))

  min=grepl("[0-9.]+min",s)
  time[min]=as.numeric(substr(s[min],1,nchar(s[min])-3))/60

  if (any(is.na(time))) stop(paste0("Time semantics cannot be used for this: ",paste(s[is.na(time)],collapse=",")))
  setNames(data.frame(time),name)
}


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
#' @seealso \link{ReadGRAND},\link{DesignSemantics},\link{Coldata}
#'
#' @examples
#' coldata <- MakeColdata(c("Mock.0h.A","Mock.0h.B","Mock.2h.A","Mock.2h.B"),
#'                                    design=c("Cell",Design$dur.4sU,Design$Replicate))
#'
#' @concept load
MakeColdata=function(names,design,semantics=DesignSemantics(),rownames=TRUE,keep.originals=TRUE) {
  coldata=data.frame(Name=factor(names,levels=unique(names)),check.names=FALSE,stringsAsFactors = FALSE)
  spl=strsplit(as.character(coldata$Name),".",fixed=TRUE)
  if (any(sapply(spl, length)!=length(design))) stop(paste0("Design parameter is incompatible with input data (e.g., ",paste(coldata$Name[which(sapply(spl, length)!=length(design))[1]]),")"))

  str2fac=function(s) if (is.character(s)) factor(s,levels=unique(s)) else s
  for (i in 1:length(design)) if (!is.na(design[i])) coldata=cbind(coldata,str2fac(type.convert(sapply(spl,function(v) v[i]),as.is=TRUE,na.strings=c("DKSALDJLKSADLKSJDLKSJDLJDA")))) # stupid function!
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
#' @param classify.genes A function that is used to add the \emph{type} column to the gene annotation table, always a call to \link{ClassifyGenes}
#' @param read.percent.conv Should the percentage of conversions also be read?
#' @param verbose Print status updates
#' @param rename.sample function that is applied to each sample name before parsing (or NULL)
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
#' @details Sometimes you might have forgotten to name all samples consistently (or you simply messed something up).
#' In this case, the rename.sample parameter can be handy (e.g. to rename a particular misnamed sample).
#'
#' @seealso \link{ReadGRAND3},\link{ClassifyGenes},\link{MakeColdata},\link{DesignSemantics}
#'
#' @examples
#' \donttest{
#' sars <- ReadGRAND("https://zenodo.org/record/5834034/files/sars.tsv.gz",
#'                       design=c("Cell",Design$dur.4sU,Design$Replicate), verbose=TRUE)
#' }
#'
#' @export
#'
#' @concept load
ReadGRAND=function(prefix,
                   design=c(Design$Condition,Design$Replicate),
                   classify.genes=ClassifyGenes(),
                   read.percent.conv=FALSE,
                   rename.sample=NULL,
                   verbose=FALSE) {
  annotations=c("Gene","Symbol","Length")
  slots=c(`Readcount`="count",MAP="ntr",alpha="alpha",beta="beta")
  if (read.percent.conv) slots=c(slots,`Conversions`="conv",Coverage="cove")
  re=read.grand.internal(prefix = prefix, design = design, slots=slots, annotations=annotations,classify.genes = classify.genes,verbose = verbose,rename.sample = rename.sample)
  if (read.percent.conv) {
    re=AddSlot(re,"percent_conv",re$data$conv/re$data$cove)
    re=DropSlot(re,"conv|cove")
  }
  re
}



ReadCounts=function(file, design=c(Design$Condition,Design$Replicate),classify.genes=ClassifyGenes(),verbose=FALSE,sep="\t") {

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

  prefix=file
  url=NULL

  if (suppressWarnings(requireNamespace("RCurl",quietly = TRUE))) {
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
  header <- strsplit(readLines(con,n=1),sep)[[1]]
  close(con)

  terms=strsplit(header[length(header)],".",fixed=TRUE)[[1]]

  if (length(terms)!=length(design)) stop(paste0("Design parameter is incompatible with input data: ",paste(terms,collapse=".")))


  if (verbose) cat("Reading file...\n")
  data=read.table(file,sep=sep,stringsAsFactors=FALSE,check.names=FALSE,header=TRUE)
  clss=sapply(data,class)
  firstnumeric=min(which(clss=="numeric"))
  if (!all(clss[firstnumeric:length(clss)]=="numeric") || firstnumeric==1) stop("Columns (except for the first n) must be numeric!")
  anno.names=colnames(data)[1:(firstnumeric-1)]

  if (anyDuplicated(data[[1]])) {
    dupp=table(data[[1]])
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    data[[1]]=make.unique(data[[1]])
  }
  if (verbose) cat("Processing...\n")

  coldata=MakeColdata(colnames(data)[firstnumeric:ncol(data)],design)

  gene.info = data.frame(Gene=as.character(data[[1]]),Symbol=as.character(data[[1]]),stringsAsFactors=FALSE)
  for (i in 1:(firstnumeric-1)) gene.info[[anno.names[i]]]=data[[i]]
  gene.info$Type=classify.genes(gene.info)

  re=list()
  re$count=tomat(data[,firstnumeric:ncol(data)],data$Gene,names(data)[firstnumeric:ncol(data)])
  re$ntr=re$count
  re$ntr[,]=NA

  coldata$no4sU=TRUE
  do.callback()

  # insert name no rep and set this as default condition, if there is no condition field!
  re=grandR(prefix=prefix,gene.info=gene.info,slots=re,coldata=coldata,metadata=list(Description="Count data"))
  DefaultSlot(re)="count"
  re
}


GetTableQC=function(data,name) {
  fn=paste0(data$prefix,".",name,".tsv.gz")
  if (!file.exists(fn)) {
    fn2=paste0(data$prefix,".",name,".tsv")
    if (!file.exists(fn2)) stop(paste0("Cannot find QC table ",fn," or ",fn2))
    fn = fn2
  }

  header = !name %in% c("clip","strandness")
  read.tsv(fn,header=header)
}



#' Read the output of GRAND-SLAM 3.0 into a grandR object.
#'
#' Metabolic labeling - nucleotide conversion RNA-seq data (such as generated by SLAM-seq,TimeLapse-seq or TUC-seq)
#' must be carefully analyzed to remove bias due to incomplete labeling. GRAND-SLAM is a software package that
#' employs a binomial mixture modeling approach to obtain precise estimates of the new-to-total RNA ratio (NTR) per gene and sample (or cell).
#' This function directly reads the output of GRAND-SLAM 3.0 into a grandR object.
#'
#' @param prefix the prefix used to call GRAND-SLAM
#' @param design Either a design vector (see details), or a data.frame providing metadata for all columns (samples/cells),
#' or a function that is called with the condition name vector and is supposed to return this data.frame. if NULL, a
#' library,sample,barcode design is used for sparse data, and a condition,replicate design for dense data
#' @param label which nucleoside analog
#' @param estimator which estimator to use (one of Binom,TbBinom,TbBinomShape)
#' @param classify.genes A function that is used to add the \emph{type} column to the gene annotation table, always a call to \link{ClassifyGenes}
#' @param read.posterior also read the posterior parameters alpha and beta? if NULL, TRUE for dense data, FALSE for sparse data
#' @param rename.sample function that is applied to each sample name before parsing (or NULL)
#' @param verbose Print status updates
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
#' @details Sometimes you might have forgotten to name all samples consistently (or you simply messed something up).
#' In this case, the rename.sample parameter can be handy (e.g. to rename a particular misnamed sample).
#'
#' @seealso \link{ReadGRAND},\link{ClassifyGenes},\link{MakeColdata},\link{DesignSemantics}
#'
#'
#' @export
#'
#' @concept load
ReadGRAND3=function(prefix,
                    design=NULL,
                    label="4sU",
                    estimator="Binom",
                    classify.genes=ClassifyGenes(),
                    read.posterior=NULL,
                    rename.sample=NULL,
                    verbose=FALSE) {

  if (length(estimator)!=1 || !estimator %in% c("Binom","TbBinom","TbBinomShape")) stop("Invalid estimator!")

  # TODO: do not use read.delim or .csv, nor file.exists (to allow using urls!)
  isSparse=function(prefix) !file.exists(paste0(prefix,".targets/data.tsv.gz")) & !file.exists(prefix)

  if (isSparse(prefix)) {
    if (is.null(design)) design=c(Design$Library,Design$Sample,Design$Barcode)
    if (is.null(read.posterior)) read.posterior=FALSE
    ReadGRAND3_sparse(prefix=prefix,design=design,label=label,estimator=estimator,classify.genes=classify.genes,read.posterior=read.posterior,rename.sample = rename.sample,verbose=verbose)
  } else {
    if (is.null(design)) design=c(Design$Condition,Design$Replicate)
    if (is.null(read.posterior)) read.posterior=TRUE
    ReadGRAND3_dense(prefix=prefix,design=design,label=label,estimator=estimator,classify.genes=classify.genes,read.posterior=read.posterior,rename.sample = rename.sample,verbose=verbose)
  }
}
ReadGRAND3_sparse=function(prefix,
                           design=c(Design$Library,Design$Sample,Design$Barcode),
                           label="4sU",
                           estimator="Binom",
                           classify.genes=ClassifyGenes(),
                           read.posterior=FALSE,
                           rename.sample=NULL,
                           verbose=FALSE) {

  cols=readLines(paste0(prefix, ".targets/barcodes.tsv.gz"))
  conds=strsplit(cols,".",fixed=TRUE)[[1]]
  if (!is.null(rename.sample)) conds=sapply(conds,rename.sample) # let's use sapply, we have no idea whether the function works with vectorization

  if (length(conds)!=length(design)) stop(paste0("Design parameter is incompatible with input data: ",paste(conds,collapse=".")))


  if (verbose) cat("Reading count matrix...\n")
  count=methods::as(Matrix::readMM(paste0(prefix, ".targets/matrix.mtx.gz")),Class = "dgCMatrix")
  gene.info=read.delim(paste0(prefix, ".targets/features.tsv.gz"),header = FALSE,stringsAsFactors = FALSE)
  if (ncol(gene.info)==4) gene.info$Length=1
  gene.info=setNames(gene.info,c("Gene","Symbol","Mode","Category","Length"))

  if (anyDuplicated(gene.info$Gene)) {
    dupp=table(gene.info$Gene)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Gene=make.unique(gene.info$Gene)
  }
  if (any(is.na(gene.info$Symbol) | gene.info$Symbol=="")) {
    warning(sprintf("No gene symbols (e.g. %s) present, replacing by gene ids!",paste(head(gene.info$Gene[is.na(gene.info$Symbol) | gene.info$Symbol==""]),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Symbol[is.na(gene.info$Symbol) | gene.info$Symbol==""]=gene.info$Gene[is.na(gene.info$Symbol) | gene.info$Symbol==""]
  }
  if (anyDuplicated(gene.info$Symbol)) {
    dupp=table(gene.info$Symbol)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene symbols (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    gene.info$Symbol=make.unique(gene.info$Symbol)
  }
  if (any(gene.info$Gene=="")) stop("Gene name may not be empty!")

#  if (any(make.names(gene.info$Symbol)!=gene.info$Symbol)) {
#    ill=gene.info$Symbol[which(make.names(gene.info$Symbol)!=gene.info$Symbol)]
#    warning(sprintf("Illegal identifiers among the gene symbols (e.g. %s), making legal!",paste(head(ill),collapse=",")),call. = FALSE,immediate. = TRUE)
#    gene.info$Symbol=make.names(gene.info$Symbol)
#  }

  colnames(count)=cols
  rownames(count)=gene.info$Gene
  re=list()
  re$count=count

  if (verbose) cat("Reading NTRs...\n")
  ntr=methods::as(Matrix::readMM(sprintf("%s.targets/%s.%s.ntr.mtx.gz",prefix,label,estimator)),Class = "dgCMatrix")
  colnames(ntr)=cols
  rownames(ntr)=gene.info$Gene
  re$ntr=ntr

  if (read.posterior && file.exists(sprintf("%s.targets/%s.%s.alpha.mtx.gz",prefix,label,estimator)) && file.exists(sprintf("%s.targets/%s.%s.beta.mtx.gz",prefix,label,estimator))) {
    if (verbose) cat("Reading posterior beta parameters...\n")
    alpha=methods::as(Matrix::readMM(sprintf("%s.targets/%s.%s.alpha.mtx.gz",prefix,label,estimator)),Class = "dgCMatrix")
    colnames(alpha)=cols
    rownames(alpha)=gene.info$Gene
    re$alpha=alpha

    beta=methods::as(Matrix::readMM(sprintf("%s.targets/%s.%s.beta.mtx.gz",prefix,label,estimator)),Class = "dgCMatrix")
    colnames(beta)=cols
    rownames(beta)=gene.info$Gene
    re$beta=beta
  }
  gene.info$Mode=gsub(".*\\(","",gsub(")","",gene.info$Category,fixed=TRUE))
  gene.info$Mode=factor(gene.info$Mode,levels=unique(gene.info$Mode))

  gene.info$Type=classify.genes(gene.info)


  # use make coldata instead. add no4sU column!
  coldata=data.frame(Name=cols)
  spl=strsplit(as.character(coldata$Name),".",fixed=TRUE)
  for (i in 1:length(design)) coldata=cbind(coldata,factor(sapply(spl,function(v) v[i]),levels=unique(sapply(spl,function(v) v[i]))))
  names(coldata)[-1]=design
  rownames(coldata)=coldata$Name

  coldata$no4sU=Matrix::colSums(ntr)==0

  re=grandR(prefix=prefix,gene.info=gene.info,slots=re,coldata=coldata,metadata=list(Description="GRAND-SLAM 3.0 sparse data"))
  DefaultSlot(re)="count"
  re
}


ReadGRAND3_dense=function(prefix,
                          design=c(Design$Condition,Design$Replicate),
                          label="4sU",
                          estimator="Binom",
                          classify.genes=ClassifyGenes(),
                          read.posterior=TRUE,
                          rename.sample=NULL,
                          verbose=FALSE) {
  annotations=c("Gene","Symbol","Category","Length")
  slots=c("count","ntr","alpha","beta")
  names(slots)=c("Read count",sprintf("%s %s %s",label,estimator,c("NTR MAP","alpha","beta")))
  if (!read.posterior) slots=slots[1:2]
  if (estimator=="TbBinomShape") slots=c(slots,Shape="shape",LLR="llr")

  re=read.grand.internal(description="GRAND-SLAM 3.0 dense data",prefix = prefix, design = design, slots=slots, annotations=annotations,classify.genes = classify.genes,rename.sample = rename.sample,verbose = verbose)
  re
}

ReadNewTotal=function(genes, cells, new.matrix, total.matrix, detection.rate=1,verbose=FALSE) {

  cols=read.csv(cells,check.names = FALSE,stringsAsFactors = FALSE)
  gene.info=setNames(read.csv(genes,check.names = FALSE,stringsAsFactors = FALSE),c("Gene","Biotype","Symbol"))

  if (verbose) cat("Reading total count matrix...\n")
  count=methods::as(Matrix::readMM(total.matrix),Class = "dgCMatrix")

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
  new=methods::as(Matrix::readMM(new.matrix),Class = "dgCMatrix")

  if (verbose) cat("Computing NTRs...\n")
  new=new/detection.rate
  ntr=new/count
  ntr@x[ntr@x>1]=1
  ntr@x[is.nan(ntr@x)]=0
  ntr=methods::as(ntr,Class = "dgCMatrix")

  colnames(ntr)=colnames(count)
  rownames(ntr)=rownames(count)

  re=list()
  re$count=count
  re$ntr=ntr

  grandR("",gene.info,slots=re,cols)
}


read.grand.internal=function(prefix, design=c(Design$Condition,Design$Replicate),
                             slots,
                             annotations,
                             classify.genes=ClassifyGenes(),
                             rename.sample=NULL, verbose=FALSE, description="") {

  if (!all(c("count","ntr") %in% slots) || !all(c("Gene","Symbol") %in% annotations)) stop("Invalid call to read.grand.internal!")

  tomat=function(m,names,cnames){
    m=as.matrix(m)
    m[is.na(m)]=0
    colnames(m)=gsub(" .*","",cnames)
    if (!is.null(rename.sample)) colnames(m)=sapply(colnames(m),rename.sample)
    rownames(m)=names
    m
  }

  do.callback=function() {}

  url=NULL
  if (suppressWarnings(requireNamespace("RCurl",quietly = TRUE))) {
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
    if (file.exists(paste0(prefix,".targets/data.tsv.gz"))) {
      file=paste0(prefix,".targets/data.tsv.gz")
    } else {
      file=if (file.exists(prefix)) prefix else paste0(prefix,".tsv")
      if (!file.exists(file) && file.exists(paste0(file,".gz"))) file = paste0(file,".gz")
      prefix=gsub(".tsv(.gz)?$","",file)
    }
  }

  if (!file.exists(file)) stop("File not found; If you want to access non-local files directly, please install the RCurl package!")


  if (verbose) cat("Checking file...\n")
  con <- file(file,"r")
  header <- strsplit(readLines(con,n=1),"\t")[[1]]
  close(con)

  count.name=names(slots)[slots=="count"]
  if (length(count.name)!=1) stop("Invalid call to read.grand.internal!")

  if (header[1]!="Gene" || header[2]!="Symbol" || !grepl(paste0(count.name,"$"),header[!header %in% annotations][1])) stop("File is not a GRAND-SLAM output file!")
  conds=gsub(paste0(" ",count.name),"",header[grepl(paste0(" ",count.name),header)])

  if (!is.null(rename.sample)) conds=sapply(conds,rename.sample) # let's use sapply, we have no idea whether the function works with vectorization

  terms=strsplit(conds[1],".",fixed=TRUE)[[1]]

  if (is.data.frame(design)) {
    if (length(conds)!=nrow(design)) stop(paste0("Design parameter (table) is incompatible with input data: ",paste(conds,collapse=", ")))
    if (is.null(design$Name) || !all(design$Name==conds)) stop(paste0("Design parameter (table) must contain a Name column corresponding to the sample names!"))
    coldata=design
    rownames(coldata)=coldata$Name
  } else if (is.function(design)) {
    coldata=design(conds)
    if (length(conds)!=nrow(coldata)) stop(paste0("Design parameter (function) is incompatible with input data: ",paste(conds,collapse=", ")))
  } else {
    if (length(terms)!=length(design)) stop(paste0("Design parameter is incompatible with input data: ",paste(conds,collapse=", ")))
    coldata=MakeColdata(conds,design)
  }


  if (verbose) cat("Reading files...\n")
  data=read.delim(file,stringsAsFactors=FALSE,check.names=FALSE)
  if (!is.null(rename.sample)) colnames(data)=sapply(colnames(data),rename.sample)

  if (anyDuplicated(data$Gene)) {
    dupp=table(data$Gene)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    data$Gene=make.unique(data$Gene)
  }
  if (any(data$Gene=="")) stop("Gene name may not be empty!")
  if (any(is.na(data$Symbol) | data$Symbol=="")) {
    warning(sprintf("No gene symbols (e.g. %s) present, replacing by gene ids!",paste(head(data$Gene[is.na(data$Symbol) | data$Symbol==""]),collapse=",")),call. = FALSE,immediate. = TRUE)
    data$Symbol[is.na(data$Symbol) | data$Symbol==""]=data$Gene[is.na(data$Symbol) | data$Symbol==""]
  }
  if (anyDuplicated(data$Symbol)) {
    dupp=table(data$Symbol)
    dupp=names(dupp)[which(dupp>1)]
    warning(sprintf("Duplicate gene symbols (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
    data$Symbol=make.unique(data$Symbol)
  }
  if (verbose) cat("Processing...\n")

  #ntr.mean = grepl("Mean",names(data))
  #ntr = grepl("MAP",names(data))
  #alpha = grepl("alpha",names(data))
  #beta = grepl("beta",names(data))
  #count = grepl("Readcount",names(data))
  #conv = grepl("Conversions",names(data))
  #cove = grepl("Coverage",names(data)) & !grepl("Double-Hit Coverage",names(data))


  gene.info = data.frame(Gene=as.character(data$Gene),Symbol=as.character(data$Symbol),stringsAsFactors=FALSE)
  for (a in setdiff(annotations,c("Gene","Symbol"))) gene.info[[a]]=data[[a]]
  gene.info$Type=classify.genes(gene.info)


  re=list()
  for (n in names(slots)) {
    cols=intersect(paste0(conds," ",n),names(data))
    re[[slots[n]]]=tomat(data[,cols],data$Gene,cols)
  }
  #re$count=tomat(data[,count],data$Gene,names(data)[count])
  ##	re$tpm=comp.tpm(re$count,gene.info$Length)
  ##	checknames(re$count,re$tpm)
  ##re$ntr.mean=tomat(data[,alpha]/(data[,alpha]+data[,beta]),data$Gene,names(data)[ntr.mean])
  #re$ntr=tomat(data[,ntr],data$Gene,names(data)[ntr])
  #re$alpha=tomat(data[,alpha],data$Gene,names(data)[alpha])
  #re$beta=tomat(data[,beta],data$Gene,names(data)[beta])


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

  checknames=function(a){
    if (!all(colnames(a)==coldata$Name)) stop("Column names do not match!")
    if (!all(rownames(a)==gene.info$Gene)) stop("Row names do not match!")
  }

  for (slot in slots) checknames(re[[slot]])

  coldata$no4sU=no4sU.cols

  do.callback()

  # insert name no rep and set this as default condition, if there is no condition field!
  re=grandR(prefix=prefix,gene.info=gene.info,slots=re,coldata=coldata,metadata=list(Description=description))
  DefaultSlot(re)="count"
  re
}
