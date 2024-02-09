




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


check.and.make.unique = function(v,ref=NULL,label="entries",ref.label="reference",do.error=FALSE) {
  fun=if(do.error) stop else warning
  if (!is.null(ref) && (any(is.na(ref) | ref==""))) stop("References must be unique")
  if (!is.null(ref) && (length(ref)!=length(v))) stop("References have invalid length")

  if (any(is.na(v) | v=="")) {
    if (is.null(ref)) stop(sprintf("When %s are empty, %s must be provided!",label,ref.label))
    howmany=if (sum(is.na(v)|v=="")==length(v)) "All" else "Some"
    fun(sprintf("%s %s are empty (n=%d, e.g. %s), replacing by %s ids!",howmany,label,sum(is.na(v) | v==""),paste(sample(utils::head(ref[is.na(v) | v==""])),collapse=","),ref.label),call. = FALSE,immediate. = TRUE)
    v[is.na(v) | v==""]=ref[is.na(v) | v==""]
  }

  if (anyDuplicated(v)) {
    ext = ""
    if (is.null(ref) || anyDuplicated(ref)) ext = sprintf(" Cannot guarantee maintaining consistency for %s across reading several files, watch out if you merge grandR objects!",label)
    dupp=table(v)
    dupp=names(dupp)[which(dupp>1)]
    fun(sprintf("Duplicate %s (n=%d, e.g. %s) present, making unique!%s",label,length(dupp),paste(utils::head(sample(dupp)),collapse=","),ext),call. = FALSE,immediate. = TRUE)

    if (!is.null(ref)) {
      df=data.frame(id=1:length(v),v=as.character(v),ref=as.character(ref),stringsAsFactors = FALSE)
      df=df[order(df$ref),]
      df$v=make.unique(df$v)
      df=df[order(df$id),]
      v=df$v
    } else {
      v=make.unique(v)
    }
  }
  v
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
#' @details Semantics.concentration is such a predefined function: Contents such as 200uM or 1mM are converted into a numerical value (in uM), and no4sU is converted into 0.
#'
#' @details By default, Semantics.time is used for the names duration.4sU and Experimental.time, and Semantics.concentration is used for concentration.4sU
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
  if (!"Time" %in% names(ll)) ll=c(ll,list(Time=Semantics.time))
  if (!"concentration.4sU" %in% names(ll)) ll=c(ll,list(concentration.4sU=Semantics.concentration))
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
  min=grepl("[0-9.]+min",s)

  if (sum(h)>0) {
    time[h]=as.numeric(gsub("([0-9.]+)h(p.)?","\\1",s[h]))
    time[min]=as.numeric(substr(s[min],1,nchar(s[min])-3))/60
  } else {
    time[min]=as.numeric(substr(s[min],1,nchar(s[min])-3))
  }

  nounit=grepl("^[0-9.]+$",s)
  time[nounit]=as.numeric(s[nounit])

  if (any(is.na(time))) stop(paste0("Time semantics cannot be used for this: ",paste(s[is.na(time)],collapse=",")))
  setNames(data.frame(time),name)
}



#' Semantics for concentration columns
#'
#' Defines additional semantics for columns representing concentrations
#'
#' @param s original column
#' @param name the column name
#'
#' @return a data frame with a single numeric column, where <x>uM from s is replaced by x, <x>mM is replaced by
#' x*1000, and no4sU is replaced by 0
#'
#' @export
#'
#' @concept load
Semantics.concentration=function(s,name) {
  conc=rep(NA,length(s))

  no4sU=c("nos4U","no4sU","-")
  conc[s %in% no4sU]=0
  s=gsub("_",".",s)

  h=grepl("[0-9.]+uM",s)
  conc[h]=as.numeric(gsub("([0-9.]+)uM","\\1",s[h]))

  min=grepl("[0-9.]+mM",s)
  conc[min]=as.numeric(gsub("([0-9.]+)mM","\\1",s[min]))*1000

  nounit=grepl("^[0-9.]+$",s)
  conc[nounit]=as.numeric(s[nounit])


  if (any(is.na(conc))) stop(paste0("Concentration semantics cannot be used for this: ",paste(s[is.na(conc)],collapse=",")))
  setNames(data.frame(as.numeric(conc)),name)
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
  for (i in 1:length(design)) if (!is.na(design[i])) coldata=cbind(coldata,str2fac(utils::type.convert(sapply(spl,function(v) v[i]),as.is=TRUE,na.strings=c("DKSALDJLKSADLKSJDLKSJDLJDA")))) # stupid function!
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
#' @param read.min2 Should the read count with at least 2 mismatches also be read?
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
                   read.min2=FALSE,
                   rename.sample=NULL,
                   verbose=FALSE) {
  annotations=c("Gene","Symbol","Length")
  slots=c(`Readcount`="count",MAP="ntr",alpha="alpha",beta="beta")
  if (read.percent.conv) slots=c(slots,`Conversions`="conv",Coverage="cove")
  if (read.min2) slots=c(slots,`min2`="min2")
  re=read.grand.internal(prefix = prefix, design = design, slots=slots, annotations=annotations,classify.genes = classify.genes,verbose = verbose,rename.sample = rename.sample)
  re$metadata=c(re$metadata,list(`GRAND-SLAM version`=2,Output="dense"))
  if (read.percent.conv) {
    re=AddSlot(re,"percent_conv",re$data$conv/re$data$cove)
    re=DropSlot(re,"^conv|cove")
  }
  re
}



#' Read a count table
#'
#' grandR can also be used to analyze standard RNA-seq data, and this function is here to read such data.
#'
#' @param file a file containing a count matrix
#' @param design Either a design vector (see details), or a data.frame providing metadata for all columns (samples/cells),
#' or a function that is called with the condition name vector and is supposed to return this data.frame.
#' @param classify.genes A function that is used to add the \emph{type} column to the gene annotation table, always a call to \link{ClassifyGenes}
#' @param rename.sample function that is applied to each sample name before parsing (or NULL)
#' @param filter.table function that is applied to the table directly after read it (or NULL)
#' @param num.samples number of sample columns containing read counts (can be NULL, see details)
#' @param verbose Print status updates
#' @param sep The column separator used in the file
#'
#' @return a grandR object
#'
#' @details The table is assumed to have read counts in the last n columns, which must be named according to sample names.
#' If num.samples is NULL this n is automatically recognized as the number of numeric columns (so make sure to either
#' specify num.samples, or that the column immediately prior to the first sample column is *not* numeric).
#'
#' @details If these columns are named systematically in a particular way, the design vector provides
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
#' @details Sometimes the table contains more than you want to read. In this case, use the filter.table parameter to preprocess it.
#' This should be a function that receives a data.frame, and returns a data.frame.
#'
#' @details If there are no columns named "Gene" or "Symbol", the first column is used!
#'
#' @export
#'
#' @concept load
ReadCounts=function(file, design=c(Design$Condition,Design$Replicate),classify.genes=ClassifyGenes(),rename.sample=NULL,filter.table=NULL, num.samples=NULL,
                    verbose=FALSE,sep="\t") {

  tomat=function(m,names,cnames){
    m=as.matrix(m)
    m[is.na(m)]=0
    colnames(m)=gsub(" .*","",cnames)
    rownames(m)=names
    m
  }

  checknames=function(a,b){
    if (nrow(a)!=nrow(b)) stop("Number of rows do not match!")
    if (ncol(a)!=ncol(b)) stop("Number of columns do not match!")
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
      utils::download.file(url, file, quiet=!verbose)
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

  conds=header[length(header)]
  if (!is.null(rename.sample)) conds=rename.sample(conds)
  terms=strsplit(conds,".",fixed=TRUE)[[1]]

  if (length(terms)!=length(design)) stop(paste0("Design parameter is incompatible with input data: ",paste(terms,collapse=".")))


  if (verbose) cat("Reading file...\n")
  data=utils::read.table(file,sep=sep,stringsAsFactors=FALSE,check.names=FALSE,header=TRUE)
  if (!is.null(filter.table)) data=filter.table(data)

  clss=sapply(data,class)
  if (is.null(num.samples)) firstnumeric=max(which(clss!="numeric"))+1 #min(which(clss=="numeric"))
  else firstnumeric = ((ncol(data)-num.samples+1))
  if (!all(clss[firstnumeric:length(clss)]=="numeric") || firstnumeric==1 || firstnumeric>ncol(data)) stop("Columns (except for the first n) must be numeric!")
  anno.names=colnames(data)[1:(firstnumeric-1)]


  #if (anyDuplicated(data[[1]])) {
  #  dupp=table(data[[1]])
  #  dupp=names(dupp)[which(dupp>1)]
  #  warning(sprintf("Duplicate names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
  #  data[[1]]=make.unique(data[[1]])
  #}
  if (verbose) cat("Processing...\n")

  if (!is.null(rename.sample)) colnames(data)[firstnumeric:ncol(data)]=sapply(colnames(data)[firstnumeric:ncol(data)],rename.sample) # let's use sapply, we have no idea whether the function works with vectorization
  coldata=MakeColdata(colnames(data)[firstnumeric:ncol(data)],design)

  if ("Symbol" %in% colnames(data)) {
    index <- which(colnames(data) == "Symbol")
  } else if ("Gene" %in% colnames(data)) {
    index <- which(colnames(data) == "Gene")
  } else {
    # If neither "Gene" nor "Symbol" are present, use first column as index
    index <- 1
  }
  data[[index]] = check.and.make.unique(data[[index]],label="names")
  gene.info = data.frame(Gene=as.character(data[[index]]),Symbol=as.character(data[[index]]),stringsAsFactors=FALSE)
  for (i in 1:(firstnumeric-1)) gene.info[[anno.names[i]]]=data[[i]]
  gene.info$Type=classify.genes(gene.info)

  re=list()
  re$count=tomat(data[,firstnumeric:ncol(data)],gene.info$Gene,names(data)[firstnumeric:ncol(data)])
  re$ntr=re$count
  re$ntr[,]=NA

  coldata$no4sU=TRUE
  do.callback()

    # insert name no rep and set this as default condition, if there is no condition field!
  re=grandR(prefix=prefix,gene.info=gene.info,slots=re,coldata=coldata,metadata=list(Description="Count data"))
  DefaultSlot(re)="count"
  re
}

#' Read featureCounts
#'
#' grandR can also be used to analyze standard RNA-seq data, and this function is here to read such data.
#'
#' @param file a file containing featureCounts
#' @param design Either a design vector (see details), or a data.frame providing metadata for all columns (samples/cells),
#' or a function that is called with the condition name vector and is supposed to return this data.frame.
#' @param classify.genes A function that is used to add the \emph{type} column to the gene annotation table, always a call to \link{ClassifyGenes}
#' @param rename.sample function that is applied to each sample name before parsing (or NULL)
#' @param filter.table function that is applied to the table directly after read it (or NULL)
#' @param num.samples number of sample columns containing read counts (can be NULL, see details)
#' @param verbose Print status updates
#' @param sep The column separator used in the file
#'
#' @return a grandR object
#'
#' @details The table is assumed to have read counts in the last n columns, which must be named according to sample names.
#' If num.samples is NULL this n is automatically recognized as the number of columns containing .bam (so make sure to either
#' specify num.samples, or that the count columns are called after the bam files).
#'
#' @details If these columns are named systematically in a particular way, the design vector provides
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
#' @details Sometimes the table contains more than you want to read. In this case, use the filter.table parameter to preprocess it.
#' This should be a function that receives a data.frame, and returns a data.frame.
#'
#' @details If there are no columns named "Geneid", "Gene" or "Symbol", the first column is used!
#'
#' @export
#'
#' @concept load
ReadFeatureCounts=function(file, design=c(Design$Condition,Design$Replicate),classify.genes=ClassifyGenes(),rename.sample=NULL,filter.table=NULL, num.samples=NULL,
                           verbose=FALSE,sep="\t") {

  tomat=function(m,names,cnames){
    m=as.matrix(m)
    m[is.na(m)]=0
    colnames(m)=gsub(" .*","",cnames)
    rownames(m)=names
    m
  }

  checknames=function(a,b){
    if (nrow(a)!=nrow(b)) stop("Number of rows do not match!")
    if (ncol(a)!=ncol(b)) stop("Number of columns do not match!")
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
      utils::download.file(url, file, quiet=!verbose)
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
  header <- strsplit(readLines(con,n=2),sep)[[2]]
  close(con)


  conds=header[length(header)]
  if (!is.null(rename.sample)) conds=rename.sample(conds)
  terms=strsplit(conds,".",fixed=TRUE)[[1]]

  if (length(terms)!=length(design)) stop(paste0("Design parameter is incompatible with input data: ",paste(terms,collapse=".")))


  if (verbose) cat("Reading file...\n")
  data=utils::read.table(file,sep=sep,stringsAsFactors=FALSE,check.names=FALSE,header=TRUE)
  if (!is.null(filter.table)) data=filter.table(data)



  clss=sapply(data,class)
  if (is.null(num.samples)) firstsample = which(grepl(".bam", colnames(data)), arr.ind=TRUE)[1]
  else firstsample = (ncol(data)-num.samples+1)
  anno.names=colnames(data)[firstsample:ncol(data)]


  #if (anyDuplicated(data[[1]])) {
  #  dupp=table(data[[1]])
  #  dupp=names(dupp)[which(dupp>1)]
  #  warning(sprintf("Duplicate names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
  #  data[[1]]=make.unique(data[[1]])
  #}

  if (verbose) cat("Processing...\n")

  if (!is.null(rename.sample)) colnames(data)[firstsample:ncol(data)]=sapply(colnames(data)[firstsample:ncol(data)],rename.sample) # let's use sapply, we have no idea whether the function works with vectorization
  coldata=MakeColdata(colnames(data)[firstsample:ncol(data)],design)

  if ("Symbol" %in% colnames(data)) {
    index <- which(colnames(data) == "Symbol")
  }
  if ("Geneid" %in% colnames(data)) {
    index <- which(colnames(data) == "Geneid")
  }else if ("Gene" %in% colnames(data)) {
    index <- which(colnames(data) == "Gene")
  } else {
    # If neither "Gene" nor "Symbol" are present, use first column as index
    index <- 1
  }
  data[[index]] = check.and.make.unique(data[[index]],label="names")
  gene.info = data.frame(Gene=as.character(data[[index]]),Symbol=as.character(data[[index]]),stringsAsFactors=FALSE)
  for (i in 1:(firstsample-1)) gene.info[[anno.names[i]]]=data[[i]]
  gene.info$Type=classify.genes(gene.info)

  re=list()
  re$count=tomat(data[,firstsample:ncol(data)],data$Gene,names(data)[firstsample:ncol(data)])
  re$ntr=re$count
  re$ntr[,]=NA

  coldata$no4sU=TRUE
  do.callback()

  # insert name no rep and set this as default condition, if there is no condition field!
  re=grandR(prefix=prefix,gene.info=gene.info,slots=re,coldata=coldata,metadata=list(Description="Count data"))
  DefaultSlot(re)="count"
  re
}

GetTableQC=function(data,name,stop.if.not.exist=TRUE) {
  ll=try.file(paste0(data$prefix,".",name),possible.suffixes = c(".tsv.gz",".tsv",""),stop.if.not.exist=stop.if.not.exist)
  if (is.null(ll)) {
    warning(paste0("Cannot find QC table ",name))
    return(NULL)
  }
  header = !name %in% c("clip","strandness")
  re=read.tsv(ll$file,header=header)
  ll$callback()
  return(re)

#  fn=paste0(data$prefix,".",name,".tsv.gz")
#  if (!file.exists(fn)) {
#    fn2=paste0(data$prefix,".",name,".tsv")
#    if (!file.exists(fn2)) {
#      fn3=paste0(data$prefix,".",name)
#      if (!file.exists(fn3)) {
#        if (stop.if.not.exist) stop(paste0("Cannot find QC table ",fn," or ",fn2))
#        else warning(paste0("Cannot find QC table ",fn," or ",fn2))
#      }
#      fn=fn3
#    } else {
#      fn = fn2
#    }
#  }
#  if (!file.exists(fn)) return(NULL)
#
#  header = !name %in% c("clip","strandness")
#  read.tsv(fn,header=header)
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
  count=Matrix::readMM(paste0(prefix, ".targets/matrix.mtx.gz"))
  gene.info=utils::read.delim(paste0(prefix, ".targets/features.tsv.gz"),header = FALSE,stringsAsFactors = FALSE)
  if (ncol(gene.info)==4) gene.info$Length=1
  gene.info=setNames(gene.info,c("Gene","Symbol","Mode","Category","Length"))

  gene.info$Gene = check.and.make.unique(gene.info$Gene,label="gene names")
  gene.info$Symbol = check.and.make.unique(gene.info$Symbol,ref=gene.info$Gene,label="gene symbols",ref.label = "gene names")

  #if (anyDuplicated(gene.info$Gene)) {
  #  dupp=table(gene.info$Gene)
  #  dupp=names(dupp)[which(dupp>1)]
  #  warning(sprintf("Duplicate gene names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
  #  gene.info$Gene=make.unique(gene.info$Gene)
  #}
  #if (any(is.na(gene.info$Symbol) | gene.info$Symbol=="")) {
  #  warning(sprintf("No gene symbols (e.g. %s) present, replacing by gene ids!",paste(head(gene.info$Gene[is.na(gene.info$Symbol) | gene.info$Symbol==""]),collapse=",")),call. = FALSE,immediate. = TRUE)
  #  gene.info$Symbol[is.na(gene.info$Symbol) | gene.info$Symbol==""]=gene.info$Gene[is.na(gene.info$Symbol) | gene.info$Symbol==""]
  #}
  #if (anyDuplicated(gene.info$Symbol)) {
  #  dupp=table(gene.info$Symbol)
  #  dupp=names(dupp)[which(dupp>1)]
  #  warning(sprintf("Duplicate gene symbols (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
  #  gene.info$Symbol=make.unique(gene.info$Symbol)
  #}
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
  ntr=Matrix::readMM(sprintf("%s.targets/%s.%s.ntr.mtx.gz",prefix,label,estimator))
  colnames(ntr)=cols
  rownames(ntr)=gene.info$Gene
  re$ntr=ntr

  if (read.posterior && file.exists(sprintf("%s.targets/%s.%s.alpha.mtx.gz",prefix,label,estimator)) && file.exists(sprintf("%s.targets/%s.%s.beta.mtx.gz",prefix,label,estimator))) {
    if (verbose) cat("Reading posterior beta parameters...\n")
    alpha=Matrix::readMM(sprintf("%s.targets/%s.%s.alpha.mtx.gz",prefix,label,estimator))
    colnames(alpha)=cols
    rownames(alpha)=gene.info$Gene
    re$alpha=alpha

    beta=Matrix::readMM(sprintf("%s.targets/%s.%s.beta.mtx.gz",prefix,label,estimator))
    colnames(beta)=cols
    rownames(beta)=gene.info$Gene
    re$beta=beta
  }
  if (is.null(gene.info$Mode)) {
    gene.info$Mode=gsub(".*\\(","",gsub(")","",gene.info$Category,fixed=TRUE))
    gene.info$Mode=factor(gene.info$Mode,levels=unique(gene.info$Mode))
  }

  gene.info$Type=classify.genes(gene.info)

  coldata=MakeColdata(cols,design)
  # use make coldata instead. add no4sU column!
  #coldata=data.frame(Name=cols)
  #spl=strsplit(as.character(coldata$Name),".",fixed=TRUE)
  #for (i in 1:length(design)) coldata=cbind(coldata,factor(sapply(spl,function(v) v[i]),levels=unique(sapply(spl,function(v) v[i]))))
  #names(coldata)[-1]=design
  #rownames(coldata)=coldata$Name

  coldata$no4sU=Matrix::colSums(ntr)==0

  re=grandR(prefix=prefix,gene.info=gene.info,slots=re,coldata=coldata,metadata=list(`GRAND-SLAM version`=3,Output="sparse"))
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
  re$metadata=c(re$metadata,list(`GRAND-SLAM version`=3,Output="dense"))
  re
}


#' Read sparse new/total matrices
#'
#' This function can be used to load matrix market data in case genes were quantified by
#' (i) counting all reads (for total RNA) and (ii) counting T-to-C mismatch reads (for new RNA)
#'
#' @param genes csv file (or URL) containing gene information
#' @param cells csv file (or URL) containing cell information
#' @param new.matrix Matrix market file of new counts
#' @param total.matrix Matrix market file of total counts
#' @param detection.rate the detection rate of T-to-C mismatch reads (see details)
#' @param verbose verbose output
#'
#' @details  Metabolic labeling - nucleotide conversion RNA-seq data (such as generated by SLAM-seq,TimeLapse-seq or TUC-seq)
#' must be carefully analyzed to remove bias due to incomplete labeling. We advice against counting read with and without T-to-C mismatches
#' for quantification, and encourage using a statistical method such as GRAND-SLAM that properly deals with incomplete labeling.
#'
#' @details To correct for some bias, a detection rate (as suggested by Cao et al., Nature Biotech 2020) should be provided. This detection rate
#' defines, how much new RNA is detected on average using the T-to-C mismatch reads.
#'
#' @return a grandR object
#' @export
#'
#' @useDynLib grandR, .registration = TRUE
#' @concept load
ReadNewTotal=function(genes, cells, new.matrix, total.matrix, detection.rate=1,verbose=FALSE) {

  genes=try.file(genes,possible.suffixes = "",verbose=verbose)
  cells=try.file(cells,possible.suffixes = "",verbose=verbose)
  new.matrix=try.file(new.matrix,possible.suffixes = "",verbose=verbose)
  total.matrix=try.file(total.matrix,possible.suffixes = "",verbose=verbose)

  callbacks=list(genes$callback,cells$callback,new.matrix$callback,total.matrix$callback)

  cols=utils::read.csv(cells$file,check.names = FALSE,stringsAsFactors = FALSE)
  rownames(cols)=cols[,1]
  gene.info=setNames(utils::read.csv(genes$file,check.names = FALSE,stringsAsFactors = FALSE),c("Gene","Biotype","Symbol"))

  if (verbose) cat("Reading total count matrix...\n")
  count=Matrix::readMM(total.matrix$file)

  gene.info$Gene = check.and.make.unique(gene.info$Gene,label="gene names")
  gene.info$Symbol = check.and.make.unique(gene.info$Symbol,ref=gene.info$Gene,label="gene symbols",ref.label = "gene names")

  colnames(count)=cols[,1]
  rownames(count)=gene.info$Gene

  if (verbose) cat("Reading new count matrix...\n")
  new=Matrix::readMM(new.matrix$file)

  if (verbose) cat("Computing NTRs...\n")
  #new=new/detection.rate

    sX=Matrix::summary(new)
    sY=Matrix::summary(count)
  dd=.Call('fastsparsematdiv',sX$i,sX$j,sX$x,sY$i,sY$j,sY$x,detection.rate)
  ntr=Matrix::sparseMatrix(i=sX$i, j=sX$j, x=dd,dimnames=dimnames(new))

  #ntr=new/count
  #ntr@x[ntr@x>1]=1
  #ntr@x[is.nan(ntr@x)]=0

  colnames(ntr)=colnames(count)
  rownames(ntr)=rownames(count)

  re=list()
  re$count=count
  re$ntr=ntr

  for (cb in callbacks) cb()

  re=grandR(prefix=genes$prefix,gene.info=gene.info,slots=re,coldata=cols,metadata=list(Description="description"))
  DefaultSlot(re)="count"
  re
}

try.file = function(prefix, possible.suffixes=c("",".tsv",".tsv.gz",".targets/data.tsv.gz"),verbose=FALSE, stop.if.not.exist=TRUE) {
  do.callback=function() {}

  cut.suffix=function(p) {for (s in possible.suffixes) p=gsub(paste0(s,"$"),"",p); p}

  for (suffix in possible.suffixes) {
    if (file.exists(paste0(prefix,suffix))) return(list(file=paste0(prefix,suffix),prefix=cut.suffix(prefix),callback=do.callback))
  }


  if (suppressWarnings(requireNamespace("RCurl",quietly = TRUE))) {
    url=NULL
    for (suffix in possible.suffixes) {
      if (RCurl::url.exists(paste0(prefix,suffix))) {
        url=paste0(prefix,suffix)
        break
      }
    }


    if (!is.null(url)) {
      fn=gsub(".*/","",gsub("\\?.*","",url))
      fn1 = gsub("\\..*","",fn)
      ext=substr(fn,nchar(fn1)+1,nchar(fn))
      file <- tempfile(pattern = fn1,fileext = ext)

      if (verbose) cat(sprintf("Downloading file (url: %s, destination: %s) ...\n",url,file))
      utils::download.file(url, file, quiet=!verbose)
      do.callback=function() {
        if (verbose) cat("Deleting temporary file...\n")
        unlink(file)
      }
      prefix=gsub("\\?.*","",prefix) # cut ?x=y extensions
      return(list(file=file,prefix=cut.suffix(prefix),callback=do.callback))
    }
  } else {
    if (stop.if.not.exist) stop("File not found; If you want to access non-local files directly, please install the RCurl package!")
  }
  if (stop.if.not.exist) stop("File not found!")
  return(NULL)
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
    rownames(m)=names
    m
  }

  # do.callback=function() {}
  #
  # hascurl = suppressWarnings(requireNamespace("RCurl",quietly = TRUE))
  # url=NULL
  # if (hascurl) {
  #   if (RCurl::url.exists(prefix)) {
  #     url=prefix
  #   } else if (RCurl::url.exists(paste0(prefix,".tsv"))) {
  #     url=paste0(prefix,".tsv")
  #   } else if (RCurl::url.exists(paste0(prefix,".tsv.gz"))) {
  #     url=paste0(prefix,".tsv.gz")
  #   } else {
  #     url=NULL
  #   }
  #   if (!is.null(url)) {
  #     file <- tempfile()
  #     if (verbose) cat(sprintf("Downloading file (destination: %s) ...\n",file))
  #     utils::download.file(url, file, quiet=!verbose)
  #     prefix=gsub(".tsv(.gz)?$","",url)
  #     do.callback=function() {
  #       if (verbose) cat("Deleting temporary file...\n")
  #       unlink(file)
  #     }
  #   }
  # }
  #
  # if (is.null(url)) {
  #   if (file.exists(paste0(prefix,".targets/data.tsv.gz"))) {
  #     file=paste0(prefix,".targets/data.tsv.gz")
  #   } else {
  #     file=if (file.exists(prefix)) prefix else paste0(prefix,".tsv")
  #     if (!file.exists(file) && file.exists(paste0(file,".gz"))) file = paste0(file,".gz")
  #     prefix=gsub(".tsv(.gz)?$","",file)
  #   }
  # }
  #
  # if (!file.exists(file) && !hascurl) stop("File not found; If you want to access non-local files directly, please install the RCurl package!")
  # if (!file.exists(file)) stop("File not found!")
  tfile=try.file(prefix,verbose=verbose)
  file=tfile$file
  prefix=tfile$prefix

  if (verbose) cat("Checking file...\n")
  con <- file(file,"r")
  header <- strsplit(readLines(con,n=1),"\t")[[1]]
  close(con)

  count.name=names(slots)[slots=="count"]
  if (length(count.name)!=1) stop("Invalid call to read.grand.internal!")

  if (header[1]!="Gene" || header[2]!="Symbol" || !grepl(paste0(count.name,"$"),header[!header %in% annotations][1])) stop("File is not a GRAND-SLAM output file!")
  conds=gsub(paste0(" ",count.name),"",header[grepl(paste0(" ",count.name),header)])
  ori.conds = conds

  if (!is.null(rename.sample)) conds=sapply(conds,rename.sample) # let's use sapply, we have no idea whether the function works with vectorization

  terms=strsplit(conds[1],".",fixed=TRUE)[[1]]

  if (is.data.frame(design)) {
    design=as.data.frame(design) # in case it's a tibble or similar
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
  data=read.tsv(file,stringsAsFactors=FALSE)
  if (!is.null(rename.sample)) colnames(data)=sapply(colnames(data),rename.sample)

  data$Gene = check.and.make.unique(data$Gene,label="gene names")
  data$Symbol = check.and.make.unique(data$Symbol,ref=data$Gene,label="gene symbols",ref.label = "gene names")

  #if (anyDuplicated(data$Gene)) {
  #  dupp=table(data$Gene)
  #  dupp=names(dupp)[which(dupp>1)]
  #  warning(sprintf("Duplicate gene names (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
  #  data$Gene=make.unique(data$Gene)
  #}
  #if (any(is.na(data$Symbol) | data$Symbol=="")) {
  #  warning(sprintf("No gene symbols (e.g. %s) present, replacing by gene ids!",paste(head(data$Gene[is.na(data$Symbol) | data$Symbol==""]),collapse=",")),call. = FALSE,immediate. = TRUE)
  #  data$Symbol[is.na(data$Symbol) | data$Symbol==""]=data$Gene[is.na(data$Symbol) | data$Symbol==""]
  #}
  #if (anyDuplicated(data$Symbol)) {
  #  dupp=table(data$Symbol)
  #  dupp=names(dupp)[which(dupp>1)]
  #  warning(sprintf("Duplicate gene symbols (e.g. %s) present, making unique!",paste(head(dupp),collapse=",")),call. = FALSE,immediate. = TRUE)
  #  data$Symbol=make.unique(data$Symbol)
  #}

  if (any(data$Gene=="")) stop("Gene name may not be empty!")
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
    a=tomat(data[,cols],data$Gene,cols)
    if (ncol(a)==0 && n=="Conversions") stop("Columns for conversions are not available in GRAND-SLAM output; rerun gedi -e Slam with parameter -full!")
    re[[slots[n]]]=a
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

  checknames=function(n,a){
    if (nrow(a)!=nrow(gene.info)) stop(sprintf("Number of rows do not match for %s!",n))
    if (ncol(a)!=nrow(coldata)) stop(sprintf("Number of columns do not match for %s!",n))
    if (!all(colnames(a)==rownames(coldata))) stop(sprintf("Column names do not match for %s!",n))
    if (!all(rownames(a)==gene.info$Gene)) stop(sprintf("Row names do not match for %s!",n))
  }

  for (slotname in names(re)) checknames(slotname,re[[slotname]])

  coldata$no4sU=no4sU.cols

  tfile$callback()

  # insert name no rep and set this as default condition, if there is no condition field!
  re=grandR(prefix=prefix,gene.info=gene.info,slots=re,coldata=coldata,metadata=list(Description=description))
  DefaultSlot(re)="count"
  re
}
