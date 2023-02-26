
#' @import stats
#' @import utils
#' @import ggplot2
#' @import patchwork
NULL
#> NULL



#' Create a grandR object and retrieve basic information
#'
#' The grandR object contains
#' \itemize{
#'   \item{metadata about the origin (file/url) of the GRAND-SLAM output}
#'   \item{the current state (e.g., what is the current default slot) of the grandR object}
#'   \item{a gene info table (i.e. metadata for the rows of the data matrices)}
#'   \item{a column annotation table (i.e. metadata for the columns of the data matrices)}
#'   \item{several data matrices for read counts, normalized expression values, NTRs, etc. (genes x samples or genes x cells; stored in so-called \emph{slots})}
#'   \item{potentially several analysis output tables (for kinetic modeling, differential gene expression testing)}
#' }
#' Usually, this constructor is not invoked directly (but by \code{\link{ReadGRAND}} or \code{\link{SimulateTimeCourse}}).
#'
#' @param prefix Can either be the prefix used to call GRAND-SLAM with, or the main output file ($prefix.tsv.gz);
#' if the RCurl package is installed, this can also be a URL
#' @param gene.info a data frame with metadata for all genes
#' @param slots A list of matrices representing the slots
#' @param coldata a data frame with metadata for all samples (or cells)
#' @param metadata a metadata list
#' @param analyses the analyses list
#' @param plots the plots list
#' @param parent A parent object containing default values for all other parameters (i.e. all parameters not specified are obtained from this object)
#' @param data,x a grandR object
#' @param columns which columns (i.e. samples or cells) to return (see details)
#' @param reorder reorder all factors in coldata (if columns for subset define a different order)
#' @param f The name of the annotation table according to which the object is split or the new annotation table column name denoting the origin after merging
#' @param list a list of grandR objects
#' @param column.name a new name for the Coldata table to annotate the merged objects
#' @param map named list or vector representing a lookup table (names are current column names)
#' @param fun a function that maps a vector of names to a new vector of names
#' @param s1,s2 column names
#' @param drop unused
#' @param ... further arguments to be passed to or from other methods.
#'
#' @return A grandR object containing the read counts, NTRs, information on the NTR posterior distribution (alpha,beta)
#' and potentially additional information of all genes detected by GRAND-SLAM
#'
#' @details The dimensions (nrow, ncol) of the grandR object are considered to be the dimensions of the data tables,
#' i.e. \code{nrow(data)} provides the number of genes and \code{ncol(data)} the number of samples (or cells).
#'
#' @details Currently, the object is implemented as a list of the above mentioned items. This implementation is subject to change.
#' Make sure to use accessor functions to obtain the information you want.
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment havin the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#'
#' @section Functions:
#' \describe{
#'   \item{Title}{Obtain a useful title for the project (from the prefix parameter)}
#'   \item{dim}{Obtain the dimensions (genes x samples or genes x cells)}
#'   \item{is}{Check whether it is a grandR object}
#'   \item{dimnames}{Obtain the row and column names of this object (genes x samples or genes x cells)}
#'   \item{print}{Print information on this grandR object}
#'   \item{subset}{Create a new grandR object with a subset of the columns (use \code{\link{FilterGenes}} to subset on genes)}
#'   \item{split}{Split the grandR object into a list of multiple grandR objects (according to the levels of an annotation table column)}
#'   \item{RenameColumns}{Rename the column names according to a lookup table (map) or a function (invoked on the current names)}
#'   \item{SwapColumns}{Swap the order of two columns (samples or cells)}
#'   \item{Metadata}{Obtain global metadata}
#'   \item{merge}{Merge several grandR objects into one}
#' }
#'
#'
#' @seealso \link{Slots}, \link{DefaultSlot}, \link{Genes}, \link{GeneInfo}, \link{Coldata}, \link{GetTable}, \link{GetData}, \link{Analyses}, \link{GetAnalysisTable}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#' # this is part of the corona data from Finkel et al.
#' dim(sars)
#' head(rownames(sars))
#'
#' @export
#'
#' @concept grandr
grandR=function(prefix=parent$prefix,gene.info=parent$gene.info,slots=parent$data,coldata=parent$coldata,metadata=parent$metadata,analyses=NULL,plots=NULL,parent=NULL) {

  checknames=function(a){
    if (!all(colnames(a)==rownames(coldata))) stop("Column names do not match!")
    if (!all(rownames(a)==gene.info$Gene)) stop("Row names do not match!")
  }
  for (slot in slots) checknames(slot)
  check.and.make.unique(gene.info$Gene,label = "Gene",do.error=TRUE)
  check.and.make.unique(gene.info$Symbol,label = "Symbol",do.error=TRUE)

  if (!"no4sU" %in% colnames(coldata)) {
    warning("No no4sU entry in coldata, assuming all samples/cells as 4sU treated!")
    coldata$no4sU=FALSE
  }

  if (!"Name" %in% colnames(coldata)) {
    coldata=cbind(data.frame(Name=factor(rownames(coldata),levels=rownames(coldata))),coldata)
  }

  if ("Condition" %in% colnames(coldata)) {
    coldata$Condition=as.factor(coldata$Condition)
  }


  info=list()
  info$prefix=prefix
  info$gene.info=gene.info
  info$data=slots
  info$coldata=coldata
  info$metadata=metadata
  info$analysis=analyses
  info$plots=plots
  class(info)="grandR"
  info
}

as.grandR=function(mat,slot="count",coldata=MakeColdata(colnames(mat)),gene.info=rownames(mat)) {
  if (!is.data.frame(gene.info) && !is.matrix(gene.info)) gene.info=data.frame(Gene=gene.info,Symbol=gene.info,Type="Unknown")
  gene.info = as.data.frame(gene.info)

  if (!all(c("Gene","Symbol","Type") %in% names(gene.info))) stop("Gene info table has to have columns Gene, Symbol and Type!")

  rownames(mat)=gene.info$Gene
  ntr=mat
  ntr[,]=0
  slots=list()
  slots[[slot]]=mat
  slots$ntr=ntr

  r=grandR(prefix="",gene.info = gene.info, slots = slots, coldata = coldata,metadata=list(Description="Converted from matrix",`GRAND-SLAM version`=0,Output="dense"))
  DefaultSlot(r)=slot
  r
}

#' @rdname grandR
#' @export
Title=function(data) {
  x=strsplit(data$prefix,"/")[[1]]
  x[length(x)]
}

#' @rdname grandR
#' @export
IsSparse=function(data) !is.matrix(data$data$count)


#' @rdname grandR
#' @export
dim.grandR=function(x) c(dim(x$gene.info)[1],dim(x$coldata)[1])
#' @rdname grandR
#' @export
is.grandR <- function(x) inherits(x, "grandR")
#' @rdname grandR
#' @export
dimnames.grandR=function(x) dimnames(x$data[[DefaultSlot(x)]])
#' @rdname grandR
#' @export
print.grandR=function(x,...) {
  cat(
  sprintf("grandR:\nRead from %s\n%d genes, %d samples/cells\nAvailable data slots: %s\nAvailable analyses: %s\nAvailable plots: %s\nDefault data slot: %s\n",
          x$prefix,
          nrow(x),
          ncol(x),
          paste(Slots(x),collapse=","),
          paste(Analyses(x),collapse=","),
          paste(unlist(Plots(x)),collapse=","),
          DefaultSlot(x))
)
}
#' @rdname grandR
#' @export
Metadata=function(x,...) {x$metadata}

#' Internal function to apply functions to all slots etc.
#'
#' @param data a grandR object
#' @param fun apply this function to each data slot (i.e. it receives each data matrix)
#' @param fun.gene.info apply this function to the gene info table
#' @param fun.coldata apply this function to the column annotation table
#' @param ... passed further to fun, fun.gene.info and fun.coldata
#'
#' @details The additional parameters are provided to each of the functions.
#' @return A new grandR object
#' @concept helper
data.apply=function(data,fun,fun.gene.info=NULL,fun.coldata=NULL,...) {
  re=list()
  for (l1 in names(data$data)) {
    re[[l1]]=fun(data$data[[l1]],...)
  }
  ngene.info=if (!is.null(fun.gene.info)) fun.gene.info(data$gene.info,...) else data$gene.info
  ncoldata=droplevels(if (!is.null(fun.coldata)) fun.coldata(data$coldata,...) else data$coldata)
  analysis=NULL
  if (!is.null(data$analysis)) {
    map=setNames(1:nrow(data$gene.info),data$gene.info$Gene)
    analysis=lapply(data$analysis,function(d) d[map[as.character(ngene.info$Gene)],,drop=FALSE])
  }
  grandR(parent=data,gene.info=ngene.info,slots=re,coldata=ncoldata,analyses = analysis,plots=data$plots)
}

#' @rdname grandR
#' @export
subset.grandR=function(x,columns,reorder=TRUE,...) {
  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(x) else eval(columns,Coldata(x),parent.frame())
  columns=Columns(x,columns,reorder=TRUE)

  dr=if(!reorder) droplevels else function(x) {
    ix <- vapply(x, is.factor, NA)
    x[ix] <- lapply(x[ix], function(v) factor(v,levels=as.character(unique(v))))
    x
  }
  data.apply(x,function(m) m[,columns],fun.coldata = function(t){
    dr(t[columns,])
  })
}

#' @rdname grandR
#' @export
split.grandR=function(x,f=Design$Condition,drop=FALSE,...) {
  col=as.factor(x$coldata[[f]])
  setNames(lapply(levels(col),function(c) {re=subset(x,col==c);
  re$coldata[[Design$Origin]]=c; re }),levels(col))
}

#' @rdname grandR
#' @export
RenameColumns=function(data,map=NULL,fun=NULL) {
  if (!is.null(fun)) {
    map=setNames(sapply(colnames(data),fun),colnames(data))
  }
  names=rownames(data$coldata)
  names[names %in% names(map)]=unlist(map[names[names %in% names(map)]])
  rownames(data$coldata)=names
  data$coldata$Name=factor(names,levels = names)
  data.apply(data,function(m) {colnames(m)=names; m})
}
#' @rdname grandR
#' @export
SwapColumns=function(data,s1,s2) {
  i1=if(is.numeric(s1)) s1 else which(rownames(data$coldata)==s1)
  i2=if(is.numeric(s2)) s2 else which(rownames(data$coldata)==s2)
  return(data.apply(data,function(t) {
    tmp=t[,i1]
    t[,i1]=t[,i2]
    t[,i2]=tmp
    t
  }))
}

#' @rdname grandR
#' @export merge.grandR
#' @export
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
    #if (any(colnames(add$coldata)!=colnames(re$coldata))) stop("Data sets have distinct coldata columns!")
    if (!all(names(add$data) %in% names(re$data))) stop("Data sets must have the same data tables!")

    # merge coldata paying attention to columns and factor levels
    cd=NULL
    for (common in intersect(names(re$coldata),names(add$coldata))) {
      if(is.factor(re$coldata[[common]])) {
        r = c(as.character(re$coldata[[common]]),as.character(add$coldata[[common]]))
        r=factor(r,levels=union(levels(re$coldata[[common]]),levels(add$coldata[[common]])))
      } else {
        r = c(re$coldata[[common]],add$coldata[[common]])
      }
      df=setNames(data.frame(r),common)
      cd=if (is.null(cd)) df else cbind(cd,df)
    }
    for (re.only in setdiff(names(re$coldata),names(add$coldata))) {
      r=c(re$coldata[[re.only]],rep(NA,nrow(add$coldata)))
      df=setNames(data.frame(r),re.only)
      cd=if (is.null(cd)) df else cbind(cd,df)
    }
    for (add.only in setdiff(names(add$coldata),names(re$coldata))) {
      r=c(rep(NA,nrow(re$coldata)),add$coldata[[add.only]])
      df=setNames(data.frame(r),add.only)
      cd=if (is.null(cd)) df else cbind(cd,df)
    }
    rownames(cd)=c(rownames(re$coldata),rownames(add$coldata))
    re$coldata=cd

    for (n in names(re$data)) re$data[[n]]=cbind(re$data[[n]],add$data[[n]])

    # add potential additional gene annotations
    gi=GeneInfo(add)
    for (n in names(gi)) {
      if (!n %in% names(GeneInfo(re))) {
        GeneInfo(re,column = n)=gi[[n]]
      } else if (!all(GeneInfo(re)[[n]]==gi[[n]])) {
        GeneInfo(re,column = paste0(n,".",names(list)[i]))=gi[[n]]
        if (is.null(names(list)[i])) stop("Found columns by the same name in the gene info tables with distinct content; please specify names when calling merge!")
      }
    }

    # add analyses
    for (ana in Analyses(add))
    re=AddAnalysis(re,name=ana,add$analysis[[ana]])
  }
  re
}


#' Get or set the default slot for a grandR object.
#'
#' The default slot is used by default by many functions including
#' \code{\link{GetData}},\code{\link{GetTable}} or \code{\link{FitKinetics}}
#'
#' @param data A grandR object
#' @param value the name of the new default slot
#'
#' @return Either the name of the default slot for DefaultSlot(data)
#' or the grandR data object having the new default slot
#'
#' @details The default slot can be set either by \code{data<-DefaultSlot(data,"norm")} or by \code{DefaultSlot(data)<-"norm"}.
#'
#' @seealso \link{Slots}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#'
#' DefaultSlot(sars)
#' sars <- Normalize(sars)     # default behavior is to update the default slot
#' DefaultSlot(sars)
#' DefaultSlot(sars)="count"
#'
#' @export
#'
#' @concept grandr
DefaultSlot <- function(data,value=NULL) {
  if (is.null(value)) data$metadata$default.slot else {
    DefaultSlot(data)=value
    data
  }
}

#' @rdname DefaultSlot
#' @export
`DefaultSlot<-` <- function(data, value) {
  if (!value %in% names(data$data)) stop("Invalid slot name!")
  data$metadata$default.slot=value
  data
}

#' Slot functions
#'
#' Get slot names and add or remove slots
#'
#' @param data A grandR object
#'
#' @return Either the slot names or a grandR data with added/removed slots
#'
#' @seealso \link{DefaultSlot}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#'
#' sars <- Normalize(sars)     # default behavior is to update the default slot
#' sars
#' sars <- DropSlot(sars,"norm")
#' sars                        # note that the defauls slot reverted to count
#'
#' @describeIn Slots Obtain the slot names
#' @export
#'
#' @concept grandr
Slots=function(data) names(data$data)

#' @describeIn Slots Remove one or several slots from this grandR object
#' @param pattern a regular expression matched against slot names
#' @export
DropSlot=function(data,pattern=NULL) {
  if (is.null(pattern)) {
    stop("Cannot drop all slots!")
  } else {
    data$data=data$data[!grepl(pattern,names(data$data))]
  }
  if (!DefaultSlot(data) %in% names(data$data)) DefaultSlot(data)=names(data$data)[1]
  data
}
#' @describeIn Slots Add an additional slot to this grandR object
#' @param name the slot name
#' @param matrix the data matrix for the new slot
#' @param set.to.default set the new slot as the default slot?
#' @param warn issue a warning if the slot name already exists and is overwritten
#' @export
AddSlot=function(data,name,matrix,set.to.default=FALSE,warn=TRUE) {
  if (!is.matrix(matrix)) stop("Must be a matrix!")
  if (!all(colnames(matrix)==colnames(data$data$count))) stop("Column names do not match!")

  rownames(matrix)=Genes(data,rownames(matrix),use.symbols=FALSE)
  missing=setdiff(rownames(data$data$count),rownames(matrix))
  if (length(missing>0)) {
    warning(sprintf("Could not find all genes in matrix, setting to 0 (n=%d missing, e.g. %s)!",length(missing),paste(head(missing,5),collapse=",")))
    matrix = rbind(matrix,matrix(0,nrow=length(missing),ncol=ncol(matrix),dimnames=list(missing,colnames(matrix))))
  }
  matrix=matrix[rownames(data$data$count),]
  if (!all(rownames(matrix)==rownames(data$data$count))) stop("Row names do not match!")

  if (grepl(".",name,fixed=TRUE)) stop("Name may not contain a dot!")
  if (!is.null(data$data[[name]])) warning(sprintf("Slot %s already exists, overwriting!",name))
  data$data[[name]]=matrix
  if (set.to.default) DefaultSlot(data)=name
  data
}

#' Get or set the conditions in the column annotation table.
#'
#' The conditions column from the column annotation table is used by several functions to stratify the columns (samples or cells)
#' during the analysis (e.g. to estimate separate kinetic parameters with \code{\link{FitKinetics}} or it is used as covariate for
#' \code{\link{LFC}} or \code{\link{LikelihoodRatioTest}}). For that reason there are special functions to set and get this column.
#'
#' @param data A grandR object
#' @param value Either a vector of column names from the column annotation table, or the condition names themselves
#'
#' @details If the conditions column does not exist (or has been set to NULL), all analysis functions will work without stratifying samples or cells.
#' The condition can also be set up directly when loading data, by using \emph{Condition} as one of the design vector entries (see below).
#'
#' @details The condition can be set either by \code{data<-Condition(data,names)} or by \code{Condition(data)<-names}.
#' @return Either the values of the condition column for Condition(data) or the grandR data object having the new condition column
#'
#' @seealso \link{Coldata}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#'
#' Condition(sars)
#' Condition(sars) <- c("Cell","duration.4sU.original")
#' Condition(sars)
#'
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' Condition(sars)
#'
#' @export
#'
#' @concept grandr
Condition <- function(data,value=NULL) {
  if (is.null(value)) data$coldata$Condition else {
    Condition(data)<-value
    data
  }
}
#' @rdname Condition
#' @export
`Condition<-` <- function(data, value) {
  if (is.null(value)) {
    data$coldata$Condition=NULL
  } else if (length(value)==1 && value=="") {
    stop("Empty string is not allowed as condition!")
  } else if (all(value %in% names(data$coldata))) {
    data$coldata$Condition <- interaction(data$coldata[value],drop=TRUE)
  } else {
    data$coldata$Condition <- as.factor(value)
  }
  data
}





#' Gene and sample (or cell) names
#'
#' Get the genes and sample (or cell) names for a grandR object, or add an additional gene annotation column
#'
#' @param data A grandR object
#' @param use.symbols obtain the gene symbols instead of gene names
#' @param genes which genes to use
#' @param regex treat genes as a regex, and return all that match
#' @param columns which columns (i.e. samples or cells) to return (see details)
#' @param reorder if TRUE, do not enforce the current order of columns
#'
#' @details The genes are either the (often unreadable) gene ids (e.g. Ensembl ids), or the symbols.
#'
#' @details \code{Genes(data,use.symbols=FALSE)} it the same as \code{rownames(data)}, and \code{Columns(data)} is the same as \code{colnames(data)}
#'
#' @details If both column and value are specified for \code{GeneInfo}, a new column is added to the gene annotation table
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @return Either the gene or column names of the grandR data object, or the columns of an analysis table in the grandR object
#'
#' @seealso \link{Coldata}, \link{GeneInfo}, \link{Analyses}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#'
#' all(Genes(sars,use.symbols = FALSE)==rownames(sars))
#' all(Columns(sars)==colnames(sars))
#'
#'
#' @export
#'
#' @concept grandr
Genes=function(data, genes=NULL, use.symbols=TRUE,regex=FALSE) data$gene.info[[if (use.symbols) "Symbol" else "Gene"]][ToIndex(data,genes,regex=regex)]
#' @rdname Genes
#' @export
Columns=function(data,columns=NULL, reorder=FALSE) {
  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())

  columns=unname(setNames(rownames(Coldata(data)),rownames(Coldata(data)))[columns])
  if (reorder) columns else rownames(Coldata(data))[ rownames(Coldata(data)) %in% columns]  # preserve order!
}

#' Get the gene annotation table or add additional columns to it
#'
#' The gene annotation table contains meta information for the rows of a grandR object.
#' When loaded from the GRAND-SLAM output, this this contains gene ids, gene symbols, the
#' transcript length and the type.
#'
#' @param data A grandR object
#' @param column The name of the additional annotation column
#' @param value The additional annotation per gene
#'
#' @details New columns can be added either by \code{data<-GeneInfo(data,name,values)} or by \code{GeneInfo(data,name)<-values}.
#'
#' @return Either the gene annotation table or a new grandR object having an updated gene annotation table
#'
#' @seealso \link{Genes}, \link{Coldata}, \link{ReadGRAND}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#'
#' head(GeneInfo(sars))
#' GeneInfo(sars,"LengthCategory")<-cut(GeneInfo(sars)$Length,c(0,1500,2500,Inf),
#'                                           labels=c("Short","Medium","Long"))
#' table(GeneInfo(sars)$LengthCategory)
#'
#' @export
#'
#' @concept grandr
GeneInfo=function(data,column=NULL,value=NULL) {
  if (is.null(column)) {
    data$gene.info
  } else if (is.null(value)) {
    setNames(data$gene.info[[column]],data$gene.info$Symbol)
  } else {
    data$gene.info[[column]]=value
    data
  }
}
#' @rdname GeneInfo
#' @export
`GeneInfo<-` <- function(data, column, value) {
  data$gene.info[[column]]=value
  data
}


#' Update symbols using biomaRt
#'
#' If your input files only contained ENSEMBL ids, use this to add gene symbols!
#'
#' @param data a grandR object
#' @param species the species the genes belong to (eg "Homo sapiens"); can be NULL, then the species is inferred from gene ids (see details)
#' @param current.value What it the current value in the symbols field?
#'
#' @return a grandR object with updated symbol names
#'
#' @details If no species is given, a very simple automatic inference is done, which will only work when having human or mouse ENSEMBL identifiers as gene ids.
#' If you need to specify species, it must be one of \code{biomaRt::listDatasets(biomaRt::useMart("ensembl"))$dataset}!
#'
#' @details Current.value must be one of \code{biomaRt::listAttributes(biomaRt::useMart("ensembl"))$name}!
#'
#' @export
#'
#' @concept grandr
UpdateSymbols = function(data,species=NULL,current.value="ensembl_gene_id") {

  checkPackages(c("biomaRt"))

  if (is.null(species)) {
    if (sum(grepl("ENSG0",Genes(data,use.symbols=FALSE)))>nrow(data)/2) species="hsapiens_gene_ensembl"
    if (sum(grepl("ENSMUSG0",Genes(data,use.symbols=FALSE)))>nrow(data)/2) species="mmusculus_gene_ensembl"
  }
  if (is.null(species)) stop("Cannot recognize species! Specify one of biomaRt::listDatasets(biomaRt::useMart(\"ensembl\"))$dataset")


  mart <- biomaRt::useDataset(species, biomaRt::useMart("ensembl"))
  df <- biomaRt::getBM(filters= "ensembl_gene_id", attributes= c(current.value,"hgnc_symbol"),values=Genes(data,use.symbols = TRUE),mart= mart)
  map=setNames(df$hgnc_symbol,df[[current.value]])

  GeneInfo(data,"Symbol")=check.and.make.unique(map[as.character(Genes(data,use.symbols = TRUE))],ref=as.character(Genes(data,use.symbols = TRUE)),label="symbols",ref.label=current.value)
  data
}


#' Get the column annotation table or add additional columns to it
#'
#' The columns of a grandR object are samples or cells.
#' The column annotation table contains meta information for the columns of a grandR object.
#' When loaded from the GRAND-SLAM output, this this constructed from the sample/cell names by
#' \code{\link{MakeColdata}}
#'
#' @param data A grandR object
#' @param column The name of the additional annotation column; can also be a data frame (then value is ignored and the data frame is added)
#' @param value The additional annotation per sample or cell
#'
#' @details A new column can be added either by \code{data<-Coldata(data,name,values)} or by \code{Coldata(data,name)<-values}.
#'
#' @details Several new columns can be added by \code{data<-Coldata(data,df)} where df is either a data frame or matrix.
#'
#' @details The column named \emph{Condition} has a special meaning in this table: It is used by several functions to stratify the columns
#' during the analysis (e.g. to estimate separate kinetic parameters with \code{\link{FitKinetics}} or it is used as covariate for
#' \code{\link{LFC}} or \code{\link{LikelihoodRatioTest}}). For that reason there are special functions to set and get this column.
#'
#' @return Either the column annotation table or a new grandR object having an updated column annotation table
#'
#' @seealso \link{GeneInfo}, \link{MakeColdata}, \link{Condition}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#'
#' head(GeneInfo(sars))
#' GeneInfo(sars,"LengthCategory")<-cut(GeneInfo(sars)$Length,c(0,1500,2500,Inf),
#'                                           labels=c("Short","Medium","Long"))
#' table(GeneInfo(sars)$LengthCategory)
#'
#' @export
#'
#' @concept grandr
Coldata=function(data,column=NULL,value=NULL) {
  if (is.null(column)) {
    data$coldata
  } else if (is.data.frame(column)||is.matrix(column)) {
    data$coldata=cbind(data$coldata,column)
    data
  } else if (is.null(value)) {
    setNames(data$coldata[[column]],rownames(data$coldata))
  } else {
    data$coldata[[column]]=value
    data
  }
}
#' @rdname Coldata
#' @export
`Coldata<-` <- function(data, column, value) {
  data$coldata[[column]]=value
  data
}


#' Internal functions to check for a valid analysis or slot names.
#'
#' @param data a grandR object
#' @param analyses a regex to be matched to analysis names
#' @param regex interpret as regular expression
#' @param slot a slot name
#' @param mode.slot a mode.slot
#' @param allow.ntr whether to allow for the value "ntr" (and throw an error in case)
#'
#' @details A mode.slot is a mode followed by a dot followed by a slot name, or just a slot name. A mode is either \emph{total}, \emph{new} or \emph{old}.
#'
#' @return Whether or not the given name is valid and unique for the grandR object
#'
#' @concept helper
check.analysis=function(data,analyses,regex) {
  if (!regex) return(is.logical(analyses) || is.numeric(analyses) || all(analyses %in% Analyses(data)))
  sapply(analyses,function(pattern) any(grepl(pattern,Analyses(data),fixed=!regex)))
}
#' @rdname check.analysis
check.slot=function(data,slot,allow.ntr=TRUE) {
  if (!allow.ntr && slot=="ntr") return(FALSE)
  slot %in% names(data$data)
}
#' @rdname check.analysis
check.mode.slot=function(data,mode.slot,allow.ntr=TRUE) {
  sapply(strsplit(mode.slot,".",fixed=TRUE),function(spl) {
    if (length(spl)>2 || length(spl)==0) return(FALSE)
    if (length(spl)==1) check.slot(data,spl,allow.ntr=allow.ntr) else tolower(substr(spl[1],1,1)) %in% c("t","n","o") && check.slot(data,spl[2],allow.ntr = FALSE)
  })
}

#' Internal functions to parse mode.slot strings
#'
#' @param data a grandR object
#' @param mode.slot a mode.slot
#' @param allow.ntr whether to allow for the value "ntr" (and throw an error in case)
#'
#' @details A mode.slot is a mode followed by a dot followed by a slot name, or just a slot name. A mode is either \emph{total}, \emph{new} or \emph{old}
#'
#' @return a named list with elements mode and slot (or only slot in case of \emph{ntr},\emph{alpha} or \emph{beta})
#'
#' @concept helper
get.mode.slot=function(data,mode.slot,allow.ntr=TRUE) {
  if (length(mode.slot)!=1) stop("mode.slot must be a vector of length 1")
  if (!check.mode.slot(data,mode.slot,allow.ntr = allow.ntr)) stop("Invalid mode.slot")
  tno="t"
  spl=strsplit(mode.slot,".",fixed=TRUE)[[1]]
  if (length(spl)>1) {
    tno=spl[1]
    mode.slot=spl[2]
    if (mode.slot %in% c("ntr","alpha","beta")) stop(paste0(tno," may not be used with ntr, alpha or beta"))
  }
  tno = switch(tolower(substr(tno,1,1)),t="total",n="new",o="old",stop(paste0("Mode ",tno," unknown!")))
  if (mode.slot %in% c("ntr","alpha","beta")) list(slot=mode.slot) else list(mode=tno,slot=mode.slot)
}

#' Obtain the indices of the given genes
#'
#' Genes can be referred to by their names, symbols, row numbers in the gene table, or a logical vector referring to the gene table rows.
#' This function accepts all these possibilities and returns the row number in the gene table for the given genes,
#'
#' @param data The grandR object
#' @param gene A vector of genes. Can be either numeric indices, gene names, gene symbols or a logical vector
#' @param regex Treat gene as a regex and return all that match
#'
#' @return Numeric indices corresponding to the given genes
#'
#' @seealso \link{GeneInfo}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#' ToIndex(sars,c("MYC"))
#' ToIndex(sars,GeneInfo(sars)$Symbol=="MYC")
#'
#' @export
#'
#' @concept helper
ToIndex=function(data,gene,regex=FALSE) {
  if (any(is.na(gene))) {
    warning("There were NA genes, removed!");
    gene=gene[!is.na(gene)]
  }
  if (is.factor(gene)) gene = as.character(gene) # god, I hate factors


  if (regex) gene=grepl(gene,data$gene.info$Gene)|grepl(gene,data$gene.info$Symbol)
  if (is.null(gene)) return(1:nrow(data))
  if (is.numeric(gene)) return(gene)
  if (is.logical(gene) && length(gene)==nrow(data)) return(which(gene))
  if (all(gene %in% data$gene.info$Gene)) return(setNames(1:nrow(data),data$gene.info$Gene)[gene])
  if (all(gene %in% data$gene.info$Symbol)) return(setNames(1:nrow(data),data$gene.info$Symbol)[gene])
  if (sum(gene %in% data$gene.info$Gene) > sum(gene %in% data$gene.info$Symbol)) {
    mis=setdiff(gene,data$gene.info$Gene)
    warning(sprintf("Could not find given genes (n=%d missing, e.g. %s)!",length(mis),paste(head(mis,5),collapse=",")))
    return(setNames(1:nrow(data),data$gene.info$Gene)[intersect(gene,data$gene.info$Gene)])
  }
  mis=setdiff(gene,data$gene.info$Symbol)
  warning(sprintf("Could not find given genes (n=%d missing, e.g. %s)!",length(mis),paste(head(mis,5),collapse=",")))
  return(setNames(1:nrow(data),data$gene.info$Symbol)[intersect(gene,data$gene.info$Symbol)])
}

#' Obtain a genes x values table
#'
#' This is the main function to access slot data for all genes as a large matrix. If data from a particular gene (or a small set of genes)
#' must be retrieved, use the \code{\link{GetData}} function. For analysis results, use the \code{\link{GetAnalysisTable}} function.
#'
#' @param data A grandR object
#' @param type Either a mode.slot (see details) or a regex to be matched against analysis names. Can also be a vector
#' @param columns A vector of columns (either condition/cell names if the type is a mode.slot, or names in the output table from an analysis; use \link{Columns}(data,<analysis>) to learn which columns are available); all condition/cell names if NULL
#' @param genes Restrict the output table to the given genes
#' @param ntr.na For columns representing a 4sU naive sample, should types \emph{ntr},\emph{new.count} and \emph{old.count} be 0,0 and count (ntr.na=FALSE; can be any other slot than count) or NA,NA and NA (ntr.na=TRUE)
#' @param gene.info Should the table contain the \link{GeneInfo} values as well (at the beginning)?
#' @param summarize Should replicates by summarized? see details
#' @param prefix Prepend each column in the output table (except for the gene.info columns) by the given prefix
#' @param name.by A column name of \link{Coldata}(data). This is used as the rownames of the output table
#'
#' @return A data frame containing the desired values
#'
#' @details This is a convenience wrapper for \link{GetData} (values from data slots) and \link{GetAnalysisTable} (values from analyses).
#' Types can refer to any of the two (and can be mixed). If there are types from both data and analyses, columns must be NULL.
#' Otherwise columns must either be condition/cell names (if type refers to one or several data slots), or regular expressions
#' to match against the names in the analysis tables.
#'
#' @details Columns definitions for data slots can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @details To refer to data slots via \code{type}, the mode.slot syntax can be used: Each name is either a data slot, or one of (new,old,total)
#' followed by a dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to obtain the \emph{new counts}.
#'
#' @details The summarization parameter can only be specified if columns is NULL. It is either a summarization matrix (\link{GetSummarizeMatrix}) or
#' TRUE (in which case \link{GetSummarizeMatrix}(data) is called). If there a NA values, they are imputed as the mean per group!
#'
#' @seealso \link{GetData},\link{GetAnalysisTable},\link{DefaultSlot},\link{Genes},\link{GetSummarizeMatrix}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' sars <- Normalize(FilterGenes(sars))
#'
#' head(GetTable(sars))
#' # DefaultSlot values, i.e. size factor normalized read counts for all samples
#' head(GetTable(sars,summarize=TRUE))
#' # DefaultSlot values averaged over the two conditions
#' head(GetTable(sars,type="new.count",columns=!no4sU))
#' # Estimated counts for new RNA for all samples with 4sU
#'
#' sars<-LFC(sars,contrasts=GetContrasts(sars,group = "duration.4sU"))
#' head(GetAnalysisTable(sars,columns="LFC"))
#' # Estimated fold changes SARS vs Mock for each time point
#'
#'
#'
#' @export
#'
#' @concept data
GetTable=function(data,type=DefaultSlot(data),columns=NULL,genes=Genes(data),ntr.na=TRUE,gene.info=FALSE,summarize=NULL,prefix=NULL,name.by="Symbol") {

  genes=ToIndex(data,genes)

  mode.slot=check.mode.slot(data,type)
  analysis=check.analysis(data,type,TRUE) & !mode.slot
  if (!all(analysis|mode.slot)) {
    r=NULL
    if (length(type)==1) {
      r = GetAnalysisTable(data,genes=genes,gene.info = FALSE,name.by=name.by)
      use=grepl(type,names(r))
      r=r[,use]
    }
    if (ncol(r)==0) stop(sprintf("Type %s is neither a mode.slot nor an analysis name!",paste(type[!analysis&!mode.slot],collapse=",")))
  }
  else {

    # obtain mode.slot data
    r1=NULL
    if (any(mode.slot)) {

      columns=substitute(columns)
      cols=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
      cols=Columns(data,cols)

      if (!is.null(summarize)) {
        if (is.logical(summarize) && length(summarize)==1 && !summarize) {
          summarize=NULL
        } else {
          if (is.logical(summarize) && length(summarize)==1 && summarize) summarize=GetSummarizeMatrix(data)
          summarize=summarize[cols,]
          summarize=summarize[,colSums(summarize!=0)>1,drop=FALSE]
        }
      }

      for (tt in type[mode.slot]) {
        rtt=as.data.frame(t(GetData(data,tt,columns=cols,genes,ntr.na = ntr.na,coldata=FALSE, by.rows=FALSE, name.by = name.by)))
        names(rtt)=cols
        if (!is.null(summarize)) {
          #rtt=as.data.frame(as.matrix(rtt) %*% summarize)  ## this is without imputation, which is really bad!
          mrtt=as.matrix(rtt)
          mrtt=apply(summarize,2,function(cc) {
            h=mrtt[,cc!=0,drop=FALSE]
            cc=cc[cc!=0]
            apply(h,1,function(v) { v[is.na(v)] = mean(v,na.rm = TRUE); sum(v*cc)})
          })
          if (!is.matrix(mrtt)) mrtt=matrix(mrtt,nrow=1)
          mrtt[is.nan(mrtt)]=NA
          rownames(mrtt)=rownames(rtt)
          colnames(mrtt)=colnames(summarize)
          rtt=as.data.frame(mrtt)
        }
        if (length(type[mode.slot])>1) names(rtt)=paste0(names(rtt),".",tt)
        r1=if(is.null(r1)) rtt else cbind(r1,rtt)
      }
    }

    # check that columns is only used if type is either completely analysis or mode.slot
    if (!is.null(columns) && sum(mode.slot)>0 && sum(analysis)>0) stop("Columns can only be specified if type either refers to mode.slots or analyses")

    # obtain analysis data
    r2=NULL
    if(any(analysis)) {
      r2 = GetAnalysisTable(data,type[analysis],columns = columns,genes=genes,gene.info = FALSE,name.by=name.by)
    }

    # reorder according to order in type
    if (is.null(r1)) {
      r=r2
    } else if (is.null(r2)) {
      r=r1
    } else {
      r = cbind(r1,r2)
      r[,mode.slot]=r1
      r[,analysis]=r2
    }
  }

  # add necessary stuff
  if (!is.null(prefix)) colnames(r)=paste0(prefix,colnames(r))
  if (is.logical(gene.info) && gene.info) r=cbind(GeneInfo(data)[ToIndex(data,genes),],r)
  if (is.character(gene.info)) r=cbind(GeneInfo(data)[ToIndex(data,genes),-which(!names(GeneInfo(data)) %in% gene.info)],r)

  rownames(r)=data$gene.info[[name.by]][ToIndex(data,genes)]
  r
}


#' Obtain a genes x values table as a large matrix
#'
#' This is the main function to access slot data for all genes as a (potentially sparse) matrix.
#'
#' @param data A grandR object
#' @param mode.slot Which kind of data to access (see details)
#' @param columns which columns (i.e. samples or cells) to return (see details)
#' @param genes Restrict the output table to the given genes
#' @param name.by A column name of \link{Coldata}(data). This is used as the rownames of the output table
#' @param summarize Should replicates by summarized? see details
#'
#' @return A (potentially) sparse matrix containing the desired values
#'
#' @details To refer to data slots, the mode.slot syntax can be used: It is either a data slot, or one of (new,old,total) followed by a dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to obtain the \emph{new counts}.
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment havin the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @details The summarization parameter can only be specified if columns is NULL. It is either a summarization matrix (\link{GetSummarizeMatrix}) or
#' TRUE (in which case \link{GetSummarizeMatrix}(data) is called). If there a NA values, they are imputed as the mean per group!
#'
#' @seealso \link{GetData},\link{GetAnalysisTable},\link{DefaultSlot},\link{Genes},\link{GetSummarizeMatrix}
#'
#' @export
#'
#' @useDynLib grandR, .registration = TRUE
#' @concept data
GetMatrix=function(data,mode.slot=DefaultSlot(data),columns=NULL,genes=Genes(data),name.by="Symbol",summarize=NULL) {

  if (!all(check.mode.slot(data,mode.slot))) stop(sprintf("mode.slot %s unknown!",paste(mode.slot[!check.mode.slot(data,mode.slot)],collapse=",")))
  if (length(mode.slot)!=1) stop("Specify exactly one mode.slot!")

  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)

  genes=ToIndex(data,genes)

  tno="t"
  spl=strsplit(mode.slot,".",fixed=TRUE)[[1]]
  if (length(spl)>1) {tno=spl[1]; mode.slot=spl[2];}

  re=data$data[[mode.slot]][genes,columns,drop=FALSE]
  rownames(re)=Genes(data,genes,use.symbols = name.by=="Symbol")

  conv=function(v) { mode(v)="integer"; v}

  if (is.matrix(re)) {  # plain old R matrix!
    mf = switch(tolower(substr(tno,1,1)),t=1,n=as.matrix(data$data$ntr[genes,columns,drop=FALSE]),o=1-as.matrix(data$data$ntr[genes,columns,drop=FALSE]),stop(paste0(mode.slot," unknown!")))
      mf[is.na(mf)]=if(tolower(substr(tno,1,1))=="n") 0 else 1
    re=re*mf
    if (mode.slot=="count") {
      mode(re) <- "integer"
    } else if (mode.slot=="ntr") {
      re[is.na(re)]=0
    }
  } else {
    if (tolower(substr(tno,1,1))=="n") {

      X <- Matrix::summary(re)
      Y <- Matrix::summary(data$data$ntr[genes,columns,drop=FALSE])
      R=.Call('fastsparsematcompmult',X$i,X$j,X$x,Y$i,Y$j,Y$x)

      re=(Matrix::sparseMatrix(i=R[[1]], j=R[[2]], x=conv(round(R[[3]])),dims=dim(re),
                                  dimnames=dimnames(re)))

      # all that have zero in ntr matrix will be zero, so this is fine
      #sX <- Matrix::summary(re)
      #sY <- Matrix::summary(data$data$ntr[genes,columns,drop=FALSE])
      #sRes <- merge(sX, sY, by=c("i", "j"))
      #return(Matrix::sparseMatrix(i=sRes[,1], j=sRes[,2], x=conv(sRes[,3]*sRes[,4]),dims=dim(re),
      #                            dimnames=dimnames(re)))
    } else if (tolower(substr(tno,1,1))=="o") {
      X <- Matrix::summary(re)
      Y <- Matrix::summary(data$data$ntr[genes,columns,drop=FALSE])
      R=.Call('fastsparsematcompmult1m',X$i,X$j,X$x,Y$i,Y$j,Y$x)

      re=(Matrix::sparseMatrix(i=R[[1]], j=R[[2]], x=conv(round(R[[3]])),dims=dim(re),
                                  dimnames=dimnames(re)))

      #sX <- Matrix::summary(re)
      #sY <- Matrix::summary(data$data$ntr[genes,columns,drop=FALSE])
      #sRes <- merge(sX, sY, by=c("i", "j"),all.x=TRUE)
      #sRes[is.na(sRes[,4]),4]=0
      #sRes[,4]=1-sRes[,4]
      #sRes=sRes[sRes[,4]>0,]
      #return(Matrix::sparseMatrix(i=sRes[,1], j=sRes[,2], x=conv(sRes[,3]*sRes[,4]),dims=dim(re),
      #                            dimnames=dimnames(re)))
    }
  }

  if (!is.null(summarize)) {
    if (is.logical(summarize) && length(summarize)==1 && !summarize) {
      summarize=NULL
    } else {
      if (is.logical(summarize) && length(summarize)==1 && summarize) summarize=GetSummarizeMatrix(data)
      summarize=summarize[columns,]
      summarize=summarize[,colSums(summarize!=0)>1,drop=FALSE]
    }

    re=apply(summarize,2,function(cc) {
      h=re[,cc!=0,drop=FALSE]
      cc=cc[cc!=0]
      apply(h,1,function(v) { v[is.na(v)] = mean(v,na.rm = TRUE); sum(v*cc)})
    })
    if (!is.matrix(re)) re=matrix(re,nrow=1)
    re[is.nan(re)]=NA
    rownames(re)=Genes(data,genes,use.symbols = name.by=="Symbol")
    colnames(re)=colnames(summarize)
  }

  re

}

#' Obtain a tidy table of values for a gene or a small set of genes
#'
#' This is the main function to access slot data data from a particular gene
#' (or a small set of genes) as a tidy table. If data for all genes
#' must be retrieved (as a large matrix), use the \code{\link{GetTable}}
#' function. For analysis results, use the \code{\link{GetAnalysisTable}} function.
#'
#' @param data A grandR object
#' @param mode.slot Which kind of data to access (see details)
#' @param columns A vector of columns (see details); all condition/cell names if NULL
#' @param genes Restrict the output table to the given genes (this typically is a single gene, or very few genes)
#' @param by.rows if TRUE, add rows if there are multiple genes / mode.slots; otherwise, additional columns are appended
#' @param coldata Should the table contain the \link{Coldata} values as well (at the beginning)?
#' @param ntr.na For columns representing a 4sU naive sample, should mode.slot \emph{ntr},\emph{new.count} and \emph{old.count} be 0,0 and count (ntr.na=FALSE; can be any other slot than count) or NA,NA and NA (ntr.na=TRUE)
#' @param name.by A column name of \link{Coldata}(data). This is used as the colnames of the output table
#'
#' @return A data frame containing the desired values
#'
#' @details To refer to data slots, the mode.slot syntax can be used: Each name is either a data slot, or one of (new,old,total) followed by a dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to obtain the \emph{new counts}.
#'
#' @details If only one mode.slot and one gene is given, the output table contains one column (and potentially columns from \link{Coldata}) named \emph{Value}. If one gene and multiple mode.slots are given, the columns are named according to the mode.slots. If one mode.slot and multiple genes are given, the columns are named according to the genes. If multiple genes and mode.slots are given, columns are named gene.mode.slot.
#'
#' @details If by.rows=TRUE, the table is molten such that each row contains only one value (for one of the genes and for one of the mode.slots). If only one gene and one mode.slot is given, melting does not have an effect.
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment havin the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @seealso \link{GetTable},\link{GetAnalysisTable},\link{DefaultSlot},\link{Genes}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#' GetData(sars,mode.slot="ntr",gene="MYC")
#' # one gene, one mode.slot
#' GetData(sars,mode.slot=c("count","ntr"),gene="MYC",coldata = FALSE)
#' # one gene, multiple mode.slots
#' GetData(sars,mode.slot=c("count","ntr"),gene=c("SRSF6","MYC"),by.rows=TRUE)
#' # multiple genes, multiple mode.slots, by rows
#'
#' @export
#'
#' @concept data
GetData=function(data,mode.slot=DefaultSlot(data),columns=NULL,genes=Genes(data),by.rows=FALSE,coldata=TRUE,ntr.na=TRUE,name.by="Symbol") {
  if (!all(check.mode.slot(data,mode.slot))) stop(sprintf("mode.slot %s unknown!",paste(mode.slot[!check.mode.slot(data,mode.slot)],collapse=",")))

  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)

  genes=ToIndex(data,genes)
  if (length(genes)==0) stop("Genes not found!")

  og=if (name.by %in% names(data$gene.info)) data$gene.info[[name.by]][genes] else data$gene.info[genes,1]
  uno=function(mode.slot) {
    tno="t"
    spl=strsplit(mode.slot,".",fixed=TRUE)[[1]]
    if (length(spl)>1) {tno=spl[1]; mode.slot=spl[2];}
    mf = switch(tolower(substr(tno,1,1)),t=1,n=as.matrix(data$data$ntr[genes,columns]),o=1-as.matrix(data$data$ntr[genes,columns]),stop(paste0(mode.slot," unknown!")))
    if (!ntr.na) {
      mf[is.na(mf)|is.nan(mf)]=if(tolower(substr(tno,1,1))=="n") 0 else 1
    }
    conv=if (mode.slot=="count") function(m) {mode(m) <- "integer";m} else if (mode.slot=="ntr" && !ntr.na) function(m) {m[is.na(m)]=0; m} else function(m) m

  if (!(mode.slot %in% names(data$data))) stop(paste0(mode.slot," unknown!"))
    if (length(genes)==1) data.frame(conv(as.matrix(data$data[[mode.slot]][genes,columns])*mf)) else as.data.frame(conv(t(as.matrix(data$data[[mode.slot]][genes,columns])*mf)))
  }
  re=as.data.frame(lapply(mode.slot,uno))
  if(length(mode.slot)==1 && length(genes)==1) names(re)="Value" else if (length(mode.slot)==1) names(re)=og else if (length(genes)==1) names(re)=mode.slot else names(re)=paste0(rep(og,length(mode.slot)),".",rep(mode.slot,each=length(og)))
  if (coldata) re = cbind(data$coldata[columns,],re)
  if (by.rows && (length(genes)>1 || length(mode.slot)>1)) {
    re = reshape2::melt(re,id.vars=if(coldata) names(data$coldata) else c(),value.name="Value")
    if (length(mode.slot)==1) names(re)[dim(re)[2]-1]="Gene" else if (length(genes)==1) names(re)[dim(re)[2]-1]="mode.slot" else {
      re=cbind(re[,c(1:(dim(re)[2]-2))],setNames(as.data.frame(t(as.data.frame(strsplit(as.character(re$variable),".",fixed = TRUE)))),c("Gene","mode.slot")),Value=re$Value)
      rownames(re)=NULL
    }
  }
  re
}


#' Copy the NTR slot and save under new name
#'
#' @param data the grandR object
#' @param name the name of the new slot
#'
#' @return a grandR object
#' @export
#'
#' @concept data
SaveNtrSlot=function(data,name) {
  AddSlot(data,name,data$data$ntr)
}

#' Copy the NTR slot and save under new name
#'
#' @param data the grandR object
#' @param name the name of the new slot
#'
#' @return a grandR object
#' @export
#'
#' @concept data
UseNtrSlot=function(data,name) {
  if (!check.slot(data,name,allow.ntr = FALSE)) stop("Illegal slot!")
  data$data$ntr = data$data[[name]]
  data
}


#' Significant genes
#'
#' Return significant genes for this grandR object
#'
#' @param data the grandR object
#' @param analysis the analysis to use, can be more than one and can be regexes (see details)
#' @param regex interpret analyses as regex?
#' @param criteria the criteria used to define what significant means; if NULL, Q<0.05 & abs(LFC)>=1 is used; can use the column names of the analysis table as variables,  should be a logical or numerical value per gene (see Details)
#' @param as.table return a table
#' @param use.symbols return them as symbols (gene ids otherwise)
#' @param gene.info add gene infos to the output table
#'
#' @details The analysis parameter (just like for \link{GetAnalysisTable} can be a regex (that will be matched
#' against all available analysis names). It can also be a vector (of regexes). Be careful with this, if
#' more than one table e.g. with column LFC ends up in here, only the first is used (if criteria=LFC).
#'
#' @details The criteria parameter can be used to define how analyses are performed. If criteria is a logical,
#' it obtains significant genes defined by cut-offs (e.g. on q value and LFC).
#' If it is a numerical, all genes are returned sorted (descendingly) by this value.
#' The columns of the given analysis table(s) can be used to build this expression.
#'
#' @return a vector of gene names (or symbols), or a table
#'
#' @export
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c(Design$Condition,Design$dur.4sU,Design$Replicate))
#' sars <- subset(sars,Coldata(sars,Design$dur.4sU)==2)
#' sars<-LFC(sars,mode="total",contrasts=GetContrasts(sars,contrast=c("Condition","Mock")))
#' GetSignificantGenes(sars,criteria=LFC>1)
#'
#' @concept diffexp
GetSignificantGenes=function(data,analysis=NULL,regex=TRUE,criteria=NULL,as.table=FALSE,use.symbols=TRUE,gene.info=TRUE) {
  analysis=match.analyses(data,analysis,regex)
  criteria=substitute(criteria)

    re=GeneInfo(data)
  rownames(re)=if(use.symbols) re$Symbol else re$Gene

  for (name in analysis) {
    tab=GetAnalysisTable(data,analyses=name,regex=FALSE,gene.info=FALSE,prefix.by.analysis=FALSE)
    if (is.null(criteria)) {
      use=tab$Q<0.05 & abs(tab$LFC)>=1
    } else {
      use=eval(criteria,envir=tab,enclos = parent.frame())
    }
    re[[name]]=use
    re[[name]][is.na(re[[name]])]=FALSE
  }
  if (!as.table) {
    re=re[,(ncol(GeneInfo(data))+1):ncol(re),drop=FALSE]
    cls=unique(sapply(re,class))
    if (length(cls)!=1) stop("Output table is mixed logical and numerical!")
    if (cls=="logical") {
      re=apply(re,1,any)
      re=names(re)[re]
    } else if (cls=="numeric") {
      if (ncol(re)>1) stop("Multiple numerical values, can only return as a table!")
      re=rownames(re[order(re[,1],decreasing=TRUE),,drop=FALSE])
    }
    return(re)
  }

  if (!gene.info) {
    re=re[,(ncol(GeneInfo(data))+1):ncol(re),drop=FALSE]
  }
  re=re[order(re[,ncol(re)],decreasing=TRUE),,drop=FALSE]
  re
}


#' Obtain reference columns (samples or cells) for all columns (samples or cells) in the data set
#'
#' In some situations (see examples) it is required to find a reference sample of some kind for each sample in a data set.
#' This is a convenience method to find such reference samples, and provide them as a lookup table.
#'
#' @param data A grandR object
#' @param reference Expression evaluating to a logical vector to indicate which columns are reference columns; evaluated in an environment having the columns of \link{Coldata}(data)
#' @param reference.function Function evaluating to a logical vector to indicate which columns are reference columns; called with the data frame row corresponding to the sample, and evaluated in an environment having the columns of \link{Coldata}(data)
#' @param group a vector of colnames in \link{Coldata}(data)
#' @param as.list return it as a list (names correspond to each sample, elements are the reference samples)
#' @param columns find references only for a subset of the columns (samples or cells; can be NULL)
#'
#' @return A logical matrix that contains for each sample or cell (in columns) a TRUE for the corresponding corresponding reference samples or cells in rows
#'
#' @details Without any group, the list simply contains all references for each sample/cell. With groups defined, each list entry consists of all references from the same group.
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment havin the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @seealso \link{Coldata},\link{Findno4sUPairs}, \link{Condition}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' FindReferences(sars,reference=no4sU)
#' # obtain the corresponding no4sU sample for each sample; use the Condition column
#' FindReferences(sars,Condition=="Mock",group="duration.4sU.original")
#' # obtain for each sample the corresponding sample in the Mock condition
#' FindReferences(sars,Condition=="Mock",group=c("duration.4sU.original","Replicate"))
#' # obtain for each sample the corresponding Mock sample, paying attention to replicates
#'
#' @export
#'
#' @concept snapshot
FindReferences=function(data,reference=NULL, reference.function=NULL,group=NULL, as.list=FALSE,columns=NULL) {
  if (!is.grandR(data)) stop("Data is not a grandR object!")
  if (!is.null(group) && !all(group %in% names(Coldata(data)))) stop(sprintf("No %s in Coldata!",group))

  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)


  df=Coldata(data)
  df=df[columns,]

  df$group=as.character(if(is.null(group)) "GROUP" else interaction(df[group],drop=FALSE,sep="."))
  ef=substitute(reference.function)
  if (!is.null(reference.function)) {
    setColnames=function(m,n) {colnames(m)=n; m}
    map=plyr::dlply(df,plyr::.(group),function(s) setColnames(sapply(1:nrow(s),function(row) setNames(eval(ef,s,parent.frame())(as.list(s[row,])),s$Name) ),s$Name))
    re=matrix(FALSE,nrow=nrow(df),ncol=nrow(df))
    colnames(re)=df$Name
    rownames(re)=df$Name
    for (mat in map) re[rownames(mat),colnames(mat)]=mat
    return(re)

  } else {
    e=substitute(reference)
    map=plyr::dlply(df,plyr::.(group),function(s) as.character(s$Name[eval(e,s,parent.frame())]))
    pairs=setNames(lapply(df$group,function(g) map[[g]]),df$Name)
  }
  if (as.list) return(pairs)
  mat = sapply(pairs,function(names) colnames(data) %in% names)
  rownames(mat)=colnames(data)
  mat
}


#' Analysis table functions
#'
#' Get analysis names and add or remove analyses
#'
#' @param data A grandR object
#' @param description if TRUE, also return the column names of each analysis table (i.e. a list named according to the analyses)
#' @param pattern A regular expression that is matched to analysis names
#' @param name The user-defined analysis name
#' @param table The analysis table to add
#' @param by Specify a column that contains gene names or symbols (see details)
#' @param warn.present Warn if an analysis with the same name is already present (and then overwrite)
#'
#' @return Either the analysis names or a grandR data with added/removed slots or the metatable to be used with AddAnalysis
#'
#' @details The columns in the analysis tables are defined by the analysis method (e.g. "Synthesis","Half-life" and "rmse" by \code{FitKinetics}).
#' A call to an analysis function might produce more than one table (e.g. because kinetic modeling is done for multiple \link{Condition}s). In this case,
#' AddAnalysisTable produces more than one analysis table.
#'
#' @details \code{AddAnalysis} is in most cases  not called directly by the user, but is
#' used by analysis methods to add their final result to a grandR object (e.g., \link{FitKinetics},\link{LikelihoodRatioTest},\link{LFC},\link{PairwiseDESeq2}).
#'
#' @details If it is called by the user (e.g. to add analysis results from external tools or from the literature, see pulse-chase vignette), then
#' the user must make sure that either the rownames of the given table can be recognized as genes (names or symbols), or that there is a column in the
#' table giving genes (this must be specified as the "by" parameter). The table does neither have to be sorted the same way the grandR object is, nor does
#' it have to be complete. \code{AddAnalysis} will take care or reordering and inserting NA for missing genes (and it will issue a warning in case of missing genes).
#'
#' @seealso \link{Slots}, \link{DefaultSlot}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#'
#' sars <- Normalize(sars)     # default behavior is to update the default slot; this calls AddSlot
#' Slots(sars)
#' DefaultSlot(sars)
#' sars <- DropSlot(sars,"norm")
#' sars                        # note that the default slot reverted to count
#'
#' @describeIn Analyses Obtain the analyses names
#' @export
#'
#' @concept grandr
Analyses=function(data, description=FALSE) {
  if (!description) {
    names(data$analysis)
  } else {
    setNames(lapply(names(data$analysis),function(n) colnames(data$analysis[[n]])),names(data$analysis))
  }
}

#' @describeIn Analyses Add an analysis table
#' @export
AddAnalysis=function(data,name,table,by = NULL, warn.present=TRUE) {
  if (!is.data.frame(table)) stop("Cannot add; analysis table must be a data frame!")
  #if (!equal(Genes(data,rownames(table)),Genes(data))) stop("Analysis table must contain row names corresponding to all genes!")

  if (!is.null(by)) {
    row.names(table) = table[,by]
    table <- table[, !names(table) %in% by, drop = FALSE]
  }

  if (!equal(Genes(data,rownames(table)),Genes(data))) {
      warning("Analysis table and grandR object does not have the same set of genes! Watch out for NA values!")
    table <- table[Genes(data), ]
    rownames(table) = Genes(data)
  }

  if (is.null(data$analysis)) data$analysis=list()
  if (is.null(data$analysis[[name]])) {
    data$analysis[[name]]=table
  } else {
    if (warn.present & any(names(table)%in%names(data$analysis[[name]]))) {
      ex=names(data$analysis[[name]])
      nc=names(table)
      add=paste(setdiff(nc,ex),collapse = ", ")
      over=paste(intersect(nc,ex),collapse = ", ")
      keep=paste(setdiff(ex,nc),collapse = ", ")
      if (add=="") add="<none>"
      if (over=="") over="<none>"
      if (keep=="") keep="<none>"
      warning(sprintf("Analysis %s already present! Adding: %s, Overwritting: %s, keeping: %s...",name,add,over,keep))
    }
    for (n in names(table)) data$analysis[[name]][[n]]=table[[n]]
    ana = attr(data$analysis[[name]],"analysis")
  }
  data
}


#' @describeIn Analyses Remove analyses from the grandR object
#' @export
DropAnalysis=function(data,pattern=NULL) {
  if (is.null(pattern)) {
    data$analysis=NULL
  } else {
    data$analysis=data$analysis[!grepl(pattern,names(data$analysis))]
  }
  invisible(data)
}

match.analyses=function(data,analyses=NULL,regex=TRUE) {
  if (!all(check.analysis(data,analyses,regex))) stop(sprintf("No analysis found for pattern %s!",paste(analyses[!check.analysis(data,analyses,regex)],collapse=",")))
  if (is.null(analyses)) {
    analyses=1:length(Analyses(data))
  } else if (regex) {
    analyses=unlist(lapply(analyses,function(pat) grep(pat,Analyses(data))))
  } else if (is.character(analyses)) {
    analyses=which(Analyses(data) %in% analyses)
  }
  Analyses(data)[analyses]
}

#' Obtain a table of analysis results values
#'
#' This is the main function to access analysis results. For slot data, use \code{\link{GetTable}} (as a large matrix)
#' or \code{\link{GetData}} (as tidy table).
#'
#' @param data A grandR object
#' @param analyses One or several regex to be matched against analysis names (\link{Analyses}); all analysis tables if NULL
#' @param regex Use regex for analyses (TRUE) or don't (FALSE, i.e. must specify the exact name)
#' @param columns Regular expressions to select columns from the analysis table (all have to match!); all columns if NULL
#' @param genes Restrict the output table to the given genes
#' @param by.rows if TRUE, add rows if there are multiple analyses; otherwise, additional columns are appended; TRUE also sets prefix.by.analysis to FALSE!
#' @param gene.info Should the table contain the \link{GeneInfo} values as well (at the beginning)?
#' @param name.by A column name of \link{Coldata}(data). This is used as the rownames of the output table
#' @param prefix.by.analysis Should the column names in the output prefixed by the analysis name?
#'
#' @return A data frame containing the analysis results
#'
#' @details The names for the output table are <Analysis name>.<columns name>
#'
#' @seealso \link{GetTable},\link{GetData},\link{Genes}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' sars<-LFC(sars,contrasts=GetContrasts(sars,group = "duration.4sU"))
#' head(GetAnalysisTable(sars,columns="LFC"))
#'
#' @export
#'
#' @concept data
GetAnalysisTable=function(data,analyses=NULL,regex=TRUE,columns=NULL,genes=Genes(data),by.rows=FALSE,gene.info=TRUE,name.by="Symbol",prefix.by.analysis=TRUE) {
  analyses=match.analyses(data,analyses,regex)
  genes=ToIndex(data,genes)

  if (by.rows) prefix.by.analysis=FALSE
  rbind.save=function(a,b) {
    if (length(unique(names(a)))!=length(names(a))) stop("Table names are not unique!")
    if (length(unique(names(b)))!=length(names(b))) stop("Table names are not unique!")
    for (n in setdiff(names(a),names(b))) b[[n]]=NA
    for (n in setdiff(names(b),names(a))) a[[n]]=NA
    b=b[,names(a)]
    rbind(a,b)
  }

  re=GeneInfo(data)[genes,]

  if (!is.null(name.by)) {
    rownames(re)=if (name.by %in% names(GeneInfo(data))) GeneInfo(data,column = name.by)[genes] else GeneInfo(data)[genes,1]
  }
  sintersect=function(a,b) if (is.null(b)) a else intersect(a,b)

  for (name in analyses) {
    t=data$analysis[[name]][genes,,drop=FALSE]
    if (!is.null(columns)) {
     use = rep(TRUE,ncol(t))
     for (r in columns) use = use&grepl(r,names(t))
     t=t[,use,drop=FALSE]
    }
    if (ncol(t)>0) {
      if (prefix.by.analysis) names(t)=paste0(name,".",names(t))
      if (by.rows) t=cbind(data.frame(Analysis=name),t)
      if (by.rows && name!=analyses[1]) {
        re=rbind.save(re,cbind(GeneInfo(data)[genes,],t))
      } else {
        re=cbind(re,t)
      }
    }
  }

  if (is.null(name.by)||by.rows) rownames(re)=NULL
  if (is.logical(gene.info) && !gene.info) re=re[,-(1:ncol(GeneInfo(data))),drop=FALSE]
  if (is.character(gene.info)) re=re[,-which(!names(GeneInfo(data)) %in% gene.info),drop=FALSE]
  if (by.rows) re$Analysis=factor(re$Analysis,levels=analyses)
  re
}




#' Stored plot functions
#'
#' Get plot names and add or remove plots
#'
#' @param data A grandR object
#' @param name The user-defined plot name
#' @param FUN The plotting function to add
#' @param pattern A regular expression that is matched to plot names
#' @param gene The gene to plot
#'
#' @return Either the plot names or a grandR data with added/removed plots
#'
#' @details FUN has to be a function with a single parameter for global plots (i.e., the grandR object) or two parameters for gene plots
#' (i.e., the grandR object and the gene name). Usually, it is either the name of a plotting function, such as \link{PlotGeneOldVsNew}, or, if it is
#' necessary to parametrize it, a call to \link{Defer} (which takes care of caching plots without storing an additional copy of the grandR object).
#'
#' @describeIn Plots Obtain the plot names
#' @export
#'
#' @concept grandr
Plots=function(data) {
  re=list()
  if (!is.null(data$plots$gene)) re=c(re,list(gene=names(data$plots$gene)))
  if (!is.null(data$plots$global)) re=c(re,list(gene=names(data$plots$global)))
  re
}
#' @describeIn Plots Add a gene plot to the grandR object
#' @export
AddGenePlot=function(data,name,FUN) {
  if (!is.function(FUN)) stop("Cannot add; FUN must be a function!")
  if (is.null(data$plots)) data$plots=list()
  if (is.null(data$plots$gene)) data$plots$gene=list()
  if (!is.null(data$plots$gene[[name]])) warning(sprintf("Plot %s already present! Overwriting...",name))
  data$plots$gene[[name]]=FUN
  data
}

#' @describeIn Plots Add a global plot to the the grandR object
#' @export
AddGlobalPlot=function(data,name,FUN) {
  if (!is.function(FUN)) stop("Cannot add; FUN must be a function!")
  if (is.null(data$plots)) data$plots=list()
  if (is.null(data$plots$global)) data$plots$global=list()
  if (!is.null(data$plots$global[[name]])) warning(sprintf("Plot %s already present! Overwriting...",name))
  data$plots$global[[name]]=FUN
  data
}

#' @describeIn Plots Create a gene plot
#' @export
PlotGene=function(data,name,gene) data$plots$gene[[name]](data,gene)
#' @describeIn Plots Create a global plot
#' @export
PlotGlobal=function(data,name) data$plots$global[[name]](data)

#' @describeIn Plots Remove plots from the grandR object
#' @export
DropPlots=function(data,pattern=NULL) {
  if (is.null(pattern)) {
    data$plots=NULL
  } else {
    if (!is.null(data$plots$gene)) data$plots$gene=data$plots$gene[!grepl(pattern,names(data$plots$gene))]
    if (!is.null(data$plots$global)) data$plots$global=data$plots$global[!grepl(pattern,names(data$plots$global))]
  }
  data
}
