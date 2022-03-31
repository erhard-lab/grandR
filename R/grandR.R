
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
#' Usually, this contructor is not invoked directly (but by \code{\link{ReadGRAND}} or \code{\link{SimulateTimeCourse}}).
#'
#' @param prefix Can either be the prefix used to call GRAND-SLAM with, or the main output file ($prefix.tsv.gz);
#' if the RCurl package is installed, this can also be a URL
#' @param gene.info a data frame with metadata for all genes
#' @param slots A list of matrices representing the slots
#' @param coldata a data frame with metedata for all samples (or cells)
#' @param metadata a metadata list
#' @param parent A parent object containing default values for all other parameters (i.e. all parameters not specified are obtained from this object)
#'
#' @param data a grandR object
#' @param order can either be an integer or character vector representing a permutation of the columns (samples or cells)
#' @param columns can either be a logical, integer or character vector representing a selection of the columns (samples or cells)
#' @param column.name The name of the annotation table according to which the object is split or the new annotation table column name denoting the origin after merging
#' @param map named list or vector representing a lookup table (names are current column names)
#' @param fun a function that maps a vector of names to a new vector of names
#' @param s1,s2 column names
#' @param list a list of grandR objects
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
#' @section Functions:
#' \describe{
#'   \item{Title}{Obtain a useful title for the project (from the prefix parameter)}
#'   \item{dim}{Obtain the dimensions (genes x samples or genes x cells)}
#'   \item{is}{Check whether it is a grandR object}
#'   \item{dimnames}{Obtain the row and column names of this object (genes x samples or genes x cells)}
#'   \item{print}{Print information on this grandR object}
#'   \item{reorder}{Create a new grandR object with reordered columns}
#'   \item{subset}{Create a new grandR object with a subset of the columns (use \code{\link{FilterGenes}} to subset on genes)}
#'   \item{split}{Split the grandR object into a list of multiple grandR objects (according to the levels of an annotation table column)}
#'   \item{RenameColumns}{Rename the column names according to a lookup table (map) or a function (invoked on the current names)}
#'   \item{SwapColumns}{Swap the order of two columns (samples or cells)}
#'   \item{merge}{Merge several grandR objects into one}
#' }
#'
#'
#' @seealso \link{Slots}, \link{DefaultSlot}, \link{Genes}, \link{GeneInfo}, \link{Coldata}, \link{GetTable}, \link{GetData}, \link{Analyses}, \link{GetAnalysisTable}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#'
#' dim(sars)               # this is the corona data from Finkel et al., Nature 2021, filtered for genes with >=100 reads in the SARS (3hpi) sample
#' head(rownames(sars))
#'
#' @export
#'
grandR=function(prefix=parent$prefix,gene.info=parent$gene.info,slots=parent$data,coldata=parent$coldata,metadata=parent$metadata,parent=NULL) {
  info=list()
  info$prefix=prefix
  info$gene.info=gene.info
  info$data=slots
  info$coldata=coldata
  info$metadata=metadata
  info$analysis=NULL
  class(info)="grandR"
  info
}


#' @rdname grandR
#' @export
VersionString=function(data) {
  "grandR v0.1.0"
}


#' @rdname grandR
#' @export
Title=function(data) {
  x=strsplit(data$prefix,"/")[[1]]
  x[length(x)]
}

#' @rdname grandR
#' @export
dim.grandR=function(data) c(dim(data$gene.info)[1],dim(data$coldata)[1])
#' @rdname grandR
#' @export
is.grandR <- function(x) inherits(x, "grandR")
#' @rdname grandR
#' @export
dimnames.grandR=function(data) dimnames(data$data$count)
#' @rdname grandR
#' @export
print.grandR=function(data) cat(
  sprintf("grandR: %s\nRead from %s\n%d genes, %d samples/cells\nAvailable data slots: %s\nAvailable analyses: %s\nDefault data slot: %s\n",
          data$metadata$Description,
          data$prefix,
          nrow(data),
          ncol(data),
          paste(Slots(data),collapse=","),
          paste(Analyses(data),collapse=","),
          DefaultSlot(data))
)
#' Internal function to apply functions to all slots etc.
#'
#' @param data a grandR object
#' @param fun apply this function to each data slot (i.e. it receives each data matrix)
#' @param fun.gene.info apply this function to the gene info table
#' @param fun.coldata apply this function to the column annotation table
#'
#' @details The additional parameters are provided to each of the functions.
#' @return A new grandR object
data.apply=function(data,fun,fun.gene.info=NULL,fun.coldata=NULL,...) {
  re=list()
  for (l1 in names(data$data)) {
    re[[l1]]=fun(data$data[[l1]],...)
  }
  ngene.info=if (!is.null(fun.gene.info)) fun.gene.info(data$gene.info,...) else data$gene.info
  ncoldata=if (!is.null(fun.coldata)) fun.coldata(data$coldata,...) else data$coldata
  grandR(parent=data,gene.info=ngene.info,slots=re,coldata=ncoldata)
}

#' @rdname grandR
#' @export
reorder.grandR=function(data,order) {
  # REVIEWED BY Lygeri: changed "columns" to "order"
  r=subset.grandR(data,order)
  r$coldata=type.convert(r$coldata)
  r
}
#' @rdname grandR
#' @export
subset.grandR=function(data,columns) {
  columns = substitute(columns)
  columns = eval(columns, Coldata(data), parent.frame())
  keep=rownames(data$coldata)[columns]
  data.apply(data,function(m) m[,intersect(keep,colnames(m))],fun.coldata = function(t)
    droplevels(t[columns,]))
}

#' @rdname grandR
#' @export
split.grandR=function(data,column.name=Design$Condition) {
  col=as.factor(data$coldata[[column.name]])
  setNames(lapply(levels(col),function(c) {re=subset(data,col==c);
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
  data$coldata$Name=names
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
#' @export
merge.grandR=function(...,list=NULL,column.name=Design$Origin) {
  nn=c(get.varargs.names(...),names(list))
  list=c(list(...),list)
  names(list)=nn
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

    # add potential additional gene annotations
    gi=GeneInfo(add)
    for (n in names(gi)) {
      if (!n %in% names(GeneInfo(re))) {
        GeneInfo(re,column = n)=gi[[n]]
      } else if (!all(GeneInfo(re)[[n]]==gi[[n]])) {
        GeneInfo(re,column = paste0(n,".",names(list)[i]))=gi[[n]]
      }
    }
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
DefaultSlot <- function(data,value=NULL) {
  if (is.null(value)) data$metadata$default.slot else {
    DefaultSlot(data)=value
    data
  }
}

#' @rdname DefaultSlot
#' @export
`DefaultSlot<-` <- function(data, value) {
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
#' @export
AddSlot=function(data,name,matrix,set.to.default=FALSE) {
  if (!all(colnames(matrix)==colnames(data$data$count))) stop("Column names do not match!")
  if (!all(rownames(matrix)==rownames(data$data$count))) stop("Row names do not match!")
  if (!is.matrix(matrix)) stop("Must be a matrix!")
  if (grepl(".",name,fixed=TRUE)) stop("Name may not contain a dot!")
  data$data[[name]]=matrix
  if (set.to.default) DefaultSlot(data)=name
  data
}

#' Get or set the conditions in the column annotation table.
#'
#' The conditions column from the column annotation table is used by several functions to stratify the columns (samples or cells)
#' during the analysis (e.g. to estimate separate kinetic parameters with \code{\link{FitKinetics}} or it is used as covariate for
#' \code{\link{LFC}} or \code{\link{TestGenesLRT}}). For that reason there are special functions to set and get this column.
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
  } else if (all(value %in% names(data$coldata))) {
    data$coldata$Condition <- interaction(data$coldata[value])
  } else{
    data$coldata$Condition <- as.factor(value)
  }
  data
}





#' Get the genes and sample (or cell) names for a grandR object, or add an additional gene annotation column
#'
#' The genes are either the (often unreadable) gene ids (e.g. Ensembl ids), or the symbols. The Columns function
#' can also return the columns from a given analysis table instead of the sample (or cell) names.
#'
#' @param data A grandR object
#' @param use.symbols obtain the gene symbols instead of gene names
#' @param genes which genes to use
#' @param regex treat genes as a regex, and return all that match
#' @param columns which columns to return (as numeric or logical vector)
#' @param analysis The name of an analysis table
#'
#' @details \code{Genes(data,use.symbols=FALSE)} it the same as \code{rownames(data)}, and \code{Columns(data)} is the same as \code{colnames(data)}
#'
#' @details If both column and value are specified for \code{GeneInfo}, a new column is added to the gene annotation table
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
Genes=function(data, genes=NULL, use.symbols=TRUE,regex=FALSE) data$gene.info[[if (use.symbols) "Symbol" else "Gene"]][ToIndex(data,genes,regex=regex)]
#' @rdname Genes
#' @export
Columns=function(data,columns=NULL,analysis=NULL) {
  if (is.null(analysis)) {
    if (is.null(columns)) rownames(Coldata(data)) else unname(setNames(rownames(Coldata(data)),rownames(Coldata(data)))[columns])
  } else names(data$analysis[[analysis]])
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
#' GeneInfo(sars,"LengthCategory")<-cut(GeneInfo(sars)$Length,c(0,1500,2500,Inf),labels=c("Short","Medium","Long"))
#' table(GeneInfo(sars)$LengthCategory)
#'
#' @export
#'
GeneInfo=function(data,column=NULL,value=NULL) {
  if (is.null(column)) data$gene.info else {
    data$gene.info[[colum]]=value
    data
  }
}
#' @rdname GeneInfo
#' @export
`GeneInfo<-` <- function(data, column, value) {
  data$gene.info[[column]]=value
  data
}


#' Get the column annotation table or add additional columns to it
#'
#' The colums of a grandR object are samples or cells.
#' The column annotation table contains meta information for the columns of a grandR object.
#' When loaded from the GRAND-SLAM output, this this constructed from the sample/cell names by
#' \code{\link{MakeColdata}}
#'
#' @param data A grandR object
#' @param column The name of the additional annotation column
#' @param value The additional annotation per sample or cell
#'
#' @details New columns can be added either by \code{data<-Coldata(data,name,values)} or by \code{Coldata(data,name)<-values}.
#'
#' @details The column named \emph{Condition} has a special meaning in this table: It is used by several functions to stratify the columns
#' during the analysis (e.g. to estimate separate kinetic parameters with \code{\link{FitKinetics}} or it is used as covariate for
#' \code{\link{LFC}} or \code{\link{TestGenesLRT}}). For that reason there are special functions to set and get this column.
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
#' GeneInfo(sars,"LengthCategory")<-cut(GeneInfo(sars)$Length,c(0,1500,2500,Inf),labels=c("Short","Medium","Long"))
#' table(GeneInfo(sars)$LengthCategory)
#'
#' @export
#'
Coldata=function(data,column=NULL,value=NULL) {
  if (is.null(column)) data$coldata else {
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
#' @param slot a slot name
#' @param mode.slot a mode.slot
#'
#' @details A mode.slot is a mode followed by a dot followed by a slot name, or just a slot name. A mode is either \emph{total}, \emph{new} or \emph{old}.
#'
#' @return Whether or not the given name is valid and unique for the grandR object
#'
check.analysis=function(data,analyses,regex) sapply(analyses,function(pattern) any(grepl(pattern,Analyses(data),fixed=!regex)))
#' @rdname check.analysis
check.slot=function(data,slot) slot %in% names(data$data)
#' @rdname check.analysis
check.mode.slot=function(data,mode.slot) {
  sapply(strsplit(mode.slot,".",fixed=TRUE),function(spl) {
    if (length(spl)>2 || length(spl)==0) return(FALSE)
    if (length(spl)==1) check.slot(data,spl) else tolower(substr(spl[1],1,1)) %in% c("t","n","o") && check.slot(data,spl[2])
  })
}

#' Internal functions to parse mode.slot strings
#'
#' @param data a grandR object
#' @param mode.slot a mode.slot
#'
#' @details A mode.slot is a mode followed by a dot followed by a slot name, or just a slot name. A mode is either \emph{total}, \emph{new} or \emph{old}
#'
#' @return a named list with elements mode and slot (or only slot in case of \emph{ntr},\emph{alpha} or \emph{beta})
#'
get.mode.slot=function(data,mode.slot) {
  if (length(mode.slot)!=1) stop("mode.slot must be a vector of length 1")
  if (!check.mode.slot(data,mode.slot)) stop("Invalid mode.slot")
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
ToIndex=function(data,gene,regex=FALSE) {
  if (regex) gene=grepl(gene,data$gene.info$Gene)|grepl(gene,data$gene.info$Symbol)
  if (is.null(gene)) return(1:nrow(data))
  if (is.numeric(gene)) return(gene)
  if (is.logical(gene) && length(gene)==nrow(data)) return(which(gene))
  if (all(gene %in% data$gene.info$Gene)) return(setNames(1:nrow(data),data$gene.info$Gene)[gene])
  if (all(gene %in% data$gene.info$Symbol)) return(setNames(1:nrow(data),data$gene.info$Symbol)[gene])
  stop("Could not find all genes!")
}

#' Obtain a genes x values table
#'
#' This is the main function to access slot data for all genes as a large matrix. If data from a particular gene (or a small set of genes)
#' must be retrieved, use the \code{\link{GetData}} function. For analysis results, use the \code{\link{GetAnalysisTable}} function.
#'
#' @param data A grandR object
#' @param type Either a mode.slot (see details) or a regex to be matched against analysis names. Can also be a vector; If NULL, \link{DefaultSlot}(data) is used
#' @param columns A vector of columns (either condition/cell names if the type is a mode.slot, or names in the output table from an analysis; use \link{Columns}(data,<analysis>) to learn which columns are available); all condition/cell names if NULL
#' @param genes Restrict the output table to the given genes
#' @param ntr.na For columns representing a 4sU naive sample, should types \emph{ntr},\emph{new.count} and \emph{old.count} be 0,0 and count (ntr.na=FALSE; can be any other slot than count) or NA,NA and NA (ntr.na=TRUE)
#' @param gene.info Should the table contain the \link{GeneInfo} values as well (at the beginning)?
#' @param summarize Should replicates by summarized? Can only be specified if columns is NULL; either a summarization matrix (\link{GetSummarizeMatrix}) or TRUE (in which case \link{GetSummarizeMatrix}(data) is called)
#' @param prefix Prepend each column in the output table (except for the gene.info columns) by the given prefix
#' @param name.by A column name of \link{Coldata}(data). This is used as the rownames of the output table
#'
#' @return A data frame containing the desired values
#'
#' @details This is a convenience wrapper for \link{GetData} (values from data slots) and \link{GetAnalysisTable} (values from analyses). Types can refer to any of the two (and can be mixed). If there are types from both data and analyses, columns must be NULL. Otherwise columns must either be condition/cell names (if type refers to one or several data slots), or regular expressions to match against the names in the analysis tables.
#' @details To refer to data slots, the mode.slot syntax can be used: Each name is either a data slot, or one of (new,old,total) followed by a dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to obtain the \emph{new counts}.
#'
#' @seealso \link{GetData},\link{GetAnalysisTable},\link{DefaultSlot},\link{Genes},\link{GetSummarizeMatrix}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#' sars <- Normalize(FilterGenes(sars))
#'
#' head(GetTable(sars)) # DefaultSlot values, i.e. size factor normalized read counts for all samples
#' head(GetTable(sars,summarize=TRUE)) # DefaultSlot values averaged over the two conditions
#' head(GetTable(sars,type="new.count",columns=!Coldata(sars)$no4sU)) # Estimated counts for new RNA for all samples with 4sU
#'
#' sars<-FitKinetics(sars,name = "kinetics",steady.state=list(Mock=TRUE,SARS=FALSE))
#' head(GetTable(sars,type="kinetics",columns="Half-life")) # Estimated RNA half-lives for both conditions
#'
#'
#'
#' @export
#'
GetTable=function(data,type=NULL,columns=NULL,genes=Genes(data),ntr.na=TRUE,gene.info=FALSE,summarize=NULL,prefix=NULL,name.by="Symbol") {
  if (is.null(type)) type=DefaultSlot(data)

  genes=ToIndex(data,genes)

  mode.slot=check.mode.slot(data,type)
  analysis=check.analysis(data,type,TRUE) & !mode.slot
  if (!all(analysis|mode.slot)) stop(sprintf("Type %s is neither a mode.slot nor an analysis name!",paste(type[!analysis&!mode.slot],collapse=",")))

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
      rtt=as.data.frame(t(GetData(data,tt,columns=cols,genes,ntr.na = ntr.na,coldata=FALSE, melt=FALSE, name.by = name.by)))
      names(rtt)=cols
      if (!is.null(summarize)) {
        rtt=as.data.frame(as.matrix(rtt) %*% summarize)
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

  # add necessary stuff
  if (!is.null(prefix)) colnames(r)=paste0(prefix,colnames(r))
  if (is.logical(gene.info) && gene.info) r=cbind(GeneInfo(data)[ToIndex(data,genes),],r)
  if (is.character(gene.info)) r=cbind(GeneInfo(data)[ToIndex(data,genes),-which(!names(GeneInfo(data)) %in% gene.info)],r)

  rownames(r)=data$gene.info[[name.by]][ToIndex(data,genes)]
  r
}

#' Obtain a genes x values table as a sparse matrix
#'
#' This is the main function to access slot data for all genes as a sparse matrix.
#'
#' @param data A grandR object
#' @param mode.slot Which kind of data to access (see details)
#' @param columns A vector of columns (either condition/cell names if the type is a mode.slot, or names in the output table from an analysis; use \link{Columns}(data,<analysis>) to learn which columns are available); all condition/cell names if NULL
#' @param genes Restrict the output table to the given genes
#' @param name.by A column name of \link{Coldata}(data). This is used as the rownames of the output table
#'
#' @return A sparse matrix containing the desired values
#'
#' @details To refer to data slots, the mode.slot syntax can be used: It is either a data slot, or one of (new,old,total) followed by a dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to obtain the \emph{new counts}.
#'
#' @seealso \link{GetData},\link{GetAnalysisTable},\link{DefaultSlot},\link{Genes},\link{GetSummarizeMatrix}
#'
#' @export
#'
GetSparseMatrix=function(data,mode.slot=DefaultSlot(data),columns=NULL,genes=Genes(data),name.by="Symbol") {

  if (!all(check.mode.slot(data,mode.slot))) stop(sprintf("mode.slot %s unknown!",paste(mode.slot[!check.mode.slot(data,mode.slot)],collapse=",")))
  if (length(mode.slot)!=1) stop("Specify exactly one mode.slot!")

  if (is.null(columns)) columns=colnames(data)
  genes=ToIndex(data,genes)

  tno="t"
  spl=strsplit(mode.slot,".",fixed=TRUE)[[1]]
  if (length(spl)>1) {tno=spl[1]; mode.slot=spl[2];}

  re=data$data[[mode.slot]][genes,columns,drop=FALSE]
  rownames(re)=Genes(data,genes)

  conv=function(v) { mode(v)="integer"; v}

  if (is.matrix(re)) {
    mf = switch(tolower(substr(tno,1,1)),t=1,n=as.matrix(data$data$ntr[genes,columns,drop=FALSE]),o=1-as.matrix(data$data$ntr[genes,columns,drop=FALSE]),stop(paste0(mode.slot," unknown!")))
    if (!ntr.na) {
      mf[is.na(mf)]=if(tolower(substr(tno,1,1))=="n") 0 else 1
    }
    re=re*mf
    if (mode.slot=="count") {
      mode(re) <- "integer"
    } else if (mode.slot=="ntr") {
      re[is.na(re)]=0
    }
    return(as(re,Class = Matrix::"dgCMatrix"))
  } else {
    if (tolower(substr(tno,1,1))=="t") return(re)
    if (tolower(substr(tno,1,1))=="n") {
      # all that have zero in ntr matrix will be zero, so this is fine
      sX <- Matrix::summary(re)
      sY <- Matrix::summary(data$data$ntr[genes,columns,drop=FALSE])
      sRes <- merge(sX, sY, by=c("i", "j"))
      print(dimnames(re))
      return(Matrix::sparseMatrix(i=sRes[,1], j=sRes[,2], x=conv(sRes[,3]*sRes[,4]),dims=dim(re),
                                  dimnames=dimnames(re)))
    }
    if (tolower(substr(tno,1,1))=="o") {
      sX <- Matrix::summary(re)
      sY <- Matrix::summary(data$data$ntr[genes,columns,drop=FALSE])
      sRes <- merge(sX, sY, by=c("i", "j"),all.x=TRUE)
      sRes[is.na(sRes[,4]),4]=0
      sRes[,4]=1-sRes[,4]
      sRes=sRes[sRes[,4]>0,]
      return(Matrix::sparseMatrix(i=sRes[,1], j=sRes[,2], x=conv(sRes[,3]*sRes[,4]),dims=dim(re),
                                  dimnames=dimnames(re)))
    }
    stop(paste0(mode.slot," unknown!"))
  }


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
#' @param columns A vector of columns (i.e. condition/cell names; use colnames(data) to learn which columns are available); all condition/cell names if NULL
#' @param genes Restrict the output table to the given genes (this typically is a single gene, or very few genes)
#' @param melt Should the table be melted if multiple genes / mode.slots are given
#' @param coldata Should the table contain the \link{Coldata} values as well (at the beginning)?
#' @param ntr.na For columns representing a 4sU naive sample, should mode.slot \emph{ntr},\emph{new.count} and \emph{old.count} be 0,0 and count (ntr.na=FALSE; can be any other slot than count) or NA,NA and NA (ntr.na=TRUE)
#' @param name.by A column name of \link{Coldata}(data). This is used as the colnames of the output table
#'
#' @return A data frame containing the desired values
#'
#' @details To refer to data slots, the mode.slot syntax can be used: Each name is either a data slot, or one of (new,old,total) followed by a dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to obtain the \emph{new counts}.
#' @details If only one mode.slot and one gene is given, the output table contains one column (and potentially columns from \link{Coldata}) named \emph{Value}. If one gene and multiple mode.slots are given, the columns are named according to the mode.slots. If one mode.slot and multiple genes are given, the columns are named according to the genes. If multiple genes and mode.slots are given, columns are named gene.mode.slot.
#' @details If melt=TRUE, the table is molten such that each row contains only one value (for one of the genes and for one of the mode.slots). If only one gene and one mode.slot is given, melting does not have an effect.
#'
#' @seealso \link{GetTable},\link{GetAnalysisTable},\link{DefaultSlot},\link{Genes}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#' GetData(sars,mode.slot="ntr",gene="MYC") # one gene, one mode.slot
#' GetData(sars,mode.slot=c("count","ntr"),gene="MYC",coldata = F) # one gene, multiple mode.slots
#' GetData(sars,mode.slot=c("count","ntr"),gene=c("SRSF6","MYC"),melt=TRUE) # multiple genes, multiple mode.slots, molten
#'
#' @export
#'
GetData=function(data,mode.slot=DefaultSlot(data),columns=NULL,genes=Genes(data),melt=FALSE,coldata=TRUE,ntr.na=TRUE,name.by="Symbol") {
  if (!all(check.mode.slot(data,mode.slot))) stop(sprintf("mode.slot %s unknown!",paste(mode.slot[!check.mode.slot(data,mode.slot)],collapse=",")))


  if (is.null(columns)) columns=colnames(data)
  genes=ToIndex(data,genes)
  og=if (name.by %in% names(data$gene.info)) data$gene.info[[name.by]][genes] else data$gene.info[genes,1]
  uno=function(mode.slot) {
    tno="t"
    spl=strsplit(mode.slot,".",fixed=TRUE)[[1]]
    if (length(spl)>1) {tno=spl[1]; mode.slot=spl[2];}
    mf = switch(tolower(substr(tno,1,1)),t=1,n=as.matrix(data$data$ntr[genes,columns]),o=1-as.matrix(data$data$ntr[genes,columns]),stop(paste0(mode.slot," unknown!")))
    if (!ntr.na) {
      mf[is.na(mf)]=if(tolower(substr(tno,1,1))=="n") 0 else 1
    }
    conv=if (mode.slot=="count") function(m) {mode(m) <- "integer";m} else if (mode.slot=="ntr" && !ntr.na) function(m) {m[is.na(m)]=0; m} else function(m) m

  if (!(mode.slot %in% names(data$data))) stop(paste0(mode.slot," unknown!"))
    if (length(genes)==1) data.frame(conv(as.matrix(data$data[[mode.slot]][genes,columns])*mf)) else as.data.frame(conv(t(as.matrix(data$data[[mode.slot]][genes,columns])*mf)))
  }
  re=as.data.frame(lapply(mode.slot,uno))
  if(length(mode.slot)==1 && length(genes)==1) names(re)="Value" else if (length(mode.slot)==1) names(re)=og else if (length(genes)==1) names(re)=mode.slot else names(re)=paste0(rep(og,length(mode.slot)),".",rep(mode.slot,each=length(og)))
  if (coldata) re = cbind(data$coldata[columns,],re)
  if (melt && (length(genes)>1 || length(mode.slot)>1)) {
    re = melt(re,id.vars=if(coldata) names(data$coldata) else c(),value.name="Value")
    if (length(mode.slot)==1) names(re)[dim(re)[2]-1]="Gene" else if (length(genes)==1) names(re)[dim(re)[2]-1]="Type" else {
      re=cbind(re[,c(1:(dim(re)[2]-2))],setNames(as.data.frame(t(as.data.frame(strsplit(as.character(re$variable)," ")))),c("Gene","Type")),Value=re$Value)
    }
  }
  re
}


#' Obtain reference columns (samples or cells) for all columns (samples or cells) in the data set
#'
#' In some situations (see examples) it is required to find a reference sample of some kind for each sample in a data set.
#' This is a convenience method to find such reference samples, and provide them as a lookup table.
#'
#' @param data A grandR object
#' @param reference Expression evaluating to a logical vector to indicate which columns are reference columns; evaluated in an environment having the columns of \link{Coldata}(data)
#' @param group a vector of colnames in \link{Coldata}(data)
#' @param as.list return it as a list (names correspond to each sample, elements are the reference samples)
#' @param columns find references only for a subset of the columns (samples or cells; can be NULL)
#'
#' @return A 0-1 matrix that contains for each sample or cell (in columns) a 1 for the corresponding corresponding reference samples or cells in rows
#'
#' @details Without any group, the list simply contains all references for each sample/cell. With groups defined, each list entry consists of all references from the same group.
#'
#' @seealso \link{Coldata},\link{Findno4sUPairs}, \link{Condition}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' FindReferences(sars,reference=no4sU) # obtain the corresponding no4sU sample for each sample; use the Condition column
#' FindReferences(sars,Condition=="Mock",group="duration.4sU.original") # obtain for each sample the corresponding sample in the Mock condition
#' FindReferences(sars,Condition=="Mock",group=c("duration.4sU.original","Replicate")) # obtain for each sample the corresponding sample in the Mock condition, paying attention to replicates
#'
#' @export
#'
FindReferences=function(data,reference, group="Condition", as.list=FALSE,columns=NULL) {
  if (!is.grandR(data)) stop("Data is not a grandR object!")
  if (!is.null(group) && !group %in% names(Coldata(data))) stop(sprintf("No %s in Coldata!",group))

  df=Coldata(data)
  if (!is.null(columns)) df=df[columns,]
  df$group=as.character(if(is.null(group)) 1 else interaction(df[group],drop=FALSE,sep="."))
  e=substitute(reference)
  map=dlply(df,.(group),function(s) as.character(s$Name[eval(e,s,parent.frame())]))
  pairs=setNames(lapply(df$group,function(g) map[[g]]),df$Name)
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
#' @param pattern A regular expression that is matched to analysis names
#' @param description A metatable created by MakeAnalysis
#' @param table The analysis table to add
#' @param warn.present Warn if an analysis with the same name is already present (and then overwrite)
#' @param name The user-defined analysis name
#' @param analysis The name of the analysis tool
#' @param mode An optional mode (new,old,total) on which the analysis has been run
#' @param slot An optional data slot on which the analysis has been run
#' @param columns An optional vector of columns the analysis was run on
#'
#' @return Either the analysis names or a grandR data with added/removed slots or the metatable to be used with AddAnalysis
#'
#' @details The columns in the analysis tables are defined by the analysis method (e.g. "Synthesis","Half-life" and "rmse" by \code{FitKinetics}).
#' A call to an analysis function might produce more than one table (e.g. because kinetic modeling is done for multiple \link{Condition}s). In this case,
#' AddAnalysisTable produces more than one analysis table.
#'
#' @details \code{AddAnalysis} (and therefore also \code{MakeAnalysis}) is usually not called directly by the user, but is
#' used by analysis methods to add their final result to a grandR object (e.g., \link{FitKinetics},\link{TestGenesLRT},\link{TestPairwise},\link{LFC}).
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
#' sars                        # note that the defauls slot reverted to count
#'
#' @describeIn Analyses Obtain the analyses names
#' @export
#'
Analyses=function(data) names(data$analysis)

#' @describeIn Analyses Add an analysis table
#' @export
AddAnalysis=function(data,description,table,warn.present=TRUE) {
  if (!is.data.frame(table)) stop("Cannot add; analysis table must be a data frame!")
  stopifnot(!is.null(description$name))
  description$results=names(table)
  if (is.null(data$analysis)) data$analysis=list()
  if (is.null(data$analysis[[description$name]])) {
    data$analysis[[description$name]]=table
    attr(data$analysis[[description$name]],"analysis")=list(description)
  } else {
    if (warn.present & any(names(table)%in%names(data$analysis[[description$name]]))) warning(sprintf("Analysis %s already present! Overwritting...",description$name))
    for (n in names(table)) data$analysis[[description$name]][[n]]=table[[n]]
    ana = attr(data$analysis[[description$name]],"analysis")
    attr(data$analysis[[description$name]],"analysis") = list(ana,description)
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
#' @describeIn Analyses Create a metatable for an analysis
#' @export
MakeAnalysis=function(name,analysis,mode=NULL,slot=NULL,columns=NULL) {
  list(name=name,mode=mode,analysis=analysis,slot=slot,columns=columns)
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
#' @param gene.info Should the table contain the \link{GeneInfo} values as well (at the beginning)?
#' @param name.by A column name of \link{Coldata}(data). This is used as the rownames of the output table
#'
#' @return A data frame containing the analysis results
#'
#' @details The names for the output table are <Analysis name>.<columns name>
#'
#' @seealso \link{GetTable},\link{GetData},\link{Genes}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#' SetParallel()
#' sars<-FitKinetics(sars,name = "kinetics",steady.state=list(Mock=TRUE,SARS=FALSE))
#' head(GetAnalysisTable(sars,columns="Half-life"))
#'
#' @export
#'
GetAnalysisTable=function(data,analyses=NULL,regex=TRUE,columns=NULL,genes=Genes(data),gene.info=TRUE,name.by="Symbol") {
  if (!all(check.analysis(data,analyses,regex))) stop(sprintf("No analysis found for pattern %s!",paste(analyses[!check.analysis(data,analyses,regex)],collapse=",")))

  genes=ToIndex(data,genes)


  re=data$gene.info[genes,]

  if (!is.null(name.by)) {
    rownames(re)=if (name.by %in% names(data$gene.info)) data$gene.info[[name.by]][genes] else data$gene.info[genes,1]
  }
  sintersect=function(a,b) if (is.null(b)) a else intersect(a,b)

  analyses=if (is.null(analyses)) 1:length(Analyses(data)) else unlist(lapply(analyses,function(pat) grep(pat,Analyses(data),fixed=!regex)))
  for (name in Analyses(data)[analyses]) {
    t=data$analysis[[name]][genes,,drop=FALSE]
    if (!is.null(columns)) {
     use = rep(TRUE,ncol(t))
     for (r in columns) use = use&grepl(r,names(t))
     t=t[,use,drop=FALSE]
    }
    if (ncol(t)>0) {
      names(t)=paste0(name,".",names(t))
      re=cbind(re,t)
    }
  }

  if (is.logical(gene.info) && !gene.info) re=re[,(ncol(data$gene.info)+1):ncol(re),drop=FALSE]
  if (is.character(gene.info)) re=re[,-which(!names(data$gene.info) %in% gene.info),drop=FALSE]
  re
}
