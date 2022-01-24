
grandR=function(prefix=parent$prefix,gene.info=parent$gene.info,data=parent$data,coldata=parent$coldata,metadata=parent$metadata,parent=NULL) {
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
DefaultSlot <- function(data,value=NULL) {
  if (is.null(value)) data$metadata$default.slot else {
    data$metadata$default.slot=value
    data
  }
}

Slots=function(data) names(data$data)

DropSlot=function(data,pattern=NULL) {
  if (is.null(pattern)) {
    stop("Cannot drop all slots!")
  } else {
    data$data=data$data[!grepl(pattern,names(data$data))]
  }
  invisible(data)
}
AddSlot=function(data,name,matrix) {
  if (!all(colnames(matrix)==colnames(data$data$count))) stop("Column names do not match!")
  if (!all(rownames(matrix)==rownames(data$data$count))) stop("Row names do not match!")
  if (!is.matrix(matrix)) stop("Must be a matrix!")
  if (grepl(".",name,fixed=TRUE)) stop("Name may not contain a dot!")
  data$data[[name]]=matrix
  data
}

`Condition<-` <- function(data, value) {
  data$coldata$Condition <- interaction(data$coldata[value])
  data
}
Condition <- function(data,value=NULL) {
  if (is.null(value)) data$coldata$Condition else {
    Condition(data)<-value
    data
  }
}


Genes=function(data, use.symbols=TRUE) data$gene.info[[if (use.symbols) "Symbol" else "Gene"]]


dim.grandR=function(data) c(dim(data$gene.info)[1],dim(data$coldata)[1])
is.grandR <- function(x) inherits(x, "grandR")
dimnames.grandR=function(data) dimnames(data$data$count)
print.grandR=function(data) cat(sprintf("grandR: %s\nRead from %s\n%d genes, %d samples/cells\nAvailable data slots: %s\nAvailable analyses: %s\nDefault data slot: %s\n",data$metadata$Description,data$prefix,nrow(data),ncol(data),paste(Slots(data),collapse=","),paste(Analyses(data),collapse=","),DefaultSlot(data)))
data.apply=function(data,fun,fun.gene.info=NULL,fun.coldata=NULL,...) {
  re=list()
  for (l1 in names(data$data)) {
    re[[l1]]=fun(data$data[[l1]],...)
  }
  ngene.info=if (!is.null(fun.gene.info)) fun.gene.info(data$gene.info,...) else data$gene.info
  ncoldata=if (!is.null(fun.coldata)) fun.coldata(data$coldata,...) else data$coldata
  invisible(grandR(parent=data,gene.info=ngene.info,data=re,coldata=ncoldata))
}

GeneInfo=function(data,column=NULL,value=NULL) {
  if (is.null(column)) data$gene.info else {
    data$gene.info[[colum]]=value
    data
  }
}
ColData=function(data,column=NULL,value=NULL) {
  if (is.null(column)) data$coldata else {
    data$coldata[[column]]=value
    data
  }
}
Analyses=function(data) names(data$analysis)
Columns=function(data,analysis=NULL) {
  if (is.null(analysis)) ColData(data)$Name else names(data$analysis[[analysis]])
}

`GeneInfo<-` <- function(data, column, value) {
  data$gene.info[[column]]=value
  data
}
`ColData<-` <- function(data, column, value) {
  data$coldata[[column]]=value
  data
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

check.analysis=function(data,analysis) analysis %in% names(data$analysis)
check.slot=function(data,slot) slot %in% names(data$data)
check.mode.slot=function(data,mode.slot) {
  sapply(strsplit(mode.slot,".",fixed=TRUE),function(spl) {
    if (length(spl)>2 || length(spl)==0) return(FALSE)
    if (length(spl)==1) check.slot(data,spl) else tolower(substr(spl[1],1,1)) %in% c("t","n","o") && check.slot(data,spl[2])
  })
}


#' Obtain the indices of the given genes
#'
#' @title ToIndex
#' @param data The grandR object
#' @param gene A vector of genes. Can be either numeric indices, gene names, gene symbols or a logical vector
#'
#' @return Numeric indices corresponding to the given genes
#'
#' @examples
#' sars <- ReadGRAND("https://zenodo.org/record/5834034/files/sars.tsv.gz", design=c("Cell",Design$dur.4sU,Design$Replicate), verbose=TRUE)
#' ToIndex(sars,c("MYC"))
#' ToIndex(sars,GeneInfo(sars)$Symbol=="MYC")
#'
#' @export
#'
ToIndex=function(data,gene) {
  if (is.numeric(gene)) return(gene)
  if (is.logical(gene) && length(gene)==nrow(data)) return(which(gene))
  if (all(gene %in% data$gene.info$Gene)) return(setNames(1:nrow(data),data$gene.info$Gene)[gene])
  if (all(gene %in% data$gene.info$Symbol)) return(setNames(1:nrow(data),data$gene.info$Symbol)[gene])
  stop("Could not find all genes!")
}



#' Obtain a genes x values table
#'
#' @title GetTable
#' @param data A grandR object
#' @param type Either a mode.slot (see details) or an analysis name. Can also be a vector; If NULL, \link{DefaultSlot}(data) is used
#' @param columns A vector of columns (either condition/cell names if the type is a mode.slot, or names in the output table from an analysis; use \link{Columns}(data,<analysis>) to learn which columns are available); all condition/cell names if NULL
#' @param genes Restrict the output table to the given genes
#' @param ntr.na For columns representing a 4sU naive sample, should types \emph{ntr},\emph{new.count} and \emph{old.count} be 0,0 and count (ntr.na=FALSE; can be any other slot than count) or NA,NA and NA (ntr.na=TRUE)
#' @param gene.info Should the table contain the \link{GeneInfo} values as well (at the beginning)?
#' @param summarize Should replicates by summarized? Can only be specified if columns is NULL; either a summarization matrix (\link{GetSummarizeMatrix}) or TRUE (in which case \link{GetSummarizeMatrix}(data) is called)
#' @param prefix Prepend each column in the output table (except for the gene.info columns) by the given prefix
#' @param name.by A column name of \link{ColData}(data). This is used as the rownames of the output table
#'
#' @return A data frame containing the desired values
#'
#' @details This is a convenience wrapper for \link{GetData} (values from data slots) and \link{GetAnalysisTable} (values from analyses). Types can refer to any of the two (and can be mixed). If there are types from both data and analyses, columns must be NULL. Otherwise columns must either be condition/cell names (if type refers to one or several data slots), or regular expressions to match against the names in the analysis tables.
#' @details To refer to data slots, the mode.slot syntax can be used: Each name is either a data slot, or one of (new,old,total) followed by a dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to obtain the \emph{new counts}.
#'
#' @seealso \link{GetData},\link{GetAnalysisTable},\link{DefaultSlot},\link{Genes},\link{GetSummarizeMatrix}
#'
#' @examples
#' sars <- ReadGRAND("https://zenodo.org/record/5834034/files/sars.tsv.gz", design=c("Condition",Design$dur.4sU,Design$Replicate), verbose=TRUE)
#' sars <- Normalize(FilterGenes(sars))
#'
#' head(GetTable(sars)) # DefaultSlot values, i.e. size factor normalized read counts for all samples
#' head(GetTable(sars,summarize=TRUE)) # DefaultSlot values averaged over the two conditions
#' head(GetTable(sars,type="new.count",columns=!ColData(sars)$no4sU)) # Estimated counts for new RNA for all samples with 4sU
#'
#' sars<-FitKinetics(sars,name = "kinetics",steady.state=list(Mock=TRUE,SARS=FALSE))
#' head(GetTable(sars,type="kinetics",columns="Half-life")) # Estimated RNA half-lives for both conditions
#'
#'
#'
#' @export
#'
GetTable=function(data,type=NULL,columns=NULL,genes=Genes(data),ntr.na=TRUE,gene.info=FALSE,summarize=NULL,prefix=NULL,name.by="Symbol") {
  if (!is.null(columns) && !is.null(summarize)) stop("columns and summarize may not be set simultaneously!")
  if (is.null(type)) type=DefaultSlot(data)
  if (is.null(columns)) columns=colnames(data)

  genes=ToIndex(data,genes)

  mode.slot=check.mode.slot(data,type)
  analysis=check.analysis(data,type) & !mode.slot
  if (!all(analysis|mode.slot)) stop(sprintf("Type %s is neither a mode.slot nor an analysis name!",paste(type[!analysis&!mode.slot],collapse=",")))

  # check that columns is only used if type is either completely analysis or mode.slot
  if (!is.null(columns) && sum(mode.slot)>0 && sum(analysis)>0) stop("Columns can only be specified if type either refers to mode.slots or analyses")

  # obtain mode.slot data
  r1=NULL
  if (any(mode.slot)) {
    r1=as.data.frame(t(GetData(data,type[mode.slot],columns=columns,genes,ntr.na = ntr.na,coldata=FALSE, melt=FALSE, name.by = name.by)))
    names(r1)=data$coldata[columns,"Name"]
    if (!is.null(summarize)) {
      if (is.logical(summarize) && length(summarize)==1) summarize=GetSummarizeMatrix(data)
      r1=as.data.frame(as.matrix(r1) %*% summarize)
    }
  }

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


#' Obtain a tidy table of values for a gene or a (preferentially) small set of genes
#'
#' @title GetData
#' @param data A grandR object
#' @param mode.slot Which kind of data to access (see details)
#' @param columns A vector of columns (i.e. condition/cell names; use colnames(data) to learn which columns are available); all condition/cell names if NULL
#' @param genes Restrict the output table to the given genes (this typically is a single gene, or very few genes)
#' @param melt Should the table be melted if multiple genes / mode.slots are given
#' @param coldata Should the table contain the \link{ColData} values as well (at the beginning)?
#' @param ntr.na For columns representing a 4sU naive sample, should mode.slot \emph{ntr},\emph{new.count} and \emph{old.count} be 0,0 and count (ntr.na=FALSE; can be any other slot than count) or NA,NA and NA (ntr.na=TRUE)
#' @param name.by A column name of \link{ColData}(data). This is used as the colnames of the output table
#'
#' @return A data frame containing the desired values
#'
#' @details To refer to data slots, the mode.slot syntax can be used: Each name is either a data slot, or one of (new,old,total) followed by a dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to obtain the \emph{new counts}.
#' @details If only one mode.slot and one gene is given, the output table contains one column (and potentially columns from \link{ColData}) named \emph{Value}. If one gene and multiple mode.slots are given, the columns are named according to the mode.slots. If one mode.slot and multiple genes are given, the columns are named according to the genes. If multiple genes and mode.slots are given, columns are named gene.mode.slot.
#' @details If melt=TRUE, the table is molten such that each row contains only one value (for one of the genes and for one of the mode.slots). If only one gene and one mode.slot is given, melting does not have an effect.
#'
#' @seealso \link{GetTable},\link{GetAnalysisTable},\link{DefaultSlot},\link{Genes}
#'
#' @examples
#' sars <- ReadGRAND("https://zenodo.org/record/5834034/files/sars.tsv.gz", design=c("Condition",Design$dur.4sU,Design$Replicate), verbose=TRUE)
#' GetData(sars,mode.slot="ntr",gene="MYC") # one gene, one mode.slot
#' GetData(sars,mode.slot=c("count","ntr"),gene="MYC",coldata = F) # one gene, multiple mode.slots
#' GetData(sars,mode.slot=c("count","ntr"),gene=c("SRSF6","MYC"),melt=TRUE) # multiple genes, multiple mode.slots, molten
#'
#' @export
#'
GetData=function(data,mode.slot=DefaultSlot(data),columns=NULL,genes=Genes(data),melt=FALSE,coldata=TRUE,ntr.na=TRUE,name.by="Symbol") {
  if (!all(check.mode.slot(data,mode.slot))) stop(sprintf("mode.slot %s unknown!",paste(mode.slot[!check.mode.slot(mode.slot)],collapse=",")))

  if (is.null(columns)) columns=colnames(data)
  genes=ToIndex(data,genes)
  og=if (name.by %in% names(data$gene.info)) data$gene.info[[name.by]][genes] else data$gene.info[genes,1]
  uno=function(mode.slot) {
    tno="t"
    spl=strsplit(mode.slot,".",fixed=TRUE)[[1]]
    if (length(spl)>1) {tno=spl[1]; mode.slot=spl[2];}
    mf = switch(tolower(substr(tno,1,1)),t=1,n=data$data$ntr[genes,columns],o=1-data$data$ntr[genes,columns],stop(paste0(mode.slot," unknown!")))
    if (!ntr.na) {
      mf[is.na(mf)]=if(tolower(substr(tno,1,1))=="n") 0 else 1
    }
    conv=if (mode.slot=="count") function(m) {mode(m) <- "integer";m} else if (mode.slot=="ntr" && !ntr.na) function(m) {m[is.na(m)]=0; m} else function(m) m

    if (!(mode.slot %in% names(data$data))) stop(paste0(mode.slot," unknown!"))
    if (length(genes)==1) data.frame(conv(data$data[[mode.slot]][genes,columns]*mf)) else as.data.frame(conv(t(data$data[[mode.slot]][genes,columns]*mf)))
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
  rownames(re)=NULL
  re
}


#' Obtain reference columns (samples or cells) for all columns (samples or cells) in the data set
#'
#' @title FindReferences
#' @param data A grandR object
#' @param reference Expression evaluating to a logical vector to indicate which columns are reference columns; evaluated in an environment having the columns of \link{ColData}(data)
#' @param group a vector of colnames in \link{ColData}(data)
#'
#' @return A named list for each sample or cell containing all corresponding reference columns
#'
#' @details Without any group, the list simply contains all references for each sample/cell. With groups defined, each list entry consists of all references from the same group.
#'
#' @seealso \link{ColData},\link{Findno4sUPairs}
#'
#' @examples
#' sars <- ReadGRAND("https://zenodo.org/record/5834034/files/sars.tsv.gz", design=c("Condition",Design$dur.4sU,Design$Replicate), verbose=TRUE)
#' FindReferences(sars,reference=no4sU) # obtain the corresponding no4sU sample for each sample
#' FindReferences(sars,Condition=="Mock",group="duration.4sU.original") # obtain for each sample the corresponding sample in the Mock condition
#'
#' @export
#'
FindReferences=function(data,reference, group="Condition") {
  if (!is.grandR(data)) stop("Data is not a grandR object!")

  df=ColData(data)
  df$group=if(is.null(group)) 1 else interaction(df[group],drop=FALSE,sep=".")
  e=substitute(reference)
  map=dlply(df,.(group),function(s) as.character(s$Name[eval(e,s,parent.frame())]))
  pairs=setNames(lapply(df$group,function(g) map[[g]]),df$Name)
  pairs
}




#' Remove analyses from the grandR object
#'
#' @title DropAnalysis
#' @param data A grandR object
#' @param pattern A regular expression that is matched to analysis names
#'
#' @return A new grandR object without the dropped analyses
#'
#' @seealso \link{Analyses}
#'
#' @export
#'
DropAnalysis=function(data,pattern=NULL) {
  if (is.null(pattern)) {
    data$analysis=NULL
  } else {
    data$analysis=data$analysis[!grepl(pattern,names(data$analysis))]
  }
  invisible(data)
}
#' Create a metatable for an analysis
#'
#' @title MakeAnalysis
#' @param name The user-defined analysis name
#' @param analysis The name of the analysis tool
#' @param mode An optional mode (new,old,total) on which the analysis has been run
#' @param slot An optional data slot on which the analysis has been run
#'
#' @return A metatable to be used as description parameter for \link{AddAnalysis}
#'
#' @seealso \link{AddAnalysis}
#'
#' @export
#'
MakeAnalysis=function(name,analysis,mode=NULL,slot=NULL) {
  list(name=name,mode=mode,analysis=analysis,slot=slot)
}

#' Create a new analysis within a grandR object
#'
#' @title AddAnalysis
#' @param data The grandR object
#' @param description A metatable created by \link{MakeAnalysis}
#' @param table The analysis table
#' @param warn.present Warn if an analysis with the same name is already present (and then overwrite)
#'
#' @return A new grandR object containing the given analysis
#'
#' @seealso \link{FitKinetics},\link{TestGenesLRT},\link{TestPairwise},\link{LFC}
#'
#' @export
#'
AddAnalysis=function(data,description,table,warn.present=TRUE) {
  stopifnot(!is.null(description$name))
  if (is.null(data$analysis)) data$analysis=list()
  if (is.null(data$analysis[[description$name]])) {
    data$analysis[[description$name]]=table
    attr(data$analysis[[description$name]],"analysis")=description
  } else {
    if (warn.present) warning(sprintf("Analysis %s already present! Overwritting...",description$name))
    for (n in names(table)) data$analysis[[description$name]][[n]]=table[[n]]
    ana = attr(data$analysis[[description$name]],"analysis")
    for (n in names(description)) ana[[n]]=description[[n]]
    attr(data$analysis[[description$name]],"analysis") = ana
  }
  invisible(data)
}


#' Obtain a table of analysis results values
#'
#' @title GetAnalysisTable
#' @param data A grandR object
#' @param names One or several analysis names (\link{Analyses})
#' @param columns Regular expression to select columsn from the analysis table
#' @param genes Restrict the output table to the given genes
#' @param gene.info Should the table contain the \link{GeneInfo} values as well (at the beginning)?
#' @param name.by A column name of \link{ColData}(data). This is used as the rownames of the output table
#'
#' @return A data frame containing the desired values
#'
#' @details The names for the output table are <Analysis name>.<columns name>
#'
#' @seealso \link{GetTable},\link{Genes}
#'
#' @examples
#' sars <- ReadGRAND("https://zenodo.org/record/5834034/files/sars.tsv.gz", design=c("Condition",Design$dur.4sU,Design$Replicate), verbose=TRUE)
#' GetData(sars,mode.slot="ntr",gene="MYC") # one gene, one mode.slot
#' GetData(sars,mode.slot=c("count","ntr"),gene="MYC",coldata = F) # one gene, multiple mode.slots
#' GetData(sars,mode.slot=c("count","ntr"),gene=c("SRSF6","MYC"),melt=TRUE) # multiple genes, multiple mode.slots, molten
#'
#' @export
#'
GetAnalysisTable=function(data,names=NULL,columns=NULL,genes=Genes(data),gene.info=TRUE,name.by="Symbol") {
  if (!all(check.analysis(data,names))) stop(sprintf("Analysis name %s unknown!",paste(analysis[!check.analysis(mode.slot)],collapse=",")))

  genes=ToIndex(data,genes)

  re=data$gene.info[genes,]
  if (!is.null(name.by)) {
    rownames(re)=if (name.by %in% names(data$gene.info)) data$gene.info[[name.by]][genes] else data$gene.info[genes,1]
  }
  sintersect=function(a,b) if (is.null(b)) a else intersect(a,b)
  for (name in names) {
    t=data$analysis[[name]][genes,]
    analysis = attr(data$analysis[[name]],"analysis")
    names(t)=if (is.null(analysis$mode)) paste(name,names(t),sep=".") else paste(name,analysis$mode,names(t),sep=".")
    if (!is.null(columns)) {
      use = rep(FALSE,ncol(t))
      for (r in columns) use=use|grepl(r,names(t))
      t=t[,use,drop=FALSE]
    }
    re=cbind(re,t)
  }
  if (is.logical(gene.info) && !gene.info) re=re[,(ncol(data$gene.info)+1):ncol(re),drop=FALSE]
  if (is.character(gene.info)) re=re[,-which(!names(data$gene.info) %in% gene.info),drop=FALSE]
  re
}
