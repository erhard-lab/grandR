
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
DefaultSlot <- function(data,value=NULL) {
  if (is.null(value)) data$metadata$default.slot else {
    data$metadata$default.slot=value
    data
  }
}

Slots=function(data) {
  names(data$data)
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

AddSlot=function(data,name,matrix) {
  if (!all(colnames(matrix)==colnames(data$data$count))) stop("Column names do not match!")
  if (!all(rownames(matrix)==rownames(data$data$count))) stop("Row names do not match!")
  if (!is.matrix(matrix)) stop("Must be a matrix!")
  if (grepl(".",name,fixed=TRUE)) stop("Name may not contain a dot!")
  data$data[[name]]=matrix
  data
}

dim.grandR=function(data) c(dim(data$gene.info)[1],dim(data$coldata)[1])
is.grandR <- function(x) inherits(x, "grandR")
dimnames.grandR=function(data) dimnames(data$data$count)
print.grandR=function(data) cat(sprintf("grandR: %s\nRead from %s\n%d genes, %d samples/cells\nAvailable data slots: %s\nDefault data slot: %s\nAnalyses: %s",data$metadata$Description,data$prefix,nrow(data),ncol(data),paste(names(data$data),collapse=","),DefaultSlot(sars),Analyses(sars)))
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
    if (length(spl)==1) check.slot(data,mode.slot) else tolower(substr(spl[1],1,1)) %in% c("t","n","o") && check.slot(data,spl[2])
  })
}


ToIndex=function(data,gene) {
  if (is.numeric(gene)) return(gene)
  if (is.logical(gene) && length(gene)==nrow(data)) return(which(gene))
  if (all(gene %in% data$gene.info$Gene)) return(setNames(1:nrow(data),data$gene.info$Gene)[gene])
  if (all(gene %in% data$gene.info$Symbol)) return(setNames(1:nrow(data),data$gene.info$Symbol)[gene])
  stop("Could not find all genes!")
}



GetTable=function(data,type,columns=NULL,genes=Genes(data),ntr.na=TRUE,gene.info=FALSE,summarize=NULL,prefix=NULL,name.by="Symbol") {
  if (!is.null(columns) && !is.null(summarize)) stop("columns and summarize may not be set simultaneously!")
  if (length(type)!=1) stop("Only one type for table output!")
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
  if(length(mode.slot)==1 && length(genes)==1) names(re)="Value" else if (length(mode.slot)==1) names(re)=og else if (length(genes)==1) names(re)=mode.slot else names(re)=paste0(rep(og,length(mode.slot))," ",rep(mode.slot,each=length(og)))
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


FindReferences=function(data,reference, group=NULL) {
  stopifnot(is.grandR(data))
  df=data$coldata
  df$group=if(is.null(group)) 1 else interaction(df[group],drop=FALSE,sep=".")
  e=substitute(reference)
  map=dlply(df,.(group),function(s) as.character(s$Name[eval(e,s,parent.frame())]))
  pairs=setNames(lapply(df$group,function(g) map[[g]]),df$Name)
  pairs
}




DropAnalysis=function(data,pattern=NULL) {
  if (is.null(pattern)) {
    data$analysis=NULL
  } else {
    data$analysis=data$analysis[!grepl(pattern,names(data$analysis))]
  }
  invisible(data)
}
Analysis=function(name,analysis,mode=NULL,slot=NULL) {
  list(name=name,mode=mode,analysis=analysis,slot=slot)
}
AddAnalysis=function(data,description,table,warn.present=TRUE) {
  stopifnot(!is.null(description$name))
  if (is.null(data$analysis)) data$analysis=list()
  if (is.null(data$analysis[[description$name]])) {
    data$analysis[[description$name]]=table
    attr(data$analysis[[description$name]],"analysis")=description
  } else {
    if (warn.present) warning(sprintf("Analysis %s already present!",description$name))
    for (n in names(table)) data$analysis[[description$name]][[n]]=table[[n]]
    ana = attr(data$analysis[[description$name]],"analysis")
    for (n in names(description)) ana[[n]]=description[[n]]
    attr(data$analysis[[description$name]],"analysis") = ana
  }
  invisible(data)
}


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
