
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

comp.fpkm=function(cmat,lengths,subset=NULL) {
  scale=colSums(if(!is.null(subset)) cmat[subset,] else cmat,na.rm=T)/1E6
  rpm=t(t(cmat)/scale)

  zerolen=lengths==0
  lengths[zerolen]=1
  re=rpm/(lengths/1000)
  re[zerolen,]=NA
  re
}
comp.rpm=function(cmat,subset=NULL) {
  scale=colSums(if(!is.null(subset)) cmat[subset,] else cmat,na.rm=T)/1E6
  re=t(t(cmat)/scale)
  re
}

#' Compute NTR quantiles
#'
#' Computes quantiles from the NTR posterior and puts them into a new slot
#'
#' @param data the grandR object
#' @param quantile which quantile to compute
#' @param name the name of the new slot to put quantile values in
#' @param name.lower the name of the new slot to put the lower bound of the CI in
#' @param name.upper the name of the new slot to put the upper bound of the CI in
#' @param CI.size A number between 0 and 1 representing the size of the credible interval
#'
#' @details The NTR posterior distribution can be approximated by a beta distribution.
#'
#' @details ComputeNtrPosteriorQuantile computes any quantile from this Beta approximation
#'
#' @details ComputeNtrPosteriorLower computes the (1-CI.size)/2 quantile
#'
#' @details ComputeNtrPosteriorUpper computes the 1-(1-CI.size)/2 quantile
#'
#' @details ComputeNtrCI computes both of these quantiles.
#'
#' @return a new grandR object containing an additional slot
#' @export
#'
ComputeNtrPosteriorQuantile=function(data,quantile,name) {
  a=as.matrix(GetTable(data,type="alpha",name.by = "Gene"))
  b=as.matrix(GetTable(data,type="beta",name.by = "Gene"))
  v=qbeta(quantile,a,b)
  AddSlot(data,name,v)
}
#' @rdname ComputeNtrPosteriorQuantile
#' @export
ComputeNtrCI=function(data,CI.size=0.95,name.lower="lower",name.upper="upper") {
  data=ComputeNtrPosteriorLower(data=data,CI.size=CI.size,name=name.lower)
  data=ComputeNtrPosteriorUpper(data=data,CI.size=CI.size,name=name.upper)
  data
}
#' @rdname ComputeNtrPosteriorQuantile
#' @export
ComputeNtrPosteriorLower=function(data,CI.size=0.95,name="lower") ComputeNtrPosteriorQuantile(data=data,quantile=(1-CI.size)/2,name=name)
#' @rdname ComputeNtrPosteriorQuantile
#' @export
ComputeNtrPosteriorUpper=function(data,CI.size=0.95,name="upper") ComputeNtrPosteriorQuantile(data=data,quantile=1-(1-CI.size)/2,name=name)


#' Steady state half-lives for each sample
#'
#' Transforms each NTR to a half-life value (assuming steady state gene expression) and puts them into a new slot or adds an analysis
#'
#' @param data the grandR object
#' @param time either a number indicating the labeling time, or a name of the \link{Coldata} table
#' @param name the name of the new slot/analysis to put half-life values in
#' @param max.HL all values above this will be set to this
#' @param CI.size A number between 0 and 1 representing the size of the credible interval
#' @param compute.CI it TRUE, credible intervals are computed, this also sets as.analysis to TRUE
#' @param as.analysis if TRUE add the results as analysis and not as data slot
#'
#' @details An NTR value p can be transformed into an RNA half-live using the equation
#' log(2)/(-1/t*log(1-p))
#' This is described in our GRAND-SLAM paper (Juerges et al., Bioinformatics 2018).
#'
#' @return a new grandR object with an additional slot or analysis
#' @export
#'
ComputeSteadyStateHalfLives=function(data,time=Design$dur.4sU,name, max.HL=48, CI.size=0.95, compute.CI=FALSE, as.analysis=FALSE) {
  if (is.character(time) && length(time)==1) time=data$coldata[[time]]

  ntrs=as.matrix(GetTable(data,type="ntr",name.by = "Gene"))
  if (length(time)==1) time=rep(time,ncol(ntrs))
  stopifnot(ncol(ntrs)==length(time))

  if (compute.CI) {
    as.analysis=TRUE
    data=ComputeNtrCI(data,CI.size = CI.size)
    lower=as.matrix(GetTable(data,type="lower",name.by = "Gene"))
    upper=as.matrix(GetTable(data,type="upper",name.by = "Gene"))
    hls = do.call(cbind,lapply(1:length(time),function(i) cbind(
      pmin(comp.hl(p = upper[,i],time = time[i]),max.HL),
      pmin(comp.hl(p = ntrs[,i],time = time[i]),max.HL),
      pmin(comp.hl(p = lower[,i],time = time[i]),max.HL)
    )))
    colnames(hls)=paste0(rep(c("CI.lower.","Half-life","CI.upper."),ncol(ntrs)),rep(colnames(ntrs),each=3))
  }
  else {
    hls = sapply(1:length(time),function(i) comp.hl(p = ntrs[,i],time = time[i]))
    colnames(hls)=colnames(ntrs)
  }

  if (as.analysis) {
    AddAnalysis(data,name = name,table = as.data.frame(hls))
  } else {
    AddSlot(data,name,hls)
  }
}


#' Compute absolute expression using ERCC spike ins
#'
#' Compute absolute expression in a grandR object and puts the normalized data into a new slot
#'
#' @param data the grandR object
#' @param dilution the dilution of the spikein transcript in the lysis reaction mix
#' @param volume 	the approximate volume of the lysis chamber (nanoliters)
#' @param slot the slot containing relative expression values
#' @param name the name of the new slot to put absolute expression values in
#'
#' @seealso \link[monocle]{relative2abs}
#' @return a new grandR object with an additional slot
#' @export
#'
ComputeAbsolute=function(data,dilution=4E4,volume=10,slot="tpm",name="absolute") {
  fd=data.frame(gene_short_name=Genes(data))
  rownames(fd)=Genes(data)
  mat=GetTable(data,type=slot)
  ercc=mat[GeneInfo(data,"Type")=='ERCC',]

  cds <- newCellDataSet(mat,
                        featureData = new ("AnnotatedDataFrame",fd),
                        lowerDetectionLimit = 0.1,
                        expressionFamily = tobit(Lower = 0.1))
  rpc_matrix <- monocle::relative2abs(cds, method = "num_genes", ERCC_controls=ercc, ERCC_annotation=monocle::spike_df,dilution=dilution, volume=volume)
  rpc_matrix[mat]=0

  colnames(rpc_matrix)=colnames(data)
  rownames(rpc_matrix)=Genes(data)
  rpc_matrix[,grepl("test",colnames(rpc_matrix))]=NA
  AddSlot(data,name,rpc_matrix)
}



#' Normalization
#'
#' Normalizes data in a grandR object and puts the normalized data into a new slot
#'
#' @param data the grandR object
#' @param genes compute the normalization w.r.t. these genes (see details)
#' @param name the name of the new slot for the normalized data
#' @param slot the name of the slot for the data to normalize
#' @param set.to.default set the new slot as the default slot
#' @param size.factors numeric vector; if not NULL, use these size factors instead of computing size factors
#' @param return.sf return the size factors and not a grandR object
#' @param tlen the transcript lengths (for FPKM and TPM)
#'
#' @details Normalize will perform DESeq2 normalization, i.e. it will use \link[DESeq2]{estimateSizeFactorsForMatrix}
#' to estimate size factors, and divide each value by this. If genes are given, size factors will be computed only w.r.t. these genes (but then all genes are normalized).
#'
#' @details NormalizeFPKM will compute fragments per kilobase and million mapped reads. If genes are given, the scaling factor
#' will only be computed w.r.t. these genes (but then all genes are normalized).
#'
#' @details NormalizeRPM will compute reads per million mapped reads. If genes are given, the scaling factor
#' will only be computed w.r.t. these genes (but then all genes are normalized).
#'
#' @details NormalizeTPM will compute transcripts per million mapped reads. If genes are given, the scaling factor
#' will only be computed w.r.t. these genes (but then all genes are normalized).
#'
#' @details Genes can be referred to by their names, symbols, row numbers in the gene table, or a logical vector referring to the gene table rows.
#'
#' @return a new grandR object with a new data slot
#'
#' @seealso \code{\link{NormalizeBaseline}}
#'
#' @export
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#'
#' sars <- Normalize(sars)
#' DefaultSlot(sars)
Normalize=function(data,genes=Genes(data),name="norm",slot="count",set.to.default=TRUE,size.factors=NULL,return.sf=FALSE) {
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,type=slot,genes=genes,ntr.na = FALSE,name.by = "Gene"))
  if (is.null(size.factors)) size.factors=DESeq2::estimateSizeFactorsForMatrix(mat)
  if (return.sf) return(size.factors)

  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  data=AddSlot(data,name,matrix = t(t(mat)/size.factors),set.to.default=set.to.default)
  data
}
#' @rdname Normalize
#' @export
NormalizeFPKM=function(data,genes=Genes(data),name="fpkm",slot="count",set.to.default=TRUE,tlen=GeneInfo(data,"Length")) {
  genes=ToIndex(data,genes)
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  data=AddSlot(data,name,comp.fpkm(mat,tlen,subset = genes),set.to.default=set.to.default)
  data
}
#' @rdname Normalize
#' @export
NormalizeRPM=function(data,genes=Genes(data),name="rpm",slot="count",set.to.default=TRUE) {
  genes=ToIndex(data,genes)
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  data=AddSlot(data,name,comp.rpm(mat,subset = genes),set.to.default=set.to.default)
  data
}
#' @rdname Normalize
#' @export
NormalizeTPM=function(data,genes=Genes(data),name="tpm",slot="count",set.to.default=TRUE,tlen=GeneInfo(data,"Length")) {
  genes=ToIndex(data,genes)
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  data=AddSlot(data,name,comp.tpm(mat,tlen,subset = genes),set.to.default=set.to.default)
  data
}

#' Normalization to a baseline
#'
#' Normalizes data in a grandR object to a baseline and puts the normalized data into a new slot
#'
#' @param data the grandR object
#' @param baseline matrix defining the corresponding baseline (row) for each column (sample or cell; see details)
#' @param name the name of the new slot for the normalized data
#' @param slot the name of the slot for the data to normalize
#' @param set.to.default set the new slot as the default slot
#' @param LFC.fun either \link[lfc]{NormLFC} or \link[lfc]{PsiLFC} from the lfc package
#' @param ... forwarded to LFC.fun
#'
#' @details Baseline normalization computes the log2 fold change for a column (i.e. sample or cell) to a baseline columns (or several baseline columns). This is by default done using the
#' \code{\link[lfc]{PsiLFC}} function from the lfc package, which, by default, also normalizes log2 fold changes by adding a constant
#' such that the median is zero.
#'
#' @details Baselines are defined by a square logical matrix, defining for each sample or cell of the grandR object, represented by the column of the matrix,
#' which samples or cells are indeed the baseline (represented by the rows). Such matrices can conveniently be obtained by \code{\link{FindReferences}}.
#'
#' @return a new grandR object with an additional slot
#'
#' @seealso \code{\link{Normalize}},\code{\link{FindReferences}}
#' @export
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#' blmat <- FindReferences(sars,reference = duration.4sU==0, group = "Cell")  # the Mock.no4sU or SARS.no4sU sample are the baselines for each sample
#' sars <- NormalizeBaseline(data,baseline=blmat)
#' head(GetTable(data,type="baseline"))
NormalizeBaseline=function(data,baseline=FindReferences(data,reference=Condition==levels(Condition)[1]),name="baseline",slot=DefaultSlot(data),set.to.default=FALSE,LFC.fun=PsiLFC,...) {
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  mat=sapply(colnames(baseline),function(n) LFC.fun(mat[,n],rowMeans(mat[,baseline[,n],drop=FALSE]),...))
  data=AddSlot(data,name,mat,set.to.default=set.to.default)
  data
}


#' Scale data
#'
#' Compute values for all genes standardized (i.e. z scores) across samples.
#'
#' @param data a grandR object
#' @param name the new slot name
#' @param slot the slot from where to take values
#' @param set.to.default set the new slot as default slot
#' @param group Perform standardization per group of columns (see details)
#' @param center Perform centering (forwarded to \link[base]{scale})
#' @param scale Perform scaling (forwarded to \link[base]{scale})
#'
#' @details Standardization can be done per group. For this, the group parameter has to be a name of the \code{\link{Coldata}} table,
#' to define groups of columns (i.e. samples or cells).
#'
#' @seealso \link[base]{scale}
#'
#' @return a new grandR object with a new slot
#' @export
#'
Scale=function(data,name="scaled",slot=DefaultSlot(data),set.to.default=FALSE,group=NULL,center=TRUE,scale=TRUE) {
  stopifnot(is.grandR(data))
  mat=as.matrix(GetTable(data,type=slot,ntr.na = FALSE,name.by = "Gene"))
  if (is.null(group)) {
    mat=t(scale(t(mat),center = center, scale = scale))
  } else {
    gr=Coldata(data)[[group]]
    for (c in unique(gr)) mat[,gr==c]=t(scale(t(mat[,gr==c]),center = center, scale = scale))
  }
  data=AddSlot(data,name,mat,set.to.default=set.to.default)
  data
}



#' Filter genes
#'
#' Return a grandR object with fewer genes than the given grandR object (usually to filter out weakly expressed genes).
#'
#' @param data the grandR object
#' @param mode.slot the mode.slot that is used for filtering (see details)
#' @param minval the minimal value for retaining a gene
#' @param mincol the minimal number of columns (i.e. samples or cells) a gene has to have a value >= minval
#' @param min.cond if not NULL, do not compare values per column, but per condition (see details)
#' @param use if not NULL, defines the genes directly that are supposed to be retained (see details)
#' @param keep if not NULL, defines genes directly, that should be kept even though they do not adhere to the filtering critera (see details)
#' @param return.genes if TRUE, return the gene names instead of a new grandR object
#'
#' @details By default genes are retained, if they have 100 read counts in at least half of the columns (i.e. samples or cells).
#'
#' @details The \code{use} parameter can be used to define genes to be retained directly. The \code{keep} parameter, in contrast, defines
#' \emph{additional} genes to be retained. For both, genes can be referred to by their names, symbols, row numbers in the gene table,
#' or a logical vector referring to the gene table rows.
#'
#' @details To refer to data slots, the mode.slot syntax can be used: Each name is either a data slot, or one of (new,old,total)
#' followed by a dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to filter by \emph{new counts}.
#'
#' @details if the \code{min.cond} parameter is given, first all columns belonging to the same \code{\link{Condition}} are summed up, and then the usual filtering
#' is performed by conditions instead of by columns.
#'
#' @return either a new grandR object (if return.genes=FALSE), or a vector containing the gene names that would be retained
#'
#' @export
#'
#' @examples
#'
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#'
#' nrow(sars)                                        # This is already filtered and has 1045 genes
#' nrow(FilterGenes(sars,minval=1000))               # There are 966 genes with at least 1000 read counts in half of the samples
#' nrow(FilterGenes(sars,minval=10000,min.cond=1))   # There are 944 genes with at least 10000 read counts in the Mock or SARS condition
#' nrow(FilterGenes(sars,use=GeneInfo(sars,"Type")!="Cellular"))
#'                                                   # These are the 11 viral genes.
#'
FilterGenes=function(data,mode.slot='count',minval=100,mincol=ncol(data)/2,min.cond=NULL,use=NULL,keep=NULL,return.genes=FALSE) {
  if (!is.null(use) & !is.null(keep)) stop("Do not specify both use and keep!")

  if (length(mode.slot)!=1) stop("One mode.slot must be given!")
  if (!all(check.mode.slot(data,mode.slot))) stop(sprintf("mode.slot %s unknown!",paste(mode.slot[!check.mode.slot(data,mode.slot)],collapse=",")))

  if (is.null(use)) {
    summi = if (!is.null(min.cond)) GetSummarizeMatrix(data,no4sU=TRUE,average=FALSE) else NULL
    mincol=if (!is.null(min.cond)) min.cond else mincol
    t=GetTable(data,type=mode.slot,summarize = summi)
    use=apply(t,1,function(v) sum(v>=minval,na.rm=TRUE)>=mincol)
    if (!is.null(keep)) use = use | rownames(t) %in% rownames(t[keep,])
  } else {
    use=ToIndex(data,use)
  }

  if (return.genes) return(Genes(data,use))
  return(data.apply(data,function(t) t[use,],fun.gene.info = function(t) t[use,]))
}


#' Expression percentage computation
#'
#' Compute the expression percentage for a particular set of genes.
#'
#' @param data the grandR object
#' @param name the new name by which this is added to the Coldata
#' @param genes define the set of genes to compute the percentage for
#' @param mode.slot which mode.slot to take the values for computing the percentage from
#' @param multiply.by.100 if TRUE, compute percentage values, otherwise fractions between 0 and 1
#'
#' @seealso \code{\link{Coldata}}
#'
#' @details Genes can be referred to by their names, symbols, row numbers in the gene table, or a logical vector referring to the gene table rows.
#'
#' @details To refer to data slots, the mode.slot syntax can be used: Each name is either a data slot, or one of (new,old,total)
#' followed by a dot followed by a slot. For new or old, the data slot value is multiplied by ntr or 1-ntr. This can be used e.g. to filter by \emph{new counts}.
#'
#' @return a new grandR object having the expression percentage in its Coldata table
#' @export
#'
#' @examples
ComputeExpressionPercentage=function(data,name,genes,mode.slot=DefaultSlot(data),multiply.by.100=TRUE) {
  gof=colSums(GetTable(data,type=mode.slot,ntr.na = FALSE,genes = genes))
  total=colSums(GetTable(data,type=mode.slot,ntr.na = FALSE))
  percentage=gof/total
  if (multiply.by.100) percentage=percentage*100
  Coldata(data,name)=percentage
  data
}
