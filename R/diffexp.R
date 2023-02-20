

#' Compute the Bayesian information criterion (BIC)
#'
#' Compute the delta BIC for a list of potential models
#'
#' @param data A grandR object
#' @param name the user defined analysis name to store the results
#' @param mode either "total", "new" or "old"
#' @param normalization normalize on "total", "new", or "old" (see details)
#' @param formulas list of formulas specifying the models (you can use any column name from the \code{\link{Coldata}(data)})
#' @param no4sU Use no4sU columns (TRUE) or not (FALSE)
#' @param columns logical vector of which columns (samples or cells) to use (or NULL: use all)
#' @param verbose Print status updates
#'
#' @details DESeq2 by default performs size factor normalization. When computing differential expression of new RNA,
#' it might be sensible to normalize w.r.t. to total RNA, i.e. use the size factors computed from total RNA instead of computed from new RNA.
#' This can be accomplished by setting mode to "new", and normalization to "total"!
#'
#' @return a new grandR object including a new analysis table. The columns of the new analysis table are named as <name in list>.dBIC
#'
#' @export
#'
#' @concept diffexp
DESeq2BIC=function(data,name="BIC",mode="total",normalization=mode,formulas=list(Condition=~Condition, Background=~1),no4sU=FALSE,columns=NULL,verbose=FALSE) {

  checkPackages("DESeq2")

  mode.slot=paste0(mode,".count")
  normalization=paste0(normalization,".count")
  if (!check.mode.slot(data,mode.slot)) stop("Invalid mode")

  columns=if (is.null(columns)) no4sU | !Coldata(data)$no4sU else columns&(no4sU | !Coldata(data)$no4sU)

  colData=droplevels(Coldata(data)[columns,])
  mat=GetTable(data,type=mode.slot,columns = columns,ntr.na = FALSE)

  norm.mat=as.matrix(GetTable(data,type=normalization,columns = columns))

  colData=as.data.frame(lapply(colData,function(c) {
    if (is.factor(c)) levels(c)=make.names(levels(c),unique = TRUE)
    c
  }))

  df=sapply(names(formulas),function(name) {
    formula=formulas[[name]]
    if (verbose) cat(sprintf("Fitting %s...\n",name))
    model.mat = stats::model.matrix.default(formula, data = Coldata(data))

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnt(mat),colData=colData,design = formula)
    DESeq2::sizeFactors(dds)<-DESeq2::estimateSizeFactorsForMatrix(norm.mat)
    dds <- DESeq2::DESeq(dds,quiet=!verbose)
    #data.frame(
    #  deviance = S4Vectors::mcols(dds)$deviance,
    #  AIC = S4Vectors::mcols(dds)$deviance + ncol(model.mat)*2,
    #  BIC = S4Vectors::mcols(dds)$deviance + ncol(model.mat)*ncol(mat)
    #)
    S4Vectors::mcols(dds)$deviance + ncol(model.mat)*ncol(mat)
  })
  df=as.data.frame(df-apply(df,1,min))
  colnames(df)=paste0(colnames(df),'.dBIC')
  rownames(df)=rownames(mat)
  AddAnalysis(data,name = name,table = df)
}



#' Compute a likelihood ratio test.
#'
#' The test is computed on any of total/old/new counts using DESeq2 based on two nested models
#' specified using formulas.
#'
#' @param data A grandR object
#' @param name the user defined analysis name to store the results
#' @param mode either "total", "new" or "old"
#' @param slot which slot to use (should be a count slot, not normalized values)
#' @param normalization normalize on "total", "new", or "old" (see details)
#' @param target formula specifying the target model (you can use any column name from the \code{\link{Coldata}(data)})
#' @param background formula specifying the background model (you can use any column name from the \code{\link{Coldata}(data)})
#' @param columns logical vector of which columns (samples or cells) to use (or NULL: use all)
#' @param logFC compute and add the log2 fold change as well
#' @param verbose Print status updates
#'
#' @details This is a convenience wrapper around the likelihood ratio test implemented in DESeq2.
#'
#' @details DESeq2 by default performs size factor normalization. When computing differential expression of new RNA,
#' it might be sensible to normalize w.r.t. to total RNA, i.e. use the size factors computed from total RNA instead of computed from new RNA.
#' This can be accomplished by setting mode to "new", and normalization to "total"!
#'
#' @return a new grandR object including a new analysis table. The columns of the new analysis table are
#' \itemize{
#'  \item{"M"}{the base mean}
#'  \item{"S"}{the difference in deviance between the reduced model and the full model}
#'  \item{"P"}{the likelihood ratio test P value}
#'  \item{"Q"}{same as P but Benjamini-Hochberg multiple testing corrected}
#'  \item{"LFC"}{the log2 fold change for the target model (only with the logFC parameter set to TRUE)}
#' }
#'
#' @export
#'
#' @concept diffexp
LikelihoodRatioTest=function(data,name="LRT",mode="total",slot="count",normalization=mode,target=~Condition,background=~1,columns=NULL, logFC=FALSE, verbose=FALSE) {
  checkPackages("DESeq2")

  mode.slot=paste0(mode,".",slot)
  if (!check.mode.slot(data,mode.slot)) stop("Invalid mode")

  if (!is.null(normalization)) {
    if (!check.mode.slot(data,normalization)) normalization=paste0(normalization,".",slot)
  } else {
    normalization = mode.slot
  }
  if (!check.mode.slot(data,normalization)) stop("Invalid normalization!")


  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)


	colData=droplevels(Coldata(data)[columns,])
	mat=GetTable(data,type=mode.slot,columns = columns,ntr.na = FALSE)

	norm.mat=as.matrix(GetTable(data,type=normalization,columns = columns))

	colData=as.data.frame(lapply(colData,function(c) {
	  if (is.factor(c)) levels(c)=make.names(levels(c),unique = TRUE)
	  c
	}))


	dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnt(mat),colData=colData,design = target)
	DESeq2::sizeFactors(dds)<-DESeq2::estimateSizeFactorsForMatrix(norm.mat)
	dds <- DESeq2::DESeq(dds,test="LRT", reduced=background,quiet=!verbose)
	res <- DESeq2::results(dds)

	df=data.frame(
    M=res$baseMean,
    S=res$stat,
    P=res$pvalue,
    Q=p.adjust(res$pvalue,method="BH")
  )
	if (logFC) df$LFC=res$log2FoldChange

  rownames(df)=rownames(mat)
  AddAnalysis(data,name = name,table = df)
}


#' Apply a function over contrasts
#'
#' Helper function to run many pairwise comparisons using a contrast matrix
#'
#' @param data the grandR object
#' @param analysis a plain name, only used for status messages
#' @param name.prefix the prefix for the new analysis name; a dot and the column names of the contrast matrix are appended; can be NULL (then only the contrast matrix names are used)
#' @param contrasts contrast matrix that defines all pairwise comparisons, generated using \link{GetContrasts}
#' @param mode.slot which slot to take expression values from
#' @param FUN a function taking 1. the data matrix, 2. a logical vector indicating condition A and 3. a logical vector indicating condition B
#' @param verbose print status messages?
#' @param ... further parameters forward to FUN
#'
#' @details To implement most pairwise analyses, you only have to define FUN; see the source code of \link{LFC} for an example!
#'
#' @return a new grandR object with added analysis tables (that were returned by FUN)
#'
#' @seealso \link{LFC},\link{PairwiseDESeq2},\link{GetContrasts}
#'
#' @export
#'
#' @concept helper
ApplyContrasts=function(data,analysis,name.prefix,contrasts,mode.slot=NULL,verbose=FALSE,FUN,...) {
  if (is.null(mode.slot)) stop("Need to specify mode.slot!")

  mat=as.matrix(GetTable(data,type=mode.slot,ntr.na=FALSE))
  mode.slot=get.mode.slot(data,mode.slot)
  for (n in names(contrasts)) {
    if (verbose) cat(sprintf("Computing %s for %s...\n",analysis,n))
    re.df=FUN(mat,contrasts[[n]]==1,contrasts[[n]]==-1,...)
    data=AddAnalysis(data,name = if (is.null(name.prefix)) n else paste0(name.prefix,".",n),table = re.df)
  }
  data
}


#' Log2 fold changes and Wald tests for differential expression
#'
#' This function is a shortcut for first calling \link{PairwiseDESeq2} and then \link{LFC}.
#'
#' @param data the grandR object
#' @param name.prefix the prefix for the new analysis name; a dot and the column names of the contrast matrix are appended; can be NULL (then only the contrast matrix names are used)
#' @param contrasts contrast matrix that defines all pairwise comparisons, generated using \link{GetContrasts}
#' @param LFC.fun function to compute log fold changes (default: \link[lfc]{PsiLFC}, other viable option: \link[lfc]{NormLFC})
#' @param slot the slot of the grandR object to take the data from; should contain counts!
#' @param mode compute LFCs for "total", "new", or "old" RNA
#' @param normalization normalize on "total", "new", or "old" (see details)
#' @param verbose print status messages?
#'
#' @details Both \link[lfc]{PsiLFC} and  \link[lfc]{NormLFC}) by default perform normalization by subtracting the median log2 fold change from all log2 fold changes.
#' When computing LFCs of new RNA, it might be sensible to normalize w.r.t. to total RNA, i.e. subtract the median log2 fold change of total RNA from all the log2 fold change of new RNA.
#' This can be accomplished by setting mode to "new", and normalization to "total"!

#' @details Normalization can also be a mode.slot! Importantly, do not specify a slot containing normalized values, but specify a slot of unnormalized values
#' (which are used to compute the size factors for normalization!)
#'
#' @return a new grandR object including a new analysis table. The columns of the new analysis table are
#' \itemize{
#'  \item{"M"}{the base mean}
#'  \item{"S"}{the log2FoldChange divided by lfcSE}
#'  \item{"P"}{the Wald test P value}
#'  \item{"Q"}{same as P but Benjamini-Hochberg multiple testing corrected}
#'  \item{"LFC"}{the log2 fold change}
#' }
#' @seealso \link{PairwiseDESeq2},\link{GetContrasts}
#' @export
#'
#' @concept diffexp
Pairwise=function(data, name.prefix = mode, contrasts, LFC.fun=lfc::PsiLFC, slot="count", mode="total",
               normalization=mode,
               verbose=FALSE) {

  contrasts = contrasts[,apply(contrasts,2,function(v) all(c(-1,1) %in% v)),drop=FALSE]
  if (ncol(contrasts)==0) stop("Contrasts do not define any comparison!")

  data=LFC(data,name.prefix = name.prefix, contrasts = contrasts, LFC.fun = LFC.fun, slot=slot, mode=mode,normalization = normalization,verbose=verbose)
  data=PairwiseDESeq2(data,name.prefix = name.prefix, contrasts = contrasts, slot=slot, mode=mode,normalization = normalization,verbose=verbose)
  data

}

#' Estimation of log2 fold changes
#'
#' Estimate the log fold changes based on a contrast matrix, requires the LFC package.
#'
#' @param data the grandR object
#' @param name.prefix the prefix for the new analysis name; a dot and the column names of the contrast matrix are appended; can be NULL (then only the contrast matrix names are used)
#' @param contrasts contrast matrix that defines all pairwise comparisons, generated using \link{GetContrasts}
#' @param slot the slot of the grandR object to take the data from; for \link[lfc]{PsiLFC}, this really should be "count"!
#' @param LFC.fun function to compute log fold changes (default: \link[lfc]{PsiLFC}, other viable option: \link[lfc]{NormLFC})
#' @param mode compute LFCs for "total", "new", or "old" RNA
#' @param normalization normalize on "total", "new", or "old" (see details)
#' @param verbose print status messages?
#' @param ... further arguments forwarded to LFC.fun
#'
#' @details Both \link[lfc]{PsiLFC} and  \link[lfc]{NormLFC}) by default perform normalization by subtracting the median log2 fold change from all log2 fold changes.
#' When computing LFCs of new RNA, it might be sensible to normalize w.r.t. to total RNA, i.e. subtract the median log2 fold change of total RNA from all the log2 fold change of new RNA.
#' This can be accomplished by setting mode to "new", and normalization to "total"!
#'
#' @details Normalization can also be a mode.slot! Importantly, do not specify a slot containing normalized values, but specify a slot of unnormalized values
#' (which are used to compute the size factors for normalization!)
#'
#' @return a new grandR object including a new analysis table. The columns of the new analysis table are
#' \itemize{
#'  \item{"LFC"}{the log2 fold change}
#' }
#'
#' @seealso \link{PairwiseDESeq2},\link{GetContrasts}
#' @export
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c(Design$Condition,Design$dur.4sU,Design$Replicate))
#' sars <- subset(sars,Coldata(sars,Design$dur.4sU)==2)
#' sars<-LFC(sars,mode="total",contrasts=GetContrasts(sars,contrast=c("Condition","Mock")))
#' sars<-LFC(sars,mode="new",normalization="total",
#'                             contrasts=GetContrasts(sars,contrast=c("Condition","Mock")))
#' head(GetAnalysisTable(sars))
#'
#' @concept diffexp
LFC=function(data, name.prefix = mode, contrasts, slot="count",LFC.fun=lfc::PsiLFC, mode="total",
             normalization=NULL,
             verbose=FALSE,...) {

  contrasts = contrasts[,apply(contrasts,2,function(v) all(c(-1,1) %in% v)),drop=FALSE]
  if (ncol(contrasts)==0) stop("Contrasts do not define any comparison!")

  mode.slot=paste0(mode,".",slot)
  if (!is.null(normalization)) {
    if (!check.mode.slot(data,normalization)) normalization=paste0(normalization,".",slot)
    if (!check.mode.slot(data,normalization)) stop("Invalid normalization!")
  }
  if (!check.mode.slot(data,mode.slot)) stop("Invalid mode")
  ApplyContrasts(data,name.prefix=name.prefix,contrasts=contrasts,mode.slot=mode.slot,verbose=verbose,analysis="LFC",FUN=function(mat,A,B) {
    if (!is.null(normalization)) {
      norm.mat=as.matrix(GetTable(data,type=normalization,ntr.na=FALSE))
      nlfcs=LFC.fun(rowSums(norm.mat[,A,drop=FALSE]),rowSums(norm.mat[,B,drop=FALSE]),normalizeFun=function(i) i)
      med.element=median(nlfcs)
      lfcs=LFC.fun(rowSums(mat[,A,drop=FALSE]),rowSums(mat[,B,drop=FALSE]),normalizeFun=function(i) i-med.element,...)
    } else {
      lfcs=LFC.fun(rowSums(mat[,A,drop=FALSE]),rowSums(mat[,B,drop=FALSE]),...)
    }
    if (is.data.frame(lfcs)) lfcs else data.frame(LFC=lfcs)
  })
}

#' Perform Wald tests for differential expression
#'
#' Apply DESeq2 for comparisons defined in a contrast matrix, requires the DESeq2 package.
#'
#' @param data the grandR object
#' @param name.prefix the prefix for the new analysis name; a dot and the column names of the contrast matrix are appended; can be NULL (then only the contrast matrix names are used)
#' @param contrasts contrast matrix that defines all pairwise comparisons, generated using \link{GetContrasts}
#' @param separate model overdispersion separately for all pairwise comparison (TRUE), or fit a single model per gene, and extract contrasts (FALSE)
#' @param mode compute LFCs for "total", "new", or "old" RNA
#' @param slot which slot to use (should be a count slot, not normalized values)
#' @param normalization normalize on "total", "new", or "old" (see details)
#' @param logFC compute and add the log2 fold change as well
#' @param verbose print status messages?
#'
#' @details DESeq2 by default performs size factor normalization. When computing differential expression of new RNA,
#' it might be sensible to normalize w.r.t. to total RNA, i.e. use the size factors computed from total RNA instead of computed from new RNA.
#' This can be accomplished by setting mode to "new", and normalization to "total"!
#'
#' @details Normalization can also be a mode.slot! Importantly, do not specify a slot containing normalized values, but specify a slot of unnormalized values
#' (which are used to compute the size factors for normalization!)
#'
#' @return a new grandR object including a new analysis table. The columns of the new analysis table are
#' \itemize{
#'  \item{"M"}{the base mean}
#'  \item{"S"}{the log2FoldChange divided by lfcSE}
#'  \item{"P"}{the Wald test P value}
#'  \item{"Q"}{same as P but Benjamini-Hochberg multiple testing corrected}
#'  \item{"LFC"}{the log2 fold change (only with the logFC parameter set to TRUE)}
#' }
#'
#' @seealso \link{LFC},\link{GetContrasts}
#' @export
#'
#' @examples
#' \donttest{
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c(Design$Condition,Design$dur.4sU,Design$Replicate))
#' sars <- subset(sars,Coldata(sars,Design$dur.4sU)==2)
#' sars<-PairwiseDESeq2(sars,mode="total",
#'                               contrasts=GetContrasts(sars,contrast=c("Condition","Mock")))
#' sars<-PairwiseDESeq2(sars,mode="new",normalization="total",
#'                               contrasts=GetContrasts(sars,contrast=c("Condition","Mock")))
#' head(GetAnalysisTable(sars,column="Q"))
#' }
#'
#' @concept diffexp
PairwiseDESeq2=function(data, name.prefix=mode, contrasts, separate=FALSE, mode="total",
                        slot="count",
                        normalization=NULL,
                        logFC=FALSE, verbose=FALSE) {
  checkPackages("DESeq2")

  contrasts = contrasts[,apply(contrasts,2,function(v) all(c(-1,1) %in% v)),drop=FALSE]
  if (ncol(contrasts)==0) stop("Contrasts do not define any comparison!")

  mode.slot=paste0(mode,".",slot)
  if (!is.null(normalization)) {
    if (!check.mode.slot(data,normalization)) normalization=paste0(normalization,".",slot)
  } else {
    normalization = mode.slot
  }
  if (!check.mode.slot(data,normalization)) stop("Invalid normalization!")

  if (!check.mode.slot(data,mode.slot)) stop("Invalid mode")

  if (separate) {
    ApplyContrasts(data,name.prefix=name.prefix,contrasts=contrasts,mode.slot=mode.slot,verbose=verbose,analysis="DESeq2.Wald",FUN=function(mat,A,B) {
      norm.mat=as.matrix(GetTable(data,type=normalization,ntr.na=FALSE))
      norm.mat=cbind(norm.mat[,A,drop=FALSE],norm.mat[,B,drop=FALSE])
      A=mat[,A,drop=FALSE]
      B=mat[,B,drop=FALSE]
      coldata=data.frame(comparison=c(rep("A",ncol(A)),rep("B",ncol(B))))

      dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnt(cbind(A,B)),
                                    colData = coldata,
                                    design= ~ comparison-1)
      DESeq2::sizeFactors(dds)<-DESeq2::estimateSizeFactorsForMatrix(norm.mat)
      l=DESeq2::results(DESeq2::DESeq(dds,quiet=!verbose))
      re=data.frame(
        M=l$baseMean,
        S=l$stat,
        P=l$pvalue,
        Q=l$padj
      )
      if (logFC) re$LFC=l$log2FoldChange
      rownames(re)<-rownames(data)
      re
    })
  }
  else {
    groups=list()
    find.or.add=function(c) {
      if (length(groups)>0) for (i in 1:length(groups)) {
        ex=groups[[i]]
        if (all(c==ex)) { # it's already there
          return(as.character(i))
        } else if (any(c&ex)) {
          stop("Illegal intersection of contrasts for joint estimation of variance!")
        }
      }
      groups<<-c(groups,list(c))
      return(as.character(length(groups)))
    }
    dds.contrasts=list()
    cond.vec=rep(NA,nrow(contrasts))
    for (c in contrasts) {
      A=find.or.add(c==1)
      B=find.or.add(c==-1)
      dds.contrasts=c(dds.contrasts,list(c(A,B)))
      cond.vec[c==1]=A
      cond.vec[c==-1]=B
    }
    mat=as.matrix(GetTable(data,type=mode.slot,ntr.na=FALSE))
    mat=mat[,!is.na(cond.vec)]
    norm.mat=as.matrix(GetTable(data,type=normalization,ntr.na=FALSE))
    norm.mat=norm.mat[,!is.na(cond.vec)]

    cond.vec=cond.vec[!is.na(cond.vec)]
    coldata=data.frame(comparisons=factor(cond.vec,levels=as.character(1:length(groups))))
    mode.slot=get.mode.slot(data,mode.slot)

    dds <- DESeq2::DESeqDataSetFromMatrix(countData =cnt(mat),
                                  colData = coldata,
                                  design= ~ comparisons-1)

    DESeq2::sizeFactors(dds)<-DESeq2::estimateSizeFactorsForMatrix(norm.mat)
    dds <- DESeq2::DESeq(dds,quiet=!verbose)

    for (i in 1:length(contrasts)) {
      n=names(contrasts)[i]
      l=DESeq2::results(dds,contrast=c("comparisons",dds.contrasts[[i]][1],dds.contrasts[[i]][2]))
      re.df=data.frame(
        M=l$baseMean,
        S=l$stat,
        P=l$pvalue,
        Q=l$padj
      )
      if (logFC) re.df$LFC=l$log2FoldChange
      rownames(re.df)<-rownames(data)
      data=AddAnalysis(data,name = if (is.null(name.prefix)) n else paste0(name.prefix,".",n),table = re.df)
    }
    return(data)
  }
}


#' Estimate regulation from snapshot experiments
#'
#' Compute the posterior log2 fold change distributions of RNA synthesis and degradation
#'
#' @param data the grandR object
#' @param name.prefix the prefix for the new analysis name; a dot and the column names of the contrast matrix are appended; can be NULL (then only the contrast matrix names are used)
#' @param contrasts contrast matrix that defines all pairwise comparisons, generated using \link{GetContrasts}
#' @param reference.columns a reference matrix usually generated by \link{FindReferences} to define reference samples for each sample; can be NULL if all conditions are at steady state (see details)
#' @param slot the data slot to take f0 and totals from
#' @param time.labeling the column in the \link{Coldata} table denoting the labeling duration, or the numeric labeling duration itself
#' @param time.experiment the column in the \link{Coldata} table denoting the experimental time point (can be NULL, see details)
#' @param ROPE.max.log2FC the region of practical equivalence is [-ROPE.max.log2FC,ROPE.max.log2FC] in log2 fold change space
#' @param sample.f0.in.ss whether or not to sample f0 under steady state conditions
#' @param N the sample size
#' @param N.max the maximal number of samples (necessary if old RNA > f0); if more are necessary, a warning is generated
#' @param CI.size A number between 0 and 1 representing the size of the credible interval
#' @param seed Seed for the random number generator
#' @param dispersion overdispersion parameter for each gene; if NULL this is estimated from data
#' @param sample.level Define how the NTR is sampled from the hierarchical Bayesian model (must be 0,1, or 2; see details)
#' @param correct.labeling Labeling times have to be unique; usually execution is aborted, if this is not the case; if this is set to true, the median labeling time is assumed
#' @param verbose Print status messages
#'
#' @details The kinetic parameters s and d are computed using \link{TransformSnapshot}. For that, the sample either must be in steady state
#' (this is the case if defined in the reference.columns matrix), or if the levels at an earlier time point are known from separate samples,
#' so called temporal reference samples. Thus, if s and d are estimated for a set of samples x_1,...,x_k (that must be from the same time point t),
#' we need to find (i) the corresponding temporal reference samples from time t0, and (ii) the time difference between t and t0.
#'
#' @details The temporal reference samples are identified by the reference.columns matrix. This is a square matrix of logicals, rows and columns correspond to all samples
#' and TRUE indicates that the row sample is a temporal reference of the columns sample. This time point is defined by \code{time.experiment}. If \code{time.experiment}
#' is NULL, then the labeling time of the A or B samples is used (e.g. useful if labeling was started concomitantly with the perturbation, and the steady state samples
#' are unperturbed samples).
#'
#' @details By default, the hierarchical Bayesian model is estimated. If sample.level = 0, the NTRs are sampled from a beta distribution
#' that approximates the mixture of betas from the replicate samples. If sample.level = 1, only the first level from the hierarchical model
#' is sampled (corresponding to the uncertainty of estimating the biological variability). If sample.level = 2, the first and second levels
#' are estimated (corresponding to the full hierarchical model).
#'
#' @details if N is set to 0, then no sampling from the posterior is performed, but the transformed MAP estimates are returned
#'
#' @return a new grandR object including a new analysis table. The columns of the new analysis table are
#' \itemize{
#'  \item{"s.A"}{the posterior mean synthesis rate for sample A from the comparison}
#'  \item{"s.B"}{the posterior mean synthesis rate for sample B from the comparison}
#'  \item{"HL.A"}{the posterior mean RNA half-life for sample A from the comparison}
#'  \item{"HL.B"}{the posterior mean RNA half-life for sample B from the comparison}
#'  \item{"s.log2FC"}{the posterior mean synthesis rate log2 fold change}
#'  \item{"s.cred.lower"}{the lower CI boundary of the synthesis rate log2 fold change}
#'  \item{"s.cred.upper"}{the upper CI boundary of the synthesis rate log2 fold change}
#'  \item{"s.ROPE"}{the signed ROPE probability (negative means downregulation) for the synthesis rate fold change}
#'  \item{"HL.log2FC"}{the posterior mean half-life log2 fold change}
#'  \item{"HL.cred.lower"}{the lower CI boundary of the half-life log2 fold change}
#'  \item{"HL.cred.upper"}{the upper CI boundary of the half-life log2 fold change}
#'  \item{"HL.ROPE"}{the signed ROPE probability (negative means downregulation) for the half-life fold change}
#' }
#'
#'
#' @seealso \link{FitKineticsGeneSnapshot},\link{FitKineticsSnapshot}
#' @export
#'
#' @examples
#' banp <- ReadGRAND(system.file("extdata", "BANP.tsv.gz", package = "grandR"),
#'           design=c("Cell","Experimental.time","Genotype",
#'                        Design$dur.4sU,Design$has.4sU,Design$Replicate))
#' contrasts <- GetContrasts(banp,contrast=c("Experimental.time.original","0h"),name.format="$A")
#' reference.columns <- FindReferences(banp,reference= Experimental.time==0)
#' banp <- EstimateRegulation(banp,"Regulation",
#'                              contrasts=contrasts,
#'                              reference.columns=reference.columns,
#'                              verbose=TRUE,
#'                              time.experiment = "Experimental.time",
#'                              N=0,               # don't sample in the example
#'                              dispersion=0.1)    # don't estimate dispersion in the example
#' head(GetAnalysisTable(banp))
#'
#' @concept diffexp
EstimateRegulation=function(data,name.prefix="Regulation",
                            contrasts,reference.columns = NULL,
                            slot=DefaultSlot(data),
                            time.labeling=Design$dur.4sU,
                            time.experiment=NULL,
                            ROPE.max.log2FC=0.25,
                            sample.f0.in.ss=TRUE,
                            N=10000,
                            N.max=N*10,
                            CI.size=0.95,
                            seed=1337,
                            dispersion=NULL,
                            sample.level=2,
                            correct.labeling=FALSE,
                            verbose=FALSE) {
  if (!check.slot(data,slot)) stop("Illegal slot definition!")
  if(!is.null(seed)) set.seed(seed)
  if (!is.null(reference.columns) && (!is.matrix(reference.columns) || nrow(reference.columns)!=ncol(data) || ncol(reference.columns)!=ncol(data))) stop("Illegal reference.columns parameter!")

  for (n in names(contrasts)) {
    if (verbose) cat(sprintf("Computing Regulation for %s...\n",n))
    A=contrasts[[n]]==1
    B=contrasts[[n]]==-1

    if (is.null(reference.columns)) {
      ss=list(A=A,B=B)
    } else {
      ss=list(A=apply(reference.columns[,Columns(data,A),drop=FALSE]==1,1,any),B=apply(reference.columns[,Columns(data,B),drop=FALSE]==1,1,any))
    }
    if (length(Columns(data,ss$A))==0 || length(Columns(data,ss$B))==0) stop("No reference columns found; check your reference.columns parameter!")
    dispersion.A = if (sum(ss$A)==1) rep(0.1,nrow(data)) else if (!is.null(dispersion)) rep(dispersion,length.out=nrow(data)) else estimate.dispersion(GetTable(data,type="count",columns = ss$A))
    dispersion.B = if (sum(ss$B)==1) rep(0.1,nrow(data)) else if (!is.null(dispersion)) rep(dispersion,length.out=nrow(data)) else estimate.dispersion(GetTable(data,type="count",columns = ss$B))
    if (verbose) {
      if (any(ss$A & A)) {
        cat(sprintf("Sampling from steady state for %s...\n",paste(colnames(data)[A],collapse = ",")))
      } else {
        cat(sprintf("Sampling from non-steady state for %s (reference: %s)...\n",paste(colnames(data)[A],collapse = ","),paste(colnames(data)[ss$A],collapse = ",")))
      }
      if (any(ss$B & B)) {
        cat(sprintf("Sampling from steady state for %s...\n",paste(colnames(data)[B],collapse = ",")))
      } else {
        cat(sprintf("Sampling from non-steady state for %s (reference: %s)...\n",paste(colnames(data)[B],collapse = ","),paste(colnames(data)[ss$B],collapse = ",")))
      }
    }

    # obtain prior from expression values
    get.beta.prior=function(columns,dispersion) {
      ex=rowMeans(GetTable(data,type = slot, columns = columns))
      ex=1/dispersion/(1/dispersion+ex)
      E.ex=mean(ex)
      V.ex=var(ex)
      c(
        shape1=(E.ex*(1-E.ex)/V.ex-1)*E.ex,
        shape2=(E.ex*(1-E.ex)/V.ex-1)*(1-E.ex)
        )
    }
    beta.prior.A=get.beta.prior(A,dispersion.A)
    beta.prior.B=get.beta.prior(B,dispersion.B)
    if (verbose) cat(sprintf("Beta prior for %s: a=%.3f, b=%.3f\n",paste(colnames(data)[A],collapse = ","),beta.prior.A[1],beta.prior.A[2]))
    if (verbose) cat(sprintf("Beta prior for %s: a=%.3f, b=%.3f\n",paste(colnames(data)[B],collapse = ","),beta.prior.B[1],beta.prior.B[2]))


    re=plapply(1:nrow(data),function(i) {
    #for (i in 1:nrow(data)) { print (i);
      fit.A=FitKineticsGeneSnapshot(data=data,gene=i,columns=A,dispersion=dispersion.A[i],reference.columns=reference.columns,slot=slot,time.labeling=time.labeling,time.experiment=time.experiment,sample.f0.in.ss=sample.f0.in.ss,sample.level=sample.level,beta.prior=beta.prior.A,return.samples=TRUE,N=N,N.max=N.max,CI.size=CI.size,correct.labeling=correct.labeling)
      fit.B=FitKineticsGeneSnapshot(data=data,gene=i,columns=B,dispersion=dispersion.B[i],reference.columns=reference.columns,slot=slot,time.labeling=time.labeling,time.experiment=time.experiment,sample.f0.in.ss=sample.f0.in.ss,sample.level=sample.level,beta.prior=beta.prior.B,return.samples=TRUE,N=N,N.max=N.max,CI.size=CI.size,correct.labeling=correct.labeling)
      samp.a=fit.A$samples
      samp.b=fit.B$samples

      N=min(nrow(samp.a),nrow(samp.b))
      if (N==0 || is.null(fit.A$samples)) return(c(
        s.A=unname(fit.A$s),
        s.B=unname(fit.B$s),
        HL.A=log(2)/unname(fit.A$d),
        HL.B=log(2)/unname(fit.A$d),
        s.log2FC=unname(log2(fit.A$s/fit.B$s)),
        s.cred.lower=-Inf,
        s.cred.upper=-Inf,
        s.ROPE=NA,
        HL.log2FC=unname(log2(fit.B$d/fit.A$d)),
        HL.cred.lower=-Inf,
        HL.cred.upper=Inf,
        HL.ROPE=NA
      ))
      samp.a=samp.a[1:N,,drop=FALSE]
      samp.b=samp.b[1:N,,drop=FALSE]

      savelfc=function(a,b) ifelse((is.infinite(a) & is.infinite(b)) | (a==0&b==0),0,log2(a/b))
      lfc.s=savelfc(samp.a[,'s'],samp.b[,'s'])
      lfc.HL=savelfc(samp.b[,'d'],samp.a[,'d'])

      sc=quantile(lfc.s,c(0.5-CI.size/2,0.5+CI.size/2))
      hc=quantile(lfc.HL,c(0.5-CI.size/2,0.5+CI.size/2))

      return(
        c(
        s.A=mean(samp.a[,'s']),
        s.B=mean(samp.b[,'s']),
        HL.A=log(2)/mean(samp.a[,'d']),
        HL.B=log(2)/mean(samp.b[,'d']),
        s.log2FC=mean(lfc.s),
        s.cred.lower=unname(sc[1]),
        s.cred.upper=unname(sc[2]),
        s.ROPE=ROPE.LFC(lfc.s,ROPE.max.log2FC),
        HL.log2FC=mean(lfc.HL),
        HL.cred.lower=unname(hc[1]),
        HL.cred.upper=unname(hc[2]),
        HL.ROPE=ROPE.LFC(lfc.HL,ROPE.max.log2FC)
        )
        )
    },seed=seed)

    re.df=as.data.frame(t(simplify2array(re)))
    rownames(re.df)=Genes(data)

    data=AddAnalysis(data,name = if (is.null(name.prefix)) n else paste0(name.prefix,".",n),table = re.df)
  }
  data
}

ROPE.LFC=function(lfc,ROPE.max.log2FC) {
  na=sum(is.na(lfc))
  lfc=lfc[!is.na(lfc)]

  l=sum(lfc< -ROPE.max.log2FC)
  r=sum(lfc> ROPE.max.log2FC)
  m=sum(lfc<=ROPE.max.log2FC & lfc>=-ROPE.max.log2FC)+na

  if (l>r) -l/(m+r+l) else r/(m+l+r)
}

beta.approximate.mixture=function(a,b) {
  mix.mean=mean(bmean(a,b)) # mean of mixture mode is (weighted) average of means
  mix.var=mean(bvar(a,b)+bmean(a,b)^2)-mix.mean^2

  list(a=(mix.mean*(1-mix.mean)/mix.var-1)*mix.mean,b=(mix.mean*(1-mix.mean)/mix.var-1)*(1-mix.mean))
}

bvar=function(a,b) a*b/(a+b)^2/(a+b+1)
bmean=function(a,b) a/(a+b)

balpha=function(mu,var) mu*(mu*(1-mu)/var-1)
bbeta=function(mu,var) (1-mu)*(mu*(1-mu)/var-1)

# the hyperprior for the variance of the prior beta is determined by var.prior.min.sd;
# if the data beta posteriors (with uniform prior) overlap, the variance of the prior
# is largely determined by this hyperprior; the parameter is the minimal sd to obtain
hierarchical.beta.posterior=function(a,b,
                  compute.marginal.likelihood=FALSE,
                  compute.grid=FALSE,
                  fak.below.max=1000,
                  res=100,
                  plot=FALSE,plot.prior=FALSE) {
  if (length(a)!=length(b)) stop("Unequal length of alpha and beta!")
 # if (length(a)<2) stop("<2 observations!")
  mu=a/(a+b)

   var.prior.min.sd=sd(a/(a+b))#sqrt(min(bvar(a,b)))
   if (is.na(var.prior.min.sd) || var.prior.min.sd==0) var.prior.min.sd=sqrt(min(bvar(a,b)))
  max.size=mean(mu*(1-mu))/var.prior.min.sd^2
  f=function(x,o=max.size,s=max.size/100) log(1/(1+exp((x-o)/s))/s/log1p(exp(o/s)))
  if (plot.prior) {
    x=seq(0,max.size*1.3,length.out=1000)
    plot(x,exp(f(x)),type='l',main=sd(a/(a+b)))
  }
  lprior=function(pa,pb) f(pa+pb)
  lmarg.posterior=function(pa,pb) {
   # print(c(pa,pb,bmean(pa,pb),sqrt(bvar(pa,pb)),lprior(pa,pb),sum(lbeta(pa+a,pb+b)-lbeta(pa,pb))))
    lprior(pa,pb)+sum(lbeta(pa+a,pb+b)-lbeta(pa,pb))
  }
  ltrans.marg.posterior=function(logitmu,logsize) {
    #print(c(logitmu,logsize))
    pa=unname(exp(logitmu+logsize)/(exp(logitmu)+1))
    pb=unname(exp(logsize)/(exp(logitmu)+1))
    lJ=(logitmu+2*logsize)-2*log1p(exp(logitmu))
    lmarg.posterior(pa,pb)+lJ
  }
  to.ab=function(logitmu,logsize) list(a=unname(exp(logitmu+logsize)/(exp(logitmu)+1)),
                                       b=unname(exp(logsize)/(exp(logitmu)+1)))

  # find MAP
  size.point=if (!is.null(var.prior.min.sd)) log(mean(mu*(1-mu))/var.prior.min.sd^2) else log(min(mean(mu)*(1-mean(mu))/var(mu)-1,min(a+b)))
  opt=optim(c(logitmu=log(mean(a/b)),logsize=size.point),function(v) ltrans.marg.posterior(v[1],v[2]),control = list(fnscale=-1))

  logMAP=opt$par
  re=c(to.ab(logMAP[1],logMAP[2]),MAP=opt$value)
  re$prior=lprior(re$a,re$b)

  if (compute.marginal.likelihood || compute.grid) {
    mu.low=uniroot(function(x) ltrans.marg.posterior(x,logMAP[2])-opt$value+log(fak.below.max),interval = c(logMAP[1]-4,logMAP[1]),extendInt = "upX")$root
    mu.high=uniroot(function(x) ltrans.marg.posterior(x,logMAP[2])-opt$value+log(fak.below.max),interval = c(logMAP[1],logMAP[1]+4),extendInt = "downX")$root

    size.low=uniroot(function(x) ltrans.marg.posterior(logMAP[1],x)-opt$value+log(fak.below.max),interval = c(logMAP[2]-4,logMAP[2]),extendInt = "upX")$root
    size.high=uniroot(function(x) ltrans.marg.posterior(logMAP[1],x)-opt$value+log(fak.below.max),interval = c(logMAP[2],logMAP[2]+4),extendInt = "downX")$root

    if (compute.marginal.likelihood) {
      checkPackages("cubature")
      ml=log(cubature::adaptIntegrate(function(v) exp(ltrans.marg.posterior(v[1],v[2])-opt$value),lowerLimit = c(mu.low,size.low),upperLimit = c(mu.high,size.high))$integral)+opt$value
      re$mar.loglik=ml
    }

    if (compute.grid) {
      mu.grid=seq(mu.low,mu.high,length.out=res)
      size.grid=seq(size.low,size.high,length.out=res)
      grid=sapply(mu.grid,function(mu) sapply(size.grid, function(size)  exp(ltrans.marg.posterior(mu,size)-opt$value) ))
      dimnames(grid)=list("log size"=size.grid,"logit mu"=mu.grid)
      re$grid=grid
      #image(t(grid))
      #contour(t(grid),levels = (seq(0.05,0.95,by=0.05)))
      marginal = colSums(grid)
      re$sample.param=function(N) {
        ilogitmu=sample.int(length(marginal),N,replace = TRUE,prob = marginal)
        ilogsize=sapply(1:N,function(i) sample.int(nrow(grid),1,replace = TRUE,prob = grid[,ilogitmu[i]]))
        logitmu=mu.grid[ilogitmu]+runif(N,-(mu.grid[2]-mu.grid[1])/2,(mu.grid[2]-mu.grid[1])/2)
        logsize=size.grid[ilogsize]+runif(N,-(size.grid[2]-size.grid[1])/2,(size.grid[2]-size.grid[1])/2)
        cbind(a=unname(exp(logitmu+logsize)/(exp(logitmu)+1)),
             b=unname(exp(logsize)/(exp(logitmu)+1)))
      }
      re$sample.mu=function(N) {
        ilogitmu=sample.int(length(marginal),N,replace = TRUE,prob = marginal)
        logitmu=mu.grid[ilogitmu]+runif(N,-(mu.grid[2]-mu.grid[1])/2,(mu.grid[2]-mu.grid[1])/2)
        exp(logitmu)/(exp(logitmu)+1)
      }
      re$sample=function(N) {
        para=re$sample.param(N)
        rbeta(N,shape1=para[,1],shape2=para[,2])
      }

    }
  }

  if (plot) {
    checkPackages("graphics")
    x=seq(min(qbeta(0.001,a,b)),max(qbeta(0.999,a,b)),length.out=1000)
    plot(x,pbeta(x,a[1],b[1]),type='l',xlim=range(x),ylim=c(0,1))
    for (i in 2:length(a)) graphics::lines(x,pbeta(x,a[i],b[i]))
    graphics::lines(x,pbeta(x,re$a,re$b),col='black',lwd=2)
    if (!is.null(re$sample.mu)) graphics::lines(ecdf(re$sample.mu(10000)),col='red')
    if (!is.null(re$sample)) graphics::lines(ecdf(re$sample(10000)),col='orange')
  }
  re
}



#' List available gene sets
#'
#' Helper function to return a table with all available gene sets for \link{AnalyzeGeneSets}.
#'
#' @details This is a convenience wrapper for \link[msigdbr]{msigdbr_collections}.
#'
#' @return the gene set table; use the values in the category and subcategory columns for the corresponding parameters of \link{AnalyzeGeneSets}
#'
#'
#'
#' @seealso \link{AnalyzeGeneSets}
#' @export
#'
#' @concept genesets
ListGeneSets=function() {
  descr=c(
    H="Hallmark gene sets  are coherently expressed signatures derived by aggregating many MSigDB gene sets to represent well-defined biological states or processes.",
    C1="Positional gene sets  for each human chromosome and cytogenetic band.",
    C2="Curated gene sets  from online pathway databases, publications in PubMed, and knowledge of domain experts.",
    C3="Regulatory target gene sets  based on gene target predictions for microRNA seed sequences and predicted transcription factor binding sites.",
    C4="Computational gene sets  defined by mining large collections of cancer-oriented microarray data.",
    C5="Ontology gene sets  consist of genes annotated by the same ontology term.",
    C6="Oncogenic signature gene sets  defined directly from microarray gene expression data from cancer gene perturbations.",
    C7="Immunologic signature gene sets  represent cell states and perturbations within the immune system.",
    C8="Cell type signature gene sets  curated from cluster markers identified in single-cell sequencing studies of human tissue."
  )
  checkPackages("msigdbr")
  re=setNames(as.data.frame(msigdbr::msigdbr_collections()),c("category","subcategory","n"))
  re$`category description` = descr[as.character(re$category)]
  re
}


#' Gene set analysis
#'
#' Perform gene-set enrichment and overrepresentation analysis (GSEA/ORA) for a specified
#' set of genes
#'
#' @param data the grandR object that contains the data to analyze
#' @param analysis the analysis to use, can be more than one and can be regexes (see details)
#' @param criteria an expression to define criteria for GSEA/ORA (see details)
#' @param species the species the genes belong to (eg "Homo sapiens"); can be NULL, then the species is inferred from gene ids (see details)
#' @param category the category defining gene sets (see \link{ListGeneSets})
#' @param subcategory the category defining gene sets (see \link{ListGeneSets})
#' @param verbose Print status messages
#' @param minSize The minimal size of a gene set to be considered
#' @param maxSize The maximal size of a gene set to be considered
#' @param process.genesets a function to process geneset names; can be NULL (see details)
#'
#' @details The analysis parameter (just like for \link{GetAnalysisTable} can be a regex (that will be matched
#' against all available analysis names). It can also be a vector (of regexes). Be careful with this, if
#' more than one table e.g. with column LFC ends up in here, only the first is used (if criteria=LFC).
#'
#' @details The criteria parameter can be used to define how analyses are performed. The criteria must be an expression
#' that either evaluates into a numeric or logical vector. In the first case, GSEA is performed, in the latter it is ORA.
#' The columns of the given analysis table(s) can be used to build this expression.
#'
#' @details If no species is given, a very simple automatic inference is done, which will only work when having human or mouse ENSEMBL identifiers as gene ids.
#'
#' @details The process.genesets parameters can be function that takes the character vector representing the names of all gene sets. The original names are replaced
#' by the return value of this function.
#'
#' @return the clusterprofile object representing the analysis results.
#'
#' @seealso \link[clusterProfiler]{GSEA},\link[clusterProfiler]{enricher},\link[msigdbr]{msigdbr}
#'
#' @examples
#' # See the differential-expression vignette!
#'
#' @export
#' @concept genesets
AnalyzeGeneSets=function(data, analysis=Analyses(data)[1], criteria=LFC,
                         species = NULL, category = NULL, subcategory = NULL,
                         verbose=TRUE, minSize=10, maxSize=500,
                         process.genesets=NULL) {

  checkPackages(c("clusterProfiler","msigdbr"))

  if (is.null(species)) {
    if (sum(grepl("ENSG0",Genes(data,use.symbols=FALSE)))>nrow(data)/2) species="Homo sapiens"
    if (sum(grepl("ENSMUSG0",Genes(data,use.symbols=FALSE)))>nrow(data)/2) species="Mus musculus"
  }
  if (is.null(species)) stop("Cannot recognize species! Specify one of msigdbr::msigdbr_species()$species_name")

  if (is.null(category)) {
    stop("You cannot perform gene set analyses on all available collections; use ListGeneSets() to see all available options for the category and subcategory parameter!")
  }

  if (verbose) cat(sprintf("Querying msigdb for %s (%s) in %s ...\n",category,subcategory,species))
  gs = unique(as.data.frame(msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)[,c('gs_name','ensembl_gene')]))
  if (verbose) cat(sprintf("Found %d gene sets!\n",length(unique(gs$gs_name))))

  if (!is.null(process.genesets)) {
    ugs=unique(gs$gs_name)
    map=setNames(process.genesets(ugs),ugs)
    gs$gs_name=map[as.character(gs$gs_name)]
  }

  genes=eval(substitute(GetSignificantGenes(data,analysis=analysis,criteria=criteria,as.table=TRUE,use.symbols=FALSE,gene.info=FALSE)),enclos = parent.frame()) # this is necessary to call the eval subs function!
  if (mode(genes[,1])=="numeric") {
    checkPackages("fgsea")

    if (verbose) cat("Performing GSEA (using fgsea)...\n")
    clusterProfiler::GSEA(setNames(genes[,1],rownames(genes)),TERM2GENE = gs,minGSSize=minSize,maxGSSize = maxSize)
  }
  else {
    if (verbose) cat(sprintf("Performing ORA (for n=%d/%d genes)...",sum(genes[,1]),nrow(data)))
    clusterProfiler::enricher(gene = rownames(genes), universe = Genes(data,use.symbols=FALSE), TERM2GENE = gs,minGSSize=minSize,maxGSSize = maxSize)
  }
}

#' Create a summarize matrix
#'
#' If this matrix is multiplied with a count table (e.g. obtained by \code{\link{GetTable}}),
#' either the average (average=TRUE) or the sum (average=FALSE) of all columns (samples or cells)
#' belonging to the same \code{\link{Condition}} is computed.
#'
#' @param x A grandR object or a named vector (the names indicate the sample names, the value the conditions to be summarized)
#' @param no4sU Use no4sU columns (TRUE) or not (FALSE)
#' @param columns which columns (i.e. samples or cells) to return (see details)
#' @param average matrix to compute the average (TRUE) or the sum (FALSE)
#' @param subset logical vector of which elements of the vector v to use (or NULL: use all)
#' @param ... further arguments to be passed to or from other methods.
#'
#' @return A matrix to be multiplied with a count table
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition="x").
#'
#' @details The method for grandR object simply calls the general method
#'
#' @seealso \link{GetTable}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#'
#' GetSummarizeMatrix(sars)
#' head(as.matrix(GetTable(sars)) %*% GetSummarizeMatrix(sars))   # average by matrix multiplication
#' head(GetTable(sars,summarize = TRUE))                          # shortcut, does the same
#'
#' # See the data-matrices-and-analysis-results vignette for more examples!
#'
#' @export
#'
#' @concept data
GetSummarizeMatrix <- function (x, ...) {
  UseMethod("GetSummarizeMatrix", x)
}
#' @rdname GetSummarizeMatrix
#' @export
GetSummarizeMatrix.grandR=function(x,no4sU=FALSE,columns=NULL,average=TRUE,...) {
  data=x
  if (is.null(Condition(data))) stop("Does not have conditions!")

  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)
  if (!no4sU) columns=setdiff(columns,Columns(data,Coldata(data,"no4sU")))

  GetSummarizeMatrix.default(setNames(Condition(data),colnames(data)),subset=columns,average=average)
}
#' @rdname GetSummarizeMatrix
#' @export
GetSummarizeMatrix.default=function(x,subset=NULL,average=TRUE,...) {
  v=x
	re=NULL
	for (e in unique(v)) re=cbind(re,ifelse(v==e,1,0))
	rownames(re)=names(v)
	colnames(re)=unique(v)
	if (!is.null(subset)) {
  	save=re[subset,]
  	re[,]=0
  	re[subset,]=save
	}
	re=re[,colSums(re)>0]
	if (average) re=t(t(re)/colSums(re))
	re
}

#' Create a contrast matrix
#'
#' Each column of a contrast matrix represents a pairwise comparison of all samples or cells of
#' a grandR object (or a column annotation table). Elements being 1 are contrasted vs. elements being -1
#' (and all 0 are irrelevant for this comparison).
#'
#' @param x A grandR object or a column annotation table
#' @param contrast A vector describing what should be contrasted
#' @param no4sU Use no4sU columns (TRUE) or not (FALSE)
#' @param columns logical vector of which columns (samples or cells) to use (or NULL: use all); for grandR objects, see details
#' @param group Split the samples or cells according to this column of the column annotation table (and adapt the of the output table)
#' @param name.format Format string for generating the column from the contrast vector (see details)
#' @param ... further arguments to be passed to or from other methods.
#'
#' @details To compare one specific factor level \emph{A} against another level \emph{B} in
#' a particular column \emph{COL} of the column annotation table, specify contrast=c("COL","A","B")
#'
#' @details To compare all levels against a specific level \emph{A} in
#' a particular column \emph{COL} of the column annotation table, specify contrast=c("COL","A")
#'
#' @details To perform all pairwise comparisons of all levels from
#' a particular column \emph{COL} of the column annotation table, specify contrast=c("COL")
#'
#' @details If the column \emph{COL} only has two levels, all three are equivalent.
#'
#' @details In all cases, if groups is not NULL, the columns annotation table is first split and contrasts are applied within all samples or cells
#' with the same \emph{group} factor level.
#'
#' @details The format string specifies the column name in the generated contrast matrix (which is used as the \emph{Analysis} name when calling
#' \code{\link{ApplyContrasts}}, \code{\link{LFC}}, \code{\link{PairwiseDESeq2}}, etc.). The keywords \emph{$GRP}, \emph{$COL}, \emph{$A} and \emph{$B} are substituted
#' by the respective elements of the contrast vector or the group this comparison refers to. By default, it is "$A vs $B" if group is NULL, and "$A vs $B.$GRP" otherwise.
#'
#' @details The method for grandR objects simply calls the general method
#'
#' @details For grandR objects, columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition="x").
#'
#' @return A data frame representig a contrast matrix to be used in \code{\link{ApplyContrasts}}, \code{\link{LFC}}, \code{\link{PairwiseDESeq2}}
#'
#' @seealso \code{\link{ApplyContrasts}}, \code{\link{LFC}}, \code{\link{PairwiseDESeq2}}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition","Time",Design$Replicate))
#'
#' GetContrasts(sars,contrast="Condition")
#' # Compare all Mock vs. all SARS
#' GetContrasts(sars,contrast=c("Condition","SARS","Mock"))
#' # This direction of the comparison is more reasonable
#' GetContrasts(sars,contrast=c("Condition","SARS","Mock"),group="Time")
#' # Compare SARS vs Mock per time point
#' GetContrasts(sars,contrast=c("Time.original","no4sU"), group="Condition",no4sU=TRUE,
#'                                                 name.format="$A vs $B ($GRP)")
#' # Compare each sample against the respective no4sU sample
#'
#' # See the differential-expression vignette for more examples!
#' @export
#'
#' @concept diffexp
GetContrasts <- function (x, ...) {
  UseMethod("GetContrasts", x)
}
#' @rdname GetContrasts
#' @export
GetContrasts.grandR=function(x,contrast="Condition",no4sU=FALSE,columns=NULL,group=NULL,name.format=NULL,...) {

  columns=substitute(columns)
  columns=if (is.null(columns)) Columns(x) else eval(columns,Coldata(x),parent.frame())
  columns=Columns(x,columns)
  if (!no4sU) columns=setdiff(columns,Columns(x,Coldata(x,"no4sU")))
  columns=Columns(x) %in% columns

  GetContrasts.default(x=Coldata(x),contrast=contrast,group=group,columns=columns,name.format=name.format)
}
#' @rdname GetContrasts
#' @export
GetContrasts.default=function(x,contrast,columns=NULL,group=NULL,name.format=NULL,...) {
  coldata=x
  if (!(length(contrast) %in% 1:3) || !contrast[1]%in%names(coldata) || (length(contrast)>1 && !all(contrast[2:length(contrast)] %in% coldata[,contrast[1]]))) stop("Illegal contrasts (either a name from design (all pairwise comparisons), a name and a reference level (all comparisons vs. the reference), or a name and two levels (exactly this comparison))")

  if (is.null(name.format)) {
    name.format = if (is.null(group)) "$A vs $B" else "$A vs $B.$GRP"
  }

  make.name=function(contrast,group) gsub("$GRP",group,gsub("$COL",contrast[1], gsub("$A",contrast[2], gsub("$B",contrast[3], name.format,fixed=TRUE),fixed=TRUE),fixed=TRUE),fixed=TRUE)
  make.col=function(contrast,group,use=TRUE) {
    re=rep(0,nrow(coldata))
    re[coldata[,contrast[1]]==contrast[2] & use]=1
    re[coldata[,contrast[1]]==contrast[3] & use]=-1
    if (!is.null(columns)) {
      save=re
      re=rep(0,nrow(coldata))
      re[columns]=save[columns]
    }
    setNames(data.frame(re,check.names=FALSE),make.name(contrast,group))
  }

  contr=if (length(contrast)==3) {
    function(use=TRUE,group="") make.col(contrast,group,use)
  } else if (length(contrast)==2) {
    function(use=TRUE,group="") {
      if (is.null(columns)) columns=TRUE
      ll=if (is.factor(coldata[,contrast[1]])) levels(droplevels(coldata[columns,contrast[1]])) else unique(coldata[columns,contrast[1]])
      ll=setdiff(ll,contrast[2])
      as.data.frame(lapply(ll,function(l) make.col(c(contrast[1],l,contrast[2]),group,use)),check.names=FALSE)
    }
  } else {
    function(use=TRUE,group="") {
      if (is.null(columns)) columns=TRUE
      ll=if (is.factor(coldata[,contrast[1]])) levels(droplevels(coldata[columns,contrast[1]])) else unique(coldata[columns,contrast[1]])
      if (length(ll)<2) stop("Less than 2 levels in contrast: ")
      re=combn(ll,2,FUN=function(v) make.col(c(contrast,v),group,use)[,1])
      colnames(re)=combn(ll,2,FUN=function(v) names(make.col(c(contrast,v),group,use))[1])
      as.data.frame(re,check.names=FALSE)
    }
  }
  re=if (is.null(group)) {
    contr()
  } else {
    covvec=interaction(coldata[group],drop=FALSE,sep=".")
    names=levels(covvec)
    as.data.frame(lapply(levels(covvec),function(bin) {
      r=contr(use=covvec==bin,group=bin)
    }),check.names=FALSE)
  }
  re=re[,!apply(re==0,2,all),drop=FALSE]
  rownames(re)=rownames(coldata)

  remove=apply(re>=0,2,all) | apply(re<=0,2,all)
  if (sum(remove)>0) {
    warning(sprintf("Removed columns without matching experiment: %s",paste(colnames(re)[remove],collapse = ",")))
    re=re[,!remove]
  }
  re
}

# Normalization of NTRs such that: median logFC new RNA vs. new RNA is 0, there is no correlation of this logFC vs the NTR
NormalizeEffectiveLabeling=function(data,reference=colnames(data),slot="norm",verbose=FALSE) {
  checkPackages("quantreg")
  w=rowMeans(GetTable(data,type=slot,columns=reference),na.rm=TRUE)
  ntr=rowMeans(GetTable(data,type="ntr",columns=reference),na.rm=TRUE)
  w=w*ntr
  use=!is.na(w) && w>0
  w=w[use]

  for (s in colnames(data)[!data$coldata$no4sU])  {
    if (verbose) cat(sprintf("Fitting model for %s...\n",s))
    d.ntr=data$data$ntr[use,s]
    d=GetTable(data,type=slot,columns=s)[use,1]
    df=data.frame(ntr=ntr,lfc=log2(d.ntr*d/w))
    df=df[!is.na(df$lfc) & !is.infinite(df$lfc),]
    fit=quantreg::lprq(df$ntr,df$lfc,h = 0.05)
    fn=splinefun(fit$xx,fit$fv)
    data$data$ntr[,s]=pmax(pmin(1,data$data$ntr[,s]/2^fn(ntr)),0)
  }

  data
}

