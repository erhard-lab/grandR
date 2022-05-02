
#' Compute a likelihood ratio test.
#'
#' The test is computed on any of total/old/new counts using DESeq2 based on two models
#' specified using formulas.
#'
#' @param data A grandR object
#' @param name the user defined analysis name to store the results
#' @param mode either "total", "new" or "old"
#' @param target formula specifying the target model (you can use any column name from the \code{\link{Coldata}(data)))
#' @param background formula specifying the background model (you can use any column name from the \code{\link{Coldata}(data)))
#' @param no4sU Use no4sU columns (TRUE) or not (FALSE)
#' @param columns logical vector of which columns (samples or cells) to use (or NULL: use all)
#' @param verbose Print status updates
#'
#' @return a new grandR object including a new analysis table
#'
#' @export
#'
LikelihoodRatioTest=function(data,name="LRT",mode="total",normalization=mode,target=~Condition,background=~1,no4sU=FALSE,columns=NULL,verbose=FALSE) {
  mode.slot=paste0(mode,".count")
  normalization=paste0(normalization,".count")
  if (!check.mode.slot(data,mode.slot)) stop("Invalid mode")

  columns=if (is.null(columns)) no4sU | !Coldata(data)$no4sU else columns&(no4sU | !Coldata(data)$no4sU)

	colData=droplevels(Coldata(data)[columns,])
	mat=GetTable(data,type=mode.slot,columns = columns)


	norm.mat=as.matrix(GetTable(data,type=normalization))

	colData=as.data.frame(lapply(colData,function(c) {
	  if (is.factor(c)) levels(c)=make.names(levels(c),unique = TRUE)
	  c
	}))


	dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnt(mat),colData=colData,design = target)
	DESeq2::sizeFactors(dds)<-DESeq2::estimateSizeFactorsForMatrix(norm.mat)
	dds.tot <- DESeq2::DESeq(dds,test="LRT", reduced=background,quiet=!verbose)
	res.tot <- DESeq2::results(dds.tot)

	df=data.frame(
    M=res.tot$baseMean,
    S=res.tot$stat,
    P=res.tot$pvalue,
    Q=p.adjust(res.tot$pvalue,method="BH")
  )

  AddAnalysis(data,name = paste0(name),table = df)
}



#' Helper function to apply apply analysis methods to pairs of conditions
#'
#' @param data
#' @param analysis
#' @param name
#' @param contrasts
#' @param mode.slot
#' @param FUN
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ApplyContrasts=function(data,analysis,name,contrasts,mode.slot="count",verbose=FALSE,FUN,...) {
  mat=as.matrix(GetTable(data,type=mode.slot,ntr.na=FALSE))
  mode.slot=get.mode.slot(data,mode.slot)
  for (n in names(contrasts)) {
    if (verbose) cat(sprintf("Computing %s for %s...\n",analysis,n))
    re.df=FUN(mat,contrasts[[n]]==1,contrasts[[n]]==-1,...)
    data=AddAnalysis(data,name = paste0(name,".",n),table = re.df)
  }
  data
}

LFC=function(data,name=mode,contrasts,LFC.fun=PsiLFC,mode="total",slot="count",normalization=NULL,verbose=FALSE,...) {
  mode.slot=paste0(mode,".",slot)
  if (!is.null(normalization))normalization=paste0(normalization,".",slot)
  if (!check.mode.slot(data,mode.slot)) stop("Invalid mode")
  ApplyContrasts(data,name=name,contrasts=contrasts,mode.slot=mode.slot,verbose=verbose,analysis="LFC",FUN=function(mat,A,B) {
    if (!is.null(normalization)) {
      norm.mat=as.matrix(GetTable(data,type=normalization))
      nlfcs=LFC.fun(rowSums(norm.mat[,A,drop=FALSE]),rowSums(norm.mat[,B,drop=FALSE]),normalizeFun=function(i) i)
      med.element=median(nlfcs)
      lfcs=LFC.fun(rowSums(mat[,A,drop=FALSE]),rowSums(mat[,B,drop=FALSE]),normalizeFun=function(i) i-med.element,...)
    } else {
      lfcs=LFC.fun(rowSums(mat[,A,drop=FALSE]),rowSums(mat[,B,drop=FALSE]),...)
    }
    if (is.data.frame(lfcs)) lfcs else data.frame(LFC=lfcs)
  })
}

PairwiseDESeq2=function(data,name=mode,contrasts,separate=FALSE,mode="total",normalization=mode,logFC=FALSE,verbose=FALSE) {
  mode.slot=paste0(mode,".count")
  normalization=paste0(normalization,".count")
  if (!check.mode.slot(data,mode.slot)) stop("Invalid mode")

  if (separate) {
    ApplyContrasts(data,name=name,contrasts=contrasts,mode.slot=mode.slot,verbose=verbose,analysis="DESeq2.Wald",FUN=function(mat,A,B) {
      norm.mat=as.matrix(GetTable(data,type=normalization))
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
    norm.mat=as.matrix(GetTable(data,type=normalization))
    norm.mat=norm.mat[,!is.na(cond.vec)]

    cond.vec=cond.vec[!is.na(cond.vec)]
    coldata=data.frame(comparisons=cond.vec)
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
      data=AddAnalysis(data,name = paste0(name,".",n),table = re.df)
    }
    return(data)
  }
}


#' Compute the posterior logFC distributions of RNA synthesis and degradation
#'
#' @param data the grandR object
#' @param name the name of the analysis added to the grandR object
#' @param contrasts A contrast matrix defining pairwise comparisons
#' @param reference.columns either a reference matrix (\code{\link{FindReferences}}) to define reference samples for each sample in A and B, or a
#' name list with names A and B containing logical vectors denoting the reference samples
#' @param slot the data slot to take f0 and totals from
#' @param time.labeling the column in the column annotation table denoting the labeling duration or the labeling duration itself
#' @param time.experiment the column in the column annotation table denoting the experimental time point (can be NULL, see details)
#' @param ROPE.max.log2FC the region of practical equivalence is [-ROPE.max.log2FC,ROPE.max.log2FC] in log2 fold change space
#' @param sample.f0.in.ss whether or not to sample f0 under steady state conditions
#' @param N the sample size
#' @param N.max the maximal number of samples (necessary if old RNA > f0); if more are necessary, a warning is generated
#' @param conf.int A number between 0 and 1 representing the size of the credible interval
#' @param seed Seed for the random number generator
#' @param hierarchical Take the NTR from the hierarchical bayesian model (see details)
#' @param verbose Vebose output
#'
#' @details The kinetic parameters s and d are computed using \link{TransformSnapshot}. For that, the sample either must be in steady state
#' (this is the case if defined in the reference.columns matrix), or if the levels at a specific time point are known. This time point is
#' defined by \code{time.experiment} (i.e. the difference between the steady state samples and the A or B samples themselves). If
#' \code{time.experiment} is NULL, then the labeling time of the A or B samples is used (e.g. usefull if labeling was started concomitantly with
#' the perturbation, and the steady state samples are unperturbed samples).
#'
#' @details By default, the hierarchical bayesian model is estimated. If hierarchical = FALSE, the NTRs are sampled from a beta distribution
#' that approximates the mixture of betas from the replicate samples. (see \link{beta.approximate.mixture}).
#'
#' @return a new grandR object containing an additional analysis
#'
#'
#' @export
#'
EstimateRegulation=function(data,name,contrasts,reference.columns,slot=DefaultSlot(data),time.labeling=Design$dur.4sU,time.experiment=NULL, ROPE.max.log2FC=0.25,sample.f0.in.ss=TRUE,N=10000,N.max=N*10,conf.int=0.95,seed=1337, hierarchical=TRUE, verbose=FALSE) {
  if (!check.slot(data,slot)) stop("Illegal slot definition!")
  if(!is.null(seed)) set.seed(seed)

  for (n in names(contrasts)) {
    if (verbose) cat(sprintf("Computing Regulation for %s...\n",n))
    A=contrasts[[n]]==1
    B=contrasts[[n]]==-1

    ss=if (is.matrix(reference.columns)) list(A=apply(reference.columns[,Columns(data,A)]==1,1,any),B=apply(reference.columns[,Columns(data,B)]==1,1,any))
    dispersion.A = if (sum(ss$A)==1) rep(0.1,nrow(data)) else estimate.dispersion(GetTable(data,type="count",columns = ss$A))
    dispersion.B = if (sum(ss$B)==1) rep(0.1,nrow(data)) else estimate.dispersion(GetTable(data,type="count",columns = ss$B))

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

    re=plapply(1:nrow(data),function(i) {
    #for (i in 1:nrow(data)) { print (i);
      fit.A=FitKineticsSnapshot(data=data,gene=i,columns=A,dispersion=dispersion.A[i],reference.columns=reference.columns,slot=slot,time.labeling=time.labeling,time.experiment=time.experiment,sample.f0.in.ss=sample.f0.in.ss,hierarchical=hierarchical,return.samples=TRUE,N=N,N.max=N.max,conf.int=conf.int)
      fit.B=FitKineticsSnapshot(data=data,gene=i,columns=B,dispersion=dispersion.B[i],reference.columns=reference.columns,slot=slot,time.labeling=time.labeling,time.experiment=time.experiment,sample.f0.in.ss=sample.f0.in.ss,hierarchical=hierarchical,return.samples=TRUE,N=N,N.max=N.max,conf.int=conf.int)
      samp.a=fit.A$samples
      samp.b=fit.B$samples

      N=min(nrow(samp.a),nrow(samp.b))
      if (N==0) return(c(
        s.A=NA,
        s.B=NA,
        HL.A=NA,
        HL.B=NA,
        s.log2FC=NA,
        s.cred.lower=-Inf,
        s.cred.upper=-Inf,
        s.ROPE=NA,
        HL.log2FC=NA,
        HL.cred.lower=-Inf,
        HL.cred.upper=Inf,
        HL.ROPE=NA
      ))
      samp.a=samp.a[1:N,,drop=FALSE]
      samp.b=samp.b[1:N,,drop=FALSE]

      savelfc=function(a,b) ifelse((is.infinite(a) & is.infinite(b)) | (a==0&b==0),0,log2(a/b))
      lfc.s=savelfc(samp.a[,'s'],samp.b[,'s'])
      lfc.HL=savelfc(samp.b[,'d'],samp.a[,'d'])

      sc=quantile(lfc.s,c(0.5-conf.int/2,0.5+conf.int/2))
      hc=quantile(lfc.HL,c(0.5-conf.int/2,0.5+conf.int/2))

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
        HL.conf.upper=unname(hc[1]),
        HL.conf.lower=unname(hc[2]),
        HL.ROPE=ROPE.LFC(lfc.HL,ROPE.max.log2FC)
        )
        )
    },seed=seed)

    re.df=as.data.frame(t(simplify2array(re)))

    data=AddAnalysis(data,name = paste0(name,".",n),table = re.df)
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
                  plot=FALSE) {
  if (length(a)!=length(b)) stop("Unequal length of alpha and beta!")
 # if (length(a)<2) stop("<2 observations!")
  mu=a/(a+b)

   var.prior.min.sd=sd(a/(a+b))#sqrt(min(bvar(a,b)))
   if (is.na(var.prior.min.sd) || var.prior.min.sd==0) var.prior.min.sd=sqrt(min(bvar(a,b)))
  max.size=mean(mu*(1-mu))/var.prior.min.sd^2
  f=function(x,o=max.size,s=max.size/100) log(1/(1+exp((x-o)/s))/s/log1p(exp(o/s)))
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
      re$sample=function(N) {
        ilogitmu=sample.int(length(marginal),N,replace = TRUE,prob = marginal)
        ilogsize=sapply(1:N,function(i) sample.int(nrow(grid),1,replace = TRUE,prob = grid[,i]))
        logitmu=mu.grid[ilogitmu]+runif(N,-(mu.grid[2]-mu.grid[1])/2,(mu.grid[2]-mu.grid[1])/2)
        logsize=size.grid[ilogsize]+runif(N,-(size.grid[2]-size.grid[1])/2,(size.grid[2]-size.grid[1])/2)
        to.ab(logitmu,logsize)
      }
      re$sample.mu=function(N) {
        ilogitmu=sample.int(length(marginal),N,replace = TRUE,prob = marginal)
        logitmu=mu.grid[ilogitmu]+runif(N,-(mu.grid[2]-mu.grid[1])/2,(mu.grid[2]-mu.grid[1])/2)
        exp(logitmu)/(exp(logitmu)+1)
      }
    }
  }

  if (plot) {
    x=seq(min(qbeta(0.001,a,b)),max(qbeta(0.999,a,b)),length.out=1000)
    plot(x,pbeta(x,a[1],b[1]),type='l',xlim=range(x),ylim=c(0,1))
    for (i in 2:length(a)) lines(x,pbeta(x,a[i],b[i]))
    lines(x,pbeta(x,re$a,re$b),col='black',lwd=2)
    if (!is.null(re$sample.mu)) lines(ecdf(re$sample.mu(10000)),col='red')
  }
  re
}


AnalyzeGeneSets=function(data,analysis=Analyses(data)[1],criteria=LFC,species = NULL, category = NULL, subcategory = NULL, verbose=TRUE, minSize=10, maxSize=500, process.genesets=NULL) {
  if (is.null(species)) {
    if (sum(grepl("ENSG0",Genes(data,use.symbols=FALSE)))>nrow(data)/2) species="Homo sapiens"
    if (sum(grepl("ENSMUSG0",Genes(data,use.symbols=FALSE)))>nrow(data)/2) species="Mus musculus"
  }
  if (is.null(species)) stop("Cannot recognize species! Specify one of msigdbr::msigdbr_species()$species_name")

  if (verbose) cat(sprintf("Querying msigdb for %s (%s) in %s ...\n",category,subcategory,species))
  gs = unique(as.data.frame(msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)[,c('gs_name','ensembl_gene')]))
  if (verbose) cat(sprintf("Found %d gene sets!\n",length(unique(gs$gs_name))))

  if (!is.null(process.genesets)) {
    ugs=unique(gs$gs_name)
    map=setNames(process.genesets(ugs),ugs)
    gs$gs_name=map[as.character(gs$gs_name)]
  }

  genes=eval(substitute(GetSignificantGenes(data,analysis=analysis,criteria=criteria,return.values=TRUE,use.symbols=FALSE)),enclos = parent.frame()) # this is necessary to call the eval subs function!

  if (mode(genes)=="numeric") {
    if (verbose) cat("Performing GSEA (using fgsea)...\n")
    clusterProfiler::GSEA(genes,TERM2GENE = gs,minGSSize=minSize,maxGSSize = maxSize)
  }
  else {
    if (verbose) cat(sprintf("Performing ORA (for n=%d/%d genes)...",length(genes),nrow(d)))
    clusterProfiler::enricher(gene = genes, universe = Genes(data,use.symbols=FALSE), TERM2GENE = gs,minGSSize=minSize,maxGSSize = maxSize)
  }
}



#' Create a summarize matrix
#'
#' If this matrix is multiplied with a count table (e.g. obtained by \code{\link{GetTable}}),
#' either the average (average=TRUE) or the sum (average=FALSE) of all columns (samples or cells)
#' belonging to the same \code{\link{Condition}} is computed.
#'
#' @param data A grandR object
#' @param no4sU Use no4sU columns (TRUE) or not (FALSE)
#' @param columns logical vector of which columns (samples or cells) to use (or NULL: use all)
#' @param average matrix to compute the average (TRUE) or the sum (FALSE)
#' @param v a named vector (the names indicate the sample names, the value the conditions to be summarized)
#' @param subset logical vector of which elements of the vector v to use (or NULL: use all)
#'
#' @return A matrix to be multiplied with a count table
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
#' GetSummarizeMatrix(c(A=1,B=1,C=2,D=3))
#'
#' @export
#'
GetSummarizeMatrix <- function (x, ...) {
  UseMethod("GetSummarizeMatrix", x)
}
#' @rdname GetSummarizeMatrix
#' @export
GetSummarizeMatrix.grandR=function(data,no4sU=FALSE,columns=NULL,average=TRUE) {
  if (is.null(Condition(data))) stop("Does not have conditions!")
  columns=if (is.null(columns)) no4sU | !Coldata(data)$no4sU else columns&(no4sU | !Coldata(data)$no4sU)
  GetSummarizeMatrix.default(setNames(Condition(data),colnames(data)),subset=columns,average=average)
}
#' @rdname GetSummarizeMatrix
#' @export
GetSummarizeMatrix.default=function(v,subset=NULL,average=TRUE) {
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
#' @param data A grandR object
#' @param coldata A column annotation table
#' @param contrast A vector describing what should be contrasted
#' @param no4sU Use no4sU columns (TRUE) or not (FALSE)
#' @param columns logical vector of which columns (samples or cells) to use (or NULL: use all)
#' @param group Split the samples or cells according to this column of the column annotation table (and adapt the of the output table)
#' @param name.format Format string for generating the column from the contrast vector (see details)
#'
#' @return A contrast matrix to be used in \code{\link{ApplyContrasts}}, \code{\link{LFC}}, \code{\link{TestPairwise}}
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
#' \code{\link{ApplyContrasts}}, \code{\link{LFC}}, \code{\link{TestPairwise}}, etc.). The keywords \emph{$COL}, \emph{$A} and \emph{$B} are substituted
#' by the respective elements of the contrast vector.
#'
#' @details The method for grandR objects simply calls the general method
#'
#' @seealso \code{\link{ApplyContrasts}}, \code{\link{LFC}}, \code{\link{TestPairwise}}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition","Time",Design$Replicate))
#'
#' GetContrasts(sars,contrast="Condition")                    # Compare all Mock vs. all SARS
#' GetContrasts(sars,contrast=c("Condition","SARS","Mock"))   # This direction of the comparison is more reasonable
#' GetContrasts(sars,contrast=c("Condition","SARS","Mock"),group="Time")   # Compare SARS vs Mock per time point
#' GetContrasts(sars,contrast=c("Time","no4sU"), group="Condition",name.format="$COL ($A vs $B)")   # Compare each sample against the respective no4sU sample
#'
#' @export
#'
GetContrasts <- function (x, ...) {
  UseMethod("GetContrasts", x)
}
#' @rdname GetContrasts
#' @export
GetContrasts.grandR=function(data,contrast="Condition",no4sU=FALSE,columns=NULL,group=NULL,name.format="$A vs $B") {
  columns=if (is.null(columns)) no4sU | !Coldata(data)$no4sU else columns&(no4sU | !Coldata(data)$no4sU)
  GetContrasts.default(coldata=Coldata(data),contrast=contrast,group=group,columns=columns,name.format=name.format)
}
#' @rdname GetContrasts
#' @export
GetContrasts.default=function(coldata,contrast,columns=NULL,group=NULL,name.format="$A vs $B") {
  if (!(length(contrast) %in% 1:3) || !contrast[1]%in%names(coldata) || (length(contrast)>1 && !all(contrast[2:length(contrast)] %in% coldata[,contrast[1]]))) stop("Illegal contrasts (either a name from design (all pairwise comparisons), a name and a reference level (all comparisons vs. the reference), or a name and two levels (exactly this comparison))")

  make.name=function(contrast) gsub("$COL",contrast[1], gsub("$A",contrast[2], gsub("$B",contrast[3], name.format,fixed=TRUE),fixed=TRUE),fixed=TRUE)
  make.col=function(contrast,use=TRUE) {
    re=rep(0,nrow(coldata))
    re[coldata[,contrast[1]]==contrast[2] & use]=1
    re[coldata[,contrast[1]]==contrast[3] & use]=-1
    if (!is.null(columns)) {
      save=re
      re=rep(0,nrow(coldata))
      re[columns]=save[columns]
    }
    setNames(data.frame(re,check.names=FALSE),make.name(contrast))
  }

  contr=if (length(contrast)==3) {
    function(use=TRUE) make.col(contrast,use)
  } else if (length(contrast)==2) {
    function(use=TRUE) {
      if (is.null(columns)) columns=TRUE
      ll=if (is.factor(coldata[,contrast[1]])) levels(droplevels(coldata[columns,contrast[1]])) else unique(coldata[columns,contrast[1]])
      ll=setdiff(ll,contrast[2])
      as.data.frame(lapply(ll,function(l) make.col(c(contrast[1],l,contrast[2]),use)),check.names=FALSE)
    }
  } else {
    function(use=TRUE) {
      if (is.null(columns)) columns=TRUE
      ll=if (is.factor(coldata[,contrast[1]])) levels(droplevels(coldata[columns,contrast[1]])) else unique(coldata[columns,contrast[1]])
      if (length(ll)<2) stop("Less than 2 levels in contrast: ")
      re=combn(ll,2,FUN=function(v) make.col(c(contrast,v),use)[,1])
      colnames(re)=combn(ll,2,FUN=function(v) names(make.col(c(contrast,v),use))[1])
      as.data.frame(re,check.names=FALSE)
    }
  }
  re=if (is.null(group)) {
    contr()
  } else {
    covvec=interaction(coldata[group],drop=FALSE,sep=".")
    names=levels(covvec)
    as.data.frame(lapply(levels(covvec),function(bin) {
      r=contr(use=covvec==bin)
      setNames(r,paste0(names(r),".",bin))
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
