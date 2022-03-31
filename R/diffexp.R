
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
LikelihoodRatioTest=function(data,name="LRT",mode="total",target=~Condition,background=~1,no4sU=FALSE,columns=NULL,verbose=FALSE) {
  mode.slot=paste0(mode,".count")
  if (!check.mode.slot(data,mode.slot)) stop("Invalid mode")

  columns=if (is.null(columns)) no4sU | !Coldata(data)$no4sU else columns&(no4sU | !Coldata(data)$no4sU)

	colData=droplevels(Coldata(data)[columns,])
	mat=GetTable(data,type=mode.slot,columns = columns)

	colData=as.data.frame(lapply(colData,function(c) {
	  if (is.factor(c)) levels(c)=make.names(levels(c),unique = TRUE)
	  c
	}))

	dds.tot <- DESeq2::DESeq(DESeq2::DESeqDataSetFromMatrix(countData = cnt(mat),colData=colData,design = target),test="LRT", reduced=background,quiet=!verbose)
	res.tot <- DESeq2::results(dds.tot)

	df=data.frame(
    M=res.tot$baseMean,
    S=res.tot$stat,
    P=res.tot$pvalue,
    Q=p.adjust(res.tot$pvalue,method="BH")
  )

  AddAnalysis(data,description=MakeAnalysis(name = paste0(name),analysis = "DESeq2.LRT",mode = mode.slot$mode,slot=mode.slot$slot,columns = colnames(mat)),table = df)
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
    data=AddAnalysis(data,description=MakeAnalysis(name = paste0(name,".",n),analysis = analysis,mode = mode.slot$mode,slot=mode.slot$slot,columns = colnames(mat)[contrasts[[n]]==1|contrasts[[n]]==-1]),table = re.df)
  }
  data
}

LFC=function(data,name=mode,contrasts,LFC.fun=PsiLFC,mode="total",slot="count",verbose=FALSE,...) {
  mode.slot=paste0(mode,".",slot)
  if (!check.mode.slot(data,mode.slot)) stop("Invalid mode")
  ApplyContrasts(data,name=name,contrasts=contrasts,mode.slot=mode.slot,verbose=verbose,analysis="LFC",FUN=function(mat,A,B) {
    lfcs=LFC.fun(rowSums(mat[,A,drop=FALSE]),rowSums(mat[,B,drop=FALSE]),...)
    if (is.data.frame(lfcs)) lfcs else data.frame(LFC=lfcs)
  })
}

PairwiseDESeq2=function(data,name=mode,contrasts,separate=FALSE,mode="total",verbose=FALSE) {
  mode.slot=paste0(mode,".count")
  if (!check.mode.slot(data,mode.slot)) stop("Invalid mode")

  if (separate) {
    ApplyContrasts(data,name=name,contrasts=contrasts,mode.slot=mode.slot,verbose=verbose,analysis="DESeq2.Wald",FUN=function(mat,A,B) {
      A=mat[,A,drop=FALSE]
      B=mat[,B,drop=FALSE]
      coldata=data.frame(comparison=c(rep("A",ncol(A)),rep("B",ncol(B))))
      dds <- DESeq2::DESeqDataSetFromMatrix(countData = cnt(cbind(A,B)),
                                    colData = coldata,
                                    design= ~ comparison-1)
      l=DESeq2::results(DESeq2::DESeq(dds,quiet=!verbose))
      data.frame(
        M=l$baseMean,
        S=l$stat,
        P=l$pvalue,
        Q=l$padj
      )
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
    cond.vec=cond.vec[!is.na(cond.vec)]
    coldata=data.frame(comparisons=cond.vec)
    mode.slot=get.mode.slot(data,mode.slot)

    dds <- DESeq2::DESeqDataSetFromMatrix(countData =cnt(mat),
                                  colData = coldata,
                                  design= ~ comparisons-1)
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
      data=AddAnalysis(data,description=MakeAnalysis(name = paste0(name,".",n),analysis = "DESeq2.Wald",mode = mode.slot$mode,slot=mode.slot$slot,columns = colnames(mat)[contrasts[[n]]==1|contrasts[[n]]==-1]),table = re.df)
    }
    return(data)
  }
}


PairwiseNtrTest=function(data,name,contrasts,verbose=FALSE) {

  alpha=as.matrix(GetTable(data,type="alpha",ntr.na=TRUE))
  beta=as.matrix(GetTable(data,type="beta",ntr.na=TRUE))
  bvar=function(a,b) a*b/(a+b)^2/(a+b+1)

  for (n in names(contrasts)) {
    if (verbose) cat(sprintf("Computing NTR test for %s...\n",n))

    A=contrasts[[n]]==1
    B=contrasts[[n]]==-1

    go=function(i,plot1=F,plot2=F) {
      alpha.A=alpha[i,A][alpha[i,A]>0&beta[i,A]>0]
      beta.A=beta[i,A][alpha[i,A]>0&beta[i,A]>0]
      alpha.B=alpha[i,B][alpha[i,B]>0&beta[i,B]>0]
      beta.B=beta[i,B][alpha[i,B]>0&beta[i,B]>0]

      bayesian.beta.test(alpha.A,beta.A,alpha.B,beta.B,plot1 = plot1,plot2 = plot2) #,plot=c(0.05,0.35))
    }

    re=psapply(1:nrow(data),go)

    re.df=data.frame(P=re,Q=p.adjust(re,method="BH"))

    data=AddAnalysis(data,description=MakeAnalysis(name = paste0(name,".",n),analysis = "NTR-test",slot="ntr",columns = colnames(data)[contrasts[[n]]==1|contrasts[[n]]==-1]),table = re.df)
  }
  data
}

compute.posterior.quantiles=function(a,b,sd) {
  mu=a/(a+b)

  eps=1E-3
  rr=range(mu)
  e.rr=c(0,1)
  if (0.25>=sd^2) e.rr=c(0.5-sqrt(0.25-sd^2)+eps,0.5+sqrt(0.25-sd^2)-eps)
  rr[1]=max(rr[1],e.rr[1])
  rr[2]=min(rr[2],e.rr[2])

  optfun=function(mu) {
    pa=mu*(mu*(1-mu)/sd^2-1)
    pb=(1-mu)*(mu*(1-mu)/sd^2-1)
    qs=sapply(seq_along(a),function(i) integrate(function(x) dbeta(x,pa,pb)*pbeta(x,a[i],b[i]),lower=qbeta(0.001,pa,pb),upper=qbeta(0.999,pa,pb))$value)
    mean(qs)-0.5
  }
  if (optfun(rr[2])<0) rr=c(rr[2],e.rr[2]) else if (optfun(rr[1])>0) rr=c(e.rr[1],rr[1])
  opt=uniroot(optfun,interval = rr)

  mu=opt$root
  pa=mu*(mu*(1-mu)/sd^2-1)
  pb=(1-mu)*(mu*(1-mu)/sd^2-1)

  sapply(seq_along(a),function(i) integrate(function(x) dbeta(x,pa,pb)*pbeta(x,a[i],b[i]),lower=qbeta(0.001,pa,pb),upper=qbeta(0.999,pa,pb))$value)
}

bayesian.beta.test=function(alpha.A,beta.A,alpha.B,beta.B,plot1=FALSE,plot2=FALSE) {
  if (length(alpha.A)<2 || length(alpha.B)<2) return(NA)
  fit=hierarchical.beta.posterior(a=c(alpha.A,alpha.B),b=c(beta.A,beta.B))
  fit1=hierarchical.beta.posterior(a=alpha.A,b=beta.A)
  #print(sqrt(bvar(fit1$a,fit1$b)))
  fit2=hierarchical.beta.posterior(a=alpha.B,b=beta.B)

  #fit=hierarchical.beta.posterior2(a=c(alpha.A,alpha.B),b=c(beta.A,beta.B))
  #fit1=hierarchical.beta.posterior2(a=alpha.A,b=beta.A)
  #fit2=hierarchical.beta.posterior2(a=alpha.B,b=beta.B)

  #fit=hierarchical.beta.posterior(a=c(alpha.A,alpha.B),b=c(beta.A,beta.B),compute.marginal.likelihood = TRUE,compute.grid=T,var.prior.min.sd=0.026)
  #fit1=hierarchical.beta.posterior(a=alpha.A,b=beta.A,compute.marginal.likelihood = TRUE,compute.grid=T,var.prior.min.sd=0.026)
  #fit2=hierarchical.beta.posterior(a=alpha.B,b=beta.B,compute.marginal.likelihood = TRUE,compute.grid=T,var.prior.min.sd=0.026)
  bvar=function(a,b) a*b/(a+b)^2/(a+b+1)
  bmean=function(a,b) a/(a+b)
  bdiff.mean=bmean(fit1$a,fit1$b)-bmean(fit2$a,fit2$b)
  bdiff.sd=sqrt(bvar(fit1$a,fit1$b)+bvar(fit2$a,fit2$b))

  if (plot1) {
    x=seq(min(qbeta(0.001,c(alpha.A,alpha.B),c(beta.A,beta.B))),max(qbeta(0.999,c(alpha.A,alpha.B),c(beta.A,beta.B))),length.out=1000)
    plot(x,pbeta(x,alpha.A[1],beta.A[1]),type='l',xlim=range(x),ylim=c(0,1)); for (i in 2:length(alpha.A)) lines(x,pbeta(x,alpha.A[i],beta.A[i]));for (i in 1:length(alpha.B)) lines(x,pbeta(x,alpha.B[i],beta.B[i]),col='red'); lines(x,pbeta(x,fit2$a,fit2$b),col='red',lwd=2); lines(x,pbeta(x,fit1$a,fit1$b),col='black',lwd=2); lines(x,pbeta(x,fit$a,fit$b),col='gray',lwd=2)
    efit=empirical.bayes.mixture(a=c(alpha.A,alpha.B),b=c(beta.A,beta.B))
    efit1=empirical.bayes.mixture(a=alpha.A,b=beta.A)
    efit2=empirical.bayes.mixture(a=alpha.B,b=beta.B)
    lines(x,pbeta(x,efit$a,efit$b),col='gray',lwd=2,lty=2)
    lines(x,pbeta(x,efit2$a,efit2$b),col='red',lwd=2,lty=2)
    lines(x,pbeta(x,efit1$a,efit1$b),col='black',lwd=2,lty=2)
    #lines(x,pnorm(x,bmean(fit1$a,fit1$b),sqrt(bvar(fit1$a,fit1$b))),col='black',lwd=2,lty=2)
    #lines(x,pnorm(x,bmean(fit2$a,fit2$b),sqrt(bvar(fit2$a,fit2$b))),col='red',lwd=2,lty=2)

  }

  if (plot2) {
    x=seq(qnorm(0.001,bdiff.mean,bdiff.sd),qnorm(0.999,bdiff.mean,bdiff.sd),length.out=1000)
    plot(x,pnorm(x,bdiff.mean,bdiff.sd),type='l',xlim=range(x),ylim=c(0,1));
    abline(v=0,lt=2)
  }
  p=pnorm(0,bdiff.mean,bdiff.sd)

  # two-sided p value
  if (p>0.5) 2*(1-p) else 2*p
}

#' Compute the posterior logFC distributions of RNA synthesis and degradation
#'
#' @param data the grandR object
#' @param name the name of the analysis added to the grandR object
#' @param contrasts A contrast matrix defining pairwise comparisons
#' @param steady.state either a reference matrix (\code{\link{FindReferences}}) to define steady state samples for each sample in A and B, or a
#' name list with names A and B containing logical vectors denoting the steady state samples
#' @param slot the data slot to take f0 and totals from
#' @param time.labeling the column in the column annotation table denoting the labeling duration or the labeling duration itself
#' @param time.experiment the column in the column annotation table denoting the experimental time point (can be NULL, see details)
#' @param ROPE.max.log2FC the region of practical equivalence is [-ROPE.max.log2FC,ROPE.max.log2FC] in log2 fold change space
#' @param sample.f0.in.ss whether or not to sample f0 under steady state conditions
#' @param N the sample size
#' @param conf.int A number between 0 and 1 representing the size of the credible interval
#' @param seed Seed for the random number generator
#' @param verbose Vebose output
#'
#' @details The kinetic parameters s and d are computed using \link{TransformOneShot}. For that, the sample either must be in steady state
#' (this is the case if defined in the steady.state matrix), or if the levels at a specific time point are known. This time point is
#' defined by \code{time.experiment} (i.e. the difference between the steady state samples and the A or B samples themselves). If
#' \code{time.experiment} is NULL, then the labeling time of the A or B samples is used (e.g. usefull if labeling was started concomitantly with
#' the perturbation, and the steady state samples are unperturbed samples).
#'
#' @return a new grandR object containing an additional analysis
#'
#'
#' @export
#'
EstimateRegulation=function(data,name,contrasts,steady.state,slot=DefaultSlot(data),time.labeling=Design$dur.4sU,time.experiment=NULL, ROPE.max.log2FC=0.25,sample.f0.in.ss=FALSE,N=10000,conf.int=0.95,seed=1337, hierarchical=FALSE,verbose=FALSE) {
  if (!check.slot(data,slot)) stop("Illegal slot definition!")
  if(!is.null(seed)) set.seed(seed)

  for (n in names(contrasts)) {
    if (verbose) cat(sprintf("Computing Regulation for %s...\n",n))
    A=contrasts[[n]]==1
    B=contrasts[[n]]==-1

    ss=if (is.matrix(steady.state)) list(A=apply(steady.state[,Columns(data,A)]==1,1,any),B=apply(steady.state[,Columns(data,B)]==1,1,any))
    dispersion.A = if (sum(ss$A)==1) rep(0.1,nrow(data)) else estimate.dispersion(GetTable(data,type="count",columns = ss$A))
    dispersion.B = if (sum(ss$B)==1) rep(0.1,nrow(data)) else estimate.dispersion(GetTable(data,type="count",columns = ss$B))

    if (verbose) {
      if (any(ss$A & A)) {
        cat(sprintf("Sampling from steady state for %s...\n",paste(colnames(data)[A],collapse = ",")))
      } else {
        cat(sprintf("Sampling from non-steady state for %s (steady-state: %s)...\n",paste(colnames(data)[A],collapse = ","),paste(colnames(data)[ss$A],collapse = ",")))
      }
      if (any(ss$B & B)) {
        cat(sprintf("Sampling from steady state for %s...\n",paste(colnames(data)[B],collapse = ",")))
      } else {
        cat(sprintf("Sampling from non-steady state for %s (steady-state: %s)...\n",paste(colnames(data)[B],collapse = ","),paste(colnames(data)[ss$B],collapse = ",")))
      }
    }

    re=plapply(1:nrow(data),function(i) {
    #for (i in 1:nrow(data)) { print (i);
      re=EstimateGeneRegulation(data=data,gene=i,A=A,B=B,dispersion.A =dispersion.A[i], dispersion.B = dispersion.B[i],steady.state = steady.state,slot=slot,time.labeling=time.labeling,time.experiment=time.experiment,ROPE.max.log2FC=ROPE.max.log2FC,sample.f0.in.ss = sample.f0.in.ss,return.samples = FALSE,N=N,conf.int = conf.int, hierarchical = hierarchical)
      unlist(re[c("s.A","s.B","HL.A","HL.B","s.log2FC","s.ROPE","HL.log2FC","HL.ROPE")])
    },seed=seed)

    re.df=as.data.frame(t(simplify2array(re)))

    data=AddAnalysis(data,description=MakeAnalysis(name = paste0(name,".",n),analysis = "Regulation",mode = "total",slot=slot,columns = colnames(data)[contrasts[[n]]==1|contrasts[[n]]==-1]),table = re.df)
  }
  data
}



#' Compute the posterior logFC distributions of RNA synthesis and degradation for a particular gene
#'
#' @param data the grandR object
#' @param gene a gene name or symbol or index
#' @param A columns for condition A (must refer to a unique labeling duration)
#' @param B columns for condition B (must refer to a unique labeling duration)
#' @param dispersion.A dispersion parameter for condition A (if NULL this is estimated, takes a lot of time!)
#' @param dispersion.B dispersion parameter for condition B (if NULL this is estimated, takes a lot of time!)
#' @param steady.state either a reference matrix (\code{\link{FindReferences}}) to define steady state samples for each sample in A and B, or a
#' name list with names A and B containing logical vectors denoting the steady state samples
#' @param slot the data slot to take f0 and totals from
#' @param time.labeling the column in the column annotation table denoting the labeling duration or the labeling duration itself
#' @param time.experiment the column in the column annotation table denoting the experimental time point (can be NULL, see details)
#' @param ROPE.max.log2FC the region of practical equivalence is [-ROPE.max.log2FC,ROPE.max.log2FC] in log2 fold change space
#' @param sample.f0.in.ss whether or not to sample f0 under steady state conditions
#' @param return.samples return the sampled logFCs for s and HL?
#' @param N the sample size
#' @param conf.int A number between 0 and 1 representing the size of the credible interval
#'
#' @details The kinetic parameters s and d are computed using \link{TransformOneShot}. For that, the sample either must be in steady state
#' (this is the case if defined in the steady.state matrix), or if the levels at a specific time point are known. This time point is
#' defined by \code{time.experiment} (i.e. the difference between the steady state samples and the A or B samples themselves). If
#' \code{time.experiment} is NULL, then the labeling time of the A or B samples is used (e.g. usefull if labeling was started concomitantly with
#' the perturbation, and the steady state samples are unperturbed samples).
#'
#' @return a list containing the posterior median, the credible interval and whether
#' 0 is within the credible interval for both synthesis rate (s) and RNA half-life (HL).
#' If return.samples=TRUE  the list also contains a data frame containing all samples
#'
#'
#' @export
#'
EstimateGeneRegulation=function(data,gene,A,B,dispersion.A=NULL,dispersion.B=NULL,steady.state,slot=DefaultSlot(data),time.labeling=Design$dur.4sU,time.experiment=NULL, ROPE.max.log2FC=0.25, sample.f0.in.ss=FALSE, return.samples=FALSE,N=10000,conf.int=0.95, hierarchical=FALSE) {

  if (is.matrix(steady.state)) steady.state=list(A=apply(steady.state[,Columns(data,A)]==1,1,any),B=apply(steady.state[,Columns(data,B)]==1,1,any))

  #ntr=GetData(data,mode.slot=c("ntr",slot,"new.norm","old.norm"),genes=gene,columns = A|B)
  alpha.A=GetData(data,mode.slot="alpha",genes=gene,columns = A)
  alpha.B=GetData(data,mode.slot="alpha",genes=gene,columns = B)
  if (!is.numeric(time.labeling) && length(unique(alpha.A[[time.labeling]]))!=1) stop("A has to refer to a unique labeling duration!")
  if (!is.numeric(time.labeling) && length(unique(alpha.B[[time.labeling]]))!=1) stop("B has to refer to a unique labeling duration!")
  if (!is.null(time.experiment) && length(unique(alpha.A[[time.experiment]]))!=1) stop("A has to refer to a unique experimental time!")
  if (!is.null(time.experiment) && length(unique(alpha.B[[time.experiment]]))!=1) stop("B has to refer to a unique experimental time!")
  beta.A=GetData(data,mode.slot="beta",genes=gene,columns = A)
  beta.B=GetData(data,mode.slot="beta",genes=gene,columns = B)
  total.A=GetData(data,mode.slot=slot,genes=gene,columns = A)
  total.B=GetData(data,mode.slot=slot,genes=gene,columns = B)
  ss.A=GetData(data,mode.slot=slot,genes=gene,columns = steady.state$A)
  ss.B=GetData(data,mode.slot=slot,genes=gene,columns = steady.state$B)
  if (!is.null(time.experiment) && length(unique(ss.A[[time.experiment]]))!=1) stop("Steady state for A has to refer to a unique experimental time!")
  if (!is.null(time.experiment) && length(unique(ss.B[[time.experiment]]))!=1) stop("Steady state for B has to refer to a unique experimental time!")

  if (is.null(dispersion.A)) dispersion.A=estimate.dispersion(as.matrix(GetTable(data,type="count",columns=colnames(data)[A],gene.info = F)))[ToIndex(data,gene)]
  if (is.null(dispersion.B)) dispersion.B=estimate.dispersion(as.matrix(GetTable(data,type="count",columns=colnames(data)[B],gene.info = F)))[ToIndex(data,gene)]


  use.A=total.A$Value>0
  use.B=total.B$Value>0

  if (sum(use.A)<2 || sum(use.B)<2) {
    re=list(
      s.A=NA,
      s.B=NA,
      HL.A=NA,
      HL.B=NA,
      s.log2FC=0,
      s.conf.int=c(-Inf,Inf),
      s.ROPE=1,
      HL.log2FC=0,
      HL.conf.int=c(-Inf,Inf),
      HL.ROPE=1
    )
    if (return.samples)
      re$samples=data.frame(s.log2FC=numeric(0),HL.log2FC=numeric(0),A.s=numeric(0),A.HL=numeric(0),B.s=numeric(0),B.HL=numeric(0))
    return(re)
  }


  alpha.A=alpha.A[use.A,]
  beta.A=beta.A[use.A,]
  total.A=total.A[use.A,]
  ss.A=ss.A[use.A,]
  alpha.B=alpha.B[use.B,]
  beta.B=beta.B[use.B,]
  total.B=total.B[use.B,]
  ss.B=ss.B[use.B,]

  t.A=if (is.numeric(time.labeling)) time.labeling else unique(alpha.A[[time.labeling]])
  t.B=if (is.numeric(time.labeling)) time.labeling else unique(alpha.B[[time.labeling]])
  t0.A=if (is.null(time.experiment)) t.A else unique(alpha.A[[time.experiment]])-unique(ss.A[[time.experiment]])
  t0.B=if (is.null(time.experiment)) t.B else unique(alpha.B[[time.experiment]])-unique(ss.B[[time.experiment]])


  if (N<1) {
    ntr=sum(alpha.A$Value)/(sum(alpha.A$Value)+sum(beta.A$Value))
    comp.ss=function(disp,mod,total,t) {
      TransformOneShot(ntr=ntr,total=mean(total$Value),t=t)
      #d=-1/t*log(1-ntr)
      #s=f0*d
      #cbind(s,d)
    }
    comp.non.ss=function(disp,ss,mod,total,t,t0) {
      if (t0<=0) stop("Experimental time is not properly defined (the steady state sample must be prior to each of A and B)!")
      f0=mean(ss$Value)
      TransformOneShot(ntr=ntr,total=mean(total$Value),t=t,t0=t0,f0=f0)
    }

    samp.a=if (any(steady.state$A & A)) comp.ss(dispersion.A,mod.A,total.A,t.A) else comp.non.ss(dispersion.A,ss.A,mod.A,total.A,t.A,t0.A)
    samp.b=if (any(steady.state$B & B)) comp.ss(dispersion.B,mod.B,total.B,t.B) else comp.non.ss(dispersion.B,ss.B,mod.B,total.B,t.B,t0.B)


    savelfc=function(a,b) ifelse(is.infinite(a) & is.infinite(b),0,log2(a/b))
    lfc.s=savelfc(samp.a['s'],samp.b['s'])
    lfc.HL=savelfc(samp.b['d'],samp.a['d'])

    re=list(
      s.A=(samp.a['s']),
      s.B=(samp.b['s']),
      HL.A=log(2)/(samp.a['d']),
      HL.B=log(2)/(samp.b['d']),
      s.log2FC=unname(lfc.s),
      s.conf.int=c(-Inf,Inf),
      s.ROPE=1,
      HL.log2FC=unname(lfc.HL),
      HL.conf.int=c(-Inf,Inf),
      HL.ROPE=1
    )
    if (return.samples)
      re$samples=data.frame(s.log2FC=numeric(0),HL.log2FC=numeric(0),A.s=numeric(0),A.HL=numeric(0),B.s=numeric(0),B.HL=numeric(0))
    return(re)
  }

  bayesian.beta.test(alpha.A$Value,beta.A$Value,alpha.B$Value,beta.B$Value,TRUE,FALSE)
  if (hierarchical) {
    mod.A=hierarchical.beta.posterior(alpha.A$Value,beta.A$Value,compute.marginal.likelihood = FALSE,compute.grid = TRUE,res=50)$sample.mu
    mod.B=hierarchical.beta.posterior(alpha.B$Value,beta.B$Value,compute.marginal.likelihood = FALSE,compute.grid = TRUE,res=50)$sample.mu
  } else {
    fit1=empirical.bayes.mixture(alpha.A$Value,beta.A$Value)
    fit2=empirical.bayes.mixture(alpha.B$Value,beta.B$Value)
    mod.A=function(N) rbeta(N,fit1$a,fit1$b)
    mod.B=function(N) rbeta(N,fit2$a,fit2$b)
  }

  sample.ss=function(disp,mod,total,t,samples) {
    if (sample.f0.in.ss) {
      f0=0.1+rnbinom(N,size=1/disp,mu=mean(total$Value))
      ntr=mod(N)
      total=0.1+rnbinom(N,size=1/disp,mu=mean(total$Value))
      TransformOneShot(ntr=ntr,total=total,t=t,f0=f0,t0=t)
    } else {
      ntr=mod(N)
      total=0.1+rnbinom(N,size=1/disp,mu=mean(total$Value))
      TransformOneShot(ntr=ntr,total=total,t=t)
    }
    #d=-1/t*log(1-ntr)
    #s=f0*d
    #cbind(s,d)
  }
  sample.non.ss=function(disp,ss,mod,total,t,t0,samples) {
    if (t0<=0) stop("Experimental time is not properly defined (the steady state sample must be prior to each of A and B)!")
    f0=0.1+rnbinom(N,size=1/disp,mu=mean(ss$Value))
    total=0.1+rnbinom(N,size=1/disp,mu=mean(total$Value))
    ntr=mod(N)
    TransformOneShot(ntr=ntr,total=total,t=t,t0=t0,f0=f0)
    #Fval=pmin(mean(total$Value)*(1-ntr)/f0,1)
    #d=ifelse(Fval>=1,0,-1/t*log(Fval))
    #s=-1/t*mean(total$Value)*ntr * ifelse(Fval>=1,-1,ifelse(is.infinite(Fval),0,log(Fval)/(1-Fval)))
    #cbind(s,d)
  }

  samp.a=if (any(steady.state$A & A)) sample.ss(dispersion.A,mod.A,total.A,t.A,A) else sample.non.ss(dispersion.A,ss.A,mod.A,total.A,t.A,t0.A,A)
  samp.b=if (any(steady.state$B & B)) sample.ss(dispersion.B,mod.B,total.B,t.B,B) else sample.non.ss(dispersion.B,ss.B,mod.B,total.B,t.B,t0.B,B)


  savelfc=function(a,b) ifelse((is.infinite(a) & is.infinite(b)) | (a==0&b==0),0,log2(a/b))
  lfc.s=savelfc(samp.a[,'s'],samp.b[,'s'])
  lfc.HL=savelfc(samp.b[,'d'],samp.a[,'d'])

  sc=quantile(lfc.s,c(0.5-conf.int/2,0.5+conf.int/2))
  hc=quantile(lfc.HL,c(0.5-conf.int/2,0.5+conf.int/2))

  re=list(
    s.A=median(samp.a[,'s']),
    s.B=median(samp.b[,'s']),
    HL.A=log(2)/median(samp.a[,'d']),
    HL.B=log(2)/median(samp.b[,'d']),
    s.log2FC=median(lfc.s),
#    s.log2FC.SEM=mad(lfc.s)/sqrt(length(lfc.s)),
#    s.log2FC.SEMdist=median(lfc.s)/(mad(lfc.s)/sqrt(length(lfc.s))),
    s.conf.int=sc,
    #s.ROPE=sum(lfc.s>-ROPE.max.log2FC & lfc.s<ROPE.max.log2FC)/N,
    s.ROPE=min(sum(lfc.s>-ROPE.max.log2FC),sum(lfc.s<ROPE.max.log2FC))/N,
#    s.regulated=c0<sc[1] || 0>sc[2],
    HL.log2FC=median(lfc.HL),
#    HL.log2FC.SEM=mad(lfc.HL)/sqrt(length(lfc.HL)),
#    HL.log2FC.SEMdist=median(lfc.HL)/(mad(lfc.HL)/sqrt(length(lfc.HL))),
    HL.conf.int=hc,
    #HL.ROPE=sum(lfc.HL>-ROPE.max.log2FC & lfc.HL<ROPE.max.log2FC)/N,
    HL.ROPE=min(sum(lfc.HL>-ROPE.max.log2FC),sum(lfc.HL<ROPE.max.log2FC))/N
#    HL.regulated=0<hc[1] || 0>hc[2]
  )
  if (return.samples)
    re$samples=data.frame(s.log2FC=lfc.s,HL.log2FC=lfc.HL,A.s=samp.a[,'s'],A.HL=log(2)/samp.a[,'d'],B.s=samp.b[,'s'],B.HL=log(2)/samp.b[,'d'])

  re
}


empirical.bayes.mixture=function(a,b) {
  bmean=function(a,b) a/(a+b)
  bvar=function(a,b) a*b/(a+b)^2/(a+b+1)
  mix.mean=mean(bmean(a,b)) # mean of mixture mode is (weighted) average of means
  mix.var=mean(bvar(a,b)+bmean(a,b)^2)-mix.mean^2

  list(a=(mix.mean*(1-mix.mean)/mix.var-1)*mix.mean,b=(mix.mean*(1-mix.mean)/mix.var-1)*(1-mix.mean))
}


# the hyperprior for the variance of the prior beta is determined by var.prior.min.sd;
# if the data beta posteriors (with uniform prior) overlap, the variance of the prior
# is largely determined by this hyperprior; the parameter is the minimal sd to obtain
hierarchical.beta.posterior=function(a,b,
                  var.prior.min.sd=NULL,
                  compute.marginal.likelihood=FALSE,
                  compute.grid=FALSE,
                  fak.below.max=1000,
                  res=100,N=100) {
  if (length(a)!=length(b)) stop("Unequal length of alpha and beta!")
  if (length(a)<2) stop("<2 observations!")

  mu=a/(a+b)
  if (!is.null(var.prior.min.sd)) {
    max.size=mean(mu*(1-mu))/var.prior.min.sd^2
    f=function(x,o=max.size,s=max.size/100) log(1/(1+exp((x-o)/s))/s/log1p(exp(o/s)))
    lprior=function(pa,pb) f(pa+pb)
  } else {
    #lprior=function(pa,pb) 0
    lprior=function(pa,pb) -5/2*log(pa+pb)
    #lprior=function(pa,pb) dcauchy(pa+pb,0,1,log=T)
  }
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
  mu=a/(a+b)
  size.point=log(min(mean(mu)*(1-mean(mu))/var(mu)-1,min(a+b)))
  opt=optim(c(logitmu=log(mean(a/b)),logsize=size.point),function(v) ltrans.marg.posterior(v[1],v[2]),control = list(fnscale=-1))

  #start=c(logitmu=log(mean.par[1]/(1-mean.par[1])),logsize=log(mean.par[1]*(1-mean.par[1])/sd.par[1]^2-1))
  #opt=optim(start,function(v) ltrans.marg.posterior(v[1],v[2]),control = list(fnscale=-1))
  #lprior(to.ab(start[1],start[2])$a,to.ab(start[1],start[2])$b)
  #lprior(to.ab(logMAP[1],logMAP[2])$a,to.ab(logMAP[1],logMAP[2])$b)

  logMAP=opt$par
  re=c(to.ab(logMAP[1],logMAP[2]),MAP=opt$value)
  rmu=re$a/(re$a+re$b)
  if(rmu>max(mu) || rmu<min(mu)) warning("Prior dominates and biases estimate!")
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
  re
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
