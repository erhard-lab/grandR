
#' Simulate metabolic labeling - nucleotide conversion RNA-seq data.
#'
#' This function takes a vector of \emph{true} relative abundances and NTRs, and then simulates
#' (i) read counts per gene and (ii) 4sU incorporation and conversion events. Subsequently, it
#' uses the same approach as implemented in the GRAND-SLAM 2.0 software (Juerges et al., Bioinformatics 2018)
#' to estimate the NTR from these simulated data.
#'
#' @param num.reads the total amount of reads for simulation
#' @param rel.abundance named (according to genes) vector of the true relative abundances. Is divided by its sum.
#' @param ntr vector of true NTRs
#' @param dispersion vector of dispersion parameters (should best be estimated by DESeq2)
#' @param beta.approx should the beta approximation of the NTR posterior be computed?
#' @param conversion.reads also output the number of reads with conversion
#' @param u.content the relative frequency of uridines in the reads
#' @param u.content.sd the standard deviation of the u content
#' @param read.length the read length for simulation
#' @param p.old the probability for a conversion in reads originating from old RNA
#' @param p.new the probability for a conversion in reads originating from new RNA
#' @param p.new.fit the probability for a conversion in reads originating from new RNA that is used for fitting (to simulate bias in the estimation of p.new)
#' @param enforce.parallelized should parallelization be used (NA: use state of IsParallel())
#' @param seed seed value for the random number generator (set to make it deterministic!)
#'
#' @return a matrix containing, per column, the simulated counts, the simulated NTRs,
#' (potentially the shape parameters of the beta distribution approximation,) and the true relative frequencies and ntrs
#'
#' @details The simulation proceeds as follows:
#' \enumerate{
#'   \item Draw for each gene the number of reads from a negative binomial distribution parametrized with the relative abundances x read number and the dispersion parameter
#'   \item For each gene: Draw for each read the number of uridines according to a beta binomial distribution for the given read length
#'   (the beta prior is parametrized to match the u.content and u.content.sd parameters)
#'   \item For each read: Draw the number of conversions according to the binomial mixture model of GRAND-SLAM
#'   (parametrized with p_old, p_new, the gene specific NTR and the read specific number of uridines)
#'   \item Estimate the NTR by using the GRAND-SLAM approach
#' }
#'
#'
#' @seealso \link{SimulateTimeCourse}
#' @export
#'
#' @examples
#' SimulateReadsForSample(num.reads = 10000,rel.abundance = rep(1,5),ntr=0.9)
#' SimulateReadsForSample(num.reads = 10000,rel.abundance = rep(1,5),ntr=0.9,seed=1337)
#' SimulateReadsForSample(num.reads = 10000,rel.abundance = rep(1,5),ntr=0.9,seed=1337)
#' # the second and third matrix should be equal, the first should be distinct
#'
#' @concept simulation
SimulateReadsForSample=function(num.reads=2E7,
                                rel.abundance=setNames(rlnorm(1E4,meanlog = 4.5,sdlog = 1),paste0("Gene",1:1E4)),
                                ntr=setNames(rbeta(1E4,1.5,3),paste0("Gene",1:1E4)),
                                dispersion=0.05,
                                beta.approx=FALSE,
                                conversion.reads=FALSE,
                                u.content=0.25,
                                u.content.sd=0.05,
                                read.length=75,
                                p.old=1E-4,
                                p.new=0.04,
                                p.new.fit=p.new,
                                enforce.parallelized = NA,
                                seed=NULL) {

  if(!is.null(seed)) set.seed(seed)

  ntr=pmin(pmax(ntr,0),1)
  mat=cbind(rel.abundance/sum(rel.abundance,na.rm=TRUE),ntr)
  mu=mat[,1]*num.reads
  mat=cbind(mat,rnbinom(length(mu),mu=mu,size=1/dispersion))
  #mat=cbind(mat,rmultinom(1,num.reads,rel.abundance))  # this does not model overdispersion

  shape1=u.content*(u.content*(1-u.content)/u.content.sd^2-1)
  shape2=(1-u.content)/u.content * shape1


  sim.ntr=t(psapply(1:nrow(mat),function(i)  {
    reads=unname(mat[i,3])
    ntr=unname(mat[i,2])
    if (beta.approx) {
      if (is.na(ntr)) {if(conversion.reads) return(c(ntr=NA,alpha=NA,beta=NA,conversion.reads=0)) else return(c(ntr=NA,alpha=NA,beta=NA))}
      if(reads==0) {if(conversion.reads) return(c(ntr=0,alpha=1,beta=1,conversion.reads=0)) else return(c(ntr=0,alpha=1,beta=1))}
    } else {
      if (is.na(ntr)) {if(conversion.reads) return(c(ntr=NA,conversion.reads=0)) else return(NA)}
      if(reads==0) {if(conversion.reads) return(c(ntr=0,conversion.reads=0)) else return(0)}
    }

    us = table(rbinom(reads,size=read.length,prob=rbeta(reads,shape1,shape2)))
    us=us[as.integer(names(us))>0]
    u.histo=rep(0,max(as.integer(names(us))))
    u.histo[as.integer(names(us))]=unname(us)

    para=model.par(ntr=ntr,p.err=p.old,p.conv=p.new)
    mixmat=CreateMixMatrix(n.vector = u.histo,par=para)
    para=model.par(ntr=ntr,p.err=p.old,p.conv=p.new.fit)
    fit.ntr(mixmat,para,plot=FALSE,beta.approx=beta.approx,conversion.reads=conversion.reads)
  },seed=seed,enforce = enforce.parallelized))
  if (!beta.approx & !conversion.reads) {sim.ntr=t(sim.ntr); colnames(sim.ntr)="ntr"}

  cbind(count=mat[,3],sim.ntr,true_freq=mat[,1],true_ntr=mat[,2])
}


#' Simulate a complete time course of metabolic labeling - nucleotide conversion RNA-seq data.
#'
#' This function takes a vector of \emph{true} synthesis rates and RNA half-lives, and then simulates
#' data for multiple time points and replicates. Both synthesis rate and RNA half-lives are assumed to be constant,
#' but the system might not be in steady-state.
#'
#' @param condition A user-defined condition name (which is placed into the \code{\link{Coldata}} of the final grandR object)
#' @param gene.info either a data frame containing gene annotation or a vector of gene names
#' @param s a vector of synthesis rates
#' @param d a vector of degradation rates (to get a specific half-life HL, use d=log(2)/HL)
#' @param f0 the abundance at time t=0
#' @param s.variation biological variability of s among all samples (see details)
#' @param d.variation biological variability of d among all samples (see details)
#' @param dispersion a vector of dispersion parameters (estimate from data using DESeq2, e.g. by the estimate.dispersion utility function)
#' @param num.reads a vector representing the number of reads for each sample
#' @param timepoints a vector representing the labeling duration (in h) for each sample
#' @param beta.approx should the beta approximation of the NTR posterior be computed?
#' @param conversion.reads also output the number of reads with conversion
#' @param verbose Print status updates
#' @param seed seed value for the random number generator (set to make it deterministic!)
#' @param ... provided to \code{\link{SimulateReadsForSample}}
#'
#' @details If \emph{s.variation} or \emph{d.variation} are > 1, then for each gene a random gaussian is added to s (or d)
#' such that 90% of all s (or d) are \emph{s.variation}-fold away from s (below or above; i.e. s is multiplied with 2^rnorm(n,mean=0,sd=x), with x chosen such that the 95% quantile
#' of the gaussian is log2(s.variation).
#'
#' @return a grandR object containing the simulated data in its data slots and the true parameters in the gene annotation table
#' @export
#'
#' @concept simulation
SimulateTimeCourse=function(condition,gene.info,s,d,f0=s/d,s.variation=1, d.variation=1, dispersion,num.reads=1E7,timepoints=c(0,0,0,1,1,1,2,2,2,4,4,4),beta.approx=FALSE,conversion.reads=FALSE,verbose=TRUE,seed=NULL,...) {
  # R CMD check guard for non-standard evaluation
  Name <- NULL

  if (!is.data.frame(gene.info)) gene.info=data.frame(Gene=as.character(gene.info),Symbol=as.character(gene.info))

  num.reads=rep(num.reads,length(timepoints))

  tt=gsub("^h$","no4sU",gsub("[_0]+h$","h",gsub(".","_",sprintf("%.2fh",timepoints),fixed=TRUE)))
  names=as.character(plyr::ddply(data.frame(Name=factor(tt,levels=unique(tt))),plyr::.(Name),function(s) data.frame(Name=paste(condition,s$Name,LETTERS[1:length(s$Name)],sep=".")))$Name)
  coldata=MakeColdata(names,design=c(Design$Condition,Design$dur.4sU,Design$Replicate))
  coldata$no4sU=timepoints==0

  sd.log2s = if (s.variation>1) uniroot(function(x) qnorm(0.95,sd=x)-log2(s.variation),lower=0,upper=9999)$root else 0
  sd.log2d = if (d.variation>1) uniroot(function(x) qnorm(0.95,sd=x)-log2(d.variation),lower=0,upper=9999)$root else 0


  data=list(
    count=matrix(nrow = length(s),ncol=length(timepoints)),
    ntr=matrix(nrow = length(s),ncol=length(timepoints)),
    true_count=matrix(nrow = length(s),ncol=length(timepoints)),
    true_ntr=matrix(nrow = length(s),ncol=length(timepoints))
  )
  if (beta.approx)
    data=c(data,list(
      alpha=matrix(nrow = length(s),ncol=length(timepoints)),
      beta=matrix(nrow = length(s),ncol=length(timepoints))
    ))
  if (conversion.reads)
    data=c(data,list(
      conversion_reads=matrix(nrow = length(s),ncol=length(timepoints))
    ))
  data=lapply(data,function(m) {rownames(m)=gene.info$Gene; colnames(m)=names; m})

  if (!is.null(seed)) set.seed(seed)

  for (i in seq_along(timepoints)) {
    if (verbose) cat(sprintf("Simulating %s (%d/%d)...\n",names[i],i,length(timepoints)))
    si=if (sd.log2s>0) s*2^rnorm(length(s),mean=0,sd=sd.log2s) else s
    di=if (sd.log2d>0) d*2^rnorm(length(d),mean=0,sd=sd.log2d) else d

    total=(f.new(timepoints[i],si,di)+f.old.nonequi(timepoints[i],f0,si,di))
    ntr=if (timepoints[i]==0) rep(NA,length(si)) else f.new(timepoints[i],si,di)/total
    sim=SimulateReadsForSample(num.reads=num.reads[i],rel.abundance=total,ntr=ntr,dispersion=dispersion,beta.approx = beta.approx,conversion.reads=conversion.reads,seed=if (is.null(seed)) NULL else runif(1,min=0,max=.Machine$integer.max),...)

    data$count[,i]=sim[,"count"]
    data$ntr[,i]=sim[,"ntr"]
    data$true_count[,i]=sim[,"true_freq"]*num.reads[i]
    data$true_ntr[,i]=sim[,"true_ntr"]
    if (beta.approx){
      data$alpha[,i]=sim[,"alpha"]
      data$beta[,i]=sim[,"beta"]
    }
    if (conversion.reads){
      data$conversion_reads[,i]=sim[,"conversion.reads"]
    }
  }

  gene.info$true_f0=f0
  gene.info$true_d=d
  gene.info$true_s=s

  re=grandR(prefix="Simulated",gene.info=gene.info,slots=data,coldata=coldata,metadata=list(Description="Simulated data"))
  DefaultSlot(re)="count"

  re
}


#' Simulate a complete time course of metabolic labeling - nucleotide conversion RNA-seq data.
#'
#' This function takes a vector of \emph{true} synthesis rates and RNA half-lives, and then simulates
#' data for multiple time points and replicates. Both synthesis rate and RNA half-lives are assumed to be constant,
#' but the system might not be in steady-state.
#'
#' @param condition A user-defined condition name (which is placed into the \code{\link{Coldata}} of the final grandR object)
#' @param gene.info either a data frame containing gene annotation or a vector of gene names
#' @param s a vector of synthesis rates (see details)
#' @param d a vector of degradation rates (see details)
#' @param dispersion a vector of dispersion parameters (estimate from data using DESeq2, e.g. by the estimate.dispersion utility function)
#' @param num.reads a vector representing the number of reads for each sample
#' @param t a single number denoting the time
#' @param replicates a single number denoting the number of replicates
#' @param beta.approx should the beta approximation of the NTR posterior be computed?
#' @param conversion.reads also output the number of reads with conversion
#' @param verbose Print status updates
#' @param seed seed value for the random number generator (set to make it deterministic!)
#' @param ... provided to \code{\link{SimulateReadsForSample}}
#'
#' @details Both rates can be either (i) a single number (constant rate), (ii) a data frame with names "offset",
#' "factor" and "exponent" (for linear functions, see \link{ComputeNonConstantParam}; only one row allowed) or
#' (iii) a unary function time->rate. Functions
#'
#' @return a grandR object containing the simulated data in its data slots and the true parameters in the gene annotation table
#' @export
#'
#' @seealso \link{SimulateTimeCourse}
#'
#' @concept simulation
SimulateTimeCourseNonConstant=function(condition,gene.info,s,d, dispersion,num.reads=1E7,t=2,replicates=3,beta.approx=FALSE,conversion.reads=FALSE,verbose=TRUE,seed=NULL,...) {
  # R CMD check guard for non-standard evaluation
  Name <- NULL

  N=if (is.data.frame(s)) nrow(s) else length(s)

  f0=(if (is.data.frame(s)) s$offset else sapply(s,function(FUN) FUN(0)))/(if (is.data.frame(d)) d$offset else sapply(d,function(FUN) FUN(0)))

  if (!is.data.frame(gene.info)) gene.info=data.frame(Gene=as.character(gene.info),Symbol=as.character(gene.info))

  timepoints=rep(t,replicates)
  num.reads=rep(num.reads,length(timepoints))

  tt=gsub("^h$","no4sU",gsub("[_0]+h$","h",gsub(".","_",sprintf("%.2fh",timepoints),fixed=TRUE)))
  names=as.character(plyr::ddply(data.frame(Name=factor(tt,levels=unique(tt))),plyr::.(Name),function(s) data.frame(Name=paste(condition,s$Name,LETTERS[1:length(s$Name)],sep=".")))$Name)
  coldata=MakeColdata(names,design=c(Design$Condition,Design$dur.4sU,Design$Replicate))
  coldata$no4sU=timepoints==0

  data=list(
    count=matrix(nrow = N,ncol=length(timepoints)),
    ntr=matrix(nrow = N,ncol=length(timepoints)),
    true_count=matrix(nrow = N,ncol=length(timepoints)),
    true_ntr=matrix(nrow = N,ncol=length(timepoints))
  )
  if (beta.approx)
    data=c(data,list(
      alpha=matrix(nrow = N,ncol=length(timepoints)),
      beta=matrix(nrow = N,ncol=length(timepoints))
    ))
  if (conversion.reads)
    data=c(data,list(
      conversion_reads=matrix(nrow = N,ncol=length(timepoints))
    ))
  data=lapply(data,function(m) {rownames(m)=gene.info$Gene; colnames(m)=names; m})

  if (!is.null(seed)) set.seed(seed)

  for (i in seq_along(timepoints)) {
    if (verbose) cat(sprintf("Simulating %s (%d/%d)...\n",names[i],i,length(timepoints)))

    #new=sapply(1:nrow(s),function(gi){
    #  f.nonconst.linear(t=timepoints[i],f0 = 0, so=s$offset[gi],sf=s$factor[gi],se=s$exponent[gi], do=d$offset[gi],df=d$factor[gi],de=d$exponent[gi])
    #})
    #old=sapply(1:nrow(s),function(gi){
    #  f.nonconst.linear(t=timepoints[i],f0 = s$offset[gi]/d$offset[gi], so=0,sf=0,se=1, do=d$offset[gi],df=d$factor[gi],de=d$exponent[gi])
    #})
    new=sapply(1:N,function(gi){
      f.nonconst(t=timepoints[i],f0 = 0, s=if(is.data.frame(s)) s[gi,] else s[[gi]],d=if(is.data.frame(d)) d[gi,] else d[[gi]])
    })
    old=sapply(1:N,function(gi){
      f.nonconst(t=timepoints[i],f0 = f0[gi], s=0,d=if(is.data.frame(d)) d[gi,] else d[[gi]])
    })

    total=new+old
    ntr=if (timepoints[i]==0) rep(NA,length(total)) else new/total
    sim=SimulateReadsForSample(num.reads=num.reads[i],rel.abundance=total,ntr=ntr,dispersion=dispersion,beta.approx = beta.approx,conversion.reads=conversion.reads,seed=if (is.null(seed)) NULL else runif(1,min=0,max=.Machine$integer.max),...)

    data$count[,i]=sim[,"count"]
    data$ntr[,i]=sim[,"ntr"]
    data$true_count[,i]=sim[,"true_freq"]*num.reads[i]
    data$true_ntr[,i]=sim[,"true_ntr"]
    if (beta.approx){
      data$alpha[,i]=sim[,"alpha"]
      data$beta[,i]=sim[,"beta"]
    }
    if (conversion.reads){
      data$conversion_reads[,i]=sim[,"conversion.reads"]
    }
  }

  gene.info$true_f0=f0
  gene.info$true_d=if (is.data.frame(d)) EvaluateNonConstantParam(t,d)$value else sapply(d,function(FUN) FUN(t))
  gene.info$true_s=if (is.data.frame(s)) EvaluateNonConstantParam(t,s)$value else sapply(s,function(FUN) FUN(t))


  re=grandR(prefix="Simulated",gene.info=gene.info,slots=data,coldata=coldata,metadata=list(Description="Simulated data"))
  DefaultSlot(re)="count"

  re
}
