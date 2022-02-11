
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
#' @param u.content the relative frequency of uridines in the reads
#' @param u.content.sd the standard deviation of the u content
#' @param read.length the read length for simulation
#' @param p.old the probability for a conversion in reads originating from old RNA
#' @param p.new the probability for a conversion in reads originating from new RNA
#' @param seed seed value for the random number generator (set to make it deterministic!)
#'
#' @return a matrix containing, per column, the simulated counts, the simulated NTRs,
#' (potentially the shape parameters of the beta distribution approximation,) and the true relative frequencies and ntrs
#'
#' @details The simulation proceeds as follows:
#' \enumerate{
#'   \item Draw for each gene the number of reads from a multinomia distribution parmetrized with the relative abundances
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
#' mat=SimulateReadsForSample(num.reads = 10000,rel.abundance = rep(1,100),ntr=0.5,beta.approx = TRUE)
#' plot(ecdf(pbeta(0.5,mat[,3],mat[,4])))
#' abline(0,1)
#' # this should roughly be a uniform distibution, which means that the posterior distributions are properly approximated!
#'
SimulateReadsForSample=function(num.reads=2E7,rel.abundance=setNames(rlnorm(1E4,meanlog = 4.5,sdlog = 1),paste0("Gene",1:1E4)),ntr=setNames(rbeta(1E4,1.5,3),paste0("Gene",1:1E4)),dispersion=0.05,beta.approx=FALSE,u.content=0.25,u.content.sd=0.05,read.length=75,p.old=1E-4,p.new=0.04, seed=NULL) {

  if(!is.null(seed)) set.seed(seed)

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
      if (is.na(ntr)) return(c(ntr=NA,alpha=NA,beta=NA))
      if(reads==0) return(c(ntr=0,alpha=1,beta=1))
    } else {
      if (is.na(ntr)) return(NA)
      if(reads==0) return(0)
    }

    us = table(rbinom(reads,size=read.length,prob=rbeta(reads,shape1,shape2)))
    us=us[as.integer(names(us))>0]
    u.histo=rep(0,max(as.integer(names(us))))
    u.histo[as.integer(names(us))]=unname(us)

    para=model.par(ntr=ntr,p.err=p.old,p.conv=p.new)
    mixmat=CreateMixMatrix(n.vector = u.histo,par=para)
    fit.ntr(mixmat,para,plot=FALSE,beta.approx=beta.approx)
  },seed=seed))
  if (!beta.approx) {sim.ntr=t(sim.ntr); colnames(sim.ntr)="ntr"}

  cbind(count=mat[,3],sim.ntr,true_freq=mat[,1],true_ntr=mat[,2])
}


#' Simulate a complete time course of metabolic labeling - nucleotide conversion RNA-seq data.
#'
#' This function takes a vector of \emph{true} synthesis rates and RNA half-lives, and then simulates
#' data for multiple time points and replicates. Both synthesis rate and RNA half-lives are assumed to be constant,
#' but the system might not be in steady-state.
#'
#' @param condition A user-defined condition name (which is placed into the \link{\code{Coldata}} of the final grandR object)
#' @param gene.info either a data frame containing gene annotation or a vector of gene names
#' @param s a vector of synthesis rates
#' @param HL a vector of RNA half-lives
#' @param dispersion a vector of dispersion parameters (estimate from data using DESeq2, e.g. by the estimate.dispersion utility function)
#' @param f0 the abundance at time t=0
#' @param num.reads a vector representing the number of reads for each sample
#' @param timepoints a vector representing the labeling duration (in h) for each sample
#' @param beta.approx should the beta approximation of the NTR posterior be computed?
#' @param verbose Print status updates
#' @param seed seed value for the random number generator (set to make it deterministic!)
#' @param ... provided to \link{\code{SimulateReadsForSample}}
#'
#' @return a grandR object containing the simulated data in its data slots and the true parameters in the gene annotation table
#' @export
#'
SimulateTimeCourse=function(condition,gene.info,s,HL,dispersion,f0=s/log(2)*HL,num.reads=rep(1E7,length(timepoints)),timepoints=c(0,0,0,1,1,1,2,2,2,4,4,4),beta.approx=FALSE,verbose=TRUE,seed=NULL,...) {

  if (!is.data.frame(gene.info)) gene.info=data.frame(Gene=as.character(gene.info),Symbol=as.character(gene.info))

  tt=gsub("0h","no4sU",paste0(timepoints,"h"))
  names=as.character(ddply(data.frame(Name=factor(tt,levels=unique(tt))),.(Name),function(s) data.frame(Name=paste(condition,s$Name,LETTERS[1:length(s$Name)],sep=".")))$Name)
  coldata=MakeColdata(names,design=c(Design$Condition,Design$dur.4sU,Design$Replicate))
  coldata$no4sU=timepoints==0

  d=log(2)/HL


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
  data=lapply(data,function(m) {rownames(m)=gene.info$Gene; colnames(m)=names; m})

  if (!is.null(seed)) set.seed(seed)

  for (i in seq_along(timepoints)) {
    if (verbose) cat(sprintf("Simulating %s (%d/%d)...\n",names[i],i,length(timepoints)))
    total=(f.new(timepoints[i],s,d)+f.old.nonequi(timepoints[i],f0,s,d))
    ntr=if (timepoints[i]==0) rep(NA,length(s)) else f.new(timepoints[i],s,d)/total
    sim=SimulateReadsForSample(num.reads=num.reads[i],rel.abundance=total,ntr=ntr,dispersion=dispersion,beta.approx = beta.approx,seed=if (is.null(seed)) NULL else runif(1,min=0,max=.Machine$integer.max),...)

    data$count[,i]=sim[,"count"]
    data$ntr[,i]=sim[,"ntr"]
    data$true_count[,i]=sim[,"true_freq"]*num.reads[i]
    data$true_ntr[,i]=sim[,"true_ntr"]
    if (beta.approx){
      data$alpha[,i]=sim[,"alpha"]
      data$beta[,i]=sim[,"beta"]
    }
  }

  gene.info$true_f0=f0
  gene.info$true_HL=HL
  gene.info$true_s=s

  re=grandR(prefix="Simulated",gene.info=gene.info,slots=data,coldata=coldata,metadata=list(Description="Simulated data"))
  DefaultSlot(re)="count"

  re
}
