
SimulateReadsForSample=function(num.reads=2E7,rel.freq=setNames(rlnorm(1E4,meanlog = 4.5,sdlog = 1),paste0("Gene",1:1E4)),ntr=setNames(rbeta(1E4,1.5,3),paste0("Gene",1:1E4)),beta.approx=FALSE,u.content=0.25,u.content.sd=0.05,read.length=75,p.old=1E-4,p.new=0.04) {

  mat=cbind(rel.freq/sum(rel.freq,na.rm=TRUE),ntr)
  mat=cbind(mat,rmultinom(1,num.reads,rel.freq))

  shape1=u.content*(u.content*(1-u.content)/u.content.sd^2-1)
  shape2=(1-u.content)/u.content * shape1


  sim.ntr=t(opt$sapply(1:nrow(mat),function(i)  {
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
  }))
  if (!beta.approx) {sim.ntr=t(sim.ntr); colnames(sim.ntr)="ntr"}

  cbind(count=mat[,3],sim.ntr,true_freq=mat[,1],true_ntr=mat[,2])
}


SimulateTimeCourse=function(condition,gene.info,s,HL,s0=s,HL0=HL,num.reads=sum(s*HL)/log(2),timepoints=c(0,0,0,1,1,1,2,2,2,4,4,4),beta.approx=FALSE,verbose=TRUE,...) {

  if (!is.data.frame(gene.info)) gene.info=data.frame(Gene=as.character(gene.info),Symbol=as.character(gene.info))

  tt=gsub("0h","no4sU",paste0(timepoints,"h"))
  names=as.character(ddply(data.frame(Name=factor(tt,levels=unique(tt))),.(Name),function(s) data.frame(Name=paste(condition,s$Name,LETTERS[1:length(s$Name)],sep=".")))$Name)
  coldata=MakeColdata(names,design=c(Design$Condition,Design$dur.4sU,Design$Replicate))

  # adapt s and s0 to match num.reads
  ns=s/sum(s*HL)*log(2)*num.reads
  s0=s0/s*ns
  s=ns

  d=log(2)/HL
  d0=log(2)/HL0

  rel.freq=s0/d0


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

  for (i in seq_along(timepoints)) {
    if (verbose) cat(sprintf("Simulating %s (%d/%d)...\n",names[i],i,length(timepoints)))
    ntr=if (timepoints[i]==0) rep(NA,length(s)) else f.new(timepoints[i],s,d)/(f.new(timepoints[i],s,d)+f.old.nonequi(timepoints[i],s0/d0,s,d))
    sim=SimulateReadsForSample(num.reads=num.reads,rel.freq=rel.freq,ntr=ntr,beta.approx = beta.approx,...)

    data$count[,i]=sim[,"count"]
    data$ntr[,i]=sim[,"ntr"]
    data$true_count[,i]=sim[,"true_freq"]*num.reads
    data$true_ntr[,i]=sim[,"true_ntr"]
    if (beta.approx){
      data$alpha[,i]=sim[,"alpha"]
      data$beta[,i]=sim[,"beta"]
    }
  }

  gene.info$true_s0=s0
  gene.info$true_HL=HL
  gene.info$true_s=s
  gene.info$true_HL0=HL0

  re=grandR(prefix="Simulated",gene.info=gene.info,data=data,coldata=coldata,metadata=list(Description="Simulated data"))
  DefaultSlot(re)="count"

  re
}
