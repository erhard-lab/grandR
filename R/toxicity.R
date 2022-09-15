


#' Find equivalent no4sU samples for 4sU samples
#'
#' Identify all no4sU samples in the same condition, and return everything as a list to be used in
#' \link{PlotToxicityTest}, \link{PlotToxicityTestRank}, \link{PlotToxicityTestAll}, \link{PlotToxicityTestRankAll}
#'
#' @param data a grandR object
#' @param paired.replicates pair replicates, i.e. only no4sU.A is found for 4sU.A
#' @param discard.no4sU do not report references for no4sU samples
#'
#' @return a named list containing, for each 4sU sample, a vector of equivalent no4sU samples
#' @export
#'
#' @seealso \link{PlotToxicityTest}, \link{PlotToxicityTestRank}, \link{PlotToxicityTestAll}, \link{PlotToxicityTestRankAll}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' Findno4sUPairs(sars)
#'
#' @concept toxicity
Findno4sUPairs=function(data, paired.replicates=FALSE,discard.no4sU=TRUE) {
  # R CMD check guard for non-standard evaluation
  no4sU <- NULL

  grp=NULL
  if (paired.replicates && Design$Replicate %in% names(Coldata(data))) grp=c(grp,Design$Replicate)
  if (Design$Condition %in% names(Coldata(data))) grp=c(grp,Design$Condition)

  pairs=FindReferences(data,reference=no4sU,group = grp, as.list=TRUE)

  pairs=pairs[as.character(data$coldata$Name[!data$coldata$no4sU])]
  if (any(sapply(pairs,function(a) length(a)==0))) warning("There were samples without corresponding no4sU sample!")
  pairs=pairs[sapply(pairs,function(a) length(a)>0)]
  if (length(pairs)==0) stop("No no4sU pairs found!")
  pairs
}

MakeToxicityTestTable=function(data,w4sU,no4sU=Findno4sUPairs(data)[[w4sU]],transform=rank,ntr=w4sU,LFC.fun=lfc::PsiLFC,slot="count",correction=1,TU.len=NULL) {
  w=rowMeans(GetTable(data,type=slot,columns=w4sU))
  col=no4sU
  n=if (is.numeric(no4sU)) no4sU[data$gene.info$Gene] else rowMeans(GetTable(data,type=slot,columns=col))
  ntr=apply(GetTable(data,"ntr",columns=ntr),1,mean,rm.na=TRUE)
  use=!is.na(w+n+ntr)
  w=w[use]
  n=n[use]
  ntr=ntr[use]

  f1=1/correction-1
  nw=w+f1*ntr*w
  ntr=(ntr*w+f1*ntr*w)/nw
  w=nw

  phl=transform(ntr)

  df=data.frame(`4sU`=w,`no4sU`=n,ntr=ntr,lfc=LFC.fun(w,n),covar=phl,check.names = FALSE)
  if (!is.null(TU.len)) df$tulen=data$gene.info[[TU.len]][use]
  df
}


ComputeBiasCorrectionFactors=function(data,pairs=Findno4sUPairs(data),TU.len=NULL,...) {
  if (is.null(TU.len)) {
    sapply(names(pairs),function(n) EstimateTranscriptionLoss(data,n,pairs[[n]],...))
  } else {
    sapply(names(pairs),function(n) EstimateTranscriptionLossLen(data,n,pairs[[n]],TU.len=TU.len,...))
  }
}

CorrectBiasHLFactor=function(data,pairs=Findno4sUPairs(data),factors=ComputeBiasCorrectionFactors(data,pairs=pairs,...),...) {
  n=if(is.matrix(factors)) ncol(factors) else length(factors)
  for (i in 1:n) {
    if (is.matrix(factors)) {
      TU.len=list(...)$TU.len
      tulen=pmin(80000,setNames(data$gene.info[[TU.len]],data$gene.info$Gene)[rownames(data$data$count)])
      f=factors["f",i]*2^(-tulen/1000*factors["p",i])
      f=(1-f)/f
    } else {
      f=1/factors[i]-1
    }
    count=data$data$count[,names(pairs)[i]]
    ntr=data$data$ntr[,names(pairs)[i]]
    a=data$data$alpha[,names(pairs)[i]]
    b=data$data$beta[,names(pairs)[i]]
    data$data$count[,names(pairs)[i]] = count+f*count*ntr
    # assume all other tables are expression tables!
    for (n in setdiff(names(data$data),c("count","ntr","alpha","beta"))) data$data[[n]][,names(pairs)[i]] = data$data[[n]][,names(pairs)[i]]+f*data$data[[n]][,names(pairs)[i]]*ntr

    data$data$ntr[,names(pairs)[i]] = (ntr*count+f*ntr*count)/(count+f*count*ntr)
    ntr=data$data$ntr[,names(pairs)[i]]
    data$data$alpha[,names(pairs)[i]]=ntr*(a+b-2)+1
    data$data$beta[,names(pairs)[i]]=a+b-1-ntr*(a+b-2)
  }
  data
}

CorrectBiasHLNonlinear.old=function(data,pairs=Findno4sUPairs(data)) {
  for (i in 1:length(pairs)) {
    df=MakeToxicityTestTable(data=data,w4sU=names(pairs)[i],no4sU=pairs[[i]],transform=rank)
    df$idx=1:nrow(df)
    df=df[order(df$covar),]
    fit=lowess(df$covar,df$lfc,f=0.2)
    stopifnot(all(fit$x==df$covar))

    df$f=2^(-fit$y)
    df$f=df$f/min(df$f) # scale such that we have an increase for each gene!
    df=df[order(df$idx),]

    count=data$data$count[,names(pairs)[i]]
    ntr=data$data$ntr[,names(pairs)[i]]
    a=data$data$alpha[,names(pairs)[i]]
    b=data$data$beta[,names(pairs)[i]]

    #f1=(count*df$f-count)/(count*ntr)
    data$data$count[,names(pairs)[i]] = count*df$f #==count+f1*count*ntr
    for (n in setdiff(names(data$data),c("count","ntr","alpha","beta"))) data$data[[n]][,names(pairs)[i]] = data$data[[n]][,names(pairs)[i]]*df$f
    data$data$ntr[,names(pairs)[i]] = (ntr+df$f-1) / df$f #==(ntr*count+f1*ntr*count)/(count+f1*count*ntr)
    ntr=data$data$ntr[,names(pairs)[i]]
    data$data$alpha[,names(pairs)[i]]=ntr*(a+b-2)+1
    data$data$beta[,names(pairs)[i]]=a+b-1-ntr*(a+b-2)
 }
  data
}


CorrectBiasHLNonlinear=function(data,pairs=Findno4sUPairs(data),TU.len=NULL) {
  for (i in 1:length(pairs)) {
    df=MakeToxicityTestTable(data=data,w4sU=names(pairs)[i],no4sU=pairs[[i]],transform=rank,TU.len = TU.len)
    fit=if (is.null(TU.len))  loess(lfc~covar,data=df) else loess(lfc~covar+log(tulen),data=df)

    df$f=2^(-predict(fit))
    df$f=df$f/min(df$f) # scale such that we have an increase for each gene!

    count=data$data$count[,names(pairs)[i]]
    ntr=data$data$ntr[,names(pairs)[i]]
    a=data$data$alpha[,names(pairs)[i]]
    b=data$data$beta[,names(pairs)[i]]

    #f1=(count*df$f-count)/(count*ntr)
    data$data$count[,names(pairs)[i]] = count*df$f #==count+f1*count*ntr
    for (n in setdiff(names(data$data),c("count","ntr","alpha","beta"))) data$data[[n]][,names(pairs)[i]] = data$data[[n]][,names(pairs)[i]]*df$f
    data$data$ntr[,names(pairs)[i]] = (ntr+df$f-1) / df$f #==(ntr*count+f1*ntr*count)/(count+f1*count*ntr)
    ntr=data$data$ntr[,names(pairs)[i]]
    data$data$alpha[,names(pairs)[i]]=ntr*(a+b-2)+1
    data$data$beta[,names(pairs)[i]]=a+b-1-ntr*(a+b-2)
  }
  data
}


CorrectBiasHLLen = function(data,pairs=Findno4sUPairs(data),LFC.fun=lfc::NormLFC) {
  for (i in 1:length(pairs)) {
    df=MakeToxicityTestTable(data=data,w4sU=names(pairs)[i],no4sU=pairs[[i]],transform=rank,LFC.fun=LFC.fun)
    df$tulen=pmin(80000,setNames(data$gene.info$TU.len,data$gene.info$Symbol)[rownames(df)])
    obj=function(par) {
      p=par[1]
      f=par[2]
      df=df[!is.na(df$tulen),]
      f=f*2^(-df$tulen/1000*p)
      f1=(1-f)/f
      df2=data.frame(lfc = LFC.fun(df$`4sU`+df$`4sU`*df$ntr*f1, df$`no4sU`),covar=(df$`4sU`*df$ntr+df$`4sU`*df$ntr*f1)/(df$`4sU`+df$`4sU`*df$ntr*f1))
      sum(loess(lfc~covar,data=df2)$y^2)
    }
    par=optim(c(1,0.1),obj)
    p=par$par[1]
    f=par$par[2]
    f=f*2^(-df$tulen/1000*p)
    f=(1-f)/f

    count=data$data$count[,names(pairs)[i]]
    ntr=data$data$ntr[,names(pairs)[i]]
    a=data$data$alpha[,names(pairs)[i]]
    b=data$data$beta[,names(pairs)[i]]
    data$data$count[,names(pairs)[i]] = count+f*count*ntr
    # assume all other tables are expression tables!
    for (n in setdiff(names(data$data),c("count","ntr","alpha","beta"))) data$data[[n]][,names(pairs)[i]] = data$data[[n]][,names(pairs)[i]]+f*data$data[[n]][,names(pairs)[i]]*ntr

    data$data$ntr[,names(pairs)[i]] = (ntr*count+f*ntr*count)/(count+f*count*ntr)
    ntr=data$data$ntr[,names(pairs)[i]]
    data$data$alpha[,names(pairs)[i]]=ntr*(a+b-2)+1
    data$data$beta[,names(pairs)[i]]=a+b-1-ntr*(a+b-2)
  }
  data
}


EstimateTranscriptionLossLen = function(data,w4sU,no4sU,ntr=w4sU,LFC.fun=lfc::NormLFC,TU.len) {
  df=MakeToxicityTestTable(data=data,w4sU=w4sU,no4sU=no4sU,transform=rank,ntr=ntr,LFC.fun=LFC.fun,TU.len=TU.len)
  df$tulen[is.na(df$tulen)]=median(df$tulen,na.rm=TRUE)
  obj=function(par) {
    p=par[1]
    f=par[2]
    df=df[!is.na(df$tulen),]
    f=f*2^(-df$tulen/1000*p)
    f1=(1-f)/f
    df2=data.frame(lfc = LFC.fun(df$`4sU`+df$`4sU`*df$ntr*f1, df$`no4sU`),covar=(df$`4sU`*df$ntr+df$`4sU`*df$ntr*f1)/(df$`4sU`+df$`4sU`*df$ntr*f1))
    sum(loess(lfc~covar,data=df2)$y^2)
  }
  setNames(optim(c(0.01,0.5),obj)$par,c("p","f"))
}

EstimateTranscriptionLoss = function(data,w4sU,no4sU,ntr=w4sU,LFC.fun=lfc::NormLFC, type=c("quantreg","spearman","linear","lowess"),bootstrap=FALSE) {
  df=MakeToxicityTestTable(data=data,w4sU=w4sU,no4sU=no4sU,transform=rank,ntr=ntr,LFC.fun=LFC.fun)

  if (bootstrap) df = df[sample.int(nrow(df),nrow(df),replace=TRUE),]

  if (type[1]=="spearman") {
    obj=function(f1) {
      df=data.frame(lfc = LFC.fun(df$`4sU`+df$`4sU`*df$ntr*f1, df$`no4sU`),covar=(df$`4sU`*df$ntr+df$`4sU`*df$ntr*f1)/(df$`4sU`+df$`4sU`*df$ntr*f1))
      cor(df$lfc,df$covar,method="spearman",use='c')
    }
    l=obj(0)
    r=obj(19)
    if (l*r>=0) {
      if (abs(l)<abs(r)) return(1)
      return(0)
    }
    f1=uniroot(obj,c(0,19))$root
  }
  else if (type[1]=="linear") {
    obj=function(f1) lm(lfc~covar,data=data.frame(lfc = LFC.fun(df$`4sU`+df$`4sU`*df$ntr*f1, df$`no4sU`), covar=(df$`4sU`*df$ntr+df$`4sU`*df$ntr*f1)/(df$`4sU`+df$`4sU`*df$ntr*f1)))$coeff[2]
    l=obj(0)
    r=obj(19)
    if (l*r>=0) {
      if (abs(l)<abs(r)) return(1)
      return(0)
    }
    f1=uniroot(obj,c(0,19))$root
  }
  else if (type[1]=="lowess") {
    obj=function(f1) sum(loess(lfc~covar,data=data.frame(lfc = LFC.fun(df$`4sU`+df$`4sU`*df$ntr*f1, df$`no4sU`), covar=(df$`4sU`*df$ntr+df$`4sU`*df$ntr*f1)/(df$`4sU`+df$`4sU`*df$ntr*f1)))$y^2)
    f1=optimize(obj,c(0,19))$minimum
  }
  else if (type[1]=="quantreg") {
    obj=function(f1) quantreg::rq(lfc~covar,data=data.frame(lfc = LFC.fun(df$`4sU`+df$`4sU`*df$ntr*f1, df$`no4sU`), covar=(df$`4sU`*df$ntr+df$`4sU`*df$ntr*f1)/(df$`4sU`+df$`4sU`*df$ntr*f1)))$coeff[2]
    l=obj(0)
    r=obj(19)
    if (l*r>=0) {
      if (abs(l)<abs(r)) return(1)
      return(0)
    }
    f1=uniroot(obj,c(0,19))$root
  }

  1/(f1+1)

}

PlotToxicityTestLengthAll=function(data,pairs=Findno4sUPairs(data),TU.len="TU.len",...) {
  setNames(lapply(names(pairs),function(n) PlotToxicityTestLength(data,n,pairs[[n]],TU.len = TU.len,...)),names(pairs))
}
PlotToxicityTestLengthDeferAll=function(data,pairs=NULL,TU.len="TU.len",...) {
  if (is.null(pairs)) pairs=Findno4sUPairs(data)
  rm(data)
  setNames(lapply(names(pairs),function(n) Defer(PlotToxicityTestRank,w4sU=n,no4sU=pairs[[n]],TU.len = TU.len,add=ggtitle(n),height=4,...)),names(pairs))
}
PlotToxicityTestLength=function(data,w4sU,no4sU=Findno4sUPairs(data)[[w4sU]],ntr=w4sU,ylim=NULL,LFC.fun=lfc::PsiLFC,slot="count",TU.len="TU.len") {
  # R CMD check guard for non-standard evaluation
  tulen <- lfc <- NULL

  df=MakeToxicityTestTable(data=data,w4sU=w4sU,no4sU=no4sU,transform=rank,ntr=ntr,LFC.fun=LFC.fun,slot=slot,correction=1,TU.len = TU.len)
  if (is.null(ylim)) {
    d=max(abs(quantile(df$lfc,c(0.01,0.99))))*1.5
    ylim=c(-d,d)
  }
  ggplot(df,aes(tulen,lfc,color=density2d(log(tulen), lfc, n = 100)))+
    cowplot::theme_cowplot()+
    scale_color_viridis_c(name = "Density",guide=FALSE)+
    geom_point(alpha=1)+
    scale_x_log10("TU length")+
    geom_smooth(method="loess",formula=y~x)+
    ylab("log FC 4sU/no4sU")+
    coord_cartesian(ylim=ylim)+
    ggtitle(w4sU)
}

#' Perform toxicity tests
#'
#' Testing for toxicity of a 4sU sample is performed by comparing half-lives or NTR ranks against
#' the log2 fold change of the 4sU sample vs equivalent no4sU samples.
#'
#' @param data a grandR object
#' @param pairs a no4sU pairs list as generated by \link{Findno4sUPairs}
#' @param w4sU the name of a 4sU sample
#' @param no4sU the name(s) of equivalent no4sU sample(s)
#' @param ntr the name of a sample to take NTRs from (usually equal to w4sU)
#' @param ylim y axis limits
#' @param LFC.fun function to compute log fold change (default: \link[lfc]{PsiLFC}, other viable option: \link[lfc]{NormLFC})
#' @param slot the slot of the grandR object to take the data from; for \link[lfc]{PsiLFC}, this really should be "count"!
#' @param hl.quantile the half-life quantile to cut the plot
#' @param correction correction factor
#' @param ... further arguments to be passed to or from other methods.
#'
#' @details The deferred versions are useful to be used in conjunction with \link{ServeGrandR} plot.static. Their implementation
#' make sure that they are lightweight, i.e. when saving the returned function to an Rdata file, the grandR object is not stored.
#'
#' @return either a ggplot object, a list of ggplot objects, or a list of deferred functions for plotting
#'
#' @seealso \link{Findno4sUPairs},\link{Defer}
#'
#' @name toxicity
#' @concept toxicity
NULL
#> NULL

#' @rdname toxicity
#' @export
PlotToxicityTestRankAll=function(data,pairs=Findno4sUPairs(data),...) {
  setNames(lapply(names(pairs),function(n) PlotToxicityTestRank(data,n,pairs[[n]],...)),names(pairs))
}
#' @rdname toxicity
#' @export
PlotToxicityTestAll=function(data,pairs=Findno4sUPairs(data),...) {
  setNames(lapply(names(pairs),function(n) PlotToxicityTest(data,n,pairs[[n]],...)),names(pairs))
}

#' @rdname toxicity
#' @export
PlotToxicityTestDeferAll=function(data,pairs=NULL,...) {
  if (is.null(pairs)) pairs=Findno4sUPairs(data)
  rm(data)
  setNames(lapply(names(pairs),function(n) Defer(PlotToxicityTest,w4sU=n,no4sU=pairs[[n]],add=ggtitle(n),...)),names(pairs))
}
#' @rdname toxicity
#' @export
PlotToxicityTestRankDeferAll=function(data,pairs=NULL,...) {
  if (is.null(pairs)) pairs=Findno4sUPairs(data)
  rm(data)
  setNames(lapply(names(pairs),function(n) Defer(PlotToxicityTestRank,w4sU=n,no4sU=pairs[[n]],add=ggtitle(n),...)),names(pairs))
}

#' @rdname toxicity
#' @export
PlotToxicityTestRank=function(data,w4sU,no4sU=Findno4sUPairs(data)[[w4sU]],ntr=w4sU,ylim=NULL,LFC.fun=lfc::PsiLFC,slot="count",correction=1) {
  # R CMD check guard for non-standard evaluation
  covar <- lfc <- NULL

  df=MakeToxicityTestTable(data=data,w4sU=w4sU,no4sU=no4sU,transform=function(v) rank(v),ntr=ntr,LFC.fun=LFC.fun,slot=slot,correction=correction)
  if (is.null(ylim)) {
    d=max(abs(quantile(df$lfc,c(0.01,0.99))))*1.5
    ylim=c(-d,d)
  }
  rho=round(cor(df$covar,df$lfc,method="spearman"),digits = 2)
  p=cor.test(df$covar,df$lfc,method="spearman")$p.value
  p=if (p<2.2E-16) p = bquote("<"~2.2 %*% 10^-16) else p = sprintf("= %.2g",p)
  df$lfc=ifelse(df$lfc<ylim[1],-Inf,df$lfc)
  df$lfc=ifelse(df$lfc>ylim[2],+Inf,df$lfc)
  ggplot(df,aes(covar,lfc,color=density2d(covar, lfc, n = 100,margin = 'x')))+
    cowplot::theme_cowplot()+
    scale_color_viridis_c(name = "Density",guide='none')+
    geom_point(alpha=1)+
    geom_hline(yintercept=0)+
    #geom_smooth(method="loess",formula=y~x)+
    xlab("NTR rank")+ylab("log FC 4sU/no4sU")+
    coord_cartesian(ylim=ylim)+
    ggtitle(w4sU,subtitle = bquote(rho == .(rho) ~ "," ~ p ~ .(p)))
}

#' @rdname toxicity
#' @export
PlotToxicityTest=function(data,w4sU,no4sU=Findno4sUPairs(data)[[w4sU]],ntr=w4sU,ylim=NULL,LFC.fun=lfc::PsiLFC,slot="count",hl.quantile=0.8,correction=1) {
  # R CMD check guard for non-standard evaluation
  covar <- lfc <- NULL

  time=if(Design$dur.4sU %in% names(data$coldata)) data$coldata[ntr,Design$dur.4sU] else 1
  df=MakeToxicityTestTable(data=data,w4sU=w4sU,no4sU=no4sU,transform=function(x) comp.hl(x,time=time),ntr=ntr,LFC.fun=LFC.fun,slot=slot,correction=correction)
  df=df[df$covar<quantile(df$covar[is.finite(df$covar)],hl.quantile) & df$ntr<1 & df$ntr>0,]
  if (is.null(ylim)) {
    d=max(abs(quantile(df$lfc[is.finite(df$lfc)],c(0.01,0.99))))*1.5
    ylim=c(-d,d)
  }
  ggplot(df,aes(covar,lfc,color=density2d(covar, lfc, n = 100)))+
    cowplot::theme_cowplot()+
    scale_color_viridis_c(name = "Density",guide=FALSE)+
    geom_point(alpha=1)+
    geom_hline(yintercept=0)+
   # geom_smooth(method="loess")+
    xlab("RNA half-life")+ylab("log FC 4sU/no4sU")+
    coord_cartesian(ylim=ylim)+
    ggtitle(w4sU)
}
