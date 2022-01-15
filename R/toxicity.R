


Findno4sUPairs=function(data, paired.replicates=FALSE,discard.no4sU=TRUE) {
  pairs=FindReferences(data,reference=no4sU,group = if(paired.replicates) c(Design$Replicate,Design$Condition) else Design$Condition)
  #stopifnot(is.grandR(data))
  #df=data$coldata
  #df$group=if (paired.replicates) interaction(df[[Design$Condition]],df[[Design$Replicate]]) else df[[Design$Condition]]
  #map=dlply(df,.(group),function(s) as.character(s$Name[s$no4sU]))
  #pairs=setNames(lapply(df$group,function(g) map[[g]]),df$Name)
  if (discard.no4sU && "no4sU" %in% names(data$coldata)) pairs=pairs[as.character(data$coldata$Name[!data$coldata$no4sU])]
  pairs
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


CorrectBiasHLLen = function(data,pairs=Findno4sUPairs(data),LFC.fun=NormLFC) {
  for (i in 1:length(pairs)) {
    df=MakeToxicityTestTable(data=data,w4sU=names(pairs)[i],no4sU=pairs[[i]],transform=rank,LFC.fun=LFC.fun)
    df$tulen=pmin(80000,setNames(fil$gene.info$TU.len,fil$gene.info$Symbol)[rownames(df)])
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


EstimateTranscriptionLossLen = function(data,w4sU,no4sU,ntr=w4sU,LFC.fun=NormLFC,TU.len) {
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

EstimateTranscriptionLoss = function(data,w4sU,no4sU,ntr=w4sU,LFC.fun=NormLFC, type=c("quantreg","spearman","linear","lowess"),bootstrap=FALSE) {
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
  if ("no4sU" %in% names(data$coldata)) pairs=pairs[as.character(data$coldata$Name[!data$coldata$no4sU])]
  setNames(lapply(names(pairs),function(n) PlotToxicityTestLength(data,n,pairs[[n]],TU.len = TU.len,...)+ggtitle(n)),names(pairs))
}
PlotToxicityTestRankAll=function(data,pairs=Findno4sUPairs(data),...) {
  if ("no4sU" %in% names(data$coldata)) pairs=pairs[as.character(data$coldata$Name[!data$coldata$no4sU])]
  setNames(lapply(names(pairs),function(n) PlotToxicityTestRank(data,n,pairs[[n]],...)+ggtitle(n)),names(pairs))
}
PlotToxicityTestAll=function(data,pairs=Findno4sUPairs(data),...) {
  if ("no4sU" %in% names(data$coldata)) pairs=pairs[as.character(data$coldata$Name[!data$coldata$no4sU])]
  setNames(lapply(names(pairs),function(n) PlotToxicityTest(data,n,pairs[[n]],...)+ggtitle(n)),names(pairs))
}

DPlotToxicityTestAll=function(data,pairs=Findno4sUPairs(data),...) {
  rm(data)
  setNames(lapply(names(pairs),function(n) DPlot(PlotToxicityTest,w4sU=n,no4sU=pairs[[n]],add=ggtitle(n),height=4,...)),names(pairs))
}
DPlotToxicityTestRankAll=function(data,pairs=Findno4sUPairs(data),...) {
  rm(data)
  setNames(lapply(names(pairs),function(n) DPlot(PlotToxicityTestRank,w4sU=n,no4sU=pairs[[n]],add=ggtitle(n),height=4,...)),names(pairs))
}
DPlotToxicityTestLengthAll=function(data,pairs=Findno4sUPairs(data),TU.len="TU.len",...) {
  rm(data)
  setNames(lapply(names(pairs),function(n) DPlot(PlotToxicityTestRank,w4sU=n,no4sU=pairs[[n]],TU.len = TU.len,add=ggtitle(n),height=4,...)),names(pairs))
}

MakeToxicityTestTable=function(data,w4sU,no4sU=Findno4sUPairs(data)[[w4sU]],transform=rank,ntr=w4sU,LFC.fun=PsiLFC,slot="count",correction=1,TU.len=NULL) {
  w=rowMeans(GetTable(data,type=slot,columns=w4sU))
  n=if (is.numeric(no4sU)) no4sU[data$gene.info$Gene] else rowMeans(GetTable(data,type=slot,columns=no4sU))
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

PlotToxicityTestLength=function(data,w4sU,no4sU=Findno4sUPairs(data)[[w4sU]],ntr=w4sU,ylim=NULL,LFC.fun=PsiLFC,slot="count",TU.len="TU.len") {
  df=MakeToxicityTestTable(data=data,w4sU=w4sU,no4sU=no4sU,transform=rank,ntr=ntr,LFC.fun=LFC.fun,slot=slot,correction=1,TU.len = TU.len)
  if (is.null(ylim)) {
    d=max(abs(quantile(df$lfc,c(0.01,0.99))))*1.5
    ylim=c(-d,d)
  }
  ggplot(df,aes(tulen,lfc,color=density2d(log(tulen), lfc, n = 100)))+
    scale_color_viridis_c(name = "Density",guide=FALSE)+
    geom_point(alpha=1)+
    scale_x_log10("TU length")+
    geom_smooth(method="loess",formula=y~x)+
    ylab("log FC 4sU/no4sU")+
    coord_cartesian(ylim=ylim)
}

PlotToxicityTestRank=function(data,w4sU,no4sU=Findno4sUPairs(data)[[w4sU]],ntr=w4sU,ylim=NULL,LFC.fun=PsiLFC,slot="count",correction=1) {
  df=MakeToxicityTestTable(data=data,w4sU=w4sU,no4sU=no4sU,transform=rank,ntr=ntr,LFC.fun=LFC.fun,slot=slot,correction=correction)
  if (is.null(ylim)) {
    d=max(abs(quantile(df$lfc,c(0.01,0.99))))*1.5
    ylim=c(-d,d)
  }
  ggplot(df,aes(covar,lfc,color=density2d(covar, lfc, n = 100)))+
    scale_color_viridis_c(name = "Density",guide=FALSE)+
    geom_point(alpha=1)+
    geom_hline(yintercept=0)+
    geom_smooth(method="loess",formula=y~x)+
    xlab("NTR rank")+ylab("log FC 4sU/no4sU")+
    coord_cartesian(ylim=ylim)
}

PlotToxicityTest=function(data,w4sU,no4sU=Findno4sUPairs(data)[[w4sU]],ntr=w4sU,ylim=NULL,LFC.fun=PsiLFC,slot="count",hl.quantile=0.8,correction=1) {
  time=if(Design$dur.4sU %in% names(data$coldata)) data$coldata[ntr,Design$dur.4sU] else 1
  df=MakeToxicityTestTable(data=data,w4sU=w4sU,no4sU=no4sU,transform=function(x) comp.hl(x,time=time),ntr=ntr,LFC.fun=LFC.fun,slot=slot,correction=correction)
  df=df[df$covar<quantile(df$covar[is.finite(df$covar)],hl.quantile) & df$ntr<1 & df$ntr>0,]
  if (is.null(ylim)) {
    d=max(abs(quantile(df$lfc,c(0.01,0.99))))*1.5
    ylim=c(-d,d)
  }
  ggplot(df,aes(covar,lfc,color=density2d(covar, lfc, n = 100)))+
    scale_color_viridis_c(name = "Density",guide=FALSE)+
    geom_point(alpha=1)+
    geom_hline(yintercept=0)+
    geom_smooth(method="loess")+
    xlab("RNA half-life")+ylab("log FC 4sU/no4sU")+
    coord_cartesian(ylim=ylim)
}
