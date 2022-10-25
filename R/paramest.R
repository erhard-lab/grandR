
EPS=1E-7
MAX_ERR=1E-2


model.par=function(ntr=0.1,p.err=4E-4,p.conv=0.05,p.mconv=0.05,shape=5) {
	r=list(ntr=ntr,p.err=p.err,p.conv=p.conv,p.mconv=p.mconv,shape=shape,logLik=NA)
	class(r)="grandR.model"
	r
}
model.par.empty=function(ntr=NA,p.err=NA,p.conv=NA,p.mconv=NA,shape=NA) {
	r=list(ntr=ntr,p.err=p.err,p.conv=p.conv,p.mconv=p.mconv,shape=shape,logLik=NA)
	class(r)="grandR.model"
	r
}
default.model.par=model.par()

print.grandR.model=function(m) {
	if (!is.na(m$p.conv) && !is.na(m$p.mconv)) {
		cat(sprintf("\nUntrained grandR model:\n\n ntr=%.3f\n p.err=%.3g\n p.conv=%.3g\n p.mconv=%.3g\n shape=%.2f",m$ntr,m$p.err,m$p.conv,m$p.mconv,m$shape))
	} else if (!is.na(m$p.conv)) {
		cat(sprintf("\nBinomMix grandR model:\n\n ntr=%.3f\n p.err=%.3g\n p.conv=%.3g\n",m$ntr,m$p.err,m$p.conv))
	} else if (!is.na(m$p.mconv)) {
		cat(sprintf("\nTB-BinomMix grandR model:\n\n ntr=%.3f\n p.err=%.3g\n E(p.conv)=%.3g\n\n p.mconv=%.3g\n shape=%.2f\n",m$ntr,m$p.err,etbeta(m$p.err,m$p.mconv,exp(m$shape),exp(-m$shape)),m$p.mconv,m$shape))
	}
	if (!is.na(m$logLik)) cat(sprintf("\nLog likelihood: %.1f\n",m$logLik))
	cat("\n")
}

is.binom=function(par) !is.na(par$p.conv)
is.tbbinom=function(par) is.na(par$p.conv)


logLik.grandR.model=function(par, k, size) {
	if (!is.na(par$p.conv) && !is.na(par$p.mconv)) {
		stop("Model not trained!")
	} else if (!is.na(par$p.conv)) {
		old=dbinom(k,size=size,prob=par$p.err,log=TRUE)
		new=dbinom(k,size=size,prob=par$p.conv,log=TRUE)
		mix=lse(old+log(1-par$ntr),new+log(par$ntr))
		return (cbind(old=old,new=new,mix=mix))
	} else if (!is.na(par$p.mconv)) {
		old=dbinom(k,size=size,prob=par$p.err,log=TRUE)
		new=dtbbinom(k,size=size,l=par$p.err,u=par$p.mconv,shape=par$shape,log=TRUE)
		mix=lse(old+log(1-par$ntr),new+log(par$ntr))
		return (cbind(old=old,new=new,mix=mix))
	}
	stop("Model invalid!")
}



ibeta=function(x,a,b) pbeta(x,a,b)*beta(a,b)
#libeta=function(x,a,b) pbeta(x,a,b,log.p=TRUE)+lbeta(a,b)
libeta=function(x,a,b) log(gsl::hyperg_2F1(a+b,1,a+1,x)) + a*log(x)+b*log(1-x)-log(a)

lsse=function(x){
   xmax <- which.max(x)
   log1p(sum(exp(x[-xmax]-x[xmax])))+x[xmax]
}
lse=function(u,v) {
    m = pmax(u,v);
    log(exp(u-m)+exp(v-m))+m;
}
logdiff <- function(l1, l2) { l1 + log1p(-exp(-(l1-l2))); }

etbeta=function(l=0,u=1,s1=1,s2=1) exp(logdiff(libeta(u,1+s1,s2),libeta(l,1+s1,s2))-logdiff(libeta(u,s1,s2),libeta(l,s1,s2)))
dtbeta=function(x,l=0,u=1,s1=1,s2=1,log=FALSE) {
 	r=if(log) dbeta(x,s1,s2,log=T)-logdiff(pbeta(u,s1,s2,log.p=TRUE),pbeta(l,s1,s2,log.p=TRUE)) else dbeta(x,s1,s2)/(pbeta(u,s1,s2)-pbeta(l,s1,s2))
	r[x<l|x>u]=0
	r
}

ptbeta=function(x,l=0,u=1,s1=1,s2=1,log.p=FALSE) {
	if (log.p) logdiff(pbeta(x,s1,s2,log.p=TRUE),pbeta(l,s1,s2,log.p=TRUE))-logdiff(pbeta(u,s1,s2,log.p=TRUE),pbeta(l,s1,s2,log.p=TRUE)) else (pbeta(x,s1,s2)-pbeta(l,s1,s2))/(pbeta(u,s1,s2)-pbeta(l,s1,s2))
}


qtbeta=function(p,l=0,u=1,s1=1,s2=1,log.p=FALSE) {
	if (!log.p) p=log(p)
	lf=lse(p+logdiff(pbeta(u,s1,s2,log.p=TRUE),pbeta(l,s1,s2,log.p=TRUE)),pbeta(l,s1,s2,log.p=TRUE))
	qbeta(lf,s1,s2,log.p=TRUE)
}
#dtbbinom=function(x,size,l,u,s1=1,s2=1) choose(size,x)*(ibeta(u,x+s1,size-x+s2)-ibeta(l,x+s1,size-x+s2))/(ibeta(u,s1,s2)-ibeta(l,s1,s2))
dtbbinom=function(x,size,l,u,shape,log=FALSE) {
	s1=exp(shape)
	s2=exp(-shape)
	cx=pmin(x,size)
	r=ifelse(u<=l | x>size,-Inf,lchoose(size,cx)+logdiff(libeta(u,cx+s1,size-cx+s2),libeta(l,cx+s1,size-cx+s2))-logdiff(libeta(u,s1,s2),libeta(l,s1,s2)))
	if (log) r else exp(r)
}

rtbeta=function(n,l,u,s1,s2) {
	qtbeta(runif(n),l,u,s1,s2)
#	ql=pbeta(l,s1,s2)
#	qu=pbeta(u,s1,s2)
#	qbeta(runif(n,min=ql,max=qu),s1,s2)
}

rtbbinom=function(n,size,l,u,shape) {
	s1=exp(shape)
	s2=exp(-shape)
	p=rtbeta(n,l,u,s1,s2)
	rbinom(n,size=size,prob=p)
}


rbinommix=function(n,size,par=default.model.par,ntr=par$ntr,p.err=par$p.err,p.conv=par$p.conv) {
	n1=if (ntr<=0) round(-ntr) else rbinom(1,size=n,prob=ntr)
	c(rbinom(n1,size=size,prob=p.conv),rbinom(n-n1,size=size,prob=p.err))
}


rtbbinommix=function(n,size,par=default.model.par,ntr=par$ntr,p.err=par$p.err,p.mconv=par$p.mconv,shape=par$shape) {
	n1=if (ntr<=0) round(-ntr) else rbinom(1,size=n,prob=ntr)
	c(rtbbinom(n1,size=size,p.err,p.mconv,shape),rbinom(n-n1,size=size,prob=p.err))
}
dtbbinommix=function(k,size,par=default.model.par,ntr=par$ntr,p.err=par$p.err,p.mconv=par$p.mconv,shape=par$shape,log=FALSE) {
	if (log) {
		lse(dbinom(k,size=size,prob=p.err,log=TRUE)+log(1-ntr),dtbbinom(k,size=size,l=p.err,u=p.mconv,shape=shape,log=TRUE)+log(ntr))
	} else {
		dbinom(k,size=size,prob=p.err)*(1-ntr)+dtbbinom(k,size=size,l=p.err,u=p.mconv,shape=shape)*ntr
	}
}

dbinommix=function(k,size,par=default.model.par,ntr=par$ntr,p.err=par$p.err,p.conv=par$p.conv,log=FALSE) {
	if (log) {
		lse(dbinom(k,size=size,prob=p.err,log=TRUE)+log(1-ntr),dbinom(k,size=size,prob=p.conv,log=TRUE)+log(ntr))
	} else {
		dbinom(k,size=size,prob=p.err)*(1-ntr)+dbinom(k,size=size,prob=p.conv)*ntr
	}
}


# Generate expected number vectors
etbbinommix=function(n,size,par=default.model.par,ntr=par$ntr,p.err=par$p.err,p.mconv=par$p.mconv,shape=par$shape) setNames(dtbbinommix(0:size,size,ntr=ntr,p.err=p.err,p.mconv=p.mconv,shape=shape)*n,0:size)
ebinommix=function(n,size,par=default.model.par,ntr=par$ntr,p.err=par$p.err,p.conv=par$p.conv)  setNames(dbinommix(0:size,size,ntr,p.err,p.conv)*n,0:size)
etbbinom=function(n,size,l,u,shape=0)  setNames(dtbbinom(0:size,size,l,u,exp(shape),exp(-shape))*n,0:size)
ebinom=function(n,size,prob=1E-4)  setNames(dbinom(0:size,size,prob)*n,0:size)

tabrtbbinommix=function(n,size,par=default.model.par,ntr=par$ntr,p.err=par$p.err,p.mconv=par$p.mconv,shape=par$shape) table(rtbbinommix(n,size,ntr=ntr,p.err=p.err,p.mconv=p.mconv,shape=shape))
tabrbinommix=function(n,size,par=default.model.par,ntr=par$ntr,p.err=par$p.err,p.conv=par$p.conv)  table(rbinommix(n,size,ntr=ntr,p.err=p.err,p.conv=p.conv))

ReadOldMixMatrices=function(data,types=c("binom","binomOverlap"),...) {
	r=lapply(types,function(type) {
		t=read.tsv(paste0(data$prefix,".",type,".tsv"))
		if ("Type" %in% names(t)) t=t[t$Type=="",]
		l=plyr::dlply(t,plyr::.(Condition),function(s) {
			m=reshape2::acast(d~n,data=s[,c("n","d","count")],value.var="count")
			m[is.na(m)]=0
			m=m[rowSums(m)>0,colSums(m)>0,drop=FALSE]
			m=m[,colnames(m)!="0",drop=FALSE]
			class(m)="MixMatrix"
			m
		})
	})
	names(r)=types
	r
}

ReadMixMatrices=function(data,...) {
  # R CMD check guard for non-standard evaluation
  Subread <- Label <- NULL

  t=read.tsv(paste0(data$prefix,".conversion.knmatrix.tsv.gz"))
	l=plyr::dlply(t,plyr::.(Condition,Subread,Label),function(s) {
		m=reshape2::acast(k~n,data=s[,c("n","k","Count")],value.var="Count")
		m[is.na(m)]=0
		m=m[rowSums(m)>0,colSums(m)>0,drop=FALSE]
		m=m[,colnames(m)!="0",drop=FALSE]
		class(m)="MixMatrix"
		attributes(m)=c(attributes(m),list(Condition=s$Condition[1],Subread=s$Subread[1],Label=s$Label[1]))
		m
	})
	attributes(l)=NULL
	l
}

ReadMixMatrix=function(data,condition,subread,label="4sU") {
	t=read.tsv(paste0(data$prefix,".conversion.knmatrix.tsv.gz"))
	s=t[t$Condition==condition&t$Subread==subread&t$Label==label,]
	m=reshape2::acast(k~n,data=s[,c("n","k","Count")],value.var="Count")
	m[is.na(m)]=0
	m=m[rowSums(m)>0,colSums(m)>0,drop=FALSE]
	m=m[,colnames(m)!="0",drop=FALSE]
	class(m)="MixMatrix"
	attributes(m)=c(attributes(m),list(Condition=s$Condition[1],Subread=s$Subread[1],Label=s$Label[1]))
	m
}


n.vector=function(n,c) {
	re=rep(0,max(n))
	re[n]=c
	re
}

CreateMixMatrix=function(n.vector=round(dnorm(0:50,mean=30,sd=10)*3E7),tabfun=tabrbinommix,...) {
	m=matrix(0,nrow=length(n.vector)+1,ncol=length(n.vector))
	rownames(m)=0:length(n.vector)
	colnames(m)=1:length(n.vector)

	for (n in 1:length(n.vector)) {
		if (n.vector[n]>0) {
			tab=tabfun(n.vector[n],n,...)
			m[names(tab),n]=tab
		}
	}
	m=m[rowSums(m)>0,colSums(m)>0,drop=FALSE]
	class(m)="MixMatrix"
	m
}

matrix2MixMatrix=function(m,n,k) {
	rownames(m)=k
	colnames(m)=n
	m=m[rowSums(m)>0,colSums(m)>0,drop=FALSE]
	class(m)="MixMatrix"
	m
}



ExpectedOld=function(mixmat,par=default.model.par,p.err=par$p.err) {
	ns=colSums(mixmat)
	n.vector=rep(0,max(as.integer(names(ns))))
	n.vector[as.integer(names(ns))]=ns
	CreateMixMatrix(n.vector,ebinom,prob=p.err)[GetMixMatk(mixmat),GetMixMatn(mixmat)]
}
ExpectedOldFraction=function(mixmat,par=default.model.par,p.err=par$p.err) {
	em=ExpectedOld(mixmat,par,p.err)
	em/(em+mixmat+1)
}

original.estimate=function(x, errp=5E-4) {
	comp.p=function(x) sum(apply(x,2,function(c) sum(c*as.integer(rownames(x))))) / sum(apply(x,1,function(c) sum(c*as.integer(colnames(x)))))

	exp.x=function(x,p,oc) {
		if (!all(as.integer(rownames(x))==0:(dim(x)[1]-1))) stop("Incomplete!")
		for (c in 1:dim(x)[2]) {
			occ = oc[,c]
			k=as.integer(rownames(x))
			n=as.integer(colnames(x))[c]
			x[!occ,c]=sum(x[occ,c])/sum(dbinom(k[occ],n,p))*dbinom(k[!occ],n,p)
		}
		x[is.na(x)]=0
		x
	}

	onlyconv=apply(sapply(GetMixMatn(x),function(n) dbinom(GetMixMatk(x),n,errp)*sum(x[,n]))<0.01*unclass(x),2,cummax)==1
	emat=x

	onlyconv=onlyconv&x>0
	onlyconv=onlyconv&matrix(rep(apply(onlyconv,2,sum)>1,each=dim(onlyconv)[1]),nrow=dim(onlyconv)[1])
	if (sum(onlyconv)==0) stop("Cannot estimate!")

	init.x=x
	init.x[matrix(rep(apply(onlyconv,2,sum)==0,each=dim(x)[1]),nrow=dim(x)[1])]=0


	inter=c(0,1)
	while(TRUE) {
		p=sum(inter)/2
		x=exp.x(x,p,onlyconv)
		np=comp.p(x)
		if (np<p) inter=c(inter[1],p) else inter=c(p,inter[2])
		print(c(inter,p,np))
		if (inter[2]-inter[1]<1E-12) break;
	}
	conv=sum(inter)/2

	semat=sapply(GetMixMatn(emat),function(n) {
		dbinom(GetMixMatk(emat),n,conv)*sum(emat[,n])
	})
	rmat=(semat*sum(unclass(emat)[onlyconv])/sum(unclass(semat)[onlyconv])-unclass(emat))/unclass(emat)
	rmat[!onlyconv]=NA

	#print(rmat)
	conv
}



computeMixMatrix=function(m1,m2,fun) {
	if (!("MixMatrix" %in% class(m2))) return(structure(fun(unclass(m1),m2), class = "MixMatrix"))
	if (!("MixMatrix" %in% class(m2))) return(structure(fun(m1,unclass(m2)), class = "MixMatrix"))
	rows=sort(as.integer(union(rownames(m1),rownames(m2))))
	cols=sort(as.integer(union(colnames(m1),colnames(m2))))
	r1=matrix(0,nrow=length(rows),ncol=length(cols))
	r2=matrix(0,nrow=length(rows),ncol=length(cols))
	rownames(r1)=rows
	colnames(r1)=cols
	rownames(r2)=rows
	colnames(r2)=cols
	r1[rownames(m1),colnames(m1)]=m1
	r2[rownames(m2),colnames(m2)]=m2
	fun(r1,r2)
}

GetMixMatk=function(mixmat) as.integer(rownames(mixmat))
GetMixMatn=function(mixmat) as.integer(colnames(mixmat))
as.data.frame.MixMatrix=function(x) {
	df=reshape2::melt(unclass(x),varnames=c("k","n"),value.name="Count")
	df=df[df$Count>0,]
	rownames(df)=NULL
	df
}


`[.MixMatrix`=function(x,k=0:dim(x)[1],n=1:dim(x)[2]) {
	if (is.logical(k)) k=GetMixMatk(x)[k]
	if (is.logical(n)) n=GetMixMatk(x)[n]
	k=as.character(k)
	n=as.character(n)
	if (!all(k %in% rownames(x)) || !all(n %in% colnames(x))) {
		m=matrix(0,nrow=length(k),ncol=length(n))
		rownames(m)=k
		colnames(m)=n
		class(m)="MixMatrix"
		x=x+m
	}
	structure(unclass(x)[k,n,drop=FALSE], class = "MixMatrix")
}
"+.MixMatrix"=function(m1, m2) computeMixMatrix(m1,m2,match.fun(FUN="+"))
"-.MixMatrix"=function(m1, m2) computeMixMatrix(m1,m2,match.fun(FUN="-"))
"*.MixMatrix"=function(m1, m2) computeMixMatrix(m1,m2,match.fun(FUN="*"))
"/.MixMatrix"=function(m1, m2) computeMixMatrix(m1,m2,match.fun(FUN="/"))


ComputeComponentLikelihoods=function(mixmat,par) {
	k=rep(rep(GetMixMatk(mixmat),dim(mixmat)[2]),times=as.numeric(mixmat))
	n=rep(rep(GetMixMatn(mixmat),each=dim(mixmat)[1]),times=as.numeric(mixmat))
	logLik(par,k,n)
}

ComputeCombinationLogLikelihood=function(llold,llnew) {
	a=llold
	b=llnew
	re = rep(0,length(a)+1)
	re[1]=a[1]
	re[2]=b[1]
	for (k in 2:length(a)) {
		re[k+1]=re[k]+b[k]
		re[k:2]=sapply(k:2,function(i)lse(re[i]+a[k],re[i-1]+b[k]))
		re[1]=re[1]+a[k]
	}
	re=re-lchoose(length(a),0:length(a))
	re
}
mask.MixMat=function(m,p_err=4E-4,max.frac=0.01) {
  an=GetMixMatn(m)
  ak=GetMixMatk(m)
  re=sapply(an,function(n) {
    re=unclass(m[0:n,n])[,1]
    s=sum(re)
    E=dbinom(0:n,prob = p_err,size=n)*s
    first=min(c(which(E<max.frac*re),length(re)+1))
    re[1:(first-1)]=NA
    c(re,rep(0,max(ak)))[1:max(ak+1)]
  })
  rownames(re)=0:max(ak)
  colnames(re)=an
  structure(re[ak+1,],class="MixMatrix")
}
logLik.MixMat=function(m,fun,...) {
	an=GetMixMatn(m)
	ak=GetMixMatk(m)
	re=sum(sapply(an,function(n) {
		hak=ak[ak<=n]
		sum(fun(hak,n,log=T,...)*unclass(m[hak,n]),na.rm = TRUE)
	}))
	re
}

fit.ntr=function(mixmat,par,beta.approx=FALSE,conversion.reads=FALSE,plot=FALSE) {
  start=0
  optfun=function(p) logLik.MixMat(mixmat,dbinommix,ntr=p,p.err=par$p.err,p.conv=par$p.conv)-start
  start=optfun(par$ntr)
  opt=optimize(optfun,maximum=TRUE,lower=0,upper=1)
  if (!beta.approx) {
    return(if (conversion.reads) c(ntr=opt$maximum,conversion.reads=sum(mixmat)-sum(mixmat[0,])) else opt$maximum)
  }

  start=opt$objective+start

  left=if(optfun(0)>log(1E-3)) 0 else uniroot(function(x) optfun(x)-log(1E-3),c(0,opt$maximum))$root
  right=if(optfun(1)>log(1E-3)) 1 else uniroot(function(x) optfun(x)-log(1E-3),c(opt$maximum,1))$root

  x=seq(left,right,length.out=100)
  start=0
  fs = sapply(x,optfun)-log(length(x)-1);
  fs2=rep(NA,length(fs))

  fs[1]=fs[1]-log(2);
  fs2[1] = fs[1];

  for (i in 2:length(fs)) {
    fs2[i] =  lse(fs[i-1],fs[i]-log(2));
    fs[i] = lse(fs[i-1],fs[i]);
  }
  fs2=exp(fs2-fs2[length(fs2)])
  shapes=constrOptim(c(alpha=3,beta=3),function(par) sum((pbeta(x,par[1],par[2])-fs2)^2),grad=NULL,ui=cbind(c(1,0),c(0,1)),ci=c(0,0))$par

  if (plot) {
    plot(x,fs2,xlab="ntr",ylab="Cumulative freq",type='l')
    graphics::lines(x,pbeta(x,shapes[1],shapes[2]),col='red')
    graphics::legend("topleft",legend=c("Actual distribution","Beta approximation"),fill=c("black","red"))
  }

  if (conversion.reads) c(ntr=opt$maximum,shapes,conversion.reads=sum(mixmat)-sum(mixmat[0,])) else c(ntr=opt$maximum,shapes)
}

binom.optim=function(mixmat,par,fix=c(F,F,F)) {
	fix=fix[1:3]
	pp=unlist(par[c("ntr","p.err","p.conv")])
	first=rep(0,length(fix))
	first[!fix]=1:sum(!fix)
	sel=function(p,i) if (fix[i]) pp[i] else p[first[i]]
	start=0
	optfun=function(p) logLik.MixMat(mixmat,dbinommix,ntr=sel(p,1),p.err=sel(p,2),p.conv=sel(p,3))-start
	start=optfun(pp[!fix])

	if (sum(!fix)==1) {
		l=c(0,0,0)
		u=c(1,MAX_ERR,1)
		opt=optimize(optfun,maximum=TRUE,lower=l[!fix],upper=u[!fix])
		re=list(convergence=0,par=setNames(opt$maximum,names(pp)[!fix]),value=opt$objective)
	} else {
		ui=cbind(c(1,-1,0,0,0,0),c(0,0,1,-1,0,0),c(0,0,0,0,1,-1))
		ci=c(0,-1,0,-MAX_ERR,0,-1)
		ui=ui[,!fix,drop=FALSE]
		use=apply(ui!=0,1,any)
		ui=ui[use,,drop=FALSE]
		ci=ci[use]
		re=constrOptim(pp[!fix],f=optfun,grad=NULL,ui=ui,ci=ci,control=list(fnscale=-1))
	}
	re$value=re$value+start
	re
}



tbbinom.optim=function(mixmat,par,fix=c(F,F,F,F)) {
  fix=fix[1:4]
  pp=unlist(par[c("ntr","p.err","p.mconv","shape")])
  first=rep(0,length(fix))
  first[!fix]=1:sum(!fix)
  sel=function(p,i) if (fix[i]) pp[i] else p[first[i]]

  start=0
  optfun=function(p) logLik.MixMat(mixmat,dtbbinommix,ntr=sel(p,1),p.err=sel(p,2),p.mconv=sel(p,3),shape=sel(p,4))-start
  start=optfun(pp[!fix])

  if (sum(!fix)==1) {
    l=c(0,0,0,-10)
    u=c(1,MAX_ERR,1,10)
    opt=optimize(optfun,maximum=TRUE,lower=l[!fix],upper=u[!fix])
    re=list(convergence=0,par=setNames(opt$maximum,names(pp)[!fix]),value=opt$objective)
  } else {
    ui=cbind(c(1,-1,0,0,0,0),c(0,0,1,-1,0,0),c(0,0,0,0,1,-1),c(0,0,0,0,0,0))
    ci=c(0,-1,0,-MAX_ERR,0,-1)
    ui=ui[,!fix,drop=FALSE]
    use=apply(ui!=0,1,any)
    ui=ui[use,,drop=FALSE]
    ci=ci[use]
    re=constrOptim(pp[!fix],f=optfun,grad=NULL,ui=ui,ci=ci,control=list(fnscale=-1))
  }
  re$value=re$value+start
  re
}


binom.optim2=function(mixmat,par) {
  p_err=unlist(par["p.err"])
  mixmat=mask.MixMat(mixmat,p_err=p_err)
  pp=unlist(par[c("p.conv")])

  mlsse=function(v) if (length(v)==0) 0 else lsse(v)
  ll=function(p) {
    a=logLik.MixMat(mixmat,dbinom,prob=p)
    b=sum(sapply(GetMixMatn(mixmat),function(n) {
      v=unclass(mixmat[0:n,n])[,1]
      mlsse(dbinom((0:n)[!is.na(v)],size=n,prob=p,log=TRUE))*sum(v[!is.na(v)])
     }
    ))
    a-b
  }
  start=0
  optfun=function(p) ll(p=p)-start
  start=optfun(pp)

  opt=optimize(optfun,maximum=TRUE,lower=0,upper=1)
  re=list(convergence=0,par=setNames(opt$maximum,names(pp)),value=opt$objective)

  re$value=re$value+start
  re
}


tbbinom.optim2=function(mixmat,par,fix=c(F,F)) {
  p_err=unlist(par["p.err"])
  mixmat=mask.MixMat(mixmat,p_err=p_err)
  fix=fix[1:2]
  pp=unlist(par[c("p.mconv","shape")])
  first=rep(0,length(fix))
  first[!fix]=1:sum(!fix)
  sel=function(p,i) if (fix[i]) pp[i] else p[first[i]]

  #ll=function(p,D,kmin=2) sum(dbinom(D[D>=kmin],size=20,prob=p,log = TRUE) - lsse(dbinom(kmin:20,size=20,prob=p,log=TRUE)))
  mlsse=function(v) if (length(v)==0) 0 else lsse(v)
  ll=function(p,shape) {
    a=logLik.MixMat(mixmat,dtbbinom,l=p_err,u=p,shape=shape)
    b=sum(sapply(GetMixMatn(mixmat),function(n) {
      v=unclass(mixmat[0:n,n])[,1]
      mlsse(dtbbinom((0:n)[!is.na(v)],size=n,l=p_err,u=p,shape=shape,log=TRUE))*sum(v[!is.na(v)])
    }
    ))
    a-b
  }
  start=0
  optfun=function(p) ll(p=sel(p,1),shape=sel(p,2))-start
  start=optfun(pp[!fix])

  if (sum(!fix)==1) {
    l=c(0,-10)
    u=c(1,10)
    opt=optimize(optfun,maximum=TRUE,lower=l[!fix],upper=u[!fix])
    re=list(convergence=0,par=setNames(opt$maximum,names(pp)[!fix]),value=opt$objective)
  } else {
    ui=cbind(c(1,-1),c(0,0))
    ci=c(0,-1)
    ui=ui[,!fix,drop=FALSE]
    use=apply(ui!=0,1,any)
    ui=ui[use,,drop=FALSE]
    ci=ci[use]
    re=constrOptim(pp[!fix],f=optfun,grad=NULL,ui=ui,ci=ci,control=list(fnscale=-1))
  }
  re$value=re$value+start
  re
}


fit.MixMat=function(mixmat,par=default.model.par,type=c("binom","tubinom","tbbinom"),fix=rep(FALSE,4)) {
	if (is.list(mixmat)) return(lapply(mixmat,fit.MixMat,par=par,type=type,fix=fix))

	opti=switch(type[1],
			binom=binom.optim(mixmat,par,fix),
			tubinom=tbbinom.optim(mixmat,model.par(ntr=par$ntr,p.err=par$p.err,p.mconv=par$p.mconv,0),fix=c(fix[1:3],FALSE)),
			tbbinom=tbbinom.optim(mixmat,par,fix)
		)
	if (opti$convergence!=0) stop("Did not converge!")

	opti=c(opti$par,logLik=opti$value)
	re=switch(type[1],
			binom=model.par.empty(ntr=par$ntr,p.err=par$p.err,p.conv=par$p.conv),
			tubinom=model.par.empty(ntr=par$ntr,p.err=par$p.err,p.mconv=par$p.mconv),
			tbbinom=model.par.empty(ntr=par$ntr,p.err=par$p.err,p.mconv=par$p.mconv,shape=par$shape)
		)
	for (n in names(opti)) re[[n]]=unname(opti[n])
	re
}



table.MixMatrix.Ratios=function(mm,num=2,denom=3:10,coldata=NULL) {
	if (is.null(names(mm))) names(mm)=1:length(mm)
	r=lapply(mm,function(mixmat) colSums(mixmat[num,])/colSums(mixmat[denom,]))
	re=data.frame(n=as.integer(names(r[[1]])),r[[1]])
	names(re)[2]=names(mm)[1]
	if (length(mm)>1) {
		for (i in 2:length(mm)) {
			re2=data.frame(n=as.integer(names(r[[i]])),r[[i]])
			names(re2)[2]=names(mm)[i]
			re=merge(re,re2,by='n',all=T)
		}
	}
	re[is.na(re)]=0
	df=reshape2::melt(re,id.vars="n",variable.name="Name",value.name="Ratio")
	if (!is.null(coldata)) df=merge(df,coldata,by="Name")
	df
}



table.MixMatrix.Ratio=function(mixmat,num=2,denom=3:10,model.binom=fit.MixMat(mixmat,type="binom"),model.tbbinom=fit.MixMat(mixmat,type="tbbinom")) {
	r=colSums(mixmat[num,])/colSums(mixmat[denom,])
	r=r[as.integer(names(r))>num]
	re=data.frame(n=as.integer(names(r)),Data=r)
	re[is.na(re)]=0
	re$Binom=dbinom(num,size=re$n,prob=model.binom$p.conv)/sapply(re$n,function(n) sum(dbinom(denom,size=n,prob=model.binom$p.conv)))
	re$`TB-Binom`=dtbbinom(num,size=re$n,l=model.tbbinom$p.err,u=model.tbbinom$p.mconv,shape=model.tbbinom$shape)/sapply(re$n,function(n) sum(dtbbinom(denom,size=n,l=model.tbbinom$p.err,u=model.tbbinom$p.mconv,shape=model.tbbinom$shape)))
	df=reshape2::melt(re,id.vars="n",variable.name="Name",value.name="Ratio")
	df
}

table.MixMatrix.KL=function(mixmat,percentage=0.95,models=list()) {

	ns=sort(colSums(mixmat)/sum(mixmat),decreasing=T)
	nminmax=range(as.integer(names(ns[cumsum(ns)<percentage])))

	an=GetMixMatn(mixmat)
	kl=function(fun,...) {
		sapply(an,function(n) {
			hak=0:n
			k=as.vector(mixmat[hak,n])
			P=fun(0:n,n,...)
			Q=k/sum(k)
			sum((Q*(log(Q)-log(P)))[k>0])
		})
	}

	re=data.frame(n=an)
	for (i in 1:length(models)) {
		df=data.frame(kl(fun=if(is.binom(models[[i]])) dbinommix else dtbbinommix,par=models[[i]]))
		names(df)=names(models)[i]
		re=cbind(re,df)
	}
	re=re[re$n %in% nminmax[1]:nminmax[2],]
	df=reshape2::melt(re,id.vars="n",variable.name="Name",value.name="KL")
	df
}


table.MixMatrix.compare=function(mixmat,percentage=0.95,models=list()) {

	ns=sort(colSums(mixmat)/sum(mixmat),decreasing=T)
	nminmax=range(as.integer(names(ns[cumsum(ns)<percentage])))
	mixmat.cut=mixmat[,nminmax[1]:nminmax[2]]
	mixmat.scaled=t(t(mixmat.cut)/colSums(mixmat.cut))
	df=cbind(data.frame(Type="Data"),as.data.frame(mixmat.scaled))

	for (i in 1:length(models)) {
		emixmat=CreateMixMatrix(n.vector=colSums(mixmat.cut),tabfun=if(is.binom(models[[i]])) ebinommix else etbbinommix,par=models[[i]])
		emixmat[emixmat<1]=0
		emixmat=emixmat[rowSums(emixmat)>0,]
		emixmat.scaled=t(t(emixmat)/colSums(emixmat))
		df=rbind(df,cbind(data.frame(Type=names(models)[i]),as.data.frame(mixmat.scaled)))
	}
	df=df[df$k>0,]
	df
}


shape.adjust=function(ntr,global.par) {
    # how much of the new RNA measured at <labeling.time> has been made before t
    perc=function(x,d) 1-(1-exp(-(1-x)*d))/(1-exp(-1*d))

    # cumulative distribution function for global parameters, and adjusted for given degradation rates
    x=seq(global.par$p.err,global.par$p.mconv,length.out = 100)
    pnew=function(shape) ptbeta(x,global.par$p.err,global.par$p.mconv,exp(shape),exp(-shape))
    pnew.adj=function(shape,d) perc(pnew(shape),d)

    # compute degradation rate for ntr
    ntr2d=function(ntr) -1/1*log(1-ntr)

    global.d=ntr2d(global.par$ntr)

    # find shape if there was no degradation
    opt.fun=function(shape) sum((pnew(global.par$shape)-pnew.adj(shape,global.d))^2)
    global.shape.d0=optimize(opt.fun,lower=global.par$shape-5,upper=global.par$shape+5)$minimum

    opt.fun2=function(shape,this.d) sum((pnew(shape)-pnew.adj(global.shape.d0,this.d))^2)

    sapply(ntr,function(n) {
        d=max(1E-12,ntr2d(n))
        optimize(opt.fun2,this.d=d,lower=global.par$shape-5,upper=global.par$shape+5)$minimum
    })
}

PlotShapeAdjust=function(...) {
  # R CMD check guard for non-standard evaluation
  ntr <- shape.adj <- name <- gntr <- shape <- NULL

  make.df=function(g) {
        df=data.frame(ntr=seq(0,0.99,by=0.01))
        df$shape.adj=shape.adjust(df$ntr,g)
        df$name=sprintf("ntr=%.2f shape=%.1f",g$ntr,g$shape)
        df$gntr=g$ntr
        df$shape=g$shape
        df
    }

    global.par=list(...)
    df=do.call("rbind",lapply(global.par,function(g) make.df(g)))

    ggplot(df,aes(ntr,shape.adj,color=name))+
      cowplot::theme_cowplot()+
      geom_line()+
    geom_vline(data=unique(df[,c("gntr","shape","name")]),mapping=aes(xintercept=gntr,color=name),linetype=2,show.legend = FALSE)+
    geom_hline(data=unique(df[,c("gntr","shape","name")]),mapping=aes(yintercept=shape,color=name),linetype=2,show.legend = FALSE)+
    scale_color_brewer(NULL,palette = "Dark2")
}

#PlotShapeAdjust(model.par(ntr=0.1,shape=1.1),model.par(ntr=0.3,shape=1.1),model.par(ntr=0.3,shape=1.5))
