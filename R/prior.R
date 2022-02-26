

considerations=function() {
  x=seq(0,1,length.out=1000)
  N=10000
  pa=16
  pb=64-16
  pp=rbeta(N,pa,pb)

  size=sample.int(51,N,replace = TRUE)+74
  k=rbinom(N,size=size,prob=pp)

  a=k
  b=size-k


  # this is where the pp came from
  plot(x,pbeta(x,pa,pb),type='l',xlim=c(0.1,0.5))
  # now draw from the individual posteriors (using Haldane's prior)
  ppp=rbeta(10*N,a,b)
  mean(ppp)
  lines(ecdf(ppp),col='red')
  # this does not fit, excess variance
  efit=hierarchical.beta.posterior(a,b)
  lines(x,pbeta(x,efit$a,efit$b),col='blue')
  # now draw from the individual posteriors (using the actual distribution of pp as prior)
  ppp=rbeta(N,a+pa,b+pb)
  mean(ppp)
  lines(ecdf(ppp),col='blue')
  # this fits!
}

dbetamix=function(p,a,b,log=FALSE) {
  r=sapply(p,function(x) ifelse(x<=0 | x>=1,-Inf,lsse(dbeta(x,a,b,log=T)-log(length(a)))))
  if (log) r else exp(r)
}

mm.fit.prior=function(a,b) {
  bvar=function(a,b) a*b/(a+b)^2/(a+b+1)
  bmean=function(a,b) a/(a+b)

  N=length(a)
  mu=mean(bmean(a,b)) # unbiased estimator for the mean of the prior
  paf=function(pa) {
    pb=pa/mu-pa
    post=mean(bvar(a+pa,b+pb)+bmean(a+pa,b+pb)^2)
    prior=bvar(pa,pb)+bmean(pa,pb)^2
    (post-prior)^2
    #mu*(1-mu)/(post-mu^2)-1-pa/mu
    #ff
  }
  myopt=function() {
    bound=100
    while(TRUE) {n
      opt=optimise(paf,c(0,bound),maximum = FALSE)
      cat(sprintf("b=%.1f  paf(b)=%.5g  opt=%.5g\n",bound,paf(bound),opt$objective))
      if (opt$minimum<bound*0.99) return(opt$minimum)
      bound=bound*1.5
    }
  }
  pa=myopt()
  pb=pa/mu-pa
  dd=seq(0,pa*2,length.out=100)
  plot(dd,sapply(dd,paf),log='y')
  c(a=pa,b=pb)
}


ll.fit.prior=function(a,b,plot=FALSE,prior=0) {
  dbetabinom=function(k,n,a,b,log=T) lchoose(n,k) +lbeta(k+a,n-k+b) - lbeta(a,b)
  mu=mean(bmean(a,b)) # unbiased estimator for the mean of the prior
  paf=function(pa) {pb=pa/mu-pa;sum(dbeta((a-1)/(a+b-2),a+pa,b+pb,log=T))}
  paf=function(pa) {pb=pa/mu-pa;sum(dbetabinom(a,a+b,a+pa,b+pb,log=T))}
  myopt=function() {
    bound=100
    while(TRUE) {
      opt=optimise(paf,c(0,bound),maximum = TRUE)
      if (paf(bound)<opt$objective) return(opt$maximum)
      bound=bound*1.5
    }
  }
  pa=myopt()
  pb=pa/mu-pa
  dd=seq(0,50,length.out=100)
  if (plot) plot(dd,sapply(dd,paf))
  c(a=pa,b=pb)
}

# find a distribution P~Beta(pa,pb) such that D(Q|P) is minimal for Q~BetaMix(1/N,...,1/N,k+pa,size-k+pb)
fit.prior=function(a,b,mu=sum(a)/sum(a+b)) {
  lsse=function(x){
    xmax <- which.max(x)
    log1p(sum(exp(x[-xmax]-x[xmax])))+x[xmax]
  }

  H=function(a,b) lbeta(a,b)-(a-1)*digamma(a)-(b-1)*digamma(b)+(a+b-2)*digamma(a+b)
  crossent=function(a,b) {
    samp=rbeta(1000,a,b)
    mean(dbetamix(samp,a=a+pa,b=b+pb,log=TRUE))
  }

  obj=function(nu,mu,a,b) {
    pa=mu*nu
    pb=(1-mu)*nu
    -H(pa,pb)-integrate(function(p) dbeta(p,pa,pb)*betamix(p,a+pa,b+pb,log=TRUE),lower=0,upper=1)$value
  }
  maxvar=var(a/(a+b))+mean(bvar(a,b))
  minvar=maxvar/10
  maxnu=mu*(1-mu)/minvar
  minnu=mu*(1-mu)/maxvar
  fit=optimise(obj,c(minnu,maxnu),mu=mu,a=a,b=b)
  if (fit$minimum>=maxnu*0.99) warning("Hit boundary!")
  c(a=mu*fit$minimum,b=(1-mu)*fit$minimum)
}

bb.fit.prior=function(a,b) {
  N=length(a)
  k=a
  size=a+b
  mu=sum(k)/sum(size)
  s2=(N*sum(size*(k/size-mu)^2))/((N-1)*sum(size))
  print(s2)
  nu=(mu*(1-mu)-s2)/(s2-mu*(1-mu)/N*sum(1/size))
  print(nu)
  c(a=nu*mu,b=(1-mu)*nu)
}


doit=function(N=5,total=100,mean=0.2,size=100) {
  pa=round(mean*total)
  pb=total-pa
  pp=rbeta(N,pa,pb)

  k=rbinom(N,size=size,prob=pp)

  a=k
  b=size-k
  fit=bayesian(a,b)
  fit2=hierarchical.beta.posterior(a,b,var.prior.min.sd = 0.05)
  print(unlist(c(fit,var=bvar(fit$a,fit$b),mean=bmean(fit$a,fit$b))))
  print(unlist(c(fit2,var=bvar(fit2$a,fit2$b),mean=bmean(fit2$a,fit2$b))))

  plot(x,dbeta(x,pa,pb),type='l',xlim=qbeta(c(0.0001,0.9999),pa,pb))
  rug(a/(a+b))
  #for (i in 1:N) lines(x,dbeta(x,a[i],b[i]))
  lines(x,dbetamix(x,a+fit$a,b+fit$b),col='orange')
  lines(x,dbeta(x,fit$a,fit$b),col='red',lwd=2)
  lines(x,dbetamix(x,a+fit2$a,b+fit2$b),col='blue')
  lines(x,dbeta(x,fit2$a,fit2$b),col='darkblue',lwd=2)
  #plot(x,pbeta(x,pa,pb),xlim=c(0.1,0.5),type='l',lwd=2); for (k in 1:N) lines(x,pbeta(x,a[k],b[k]))
  list(a,b)
}

bbtest=function() {
  N=10
  pa=60
  pb=200
  pp=rbeta(N,pa,pb)

  k=rbinom(N,size=400,prob=pp)

  a=k
  b=400-k

  bayesian.beta.test(a[1:5],b[1:5],a[6:10],b[6:10],plot=c(0.15,0.3))


}


bayesfac=function(N=5) {
  pa=c(rep(16,N),rep(32,N))
  pb=c(rep(64-16,N),rep(64-32,N))
  N2=N
  N=N*2
  pp=rbeta(N,pa,pb)
  size=sample.int(51,N,replace = TRUE)+74
  k=rbinom(N,size=size,prob=pp)
  a=k
  b=size-k
  fit=bayesian(a,b,compute.marginal.likelihood = TRUE)
  fit1=bayesian(a[1:N2],b[1:N2],compute.marginal.likelihood = TRUE)
  fit2=bayesian(a[(N2+1):N],b[(N2+1):N],compute.marginal.likelihood = TRUE)
  (fit1["mar.loglik"]+fit2["mar.loglik"]-fit["mar.loglik"])/log(10)
}

bvar=function(a,b) a*b/(a+b)^2/(a+b+1)
bmean=function(a,b) a/(a+b)

doit=function(N=5) {
  pa=16
  pb=64-16
  pp=rbeta(N,pa,pb)

  size=sample.int(51,N,replace = TRUE)+74
  k=rbinom(N,size=size,prob=pp)

  a=k
  b=size-k

  sum(bayesian(a,b))
}

#k=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0, 0, 0, 0, 1, 1, 1, 1, 1, 1,1, 1, 2, 2, 2, 2, 2, 2, 2, 2,2, 1, 5, 2, 5, 3, 2, 7, 7, 3,3, 2, 9, 10, 4, 4, 4, 4, 4, 4,4, 10, 4, 4, 4, 5, 11, 12, 5, 5,6, 5, 6, 6, 6, 6, 16, 15, 15, 9)
#n=c(20 ,20 ,20 ,20 ,20 ,20 ,20 ,19 ,19 ,19,19 ,18 ,18 ,17 ,20 ,20 ,20 ,20 ,19 ,19,18 ,18 ,25 ,24 ,23 ,20 ,20 ,20 ,20 ,20,20 ,10 ,49 ,19 ,46 ,27 ,17 ,49 ,47 ,20,20 ,13 ,48 ,50 ,20 ,20 ,20 ,20 ,20 ,20,20 ,48 ,19 ,19 ,19 ,22 ,46 ,49 ,20 ,20,23 ,19 ,22 ,20 ,20 ,20 ,52 ,47 ,46 ,24)
#a=k
#b=n-k

bayesian.gamma=function(a,b,
                  var.prior.k=NULL,
                  var.prior.theta=NULL,
                  var.prior.mean=NULL,
                  var.prior.sd=NULL,
                  compute.marginal.likelihood=FALSE,
                  compute.grid=FALSE,
                  fak.below.max=1000,
                  res=100) {
  if (length(a)!=length(b)) stop("Unequal length of alpha and beta!")
  if (length(a)<2) stop("<2 observations!")

  if (is.null(var.prior.k) || is.null(var.prior.theta)) {
    var.prior.k=var.prior.mean^2/var.prior.sd^2
    var.prior.theta=var.prior.sd^2/var.prior.mean
  }

  if ((is.null(var.prior.k) || length(var.prior.k)==0) || (is.null(var.prior.theta) || length(var.prior.theta)==0)) {
    lprior=function(pa,pb) (-5/2)*log(pa+pb)
  } else {
    lprior=function(pa,pb) {
      vv=pa*pb/(pa+pb)^2/(pa+pb+1)
      lJ=log(pa)+log(pb)-3*log(pa+pb)-2*log(pa+pb+1)
      dgamma(vv,shape=var.prior.k,scale=var.prior.theta,log=TRUE)+lJ
    }
  }
  lmarg.posterior=function(pa,pb) lprior(pa,pb)+sum(lbeta(pa+a,pb+b)-lbeta(pa,pb))
  ltrans.marg.posterior=function(logitmu,logsize) {
    pa=exp(logitmu+logsize)/(exp(logitmu)+1)
    pb=exp(logsize)/(exp(logitmu)+1)
    lJ=(logitmu+2*logsize)-2*log1p(exp(logitmu))
    lmarg.posterior(pa,pb)+lJ
  }
  to.ab=function(logitmu,logsize) list(a=unname(exp(logitmu+logsize)/(exp(logitmu)+1)),
                                    b=unname(exp(logsize)/(exp(logitmu)+1)))

  # find MAP
  mu=a/(a+b)
  size.point=log(min(mean(mu)*(1-mean(mu))/var(mu)-1,min(a+b)))
  opt=optim(c(logitmu=log(mean(a/b)),logsize=size.point),function(v) ltrans.marg.posterior(v[1],v[2]),control = list(fnscale=-1))
  logMAP=opt$par
  re=c(to.ab(logMAP[1],logMAP[2]),MAP=opt$value)
  if (compute.marginal.likelihood || compute.grid) {
    mu.low=uniroot(function(x) ltrans.marg.posterior(x,logMAP[2])-opt$value+log(fak.below.max),interval = c(logMAP[1]-4,logMAP[1]),extendInt = "upX")$root
    mu.high=uniroot(function(x) ltrans.marg.posterior(x,logMAP[2])-opt$value+log(fak.below.max),interval = c(logMAP[1],logMAP[1]+4),extendInt = "downX")$root

    size.low=uniroot(function(x) ltrans.marg.posterior(logMAP[1],x)-opt$value+log(fak.below.max),interval = c(logMAP[2]-4,logMAP[2]),extendInt = "upX")$root
    size.high=uniroot(function(x) ltrans.marg.posterior(logMAP[1],x)-opt$value+log(fak.below.max),interval = c(logMAP[2],logMAP[2]+4),extendInt = "downX")$root

    if (compute.marginal.likelihood) {
      ml=log(cubature::adaptIntegrate(function(v) exp(ltrans.marg.posterior(v[1],v[2])-opt$value),lowerLimit = c(mu.low,size.low),upperLimit = c(mu.high,size.high))$integral)+opt$value
      re=c(re,list(mar.loglik=ml))
    }

    if (compute.grid) {
      mu.grid=seq(mu.low,mu.high,length.out=res)
      size.grid=seq(size.low,size.high,length.out=res)
      grid=sapply(mu.grid,function(mu) sapply(size.grid, function(size)  exp(ltrans.marg.posterior(mu,size)-opt$value) ))
      colnames(grid)=mu.grid
      rownames(grid)=size.grid
      dimnames(grid)=c("logit mu","log size")
      re=c(re,list(grid=grid))
      #image(t(grid))
      #contour(t(grid),levels = (seq(0.05,0.95,by=0.05)))
    }
  }
  re
}

# the prior of the variance should ensure unimodality.
# the prior for the variance is a scaled beta distribution. instead of 0,1, it is in 0,1/3*mu(1-mu), which means that a+b>2 of the posterior (i.e. unimodal)
# the beta parameter of this scaled beta distribution must be >2 (then density and slope at max is 0)
# when it is specified via var.prior.mean and .sd,
bayesian=function(a,b,
                  var.prior.a=NULL,
                  var.prior.b=NULL,
                  var.prior.mean=NULL,
                  var.prior.sd=NULL,
                  compute.marginal.likelihood=FALSE,
                  compute.grid=FALSE,
                  fak.below.max=1000,
                  res=100) {
  if (length(a)!=length(b)) stop("Unequal length of alpha and beta!")
  if (length(a)<2) stop("<2 observations!")

  if ((is.null(var.prior.a) || is.null(var.prior.b)) && (is.null(var.prior.mean) || is.null(var.prior.sd))) {
    lprior=function(pa,pb) (-5/2)*log(pa+pb)
  } else {
    lprior=function(pa,pb) {
      mu=pa/(pa+pb)
      maxvar=1/3*mu*(1-mu)
      if (is.null(var.prior.a) || is.null(var.prior.b)) {
        var.prior.a=-var.prior.mean*(-maxvar*var.prior.mean+var.prior.mean^2+var.prior.sd^2)/(var.prior.sd^2*maxvar)
        var.prior.b=-(maxvar-var.prior.mean)*(-maxvar*var.prior.mean+var.prior.mean^2+var.prior.sd^2)/(var.prior.sd^2*maxvar)
      }
      if (var.prior.b<2) warning("Improper prior!")
      vv=pa*pb/(pa+pb)^2/(pa+pb+1)
      lJ=log(pa)+log(pb)-3*log(pa+pb)-2*log(pa+pb+1)
      dbeta(pmin(1,vv/maxvar),shape1=var.prior.a,shape2=var.prior.b,log=TRUE)+lJ
    }
  }
  lmarg.posterior=function(pa,pb) lprior(pa,pb)+sum(lbeta(pa+a,pb+b)-lbeta(pa,pb))
  ltrans.marg.posterior=function(logitmu,logsize) {
    pa=exp(logitmu+logsize)/(exp(logitmu)+1)
    pb=exp(logsize)/(exp(logitmu)+1)
    lJ=(logitmu+2*logsize)-2*log1p(exp(logitmu))
    lmarg.posterior(pa,pb)+lJ
  }
  to.ab=function(logitmu,logsize) list(a=unname(exp(logitmu+logsize)/(exp(logitmu)+1)),
                                       b=unname(exp(logsize)/(exp(logitmu)+1)))

  # find MAP
  mu=a/(a+b)
  size.point=log(min(mean(mu)*(1-mean(mu))/var(mu)-1,min(a+b)))
  opt=optim(c(logitmu=log(mean(a/b)),logsize=size.point),function(v) ltrans.marg.posterior(v[1],v[2]),control = list(fnscale=-1))
  logMAP=opt$par
  re=c(to.ab(logMAP[1],logMAP[2]),MAP=opt$value)
  if (compute.marginal.likelihood || compute.grid) {
    mu.low=uniroot(function(x) ltrans.marg.posterior(x,logMAP[2])-opt$value+log(fak.below.max),interval = c(logMAP[1]-4,logMAP[1]),extendInt = "upX")$root
    mu.high=uniroot(function(x) ltrans.marg.posterior(x,logMAP[2])-opt$value+log(fak.below.max),interval = c(logMAP[1],logMAP[1]+4),extendInt = "downX")$root

    size.low=uniroot(function(x) ltrans.marg.posterior(logMAP[1],x)-opt$value+log(fak.below.max),interval = c(logMAP[2]-4,logMAP[2]),extendInt = "upX")$root
    size.high=uniroot(function(x) ltrans.marg.posterior(logMAP[1],x)-opt$value+log(fak.below.max),interval = c(logMAP[2],logMAP[2]+4),extendInt = "downX")$root

    if (compute.marginal.likelihood) {
      ml=log(cubature::adaptIntegrate(function(v) exp(ltrans.marg.posterior(v[1],v[2])-opt$value),lowerLimit = c(mu.low,size.low),upperLimit = c(mu.high,size.high))$integral)+opt$value
      re=c(re,list(mar.loglik=ml))
    }

    if (compute.grid) {
      mu.grid=seq(mu.low,mu.high,length.out=res)
      size.grid=seq(size.low,size.high,length.out=res)
      grid=sapply(mu.grid,function(mu) sapply(size.grid, function(size)  exp(ltrans.marg.posterior(mu,size)-opt$value) ))
      colnames(grid)=mu.grid
      rownames(grid)=size.grid
      dimnames(grid)=c("logit mu","log size")
      re=c(re,list(grid=grid))
      #image(t(grid))
      #contour(t(grid),levels = (seq(0.05,0.95,by=0.05)))
    }
  }
  re
}


bayesian.prior=function(a,b,mean.var=0.004,sd.var=0.001) {
  gamma.k=mean.var^2/sd.var^2
  gamma.theta=sd.var^2/mean.var

  lprior=function(pa,pb) {
    vv=pa*pb/(pa+pb)^2/(pa+pb+1)
    lJ=log(pa)+log(pb)-3*log(pa+pb)-2*log(pa+pb+1)
    dgamma(vv,shape=gamma.k,scale=gamma.theta,log=TRUE)+lJ
  }
  lmarg.posterior=function(pa,pb) lprior(pa,pb)+sum(lbeta(pa+a,pb+b)-lbeta(pa,pb))
  ltrans.marg.posterior=function(logitmu,logsize) {
    pa=exp(logitmu+logsize)/(exp(logitmu)+1)
    pb=exp(logsize)/(exp(logitmu)+1)
    lJ=(logitmu+2*logsize)-2*log1p(exp(logitmu))
    lmarg.posterior(pa,pb)+lJ
  }
  to.ab=function(logitmu,logsize) c(a=unname(exp(logitmu+logsize)/(exp(logitmu)+1)),
                                    b=unname(exp(logsize)/(exp(logitmu)+1)))

  # find MAP
  logMAP=optim(c(logitmu=log(mean(a/b)),logsize=size.point),function(v) ltrans.marg.posterior(v[1],v[2]),control = list(fnscale=-1))$par
  to.ab(logMAP[1],logMAP[2])
}


plotgrid=function() {
  # compute grid
  mu=a/(a+b)
  size.point=log(min(mean(mu)*(1-mean(mu))/var(mu)-1,min(a+b)))
  murange=c(log(mean(mu))-log(mu.fac),log(mean(mu))+log(mu.fac))
  sizerange=c(size.point-log(size.fac),size.point+log(size.fac))



}



