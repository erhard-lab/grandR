


f.old.equi=function(t,s,d) s/d*exp(-t*d)
f.old.nonequi=function(t,f0,s,d) f0*exp(-t*d)
f.new=function(t,s,d) s/d*(1-exp(-t*d))

res.fun.equi=function(par,old,new) {
    s=par[1]
    d=par[2]
    f=function(t) f.old.equi(t,s,d)
    g=function(t) f.new(t,s,d)
    c(old$Value-f(old$time),new$Value-g(new$time))
}

res.fun.nonequi=function(par,old,new) {
    s=par[1]
    d=par[2]
    f0=par[3]
    f=function(t) f.old.nonequi(t,f0,s,d)
    g=function(t) f.new(t,s,d)
    c(old$Value-f(old$time),new$Value-g(new$time))
}
res.fun.nonequi.f0fixed=function(par,f0,old,new) {
    s=par[1]
    d=par[2]
    f=function(t) f.old.nonequi(t,f0,s,d)
    g=function(t) f.new(t,s,d)
    c(old$Value-f(old$time),new$Value-g(new$time))
}


sd.from.hessian=function(H) {
    m=-H
    keep=rep(T,nrow(m))
    while(sum(keep>0) && rcond(m[keep,keep,drop=FALSE])<.Machine$double.eps) {
        s=apply(abs(m),1,sum)
        keep=keep&s>min(s[keep])
    }
    d=rep(0,nrow(m))
    d[keep]=sqrt(diag(solve(m[keep,keep])))
    d
}
logLik.nls.lm <- function(object, REML = FALSE, ...)
{
    res <- object$fvec
    N <- length(res)
    val <-  -N * (log(2 * pi) + 1 - log(N) + log(sum(res^2)))/2
    ## the formula here corresponds to estimating sigma^2.
    attr(val, "df") <- 1L + length(coef(object))
    attr(val, "nobs") <- attr(val, "nall") <- N
    class(val) <- "logLik"
    val
}
FitKineticsGeneLeastSquares=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU,steady.state=NULL,group=Design$Condition,conf.int=0.95,use.old=TRUE,use.new=TRUE,return.vector=FALSE,return.fields=c("Synthesis","Half-life","rmse"),return.extra=function(p) NULL, maxiter=100, compute.modifier=FALSE) {
    if (length(group)>1) stop("Can only group using a single column!")
    if (is.null(group) || is.na(group)) group=paste0(names(ColData(data)),collapse="_")

    correct=function(s) {
        if (max(s$Value)==0) s$Value[s$time==1]=0.01
        s$Type="New"
        s
    }

    stopifnot(time %in% names(data$coldata))

    newdf=GetData(data,mode.slot=paste0("new.",slot),genes=gene,ntr.na = FALSE)
    newdf$use=1:nrow(newdf) %in% (1:nrow(newdf))[use.new]
    newdf$time=newdf[[time]]
    if (is.null(newdf[[group]])) {
        newdf[[group]]="Data"
        if (length(steady.state)==1) names(steady.state)="Data"
    }
    newdf=dlply(newdf,group,function(s) correct(s))

    correct=function(s) {
        if (max(s$Value)==0) s$Value=0.01
        s$Type="Old"
        s
    }
    olddf=GetData(data,mode.slot=paste0("old.",slot),genes=gene,ntr.na = FALSE)
    olddf$use=1:nrow(olddf) %in% (1:nrow(olddf))[use.old]
    olddf$time=olddf[[time]]
    if (is.null(olddf[[group]])) olddf[[group]]="Data"
    olddf=dlply(olddf,group,function(s) correct(s))


    fit.equi=function(ndf,odf) {

        tinit=min(ndf$time[ndf$Value>0])
        init.d=mean(-log(1-(0.1+ndf[ndf$time==tinit,"Value"])/(0.2+ndf[ndf$time==tinit,"Value"]+odf[odf$time==tinit,"Value"])))
        init.s=init.d*mean(odf[odf$time==tinit,"Value"])
        model.p=minpack.lm::nls.lm(c(init.s,init.d),lower=c(0,0.01),
                       fn=res.fun.equi,
                       old=odf[odf$use,],
                       new=ndf[ndf$use,],
                       control=minpack.lm::nls.lm.control(maxiter = maxiter))

        if (model.p$niter==maxiter) return(list(data=NA,modifier=NA,param=c(NA,NA),conf.lower=c(NA,NA),conf.upper=c(NA,NA),f0=NA,logLik=NA,rmse=NA, rmse.new=NA, rmse.old=NA,total=NA,type="equi"))
        conf.p=try(sd.from.hessian(-model.p$hessian)*qnorm(1-(1-conf.int)/2),silent=TRUE)
        par=setNames(model.p$par,c("s","d"))
        rmse=sqrt(model.p$deviance/(nrow(ndf)+nrow(odf)))
        fvec=res.fun.equi(par,odf,ndf)
        n=nrow(odf)
        rmse.old=sqrt(sum(fvec[1:n]^2)/n)
        rmse.new=sqrt(sum(fvec[(n+1):(n*2)]^2)/n)
        df=droplevels(rbind(odf[odf$use,],ndf[ndf$use,]))


        modifier=NA
        if (compute.modifier) {
            s=par["s"]
            d=par["d"]
            # solve f.old.equi(t)=old for t! ( and the same for)
            #time=daply(df,.(Name),function(sub) 1/d * log( 1+ sum(sub$Value[sub$Type=="New"])/sum(sub$Value[sub$Type=="Old"]) ) )
            # solve f.old.equi(t)=old for t!
            time=daply(df,.(Name),function(sub) -1/d* mean(c( log(sub$Value[sub$Type=="Old"]*d/s), log(1-sub$Value[sub$Type=="New"]*d/s)   )) )
            # new+old=s/d, so for a particular time point, multiply such that the sum is s/d
            norm.fac=daply(df,.(Name),function(sub) s/d / sum(sub$Value))
            modifier=data.frame(Name=names(norm.fac),Time=time,Norm.factor=norm.fac)
        }

        list(data=df,
             modifier=modifier,
             param=par,
             conf.lower=par-conf.p,
             conf.upper=par+conf.p,
             f0=unname(par['s']/par['d']),
             logLik=logLik.nls.lm(model.p),
             rmse=rmse, rmse.new=rmse.new, rmse.old=rmse.old,
             rmse.linear=sqrt(sum(ddply(df,c(group,"Type"),function(s) data.frame(se=s$Value-mean(s$Value))^2)$se)/nrow(df)),
             total=sum(ndf$Value)+sum(odf$Value),type="equi")
    }
    fit.nonequi=function(ndf,odf) {
        oind=union(which(odf$time==0),which(odf$use))
        nind=union(which(ndf$time==0),which(ndf$use))

        tinit=min(ndf$time[ndf$Value>0])
        init.d=mean(-log(1-(0.1+ndf[ndf$time==tinit,"Value"])/(0.2+ndf[ndf$time==tinit,"Value"]+odf[odf$time==tinit,"Value"])))
        init.s=max(ndf$Value[ndf$time>0]/ndf$time[ndf$time>0])
        f0=mean(odf[odf$time==0,"Value"])
        model.m=minpack.lm::nls.lm(c(init.s,init.d,f0),lower=c(0,0.01,0),
                       fn=res.fun.nonequi,
                       old=odf[oind,],
                       new=ndf[nind,],
                       control=minpack.lm::nls.lm.control(maxiter = maxiter))
        if (model.m$niter==maxiter) return(list(data=NA,modifier=NA,param=c(NA,NA),conf.lower=c(NA,NA),conf.upper=c(NA,NA),f0=NA,logLik=NA,rmse=NA, rmse.new=NA, rmse.old=NA, rmse.linear=NA, total=NA,type="non.equi"))
        par=setNames(model.m$par,c("s","d"))
        f0=par[3]
        conf.m=try(sd.from.hessian(-model.m$hessian)*qnorm(1-(1-conf.int)/2),silent=TRUE)
        rmse=sqrt(model.m$deviance/(nrow(ndf)+nrow(odf)))
        fvec=res.fun.nonequi.f0fixed(par,f0,odf,ndf)
        n=nrow(odf)
        rmse.old=sqrt(sum(fvec[1:n]^2)/n)
        rmse.new=sqrt(sum(fvec[(n+1):(n*2)]^2)/n)
        df=droplevels(rbind(odf[oind,],ndf[nind,]))
        rmse.linear=sqrt(sum(ddply(df,c(group,"Type"),function(s) data.frame(se=s$Value-mean(s$Value))^2)$se)/nrow(df))

        modifier=NA
        if (compute.modifier) {
            s=par["s"]
            d=par["d"]
            f0=unname(f0)
            # solve f.old.nonequi(t)/f.new(t)=old/new for t!
            #time=daply(df,.(Name),function(sub) 1/d * log( 1+ (sum(sub$Value[sub$Type=="New"])*f0*d)/(sum(sub$Value[sub$Type=="Old"])*s) ) )
            # solve f.old.equi(t)=old for t!
            time=daply(df,.(Name),function(sub) -1/d* mean(c( log(sub$Value[sub$Type=="Old"]/f0), log(1-sub$Value[sub$Type=="New"]*d/s)   )) )
            # new+old=f.old.nonequi(t)+f.new(t), so for the computed time point t, multiply such that the sum is f.old.nonequi(t)+f.new(t)
            norm.fac=daply(df,.(Name),function(sub) {
                t=time[as.character(sub$Name[1])]
                (f.old.nonequi(t,f0,s,d)+f.new(t,s,d)) / sum(sub$Value)
            })

            modifier=data.frame(Name=names(norm.fac),Time=time,Norm.factor=norm.fac)
        }

        list(data=df,
             modifier=modifier,
             param=par,
             conf.lower=par-conf.m,
             conf.upper=par+conf.m,
             f0=unname(f0),
             logLik=logLik.nls.lm(model.m),
             rmse=rmse, rmse.new=rmse.new, rmse.old=rmse.old,
             rmse.linear=rmse.linear,
             total=sum(ndf$Value)+sum(odf$Value),type="non.equi")
    }

    fits=lapply(names(newdf),function(cond) {
        equi=if (is.null(steady.state)) {
            TRUE
            } else if(length(steady.state)==1) {
            steady.state;
            } else is.na(unlist(steady.state)[cond]) || as.logical(unlist(steady.state)[cond])
        ndf=newdf[[cond]]
        odf=olddf[[cond]]
        if (equi) fit.equi(ndf,odf) else fit.nonequi(ndf,odf)
    })
    if (!is.null(data$coldata[[group]])) names(fits)=names(newdf)
    if (return.vector) fits=unlist(lapply(fits,function(p) {
        p["Synthesis"]=unname(p$param["s"])
        p["Half-life"]=unname(log(2)/p$param["d"])
        c(p[return.fields],return.extra(p))
        }))
    if (is.null(data$coldata[[group]]) && !return.vector) fits=fits[[1]]
    fits
}

FitKineticsGeneLogSpaceLinear=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU,group=Design$Condition,conf.int=0.95,return.vector=FALSE,return.fields=c("Synthesis","Half-life","rmse"),return.extra=function(p) NULL) {
    if (length(group)>1) stop("Can only group using a single column!")
    if (is.null(group) || is.na(group)) group=paste0(names(ColData(data)),collapse="_")

    correct=function(s) {
        if (max(s$Value)==0) s$Value[s$time==1]=0.01
        s$Type="New"
        s
    }

    stopifnot(time %in% names(data$coldata))

    correct=function(s) {
        if (max(s$Value)==0) s$Value=0.01
        s$Type="Old"
        s
    }
    olddf=GetData(data,mode.slot=paste0("old.",slot),genes=gene,ntr.na = FALSE)
    olddf$use=1:nrow(olddf) %in% (1:nrow(olddf))[use.old]
    olddf$time=olddf[[time]]
    if (is.null(olddf[[group]])) olddf[[group]]="Data"
    olddf=dlply(olddf,group,function(s) correct(s))


    fit.lm=function(odf) {
        odf=odf[odf$Value>0,]
        fit=lm(log(Value)~time,data=odf)
        summ=summary(fit)

        par=setNames(c(exp(coef(fit)["(Intercept)"])*-coef(fit)["time"],-coef(fit)["time"]),c("s","d"))
        conf.p=confint(fit,level=conf.int)
        conf.p=apply(conf.p,2,function(v) setNames(pmax(0,c(exp(v["(Intercept)"])*par['d'],-v["time"])),c("s","d")))

        modifier=NA

        list(data=odf,
             param=par,
             conf.lower=unname(c(conf.p[1,1],conf.p[2,2])),
             conf.upper=unname(c(conf.p[1,2],conf.p[2,1])),
             f0=unname(par['s']/par['d']),
             logLik=logLik(fit),
             adj.r.squared=summ$adj.r.squared,
             total=sum(odf$Value),type="lm")
    }

    fits=lapply(names(olddf),function(cond) {
        odf=olddf[[cond]]
        fit.lm(odf)
    })
    if (!is.null(data$coldata[[group]])) names(fits)=names(olddf)
    if (return.vector) fits=unlist(lapply(fits,function(p) {
        p["Synthesis"]=unname(p$param["s"])
        p["Half-life"]=unname(log(2)/p$param["d"])
        c(p[return.fields],return.extra(p))
    }))
    if (is.null(data$coldata[[group]]) && !return.vector) fits=fits[[1]]
    fits
}


FitKineticsGeneNtr=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU,group=Design$Condition,conf.int=0.95,return.vector=FALSE,total.fun=median,return.fields=c("Synthesis","Half-life","rmse"),return.extra=function(p) NULL) {
    if (length(group)>1) stop("Can only group using a single column!")
    if (is.null(group) || is.na(group)) group=paste0(names(ColData(data)),collapse="_")
    if (!all(c("alpha","beta") %in% Slots(data))) stop("Beta approximation data is not available in grandR object!")

    crit <- qchisq(1-conf.int, df = 2, lower.tail = FALSE) / 2

    bounds=c(log(2)/48,log(2)/0.01)

    df=GetData(data,mode.slot=paste0("total.",slot),genes=gene)
    if (is.null(df[[group]])) df[[group]]="Data"
    total=dlply(df,group,function(s) total.fun(s$Value))

    a=GetData(data,mode.slot="alpha",genes=gene)
    b=GetData(data,mode.slot="beta",genes=gene)
    ntr=GetData(data,mode.slot="ntr",genes=gene)
    if (is.null(a[[group]])) a[[group]]="Data"
    if (is.null(b[[group]])) b[[group]]="Data"
    if (is.null(ntr[[group]])) ntr[[group]]="Data"

    #sloglik=function(d,a,b,t) (a-1)*log(1-exp(-t*d))-t*d*b
    loglik=function(d,a,b,t) sum((a-1)*log(1-exp(-t*d))-t*d*b)
    #sloglik=function(d,a,b,t) (a-1)*log1p(-exp(-t*d))-t*d*(b-1)
    #loglik=function(d,a,b,t) sum((a-1)*log1p(-exp(-t*d))-t*d*(b-1))

    t=a[,time]
    use=t>0
    c=droplevels(a[[group]][use])
    a=a$Value[use]
    b=b$Value[use]
    ntr=ntr$Value[use]
    t=t[use]
    uniroot.save=function(fun,lower,upper) if (fun(lower)*fun(upper)>=0) mean(lower,upper) else uniroot(fun,lower=lower,upper=upper)$root
    fits=lapply(levels(c),function(cc) {
        ind=c==cc
        if(any(is.nan(c(a[ind],b[ind])))) return(NaN)
        ploglik=function(x) loglik(x,a=a[ind],b=b[ind],t=t[ind])
        d=optimize(ploglik,bounds,maximum=T)$maximum
        max=ploglik(d)
        lower=uniroot.save(function(x) ploglik(x)-max+crit,lower = bounds[1],upper=d)
        upper=uniroot.save(function(x) ploglik(x)-max+crit,lower = d,upper=bounds[2])

        rmse=sum(sqrt(((1-exp(-t[ind]*d)-ntr[ind]))^2))/sum(ind)

        list(param=c(s=total[[cc]]*d,d=d),
             conf.lower=c(s=total[[cc]]*lower,d=lower),
             conf.upper=c(s=total[[cc]]*upper,d=upper),
             f0=unname(total[[cc]]),
             logLik=max,
             rmse=rmse,
             total=total[[cc]],
             type="ntr.fit")
    })
    if (!is.null(data$coldata[[group]])) names(fits)=levels(c)
    if (return.vector) fits=unlist(lapply(fits,function(p) {
        p["Synthesis"]=unname(p$param["s"])
        p["Half-life"]=unname(log(2)/p$param["d"])
        c(p[return.fields],return.extra(p))
    }))
    if (is.null(data$coldata[[group]]) && !return.vector) fits=fits[[1]]
    fits

}

NormalizeKinetic=function(data,slot=DefaultSlot(data),time=Design$dur.4sU,norm.name="kinetic",time.name="norm_time",group=Design$Condition,steady.state=NULL,n.estimate=1000,set.to.default=TRUE) {

    conds=ColData(data)
    if (is.null(conds[[group]])) {
        conds[[group]]==factor("Data")
        if (length(steady.state)==1) names(steady.state)="Data"
    }
    timecol=setNames(rep(NA,nrow(conds)),conds$Name)
    norm.df=list()
    for (cond in levels(conds[[group]])) {
        sub=subset(data,columns=conds[[group]]==cond)

        # restrict to top n genes
        totals=rowSums(GetTable(sub,type=slot))
        threshold=sort(totals,decreasing = TRUE)[min(length(totals),n.estimate)]

        # fit models
        fits=opt$lapply(Genes(sub),FitKineticsGeneLeastSquares,data=sub,time=time,group=NULL,steady.state = steady.state[[cond]],compute.modifier=TRUE)
        HL=log(2)/sapply(fits,function(fit) fit$param['d'])

        # extract modifiers
        times=t(sapply(fits[Genes(sub) %in% names(which(totals>=threshold))],function(fit)fit$modifier[,"Time"]))
        colnames(times)=fits[[1]]$modifier$Name

        norm.factor=t(sapply(fits[Genes(sub) %in% names(which(totals>=threshold))],function(fit)fit$modifier[,"Norm.factor"]))
        colnames(norm.factor)=fits[[1]]$modifier$Name

        # just the median for the time!
        timecol[colnames(times)]=apply(times,2,function(v) median(v[!is.nan(v)]))

        # spline regression for the normalization factors
        for (i in 1:ncol(norm.factor)) {
            df=data.frame(x=log10(HL[Genes(sub) %in% names(which(totals>=threshold))]),y=log10(norm.factor[,i]))
            df=df[is.finite(df$y) & is.finite(df$x),]
            spl=smooth.spline(df$x,df$y,nknots = 5)
            norm.df[[colnames(norm.factor)[i]]]=GetTable(data,type=slot,columns=colnames(norm.factor)[i],name.by = "Gene")*10^predict(spl,log10(HL))$y
        }

    }

    norm.df=as.matrix(as.data.frame(norm.df)[,ColData(data)$Name])

    data=AddSlot(data,norm.name,norm.df)
    if (set.to.default) DefaultSlot(data)=norm.name
    timecol[ColData(data)[[time]]==0]=0
    data=ColData(data,time.name,timecol)

    invisible(data)
}

FitKinetics=function(data,name="kinetics",type=c("full","ntr","lm"),...) {
    slam.param=as.data.frame(t(opt$sapply(data$gene.info$Gene,
                                   switch(substr(tolower(type[1]),1,1),n=FitKineticsGeneNtr,f=FitKineticsGeneLeastSquares,l=FitKineticsGeneLogSpaceLinear),
                                   data=data,return.vector=TRUE,...)))
    AddAnalysis(data,MakeAnalysis(name=name,analysis="FitKinetics"),slam.param)
}


PlotGeneKinetics=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU,title=data$gene.info[ToIndex(data,gene),"Symbol"], group=Design$Condition, type=c("full","ntr","lm"), bare.plot=FALSE,exact.tics=TRUE,return.tables=FALSE,...) {
    if (length(ToIndex(data,gene))==0) return(NULL)

    fit=switch(substr(tolower(type[1]),1,1),
               n=FitKineticsGeneNtr(data,gene,slot=slot,time=time,group=group,...),
               f=FitKineticsGeneLeastSquares(data,gene,slot=slot,time=time,group=group,...),
               l=FitKineticsGeneLogSpaceLinear(data,gene,slot=slot,time=time,group=group,...)
    )
    if (is.null(data$coldata[[group]])) fit=setNames(list(fit),gene)
    df=rbind(
        cbind(GetData(data,mode.slot=paste0("total.",slot),genes=gene),Type="Total"),
        cbind(GetData(data,mode.slot=paste0("new.",slot),genes=gene,ntr.na = FALSE),Type="New"),
        cbind(GetData(data,mode.slot=paste0("old.",slot),genes=gene,ntr.na = FALSE),Type="Old")
    )
    if (substr(tolower(type[1]),1,1)=="n") {
        fac=unlist(lapply(as.character(df[[group]]),function(n) fit[[n]]$total))/df$Value[df$Type=="Total"]
        df$Value=df$Value*fac
    }

    df$time=df[[time]]
    if (is.data.frame(fit[[1]]$modifier)) {
        df$time=sapply(1:nrow(df),function(i) fit[[as.character(df$Condition)[i]]]$modifier[as.character(df$Name)[i],"Time"])
        df$Value=df$Value*sapply(1:nrow(df),function(i) fit[[as.character(df$Condition)[i]]]$modifier[as.character(df$Name)[i],"Norm.factor"])
    }
    df$Condition=if (group %in% names(df)) df[[group]] else gene
    tt=seq(0,max(df$time),length.out=100)
    df.median=ddply(df,c("Condition","Type","duration.4sU","time"),function(s) data.frame(Value=median(s$Value)))
    fitted=ldply(fit,function(f) data.frame(time=c(tt,tt),Value=c(f.old.nonequi(tt,f$f0,f$param["s"],f$param['d']),f.new(tt,f$param["s"],f$param['d'])),Type=rep(c("Old","New"),each=length(tt))),.id="Condition")
    breaks=if (exact.tics) sort(unique(df[[time]])) else scales::breaks_extended(5)(df[[time]])

    g=ggplot(df,aes(time,Value,color=Type))+
        geom_point()+
        geom_line(data=df.median[df.median$Type=="Total",])+
        scale_x_continuous(if (bare.plot) NULL else "4sU labeling",labels = scales::number_format(accuracy = max(0.01,my.precision(breaks)),suffix="h"),breaks=breaks)+
        scale_color_manual("RNA",values=c(Total="gray",New="red",Old="blue"),guide=if (bare.plot) "none" else "legend")+
        ylab(if (bare.plot) NULL else "Expression")+
        geom_line(data=fitted,aes(ymin=NULL,ymax=NULL),linetype=2)
    if (!is.null(data$coldata[[group]])) g=g+facet_wrap(~Condition,nrow=1)
    if (!bare.plot) g=g+ggtitle(title)
    if (return.tables) list(gg=g,df=df,df.median=df.median,fitted=fitted) else g
}




SimulateKinetics=function(s=100*d,d=log(2)/hl,hl=2,l0=s/d,times=seq(min.time,max.time,by=by),min.time=-1,max.time=10,length.out = 1000,by = ((max.time - 0)/(length.out - 1)),name=NULL,out=c("Old","New","Total","NTR")) {
    ode.new=function(t,s,d) ifelse(t<0,0,s/d*(1-exp(-t*d)))
    ode.old=function(t,f0,s,d) ifelse(t<0,s/d,f0*exp(-t*d))
    old=ode.old(times,l0,s,d)
    new=ode.new(times,s,d)
    re=data.frame(
        Time=times,
        Value=c(old,new,old+new,new/(old+new)),
        Type=factor(rep(c("Old","New","Total","NTR"),each=length(times)),levels=c("Old","New","Total","NTR"))
    )
    if (!is.null(name)) re$Name=name
    re=re[re$Type %in% out,]
    rownames(re)=1:nrow(re)
    re
}


PlotSimulation=function(sim.df,ntr=TRUE,old=TRUE,new=TRUE,total=TRUE) {
    if (!ntr) sim.df=sim.df[sim.df$Type!="NTR",]
    if (!old) sim.df=sim.df[sim.df$Type!="Old",]
    if (!new) sim.df=sim.df[sim.df$Type!="New",]
    if (!total) sim.df=sim.df[sim.df$Type!="Total",]
    ggplot(sim.df,aes(Time,Value,color=Type))+
        geom_line(size=1)+
        scale_color_manual(NULL,values=c(Old="#54668d",New="#953f36",Total="#373737",NTR="#e4c534"))+
        facet_wrap(~ifelse(Type=="NTR","NTR","Timecourse"),scales="free_y",ncol=1)+
        ylab(NULL)+
        scale_x_continuous(breaks=scales::pretty_breaks())+
        theme(
	  strip.background = element_blank(),
	  strip.text.x = element_blank()
	)
}

PlotCompareNTRs=function(...) {
    dfs=list(...)
    df=do.call("rbind",dfs)
    ggplot(df[df$Type=="NTR",],aes(Time,Value,color=Name))+
        geom_line(size=1)
}

#sim.hl8.steady=Simulate(hl=8,name="Steady-state, 8h")
#sim.hl4.steady=Simulate(hl=4,name="Steady-state, 4h")
#sim.hl2.steady=Simulate(hl=2,name="Steady-state, 2h")
#sim.hl1.steady=Simulate(hl=1,name="Steady-state, 1h")
#sim.hl8.2x=Simulate(l0=10,hl=4,s=10*log(2)/8,name="2x down, 8h")
#sim.hl4.2x=Simulate(l0=10,hl=2,s=10*log(2)/4,name="2x down, 4h")
#sim.hl2.2x=Simulate(l0=10,hl=1,s=10*log(2)/2,name="2x down, 2h")
#sim.hl1.2x=Simulate(l0=10,hl=0.5,s=10*log(2)/1,name="2x down, 1h")
#sim.hl8.10x=Simulate(l0=10,hl=0.8,s=10*log(2)/8,name="10x down, 8h")
#sim.hl4.10x=Simulate(l0=10,hl=0.4,s=10*log(2)/4,name="10x down, 4h")
#sim.hl2.10x=Simulate(l0=10,hl=0.2,s=10*log(2)/2,name="10x down, 2h")
#sim.hl1.10x=Simulate(l0=10,hl=0.1,s=10*log(2)/1,name="10x down, 1h")



#pdf("destabilized.pdf",width=10,height=10)
#plot_grid(
#    PlotSimulation(sim.hl8.steady),
#    PlotSimulation(sim.hl8.2x),
#    PlotSimulation(sim.hl8.10x),
#    PlotCompareNTRs(sim.hl8.steady,sim.hl8.2x,sim.hl8.10x),
#    ncol=2
#)


#plot_grid(
#    PlotSimulation(sim.hl4.steady),
#    PlotSimulation(sim.hl4.2x),
#    PlotSimulation(sim.hl4.10x),
#    PlotCompareNTRs(sim.hl4.steady,sim.hl4.2x,sim.hl4.10x),
#    ncol=2
#)

#plot_grid(
#    PlotSimulation(sim.hl2.steady),
#    PlotSimulation(sim.hl2.2x),
#    PlotSimulation(sim.hl2.10x),
#    PlotCompareNTRs(sim.hl2.steady,sim.hl2.2x,sim.hl2.10x),
#    ncol=2
#)

#plot_grid(
#    PlotSimulation(sim.hl1.steady),
#    PlotSimulation(sim.hl1.2x),
#    PlotSimulation(sim.hl1.10x),
#    PlotCompareNTRs(sim.hl1.steady,sim.hl1.2x,sim.hl1.10x),
#    ncol=2
#)

#dev.off()
