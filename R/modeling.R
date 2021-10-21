


f.old.equi=function(t,s,d) s/d*exp(-t*d)
f.new=function(t,s,d) s/d*(1-exp(-t*d))
f.old.nonequi=function(t,f0,s,d) f0*exp(-t*d)

res.fun.equi=function(par,old,new) {
    s=par[1]
    d=par[2]
    f=function(t) f.old.equi(t,s,d)
    g=function(t) f.new(t,s,d)
    c(old$Value-f(old$time),new$Value-g(new$time))
}

res.fun.nonequi=function(par,f0,old,new) {
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
FitKineticsGeneLeastSquares=function(data,gene,slot="norm",steady.state=NULL,conf.int=0.95,use.old=TRUE,use.new=TRUE,return.vector=FALSE,return.fields=c("rmse"),return.extra=function(p) NULL, maxiter=100) {
    correct=function(s) {
        if (max(s$Value)==0) s$Value[s$time==1]=0.01
        s$Type="New"
        s
    }
    
    stopifnot(Design$dur.4sU %in% names(data$coldata))
    
    newdf=GetData(data,type=paste0("new.",slot),gene=gene,keep.ntr.na = FALSE)
    newdf$time=newdf[[Design$dur.4sU]]
    if (is.null(newdf$Condition)) newdf$Condition="Data"
    newdf=dlply(newdf,.(Condition),function(s) correct(s))
    
    correct=function(s) {
        if (max(s$Value)==0) s$Value=0.01
        s$Type="Old"
        s
    }
    olddf=GetData(data,type=paste0("old.",slot),gene=gene,keep.ntr.na = FALSE)
    olddf$time=olddf[[Design$dur.4sU]]
    if (is.null(olddf$Condition)) olddf$Condition="Data"
    olddf=dlply(olddf,.(Condition),function(s) correct(s))
    
    
    fit.equi=function(ndf,odf) {
        oind=union(which(odf$time==0),(1:nrow(odf))[use.old])
        nind=union(which(ndf$time==0),(1:nrow(ndf))[use.new])
        
        tinit=min(ndf$time[ndf$Value>0])
        init.d=mean(-log(1-(0.1+ndf[ndf$time==tinit,"Value"])/(0.2+ndf[ndf$time==tinit,"Value"]+odf[odf$time==tinit,"Value"])))
        init.s=init.d*mean(odf[odf$time==tinit,"Value"])
        model.p=nls.lm(c(init.s,init.d),lower=c(0,0.01),
                       fn=res.fun.equi,
                       old=odf[oind,],
                       new=ndf[nind,],
                       control=nls.lm.control(maxiter = maxiter))
        
        if (model.p$niter==maxiter) return(list(data=NA,param=c(NA,NA),conf.lower=c(NA,NA),conf.upper=c(NA,NA),f0=NA,loglik=NA,rmse=NA, rmse.new=NA, rmse.old=NA,total=NA,type="equi"))
        conf.p=try(sd.from.hessian(-model.p$hessian)*qnorm(1-(1-conf.int)/2),silent=TRUE)
        par=setNames(model.p$par,c("s","d"))
        rmse=sqrt(model.p$deviance/(nrow(ndf)+nrow(odf)))
        fvec=res.fun.equi(par,odf,ndf)
        n=nrow(odf)
        rmse.old=sqrt(sum(fvec[1:n]^2)/n)
        rmse.new=sqrt(sum(fvec[(n+1):(n*2)]^2)/n)
        df=rbind(odf[oind,],ndf[nind,])
        list(data=df,
             param=par,
             conf.lower=par-conf.p,
             conf.upper=par+conf.p,
             f0=unname(par['s']/par['d']),
             loglik=logLik.nls.lm(model.p),
             rmse=rmse, rmse.new=rmse.new, rmse.old=rmse.old,
             rmse.linear=sqrt(sum(ddply(df,.(Sample,Type),function(s) data.frame(se=s$Value-mean(s$Value))^2)$se)/nrow(df)),
             total=sum(ndf$Value)+sum(odf$Value),type="equi")
    }
    fit.nonequi=function(ndf,odf) {
        oind=union(which(odf$time==0),(1:nrow(odf))[use.old])
        nind=union(which(ndf$time==0),(1:nrow(ndf))[use.new])
   
        tinit=min(ndf$time[ndf$Value>0])
        init.d=mean(-log(1-(0.1+ndf[ndf$time==tinit,"Value"])/(0.2+ndf[ndf$time==tinit,"Value"]+odf[odf$time==tinit,"Value"])))
        init.s=max(mean(ndf[ndf$time==1,"Value"]),mean(ndf[ndf$time==2,"Value"])/2,mean(ndf[ndf$time==4,"Value"])/4)
        f0=mean(odf[odf$time==0,"Value"])
        model.m=nls.lm(c(init.s,init.d),lower=c(0,0.01),
                       fn=res.fun.nonequi,
                       f0=f0,
                       old=odf[oind,],
                       new=ndf[nind,],
                       control=nls.lm.control(maxiter = maxiter))
        if (model.m$niter==maxiter) return(list(data=NA,param=c(NA,NA),conf.lower=c(NA,NA),conf.upper=c(NA,NA),f0=NA,loglik=NA,rmse=NA, rmse.new=NA, rmse.old=NA, total=NA,type="non.equi"))
        par=setNames(model.m$par,c("s","d"))
        conf.m=try(sd.from.hessian(-model.m$hessian)*qnorm(1-(1-conf.int)/2),silent=TRUE)
        rmse=sqrt(model.m$deviance/(nrow(ndf)+nrow(odf)))
        fvec=res.fun.nonequi(par,f0,odf,ndf)
        n=nrow(odf)
        rmse.old=sqrt(sum(fvec[1:n]^2)/n)
        rmse.new=sqrt(sum(fvec[(n+1):(n*2)]^2)/n)
        df=rbind(odf[oind,],ndf[nind,])
        linear.rmse=sqrt(sum(ddply(df,.(Sample,Type),function(s) data.frame(se=s$Value-mean(s$Value))^2)$se)/nrow(df))
        list(data=df,
             param=par,
             conf.lower=par-conf.m,
             conf.upper=par+conf.m,
             f0=unname(f0),
             loglik=logLik.nls.lm(model.m),
             rmse=rmse, rmse.new=rmse.new, rmse.old=rmse.old,
             total=sum(ndf$Value)+sum(odf$Value),type="non.equi")
    }
    
    fits=lapply(names(newdf),function(cond) {
        equi=is.null(steady.state) || is.na(unlist(steady.state)[cond]) || as.logical(unlist(steady.state)[cond])
        ndf=newdf[[cond]]
        odf=olddf[[cond]]
        if (equi) fit.equi(ndf,odf) else fit.nonequi(ndf,odf)
    })
    if (!is.null(data$coldata$Condition)) names(fits)=names(newdf)
    if (return.vector) fits=unlist(lapply(fits,function(p) c(Synthesis=unname(p$param["s"]),`Half-life`=unname(log(2)/p$param["d"]),p[return.fields],return.extra(p))))
    if (is.null(data$coldata$Condition) && !return.vector) fits=fits[[1]]
    fits
}


FitKineticsGeneNtr=function(data,gene,slot="norm",conf.int=0.95,return.vector=FALSE,total.fun=median,return.fields=c("rmse")) {
    
    crit <- qchisq(1-conf.int, df = 2, lower.tail = FALSE) / 2
    
    bounds=c(log(2)/48,log(2)/0.01)
    
    df=GetData(data,type=paste0("total.",slot),gene=gene)
    if (is.null(df$Condition)) df$Condition="Data"
    total=dlply(df,.(Condition),function(s) total.fun(s$Value))
    
    a=GetData(data,type="alpha",gene=gene)
    b=GetData(data,type="beta",gene=gene)
    ntr=GetData(data,type="ntr",gene=gene)
    if (is.null(a$Condition)) a$Condition="Data"
    if (is.null(b$Condition)) b$Condition="Data"
    if (is.null(ntr$Condition)) ntr$Condition="Data"
    
    #sloglik=function(d,a,b,t) (a-1)*log(1-exp(-t*d))-t*d*b
    loglik=function(d,a,b,t) sum((a-1)*log(1-exp(-t*d))-t*d*b)
    #sloglik=function(d,a,b,t) (a-1)*log1p(-exp(-t*d))-t*d*(b-1)
    #loglik=function(d,a,b,t) sum((a-1)*log1p(-exp(-t*d))-t*d*(b-1))
    
    t=a[,Design$dur.4sU]
    use=t>0
    c=droplevels(a$Condition[use])
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
             loglik=max,
             rmse=rmse,
             total=total[[cc]],
             type="ntr.fit")
    })
    if (!is.null(data$coldata$Condition)) names(fits)=levels(c)
    if (return.vector) fits=unlist(lapply(fits,function(p) c(Synthesis=unname(p$param["s"]),`Half-life`=unname(log(2)/p$param["d"]),p[return.fields])))
    if (is.null(data$coldata$Condition) && !return.vector) fits=fits[[1]]
    fits
    
}

DropModeling=function(data) {
    data$modeling=NULL
    invisible(data)
}
AddModeling=function(data,name,table) {
    if (is.null(data$modeling)) data$modeling=list()
    if (is.null(data$modeling[[name]])) {
        data$modeling[[name]]=table
    } else {
        for (n in names(table)) data$modeling[[name]][[n]]=table[[n]]
    }
    rownames(data$modeling[[name]])=data$gene.info$Gene
    invisible(data)
}


FitKinetics=function(data,name="kinetics",ntr.fit=FALSE,...) {
    slam.param=as.data.frame(t(simplify2array(mclapply(data$gene.info$Gene,
                                                       if(ntr.fit)FitKineticsGeneNtr else FitKineticsGeneLeastSquares,
                                                       data=data,mc.cores=max(1,detectCores()-2),return.vector=TRUE,...))))
    AddModeling(data,name,slam.param)
}
    

PlotGeneKinetics=function(data,gene,slot="norm",title=gene, ntr.fit=FALSE,...) {
    if (length(ToIndex(data,gene))==0) return(NULL)

    fit=if (ntr.fit) FitKineticsGeneNtr(data,gene,slot=slot,...) else FitKineticsGeneLeastSquares(data,gene,slot,...)
    df=rbind(
        cbind(GetData(data,type=paste0("total.",slot),gene=gene),Type="Total"),
        cbind(GetData(data,type=paste0("new.",slot),gene=gene,keep.ntr.na = FALSE),Type="New"),
        cbind(GetData(data,type=paste0("old.",slot),gene=gene,keep.ntr.na = FALSE),Type="Old")
    )
    if (ntr.fit) {
        fac=if (is.null(data$coldata$Condition)) fit$total/df$Value[df$Type=="Total"] else unlist(lapply(as.character(df$Condition),function(n) fit[[n]]$total))/df$Value[df$Type=="Total"]
        df$Value=df$Value*fac
    }
    
    df$time=df[[Design$dur.4sU]]
    tt=seq(0,max(df$time),length.out=100)
    df.median=ddply(df,.(Condition,Type,duration.4sU,time),function(s) data.frame(Value=median(s$Value)))
    fitted=ldply(fit,function(f) data.frame(time=c(tt,tt),Value=c(f.old.nonequi(tt,f$f0,f$param["s"],f$param['d']),f.new(tt,f$param["s"],f$param['d'])),Type=rep(c("Old","New"),each=length(tt))),.id="Condition")
    ggplot(df,aes(time,Value,color=Type))+
        geom_point()+
        geom_line(data=df.median[df.median$Type=="Total",])+
        scale_x_continuous("4sU labeling",labels = function(x) sprintf("%dh",x),breaks=sort(unique(df$time)))+
        scale_color_manual("RNA",values=c(Total="gray",New="red",Old="blue"))+
        facet_wrap(~Condition,nrow=1)+
        ylab("Expression")+
        ggtitle(title)+geom_line(data=fitted,aes(ymin=NULL,ymax=NULL),linetype=2)
}




Simulate=function(s=100*d,d=log(2)/hl,hl=2,l0=s/d,min.time=-1,max.time=10,length.out = 1000,by = ((max.time - 0)/(length.out - 1)),name=NULL) {
    ode.new=function(t,s,d) ifelse(t<0,0,s/d*(1-exp(-t*d)))
    ode.old=function(t,f0,s,d) ifelse(t<0,s/d,f0*exp(-t*d))
    t=seq(min.time,max.time,by=by)
    old=ode.old(t,l0,s,d)
    new=ode.new(t,s,d)
    re=data.frame(
        Time=t,
        Value=c(old,new,old+new,new/(old+new)),
        Type=factor(rep(c("Old","New","Total","NTR"),each=length(t)),levels=c("Old","New","Total","NTR"))
    )
    if (!is.null(name)) re$Name=name
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
