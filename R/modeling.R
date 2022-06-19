

#' Fit kinetics using pulseR
#'
#' @param data A grandR object
#' @param name the user defined analysis name to store the results
#' @param time The column in the column annotation table representing the labeling duration
#'
#' @details This is adapted code from https://github.com/dieterich-lab/ComparisonOfMetabolicLabeling
#'
#' @return
#' @export
#'
FitKineticsPulseR=function(data,name="pulseR",time=Design$dur.4sU) {


    conditions=rbind(
        cbind(Coldata(data),data.frame(fraction="labelled")),
        cbind(Coldata(data),data.frame(fraction="unlabelled"))
    )
    conditions$time=conditions[[time]]
    use=conditions$time>0 | conditions$fraction=="unlabelled"
    conditions=conditions[use,]
    if (is.null(conditions$Condition)) conditions$Condition=factor("Data")

    totals=GetTable(data,type="count")
    l=GetTable(data,type="conversion_reads")
    u=totals-l
    counts=cbind(l,u)[,use]
    colnames(counts)=conditions$Name

    norms <- pulseR:::findDeseqFactorsSingle(totals)
    norms=norms[colnames(counts)]

    re=lapply(levels(conditions$Condition),function(n) {

        norms1=norms[conditions$Condition==n]
        counts=counts[,conditions$Condition==n]
        conditions=conditions[conditions$Condition==n,]

        ## fitting options: tolerance and upper/lower bounds for parameters
        tolerance <- list(params = 0.01,
                          logLik = 0.01)
        # assume 20:1 labelling efficiency, neglecting substitution, ratio is approx. the same
        boundaries <- list(mu1  = c(log(1e-2*20), log(1e6*20)), # substituted variable
                           mu2  = c(log(1e-2), log(1e6)),
                           mu3  = c(log(1e-2), log(1e6)),
                           d    = c(1e-3, 2),
                           size = c(1, 1e3))

        ## set initial values here
        init <- function(counts) {
            fit <- list(
                mu1  = log(1e-1 + counts[,1]),
                mu2  = log(1e-1 + counts[,1]),
                mu3  = log(1e-1 + counts[,1]),
                d    = rep(5e-1, length(counts[,1])),
                size = 1e2)
            ## use prior fit
            fit
        }


        formulas=list(
            formulas = pulseR::MeanFormulas(
                unlabelled = exp(mu1) + exp(mu2) + exp(mu3) * (1 + exp(-d * time)),
                labelled = exp(mu2) + exp(mu3) * (1 - exp(-d * time))),
            formulaIndexes = list(
                unlabelled = "unlabelled",
                labelled = "labelled"))

        pd <- pulseR::PulseData(
            counts,
            conditions[c("fraction", "time")], #  data.frame; the first column corresponds to the conditions given in formulas
            formulas$formulas,
            formulas$formulaIndexes,
            groups=~fraction+time
        )
        pd$depthNormalisation <- norms1

        setOpts <- function(bounds, tolerance, cores = parallel::detectCores()-2, replicates = 5, normFactors = NULL) {
            opts <- pulseR::setFittingOptions(verbose = "verbose")
            opts$cores <- cores
            opts$replicates <- replicates
            ## if rt-conversion data, we do not fit normalisation coefficients, because they
            ## are derived from the DESeq-like normalisation
            if (is.null(normFactors)) { opts$fixedNorms <- TRUE }
            opts <- pulseR::setBoundaries(bounds, normFactors = normFactors, options = opts)
            opts <- pulseR::setTolerance(
                params = tolerance$params,
                normFactors = tolerance$normFactors,
                logLik = tolerance$logLik,
                options = opts
            )
            if (is.null(tolerance$normFactors)) {
                opts$tolerance$normFactors <- NULL
            }
            opts
        }

        opts <- setOpts(boundaries, tolerance)
        initf <- match.fun(init)
        initPars <- initf(pd$counts) # use pulseData here
        fit <- pulseR::fitModel(pd, initPars, opts)
        fit$d

 #       cis <- pulseR::ciGene("d",
 #                     par = re$fit,
 #                     geneIndexes = seq_along(re$fit$d),
 #                     pd = re$pd,
 #                     options = re$opts, confidence = conf.int)
 #       tp <- unique(re$pd$conditions$time)
    })

    names(re)=levels(conditions$Condition)

    adder=function(cond) {
        tab=data.frame(`Half-life`=log(2)/re[[cond]],check.names=FALSE)
        rownames(tab)=Genes(data,use.symbols = FALSE)
        name=if(is.null(cond)) name else paste0(name,".",cond)
        AddAnalysis(data,name=name,tab)
    }

    for (cond in levels(Condition(data))) data=adder(cond)
    data
}



#' Functions to compute the abundance of new or old RNA at time t.
#'
#' The standard mass action kinetics model of gene expression arises from the differential equation
#' \eqn{df/dt = s - d  f(t)}, with s being the constant synthesis rate, d the constant degradation rate and \eqn{f0=f(0)} (the abundance at time 0).
#'
#' @describeIn f.old.equi abundance of old RNA assuming steady state (i.e. f0=s/d)
#'
#' @param t time in h
#' @param s synthesis date in U/h (arbitrary unit U)
#' @param d degradation rate in 1/h
#' @param f0 the abundance at time t=0
#'
#' @return the RNA abundance at time t
#'
#' @examples
#'  d=log(2)/2
#'  s=10
#'
#'  f.new(2,s,d)  # Half-life 2, so after 2h the abundance should be half the steady state
#'  f.old.equi(2,s,d)
#'  s/d
#'
#'  t<-seq(0,10,length.out=100)
#'  plot(t,f.new(t,s,d),type='l',col='blue',ylim=c(0,s/d))
#'  lines(t,f.old.equi(t,s,d),col='red')
#'  abline(h=s/d,lty=2)
#'  abline(v=2,lty=2)
#'  # so old and new RNA are equal at t=HL (if it is at steady state at t=0)
#'
#'  plot(t,f.new(t,s,d),type='l',col='blue')
#'  lines(t,f.old.nonequi(t,f0=15,s,d),col='red')
#'  abline(h=s/d,lty=2)
#'  abline(v=2,lty=2)
#'  # so old and new RNA are not equal at t=HL (if it is not at steady state at t=0)
#'
#' @export
f.old.equi=function(t,s,d) s/d*exp(-t*d)

#' @describeIn f.old.equi abundance of old RNA without assuming steady state
#' @export
f.old.nonequi=function(t,f0,s,d) f0*exp(-t*d)
#' @describeIn f.old.equi abundance of new RNA (steady state does not matter)
#' @export
f.new=function(t,s,d) s/d*(1-exp(-t*d))



#' Fit a kinetic model according to non-linear least squares.
#'
#' Fit the standard mass action kinetics model of gene expression using least squares (i.e. assuming gaussian homoscedastic errors) for the given gene.
#' The fit takes both old and new RNA into account and requires proper normalization, but can be performed without assuming steady state.
#' The parameters are fit per \link{Condition}.
#'
#' @param data A grandR object
#' @param gene The gene for which to fit the model
#' @param slot The data slot to take expression values from
#' @param time The column in the column annotation table representing the labeling duration
#' @param conf.int A number between 0 and 1 representing the size of the confidence interval
#' @param steady.state either a named list of logical values representing conditions in steady state or not, or a single logical value for all conditions
#' @param use.old a logical vector to exclude old RNA from specific time points
#' @param use.new a logical vector to exclude new RNA from specific time points
#' @param maxiter the maximal number of iterations for the Levenberg-Marquardt algorithm used to minimize the least squares
#' @param compute.residuals set this to TRUE to compute the residual matrix
#'
#' @return
#' A named list containing the model fit:
#' \itemize{
#'   \item{data: a data frame containing the observed value used for fittin}
#'   \item{residuals: the computed residuals if compute.residuals=TRUE, otherwise NA}
#'   \item{Synthesis: the synthesis rate (in U/h, where U is the unit of the slot)},
#'   \item{Degradation: the degradation rate (in 1/h)}
#'   \item{Half-life: the RNA half-life (in h, always equal to log(2)/degradation-rate},
#'   \item{conf.lower: a vector containing the lower confidence bounds for Synthesis, Degradation and Half-life}
#'   \item{conf.upper: a vector containing the lower confidence bounds for Synthesis, Degradation and Half-life}
#'   \item{f0: The abundance at time 0 (in U)}
#'   \item{logLik: the log likelihood of the model}
#'   \item{rmse: the total root mean square error}
#'   \item{rmse.new: the total root mean square error for all new RNA values used for fitting}
#'   \item{rmse.old: the total root mean square error for all old RNA values used for fitting}
#'   \item{total: the total sum of all new and old RNA values used for fitting}
#'   \item{type: non-equi or equi}
#' }
#' If \code{Condition(data)} is not NULL, the return value is a named list (named according to the levels of \code{Condition(data)}), each
#' element containing such a structure.
#'
#' @details The start of labeling for all samples should be the same experimental time point. The fit gets more precise with multiple samples from multiple
#' labeling durations. In particular (but not only) without assuming steady state, also a sample without 4sU (representing time 0) is useful.
#'
#' @details The standard mass action kinetics model of gene expression arises from the following differential equation:
#'
#' @details \deqn{df/dt = s - d  f(t)}
#'
#' @details This model assumes constant synthesis and degradation rates (but not necessarily that the system is in steady state at time 0).
#' From the solution of this differential equation, it is straight forward to derive the expected abundance of old and new RNA at time t
#' for given parameters s (synthesis rate), d (degradation rate) and f0=f(0) (the abundance at time 0). These equations are implemented in
#' \code{\link{f.old.equi}} (old RNA assuming steady state gene expression, i.e. f0=s/d),
#' \code{\link{f.old.nonequi}} (old RNA without assuming steady state gene expression) and
#' \code{\link{f.new}} (new RNA; whether or not it is steady state does not matter).
#'
#' @details This function finds s and d such that the squared error between the observed values of old and new RNA and their corresponding functions
#' is minimized. For that to work, data has to be properly normalized.
#'
#' @seealso \link{FitKinetics}, \link{FitKineticsGeneLogSpaceLinear}, \link{FitKineticsGeneNtr}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' sars <- Normalize(sars)
#' SetParallel()
#' sars<-FitKinetics(sars,name = "kinetics",steady.state=list(Mock=TRUE,SARS=FALSE))   # fit per condition
#' head(GetAnalysisTable(sars,columns="Half-life"))
#'
#' @export
#'
FitKineticsGeneLeastSquares=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU,conf.int=0.95,steady.state=NULL,use.old=TRUE,use.new=TRUE, maxiter=100, compute.residuals=FALSE) {
    # residuals of the functions for usage with nls.lm
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

    correct=function(s) {
        if (max(s$Value)==0) s$Value[s$time==1]=0.01
        s$Type="New"
        s
    }

    stopifnot(time %in% names(Coldata(data)))

    newdf=GetData(data,mode.slot=paste0("new.",slot),genes=gene,ntr.na = FALSE)
    newdf$use=1:nrow(newdf) %in% (1:nrow(newdf))[use.new]
    newdf$time=newdf[[time]]
    if (is.null(newdf$Condition)) {
        newdf$Condition="Data"
        if (length(steady.state)==1) names(steady.state)="Data"
    }
    newdf=dlply(newdf,"Condition",function(s) correct(s))

    correct=function(s) {
        if (max(s$Value)==0) s$Value=0.01
        s$Type="Old"
        s
    }
    olddf=GetData(data,mode.slot=paste0("old.",slot),genes=gene,ntr.na = FALSE)
    olddf$use=1:nrow(olddf) %in% (1:nrow(olddf))[use.old]
    olddf$time=olddf[[time]]
    if (is.null(olddf$Condition)) olddf$Condition="Data"
    olddf=dlply(olddf,"Condition",function(s) correct(s))


    fit.equi=function(ndf,odf) {

        tinit=if (sum(ndf$Value>0)==0) min(ndf$time[ndf$time>0]) else min(ndf$time[ndf$Value>0])
        init.d=mean(-log(1-(0.1+ndf[ndf$time==tinit,"Value"])/(0.2+ndf[ndf$time==tinit,"Value"]+odf[odf$time==tinit,"Value"])))
        init.s=init.d*mean(odf[odf$time==tinit,"Value"])
        model.p=minpack.lm::nls.lm(c(init.s,init.d),lower=c(0,0.01),
                       fn=res.fun.equi,
                       old=odf[odf$use,],
                       new=ndf[ndf$use,],
                       control=minpack.lm::nls.lm.control(maxiter = maxiter))

        if (model.p$niter==maxiter) return(list(data=NA,residuals=if (compute.residuals) data.frame(Name=c(as.character(odf$Name[odf$use]),as.character(ndf$Name[ndf$use])),Type=c(rep("old",nrow(ndf)),rep("new",nrow(odf))),Absolute=NA,Relative=NA) else NA,Synthesis=NA,Degradation=NA,`Half-life`=NA,conf.lower=c(NA,NA),conf.upper=c(NA,NA),f0=NA,logLik=NA,rmse=NA, rmse.new=NA, rmse.old=NA,total=NA,type="equi"))
        conf.p=try(minpack.lm:::confint.nls.lm(model.p,level=conf.int),silent = TRUE)
        if (!is.matrix(conf.p)) conf.p=matrix(c(NA,NA,NA,NA),nrow=2)
        par=setNames(model.p$par,c("s","d"))
        rmse=sqrt(model.p$deviance/(nrow(ndf)+nrow(odf)))
        fvec=res.fun.equi(par,odf,ndf)
        n=nrow(odf)
        rmse.old=sqrt(sum(fvec[1:n]^2)/n)
        rmse.new=sqrt(sum(fvec[(n+1):(n*2)]^2)/n)
        df=droplevels(rbind(odf[odf$use,],ndf[ndf$use,]))


        residuals=NA
        if (compute.residuals) {
            #s=par["s"]
            #d=par["d"]
            ## solve f.old.equi(t)=old for t! ( and the same for)
            ##time=daply(df,.(Name),function(sub) 1/d * log( 1+ sum(sub$Value[sub$Type=="New"])/sum(sub$Value[sub$Type=="Old"]) ) )
            ## solve f.old.equi(t)=old for t!
            #time=daply(df,.(Name),function(sub) -1/d* mean(c( log(sub$Value[sub$Type=="Old"]*d/s), log(1-sub$Value[sub$Type=="New"]*d/s)   )) )
            ## new+old=s/d, so for a particular time point, multiply such that the sum is s/d
            #norm.fac=daply(df,.(Name),function(sub) s/d / sum(sub$Value))
            #modifier=data.frame(Name=names(norm.fac),Time=time,Norm.factor=norm.fac)

            resi=res.fun.equi(model.p$par,odf[odf$use,],ndf[ndf$use,])
            modval=c(odf[odf$use,"Value"],ndf[ndf$use,"Value"])-resi
            residuals=data.frame(Name=c(as.character(odf$Name[odf$use]),as.character(ndf$Name[ndf$use])),Type=c(rep("old",nrow(ndf)),rep("new",nrow(odf))),Absolute=resi,Relative=resi/modval)
        }
        total=sum(ndf$Value)+sum(odf$Value)
        list(data=df,
             residuals=residuals,
             Synthesis=unname(par['s']),
             Degradation=unname(par['d']),
             `Half-life`=log(2)/unname(par['d']),
             conf.lower=c(Synthesis=unname(conf.p[1,1]),Degradation=unname(conf.p[2,1]),`Half-life`=unname(log(2)/(conf.p[2,2]))),
             conf.upper=c(Synthesis=unname(conf.p[1,2]),Degradation=unname(conf.p[2,2]),`Half-life`=unname(log(2)/(conf.p[2,1]))),
             f0=unname(par['s']/par['d']),
             logLik=logLik.nls.lm(model.p),
             rmse=rmse, rmse.new=rmse.new, rmse.old=rmse.old,
             total=total,type="equi")
    }
    fit.nonequi=function(ndf,odf) {
        oind=union(which(odf$time==0),which(odf$use))
        nind=union(which(ndf$time==0),which(ndf$use))

        tinit=if (sum(ndf$Value>0)==0) min(ndf$time[ndf$time>0]) else min(ndf$time[ndf$Value>0])
        init.d=mean(-log(1-(0.1+ndf[ndf$time==tinit,"Value"])/(0.2+ndf[ndf$time==tinit,"Value"]+odf[odf$time==tinit,"Value"])))
        init.s=max(ndf$Value[ndf$time>0]/ndf$time[ndf$time>0])
        f0=mean(odf[odf$time==0,"Value"])
        model.m=minpack.lm::nls.lm(c(init.s,init.d,f0),lower=c(0,0.01,0),
                       fn=res.fun.nonequi,
                       old=odf[oind,],
                       new=ndf[nind,],
                       control=minpack.lm::nls.lm.control(maxiter = maxiter))
        if (model.m$niter==maxiter) return(list(data=NA,residuals=if (compute.residuals) data.frame(Name=c(as.character(odf$Name[odf$use]),as.character(ndf$Name[ndf$use])),Type=c(rep("old",nrow(ndf)),rep("new",nrow(odf))),Absolute=NA,Relative=NA) else NA,Synthesis=NA,Degradation=NA,`Half-life`=NA,conf.lower=c(NA,NA),conf.upper=c(NA,NA),f0=NA,logLik=NA,rmse=NA, rmse.new=NA, rmse.old=NA, total=NA,type="non.equi"))
        par=setNames(model.m$par,c("s","d"))
        f0=par[3]
        conf.m=try(minpack.lm:::confint.nls.lm(model.p,level=conf.int),silent = TRUE)
        if (!is.matrix(conf.m)) conf.m=matrix(c(NA,NA,NA,NA),nrow=2)
        rmse=sqrt(model.m$deviance/(nrow(ndf)+nrow(odf)))
        fvec=res.fun.nonequi.f0fixed(par,f0,odf,ndf)
        n=nrow(odf)
        rmse.old=sqrt(sum(fvec[1:n]^2)/n)
        rmse.new=sqrt(sum(fvec[(n+1):(n*2)]^2)/n)
        df=droplevels(rbind(odf[oind,],ndf[nind,]))

        residuals=NA
        if (compute.residuals) {
            #s=par["s"]
            #d=par["d"]
            #f0=unname(f0)
            ## solve f.old.nonequi(t)/f.new(t)=old/new for t!
            ##time=daply(df,.(Name),function(sub) 1/d * log( 1+ (sum(sub$Value[sub$Type=="New"])*f0*d)/(sum(sub$Value[sub$Type=="Old"])*s) ) )
            ## solve f.old.equi(t)=old for t!
            #time=daply(df,.(Name),function(sub) -1/d* mean(c( log(sub$Value[sub$Type=="Old"]/f0), log(1-sub$Value[sub$Type=="New"]*d/s)   )) )
            ## new+old=f.old.nonequi(t)+f.new(t), so for the computed time point t, multiply such that the sum is f.old.nonequi(t)+f.new(t)
            #norm.fac=daply(df,.(Name),function(sub) {
            #    t=time[as.character(sub$Name[1])]
            #    (f.old.nonequi(t,f0,s,d)+f.new(t,s,d)) / sum(sub$Value)
            #})
            #
            #modifier=data.frame(Name=names(norm.fac),Time=time,Norm.factor=norm.fac)
            resi=res.fun.nonequi(model.p$par,odf[odf$use,],ndf[ndf$use,])
            modval=c(odf[odf$use,"Value"],ndf[ndf$use,"Value"])-resi
            residuals=data.frame(Name=c(as.character(odf$Name[odf$use]),as.character(ndf$Name[ndf$use])),Type=c(rep("old",nrow(ndf)),rep("new",nrow(odf))),Absolute=resi,Relative=resi/modval)
        }
        total=sum(ndf$Value)+sum(odf$Value)
        list(data=df,
             residuals=residuals,
             Synthesis=unname(par['s']),
             Degradation=unname(par['d']),
             `Half-life`=log(2)/unname(par['d']),
             conf.lower=c(Synthesis=unname(conf.m[1,1]),Degradation=unname(conf.m[2,1]),`Half-life`=unname(log(2)/(conf.m[2,2]))),
             conf.upper=c(Synthesis=unname(conf.m[1,2]),Degradation=unname(conf.m[2,2]),`Half-life`=unname(log(2)/(conf.m[2,1]))),
             f0=unname(f0),
             logLik=logLik.nls.lm(model.m),
             rmse=rmse, rmse.new=rmse.new, rmse.old=rmse.old,
             total=total,type="non.equi")
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
    if (!is.null(Condition(data))) names(fits)=names(newdf) else fits=fits[[1]]
    fits
}

#' Fit a kinetic model using a linear model.
#'
#' Fit the standard mass action kinetics model of gene expression using a linear model after log-transforming the observed values
#' (i.e. assuming gaussian homoscedastic errors of the logarithmized values) for the given gene.
#' The fit takes only old RNA into account and requires proper normalization, but can be performed without assuming steady state for the degradation rate.
#' The parameters are fit per \link{Condition}.
#'
#' @param data A grandR object
#' @param gene The gene for which to fit the model
#' @param slot The data slot to take expression values from
#' @param time The column in the column annotation table representing the labeling duration
#' @param conf.int A number between 0 and 1 representing the size of the confidence interval
#'
#' @return
#' A named list containing the model fit:
#' \itemize{
#'   \item{data: a data frame containing the observed value used for fittin}
#'   \item{Synthesis: the synthesis rate (in U/h, where U is the unit of the slot)},
#'   \item{Degradation: the degradation rate (in 1/h)}
#'   \item{Half-life: the RNA half-life (in h, always equal to log(2)/degradation-rate},
#'   \item{conf.lower: a vector containing the lower confidence bounds for Synthesis, Degradation and Half-life}
#'   \item{conf.upper: a vector containing the lower confidence bounds for Synthesis, Degradation and Half-life}
#'   \item{f0: The abundance at time 0 (in U)}
#'   \item{logLik: the log likelihood of the model}
#'   \item{rmse: the total root mean square error}
#'   \item{adj.r.squared: adjusted R^2 of the linear model fit}
#'   \item{total: the total sum of all new and old RNA values used for fitting}
#'   \item{type: always "lm"}
#' }
#' If \code{Condition(data)} is not NULL, the return value is a named list (named according to the levels of \code{Condition(data)}), each
#' element containing such a structure.
#'
#' @details The start of labeling for all samples should be the same experimental time point. The fit gets more precise with multiple samples from multiple
#' labeling durations. Also a sample without 4sU (representing time 0) is useful.
#'
#' @details The standard mass action kinetics model of gene expression arises from the following differential equation:
#'
#' @details \deqn{df/dt = s - d  f(t)}
#'
#' @details This model assumes constant synthesis and degradation rates (but not necessarily that the system is in steady state at time 0).
#' From the solution of this differential equation, it is straight forward to derive the expected abundance of old and new RNA at time t
#' for given parameters s (synthesis rate), d (degradation rate) and f0=f(0) (the abundance at time 0). These equations are implemented in
#' \code{\link{f.old.equi}} (old RNA assuming steady state gene expression, i.e. f0=s/d),
#' \code{\link{f.old.nonequi}} (old RNA without assuming steady state gene expression) and
#' \code{\link{f.new}} (new RNA; whether or not it is steady state does not matter).
#'
#' @details This function primarily finds d such that the squared error between the observed values of old and new RNA and their corresponding functions
#' is minimized in log space. For that to work, data has to be properly normalized, but this is independent on any steady state assumptions. The synthesis
#' rate is computed (under the assumption of steady state) as \eqn{s=f0 \cdot d}
#'
#' @seealso \link{FitKinetics}, \link{FitKineticsGeneLeastSquares}, \link{FitKineticsGeneNtr}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' sars <- Normalize(sars)
#' SetParallel()
#' sars<-FitKinetics(sars,type='lm',name = "kinetics")    # fit per condition
#' head(GetAnalysisTable(sars,columns="Half-life"))
#'
#' @export
#'
FitKineticsGeneLogSpaceLinear=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU,conf.int=0.95) {
    correct=function(s) {
        if (max(s$Value)==0) s$Value[s$time==1]=0.01
        s$Type="New"
        s
    }

    stopifnot(time %in% names(Coldata(data)))

    correct=function(s) {
        if (max(s$Value)==0) s$Value=0.01
        s$Type="Old"
        s
    }
    olddf=GetData(data,mode.slot=paste0("old.",slot),genes=gene,ntr.na = FALSE)
    olddf$use=1:nrow(olddf) %in% (1:nrow(olddf))
    olddf$time=olddf[[time]]
    if (is.null(olddf$Condition)) olddf$Condition="Data"
    olddf=dlply(olddf,"Condition",function(s) correct(s))


    fit.lm=function(odf) {
        odf=odf[odf$Value>0,]
        fit=lm(log(Value)~time,data=odf)
        summ=summary(fit)

        par=setNames(c(exp(coef(fit)["(Intercept)"])*-coef(fit)["time"],-coef(fit)["time"]),c("s","d"))
        if (sum(abs(residuals(fit)))>0) {
            conf.p=confint(fit,level=conf.int)
            conf.p=apply(conf.p,2,function(v) setNames(pmax(0,c(exp(v["(Intercept)"])*par['d'],-v["time"])),c("s","d")))
        } else {
            conf.p=matrix(rep(NaN,4),ncol = 2)
        }

        modifier=NA
        total=sum(odf$Value)
        list(data=odf,
             Synthesis=unname(par['s']),
             Degradation=unname(par['d']),
             `Half-life`=log(2)/unname(par['d']),  # FIXME
             conf.lower=c(Synthesis=unname(conf.p[1,1]),Degradation=unname(conf.p[2,1]),`Half-life`=unname(log(2)/conf.p[2,2])),
             conf.upper=c(Synthesis=unname(conf.p[1,2]),Degradation=unname(conf.p[2,2]),`Half-life`=unname(log(2)/conf.p[2,1])),
             f0=unname(par['s']/par['d']),
             logLik=logLik(fit),
             rmse=sqrt(mean(fit$residuals^2)),
             adj.r.squared=summ$adj.r.squared,
             total=total,type="lm")
    }

    fits=lapply(names(olddf),function(cond) {
        odf=olddf[[cond]]
        fit.lm(odf)
    })
    if (!is.null(Condition(data))) names(fits)=names(olddf) else fits=fits[[1]]
    fits
}


#' Fit a kinetic model using the degradation rate transformed NTR posterior distribution.
#'
#' Fit the standard mass action kinetics model of gene expression by maximum a posteriori on a model based on the NTR posterior.
#' The fit takes only the NTRs into account and is completely independent on normalization, but it cannot be performed without assuming steady state.
#' The parameters are fit per \link{Condition}.
#'
#' @param data A grandR object
#' @param gene The gene for which to fit the model
#' @param slot The data slot to take expression values from
#' @param time The column in the column annotation table representing the labeling duration
#' @param conf.int A number between 0 and 1 representing the size of the confidence interval
#' @param unbiased Use the unbiased estimator instead of maximum likelihood on the transformed posterior
#' @param exact.ci compute exact credible intervals (see details)
#' @param total.fun use this function to summarize the expression values (only relevant for computing the synthesis rate s)
#'
#' @return
#' A named list containing the model fit:
#' \itemize{
#'   \item{data: a data frame containing the observed value used for fittin}
#'   \item{Synthesis: the synthesis rate (in U/h, where U is the unit of the slot)},
#'   \item{Degradation: the degradation rate (in 1/h)}
#'   \item{Half-life: the RNA half-life (in h, always equal to log(2)/degradation-rate},
#'   \item{conf.lower: a vector containing the lower confidence bounds for Synthesis, Degradation and Half-life}
#'   \item{conf.upper: a vector containing the lower confidence bounds for Synthesis, Degradation and Half-life}
#'   \item{f0: The abundance at time 0 (in U)}
#'   \item{logLik: the log likelihood of the model}
#'   \item{rmse: the total root mean square error}
#'   \item{total: the total sum of all new and old RNA values used for fitting}
#'   \item{type: always "ntr"}
#' }
#' If \code{Condition(data)} is not NULL, the return value is a named list (named according to the levels of \code{Condition(data)}), each
#' element containing such a structure.
#'
#' @details The start of labeling for all samples should be the same experimental time point. The fit gets more precise with multiple samples from multiple
#' labeling durations.
#'
#' @details The standard mass action kinetics model of gene expression arises from the following differential equation:
#'
#' @details \deqn{df/dt = s - d  f(t)}
#'
#' @details This model assumes constant synthesis and degradation rates. Further assuming steady state allows to derive the function transforming from
#' the NTR to the degradation rate d as \eqn{d(ntr)=-1/t log(1-ntr)}. Furthermore, if the ntr is (approximately) beta distributed, it is possible to
#' derive the distribution of the transformed random variable for the degradation rate (see Juerges et al., Bioinformatics 2018).
#'
#' @details This function primarily finds d by maximizing the degradation rate posterior distribution. For that, data does not have to be normalized,
#' but this only works under steady-state conditions. The synthesis rate is then computed (under the assumption of steady state) as \eqn{s=f0 \cdot d}
#'
#' @details The maximum-a-posteriori estimator is biased. Bias can be removed by a correction factor (which is done by default).
#'
#' @details By default the chi-squared approximation of the log-posterior function is used to compute credible intervals. If exact.ci is used, the
#' posterior is integrated numerically.
#'
#' @seealso \link{FitKinetics}, \link{FitKineticsGeneLeastSquares}, \link{FitKineticsGeneLogSpaceLinear}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#' sars <- subset(sars,Coldata(sars)$Cell=="Mock")         # we don't normalize, we are only interested in the Half-life
#' SetParallel()
#' sars<-FitKinetics(sars,type='ntr',name = "kinetics")    # note that Condition(sars)=NULL
#' head(GetAnalysisTable(sars,columns="Half-life"))
#'
#' @export
#'
FitKineticsGeneNtr=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU,conf.int=0.95,unbiased=TRUE,exact.ci=FALSE,total.fun=median) {
    if (!all(c("alpha","beta") %in% Slots(data))) stop("Beta approximation data is not available in grandR object!")
    steady.state=TRUE

    crit <- qchisq(1-conf.int, df = 2, lower.tail = FALSE) / 2

    bounds=c(log(2)/48,log(2)/0.01)

    total.df=GetData(data,mode.slot=paste0("total.",slot),genes=gene)
    if (is.null(total.df$Condition)) total.df$Condition=factor("Data")

    a=GetData(data,mode.slot="alpha",genes=gene)
    b=GetData(data,mode.slot="beta",genes=gene)
    ntr=GetData(data,mode.slot="ntr",genes=gene)
    total=GetData(data,mode.slot=paste0("total.",slot),genes=gene)
    if (is.null(a$Condition)) {
        a$Condition=factor("Data")
        b$Condition=factor("Data")
        ntr$Condition=factor("Data")
        total$Condition=factor("Data")
        if (length(steady.state)==1) names(steady.state)="Data"
    }


    #sloglik=function(d,a,b,t) (a-1)*log(1-exp(-t*d))-t*d*b
    #loglik=function(d,a,b,t) sum((a-1)*log(1-exp(-t*d))-t*d*b)
    #sloglik=function(d,a,b,t) (a-1)*log1p(-exp(-t*d))-t*d*(b-1)
    loglik=if (unbiased) function(d,a,b,t) sum((a-1)*log1p(-exp(-t*d))-t*d*(b-1)) else function(d,a,b,t) sum((a-1)*log1p(-exp(-t*d))-t*d*b)

    t=a[,time]
    use=t>0
    c=droplevels(a$Condition[use])
    a=a$Value[use]
    b=b$Value[use]
    total=total$Value[use]
    ntr=ntr$Value[use]
    t=t[use]
    #uniroot.save=function(fun,lower,upper) if (is.na(fun(lower)*fun(upper)) || fun(lower)*fun(upper)>=0) NA else uniroot(fun,lower=lower,upper=upper)$root
    uniroot.save=function(fun,lower,upper) if (fun(lower)*fun(upper)>=0) mean(lower,upper) else uniroot(fun,lower=lower,upper=upper)$root
    fits=lapply(levels(c),function(cc) {

        equi=if (is.null(steady.state)) {
            TRUE
        } else if(length(steady.state)==1) {
            steady.state;
        } else is.na(unlist(steady.state)[cc]) || as.logical(unlist(steady.state)[cc])

        # non-equi is super inaccurate: if one total count is off, the beta error model might just go crazy!
        f0=if (equi) total.fun(total.df$Value[total.df$Condition==cc]) else total.fun(total.df$Value[total.df$Condition==cc & total.df[[time]]==0])


        ind=c==cc & !is.na(a) & !is.na(b)
        if(any(is.nan(c(a[ind],b[ind])))) return(NaN)
        ploglik=if (equi) function(x) loglik(x,a=a[ind],b=b[ind],t=t[ind]) else function(x) loglik(x-log(f0/total[ind])/t[ind],a=a[ind],b=b[ind],t=t[ind])
        pbounds=if (equi) bounds else c(max(bounds[1],log(f0/total[ind])/t[ind]),bounds[2]) # ensure that: d> log(f0/total[ind])
        d=optimize(ploglik,pbounds,maximum=T)$maximum
        max=ploglik(d)
        df=data.frame(alpha=a[ind],beta=b[ind],t=t[ind])
        if (exact.ci) {
            plik=function(x) pmax(0,sapply(x,function(xx) exp(loglik(xx,a=a[ind],b=b[ind],t=t[ind])-max)))
            inte = function(u) integrate(plik,pbounds[1],u)$value
            lower1=uniroot.save(function(x) ploglik(x)-max+log(1E6),lower = pbounds[1],upper=d)
            upper1=uniroot.save(function(x) ploglik(x)-max+log(1E6),lower = d,upper=pbounds[2])

            total.inte=inte(upper1)
            lower=uniroot.save(function(x) inte(x)/total.inte-(1-conf.int)/2,lower = lower1,upper=d)
            upper=uniroot.save(function(x) inte(x)/total.inte-(1+conf.int)/2,lower = d,upper=upper1)
        } else {
            lower=uniroot.save(function(x) ploglik(x)-max+crit,lower = pbounds[1],upper=d)
            upper=uniroot.save(function(x) ploglik(x)-max+crit,lower = d,upper=pbounds[2])
        }

        rmse=sum(sqrt(((1-exp(-t[ind]*d)-ntr[ind]))^2))/sum(ind)
        if (!equi) {
            s=median(ntr/(1-exp(-t*d))*d*total)
            s.low=median(ntr/(1-exp(-t*lower))*lower*total)
            s.up=median(ntr/(1-exp(-t*upper))*upper*total)
        }
        else {
            s=f0*d
            s.low=f0*lower
            s.up=f0*upper
        }

        list(data=df,
             Synthesis=s,
             Degradation=d,
             `Half-life`=log(2)/d,
             conf.lower=c(Synthesis=s.low,Degradation=lower,`Half-life`=log(2)/upper),
             conf.upper=c(Synthesis=s.up,Degradation=upper,`Half-life`=log(2)/lower),
             f0=unname(f0),
             logLik=max,
             rmse=rmse,
             total=f0,
             type="ntr")
    })
    if (!is.null(Condition(data))) names(fits)=levels(c) else fits=fits[[1]]
    fits

}

#' Uses the kinetic model to calibrate the effective labeling time.
#'
#' The NTRs of each sample might be systematically too small (or large). This function identifies such systematic
#' deviations and computes labeling durations without systematic deviations.
#'
#' @param data A grandR object
#' @param slot The data slot to take expression values from
#' @param time The column in the column annotation table representing the labeling duration
#' @param time.name The name in the column annotation table to put the calibrated labeling durations
#' @param time.conf.name The name in the column annotation table to put the confidence values for the labeling durations (half-size of the confidence interval)
#' @param conf.int The level for confidence intervals
#' @param steady.state either a named list of logical values representing conditions in steady state or not, or a single logical value for all conditions
#' @param n.estimate the times are calibrated with the top n expressed genes
#' @param n.iter the maximal number of iterations for the numerical optimization
#' @param verbose verbose output
#' @param ... forwarded to FitKinetics
#'
#' @return
#' A new grandR object containing the calibrated durations in the column data annotation
#'
#' @details There are many reasons why the nominal (wall-clock) time of 4sU labeling might be distinct from the effective labeling time. Most
#' importantly, 4sU needs some time to enter the cells and get activated to be ready for transcription. Therefore, the 4sU concentration
#' (relative to the U concentration) rises, based on observations, over the timeframe of 1-2h. GRAND-SLAM assumes a constant 4sU incorporation rate,
#' i.e. specifically new RNA made early during the labeling is underestimated. This, especially for short labeling (<2h), the effective labeling duration
#' might be significantly less than the nominal labeling duration.
#'
#' @details It is impossible to obtain a perfect absolute calibration, i.e. all durations might be off by a factor.
#'
#' @seealso \link{FitKinetics}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' sars <- Normalize(sars)
#' SetParallel()
#' sars<-CalibrateEffectiveLabelingTime(sars,steady.state=list(Mock=TRUE,SARS=FALSE))
#' Coldata(sars)
#'
#' @export
#'
CalibrateEffectiveLabelingTime=function(data,slot=DefaultSlot(data),time=Design$dur.4sU,time.name="calibrated_time",time.conf.name="calibrated_time_conf",conf.int=0.95,steady.state=NULL,n.estimate=1000, n.iter=1000, verbose=FALSE,...) {

    conds=Coldata(data)
    if (is.null(conds$Condition)) {
        conds$Condition=factor("Data")
        if (length(steady.state)==1) names(steady.state)="Data"
    }
    re=matrix(NA,ncol=2,nrow=nrow(conds),dimnames = list(conds$Name,c(time.name,time.conf.name)))
    for (cond in levels(conds$Condition)) {
        if (verbose) cat(sprintf("Calibrating %s...\n",cond))
        sub=subset(data,columns=conds$Condition==cond)
        Condition(sub)=NULL

        # restrict to top n genes stratified by ntr
        totals=rowSums(GetTable(sub,type=slot))
        ntrs=rowMeans(GetTable(sub,type="ntr"),na.rm=TRUE)
        ntr.cat=cut(ntrs,quantile(ntrs,seq(0,1,by=0.1)),include.lowest = TRUE)
        fil=plyr::ddply(data.frame(Gene=Genes(sub),totals,ntr.cat),.(ntr.cat),function(s) {
           threshold=sort(s$totals,decreasing = TRUE)[min(length(s$totals),ceiling(n.estimate/10))]
           data.frame(Gene=s$Gene,use=s$totals>=threshold)
        })

        genes=as.character(fil$Gene[fil$use])
        sub=FilterGenes(sub,use=genes)
        sub=DropAnalysis(sub)

        opt.fun=function(times) {
            tt=init
            tt[use]=times
            Coldata(sub,Design$dur.4sU)=tt
            fit=FitKinetics(sub,return.fields='logLik',steady.state=steady.state[[cond]],...)
            #fit=FitKinetics(sub,return.fields='logLik',slot=slot,steady.state=steady.state[[cond]])
            sum(GetAnalysisTable(fit,gene.info=FALSE)[,1])
        }
        init=Coldata(sub)[[time]]
        init[init>0]=init[init>0]
        use=init>0 & init<max(init)
        fit=optim(init[use],fn=opt.fun,hessian=TRUE, control=list(fnscale=-1,maxit=n.iter))
        if (fit$convergence!=0) stop(sprintf("Did not converge!"))
        conf=try(sd.from.hessian(fit$hessian)*qnorm(1-(1-conf.int)/2),silent=TRUE)

        tt=init
        tt[use]=fit$par
        re[as.character(Coldata(sub)$Name),time.name]=tt
        tt=rep(0,length(tt))
        tt[use]=conf
        re[as.character(Coldata(sub)$Name),time.conf.name]=tt
    }
    data=Coldata(data,re)
    data
}

#' Fit kinetic models to all genes.
#'
#' Fit the standard mass action kinetics model of gene expression by different means. Depending on which method is used, either steady state must be assumed,
#' or data must be properly normalized. The parameters are fit per \link{Condition}.
#'
#' @param data A grandR object
#' @param name the user defined analysis name to store the results
#' @param type Which method to use (either one of "full","ntr","lm")
#' @param slot The data slot to take expression values from
#' @param time The column in the column annotation table representing the labeling duration
#' @param conf.int A number between 0 and 1 representing the size of the confidence interval
#' @param return.fields which statistics to return (see details)
#' @param return.extra additional statistics to return (see details)
#' @param ... forwarded to \code{\link{FitKineticsGeneNtr}}, \code{\link{FitKineticsGeneLeastSquares}} or \code{\link{FitKineticsGeneLogSpaceLinear}}
#' @return
#' A new grandR object with the fitted parameters as an analysis table
#'
#' @details The start of labeling for all samples should be the same experimental time point. The fit gets more precise with multiple samples from multiple
#' labeling durations.
#'
#' @details The standard mass action kinetics model of gene expression arises from the following differential equation:
#'
#' @details \deqn{df/dt = s - d  f(t)}
#'
#' @details This model assumes constant synthesis and degradation rates. Based on this, there are different ways for fitting the parameters:
#' \itemize{
#'   \item{\link{FitKineticsGeneLeastSquares}: non-linear least squares fit on the full model; depends on proper normalization; can work without steady state; assumption of homoscedastic gaussian errors is theoretically not justified}
#'   \item{\link{FitKineticsGeneLogSpaceLinear}: linear model fit on the old RNA; depends on proper normalization; assumes steady state for estimating the synthesis rate; assumption of homoscedastic gaussian errors in log space is problematic and theoretically not justified}
#'   \item{\link{FitKineticsGeneNtr}: maximum a posteriori fit on the NTR posterior transformed to the degradation rate; as it is based on the NTR only, it is independent on proper normalization; assumes steady state; theoretically well justified}
#' }
#'
#' @details This function is flexible in what to put in the analysis table. You can specifiy the statistics using return.fields and return.extra (see \code{\link{kinetics2vector}})
#'
#' @seealso \link{FitKineticsGeneNtr}, \link{FitKineticsGeneLeastSquares}, \link{FitKineticsGeneLogSpaceLinear}
#'
#' @examples
#' SetParallel()
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#' sars <- subset(sars,Coldata(sars)$Cell=="Mock")
#' sars<-FitKinetics(sars,name="kinetics.ntr",type='ntr')
#' sars<-Normalize(sars)
#' sars<-FitKinetics(sars,name="kinetics.nlls",type='full')
#' sars<-FitKinetics(sars,name="kinetics.lm",type='lm')
#' head(GetAnalysisTable(sars,columns="Half-life"))
#'
#' @export
#'
FitKinetics=function(data,name="kinetics",type=c("nlls","ntr","lm"),slot=DefaultSlot(data),time=Design$dur.4sU,conf.int=0.95,return.fields=c("Synthesis","Half-life","rmse"),return.extra=NULL,...) {

    fun=switch(tolower(type[1]),ntr=FitKineticsGeneNtr,nlls=FitKineticsGeneLeastSquares,lm=FitKineticsGeneLogSpaceLinear)
    funs=switch(tolower(type[1]),ntr="FitKineticsGeneNtr",nlls="FitKineticsGeneLeastSquares",lm="FitKineticsGeneLogSpaceLinear")

    if (funs=="FitKineticsGeneNtr" && !all(c("alpha","beta") %in% Slots(data))) stop("Beta approximation data is not available in grandR object!")

    if (is.null(fun)) stop(sprintf("Type %s unknown!",type))
    result=plapply(Genes(data),
                      fun,
                      data=data,
                      slot=slot,time=time,conf.int=conf.int,
                      ...)

    adder=function(cond) {
        tttr=function(re) {
            if (is.matrix(re)) as.data.frame(t(re)) else setNames(data.frame(re),names(kinetics2vector(result[[1]],condition=cond,return.fields=return.fields,return.extra=return.extra)))
        }
        slam.param=tttr(sapply(result,kinetics2vector,condition=cond,return.fields=return.fields,return.extra=return.extra))
        rownames(slam.param)=Genes(data,use.symbols = FALSE)
        name=if(is.null(cond)) name else paste0(name,".",cond)
        AddAnalysis(data,name=name,slam.param)
    }

    if (!is.null(Condition(data))) {
        for (cond in levels(Condition(data))) data=adder(cond)
    } else {
        data=adder(NULL)
    }
    data
}


#' Compute each of the three kinetic parameters given old and new abundance and one parameter.
#'
#' The standard mass action kinetics model of gene expression arises from the differential equation
#' \eqn{df/dt = s - d  f(t)}, with s being the constant synthesis rate, d the constant degradation rate and \eqn{f0=f(0)} (the abundance at time 0).
#' This model dictates the time evolution of old and new RNA abundance after metabolic labeling starting at time t=0.
#' This function computes two of the three parameters s, d or f0 given the one of them  for observed old and new RNA counts at time t.
#'
#' @param old the abundance of old RNA at time t
#' @param new the abundance of new RNA at time t
#' @param total the abundance of total RNA at time t
#' @param ntr the new-to-total RNA ratio at time t
#' @param t the time t
#' @param s the synthesis rate
#' @param d the degradation rate
#' @param f0 the total abundance at time t=0
#' @param name.prefix prepend this prefix to the column names of the output matrix (which are f0,s,d,t0)
#'
#' @return A matrix with columns f0,s,d and t0; the first three are the given or computed parameters, t0 is the time where the abundance was \eqn{f(t0)=0}
#'
#' @details Can either be parametrized using old and new, or using total and ntr
#'
#' @details Because old RNA cannot increase over time, and for \eqn{d>0} it has to decrease, it must be f0>old.
#'
#' @details \eqn{s \cdot t} would be the new RNA with d=0, i.e. it must be \eqn{s \cdot t>new}
#'
#' @details If \eqn{f0<s/d}, there must be a time t0 were \eqn{f(t0)=0}. t0 is also estimated.
#'
#' @export
#'
#' @examples
#' d=log(2)/2  # a Half-life of 2h
#' s=10*d      # a synthesis rate such that the steady-state abundance is 10
#' f0=8        # we start below the steady-state
#' t=2
#'
#' # compute old and new RNA
#' old=f.old.nonequi(t,f0,s,d)
#' new=f.new(t,s,d)
#'
#' # compare computed parameter to truth
#' rbind(
#'    c(f0=f0,s=s,d=d,t0=NA),
#'    TransformKineticParameters(old,new,t,s=s),
#'    TransformKineticParameters(old,new,t,d=d),
#'    TransformKineticParameters(old,new,t,f0=f0)
#'    )
#'
#'  # verify t0 (it has to be f.new(-t0)=f0)
#'  t0=unname(TransformKineticParameters(old,new,t,s=s)[,'t0'])
#'  f.new(-t0,s,d)
#'
TransformKineticParameters=function(old=NULL,new=NULL,total=NULL,ntr=NULL,t=2,s=NULL,d=NULL,f0=NULL,name.prefix=NULL) {
    if (is.null(s) + is.null(d) + is.null(f0)!=2) stop("Exactly two of of s,d,f0 have to be NULL!")

    total.ntr=!is.null(total) & !is.null(ntr)
    old.new=!is.null(old) & !is.null(new)

    if (total.ntr+old.new!=1) stop("Must be parametrized either with old and new, or with total and ntr!")
    if (total.ntr) {
        new=total*ntr
        old=total-new
    }

    if (!is.null(f0)) {
        if (any(f0<=old)) {
            warning("It must be f0>old!")
        }
        F=old/f0
        d=-1/t*log(F)
        s=-1/t*new * ifelse(F>=1,-1,ifelse(is.infinite(F),0,log(F)/(1-F)))
    } else if (!is.null(d)) {
        f0=old/exp(-t*d)
        s=new*d/(1-exp(-t*d))
    } else {
        if (any(s*t<=new)) {
            warning("It must be s*t>new!")
        }
        d=pmax(0,lamW::lambertW0(-exp(-s*t/new)*s*t/new)/t+s/new)
        f0=old/exp(-t*d)
    }

    t0=t-ifelse(f0<s/d,-1/d*log(1-(old+new)*d/s),NA)

    re=cbind(f0=f0,s=s,d=d,t0=t0)
    if (!is.null(name.prefix)) colnames(re)=paste0(name.prefix,".",colnames(re))
    re
}




#' Compute the posterior distributions of RNA synthesis and degradation for a particular gene
#'
#' @param data the grandR object
#' @param gene a gene name or symbol or index
#' @param columns samples or cell representing the same experimental condition (must refer to a unique labeling duration)
#' @param reference.columns either a reference matrix (\code{\link{FindReferences}}) to define reference columns for each of the samples or cell defined by the \code{columns} parameter, or a
#' logical vector denoting the reference columns; can be NULL, but then the samples or cells have to be in steady state
#' @param dispersion dispersion parameter for the given columns (if NULL, this is estimated from the data, takes a lot of time!)
#' @param slot the data slot to take f0 and totals from
#' @param time.labeling the column in the column annotation table denoting the labeling duration or the labeling duration itself
#' @param time.experiment the column in the column annotation table denoting the experimental time point (can be NULL, see details)
#' @param sample.f0.in.ss whether or not to sample f0 under steady state conditions
#' @param hierarchical Take the NTR from the hierarchical bayesian model (see details)
#' @param return.samples return the posterior samples of the parameters?
#' @param return.points return the point estimates per replicate as well?
#' @param N the posterior sample size
#' @param N.max the maximal number of posterior samples (necessary if old RNA > f0); if more are necessary, a warning is generated
#' @param conf.int A number between 0 and 1 representing the size of the credible interval
#'
#' @details The kinetic parameters s and d are computed using \link{TransformSnapshot}. For that, the sample either must be in steady state
#' (this is the case if defined in the reference.columns matrix), or if the levels of reference samples from a specific prior time point are known. This time point is
#' defined by \code{time.experiment} (i.e. the difference between the reference samples and samples themselves). If
#' \code{time.experiment} is NULL, then the labeling time of the samples is used (e.g. useful if labeling was started concomitantly with
#' the perturbation, and the reference samples are unperturbed samples).
#'
#' @details By default, the hierarchical bayesian model is estimated. If hierarchical = FALSE, the NTRs are sampled from a beta distribution
#' that approximates the mixture of betas from the replicate samples. (see \link{beta.approximate.mixture}).
#'
#' @return a list containing the posterior mean of s and s, its credible intervals and,
#' if return.samples=TRUE a data frame containing all posterior samples
#'
#'
#' @export
#'
FitKineticsSnapshot=function(data,gene,columns,
                             reference.columns=NULL,
                             dispersion=NULL,
                             slot=DefaultSlot(data),time.labeling=Design$dur.4sU,time.experiment=NULL,
                             sample.f0.in.ss=TRUE,
                             hierarchical=TRUE,
                             return.samples=FALSE,
                             return.points=FALSE,
                             N=10000,N.max=N*10,
                             conf.int=0.95) {
    if (is.matrix(reference.columns)) reference.columns=apply(reference.columns[,Columns(data,columns)]==1,1,any)

    alpha=GetData(data,mode.slot="alpha",genes=gene,columns = columns)
    if (!is.numeric(time.labeling) && length(unique(alpha[[time.labeling]]))!=1) stop("A has to refer to a unique labeling duration!")
    if (!is.null(time.experiment) && length(unique(alpha[[time.experiment]]))!=1) stop("A has to refer to a unique experimental time!")
    beta=GetData(data,mode.slot="beta",genes=gene,columns = columns)
    total=GetData(data,mode.slot=slot,genes=gene,columns = columns)
    ss=GetData(data,mode.slot=slot,genes=gene,columns = reference.columns)
    if (!is.null(time.experiment) && length(unique(ss[[time.experiment]]))!=1) stop("Steady state for A has to refer to a unique experimental time!")

    if (is.null(dispersion)) dispersion=estimate.dispersion(as.matrix(GetTable(data,type="count",columns=columns,gene.info = F)))[ToIndex(data,gene)]


    use=total$Value>0 & !is.na(alpha$Value)

    emptyres=function() {
        re=list(
            s=NA,
            d=NA,
            s.conf.int=c(-Inf,Inf),
            d.conf.int=c(-Inf,Inf)
        )
        if (return.samples)
            re$samples=data.frame(s=numeric(0),d=numeric(0))
        if (return.points)
            re$points=data.frame(s=numeric(0),d=numeric(0))
        return(re)
    }

    if (sum(use)<2 || sum(ss$Value>0)==0) return(emptyres())


    alpha=alpha[use,]
    beta=beta[use,]
    total=total[use,]
    ss=ss[ss$Value>0,]

    t=if (is.numeric(time.labeling)) time.labeling else unique(alpha[[time.labeling]])
    t0=if (is.null(time.experiment)) t else unique(alpha[[time.experiment]])-unique(ss[[time.experiment]])
    is.steady.state=any(reference.columns & columns)


    if (return.points)
        points=as.data.frame(if (is.steady.state) TransformSnapshot(ntr=alpha$Value/(alpha$Value+beta$Value),total=total$Value,t=t) else TransformSnapshot(ntr=alpha$Value/(alpha$Value+beta$Value),total=total$Value,t=t,t0=t0,f0=mean(ss$Value)))

    if (N<1) {
        ntr=sum(alpha$Value)/(sum(alpha$Value)+sum(beta$Value))
        if (is.steady.state) {
            param=TransformSnapshot(ntr=ntr,total=mean(total$Value),t=t)
        } else {
            if (t0<=0) stop("Experimental time is not properly defined (the steady state sample must be prior to each of A and B)!")
            f0=mean(ss$Value)
            param=TransformSnapshot(ntr=ntr,total=mean(total$Value),t=t,t0=t0,f0=f0)
        }

        re=list(
            s=param['s'],
            d=param['d'],
            s.conf.int=c(-Inf,Inf),
            d.conf.int=c(-Inf,Inf)
        )
        if (return.samples)
            re$samples=data.frame(s=numeric(0),d=numeric(0))
        if (return.points)
            re$points=points
        return(re)
    }

    if (hierarchical) {
        mod=hierarchical.beta.posterior(alpha$Value,beta$Value,compute.marginal.likelihood = FALSE,compute.grid = TRUE,res=50)$sample.mu
    } else {
        fit=beta.approximate.mixture(alpha$Value,beta$Value)
        mod=function(N) rbeta(N,fit$a,fit$b)
    }

    sample.ss.with.f0=function(N) {
        f0=0.1+rnbinom(N,size=1/dispersion,mu=mean(total$Value))
        ntr=mod(N)
        total=0.1+rnbinom(N,size=1/dispersion,mu=mean(total$Value))
        TransformSnapshot(ntr=ntr,total=total,t=t,f0=f0,t0=t)
    }
    sample.ss.without.f0=function(N) {
        ntr=mod(N)
        total=0.1+rnbinom(N,size=1/dispersion,mu=mean(total$Value))
        TransformSnapshot(ntr=ntr,total=total,t=t)
    }
    sample.non.ss=function(N) {
        if (t0<=0) stop("Experimental time is not properly defined (the steady state sample must be prior to each of A and B)!")
        f0=0.1+rnbinom(N,size=1/dispersion,mu=mean(ss$Value))
        total=0.1+rnbinom(N,size=1/dispersion,mu=mean(total$Value))
        ntr=mod(N)
        TransformSnapshot(ntr=ntr,total=total,t=t,t0=t0,f0=f0)
    }
    resample=function(FUN) {
        mat=FUN(N)
        inadmissible=mat[,'d']==0 | is.infinite(mat[,'d'])
        if (sum(inadmissible)>0) {
            success.prob=1-sum(inadmissible)/N
            additional.samples=sum(inadmissible)+ceiling(sum(inadmissible)*(1-success.prob)/success.prob)
            if (additional.samples>N.max-N) {
                warning("Inefficient sampling; results might be unreliable")
                additional.samples=N.max-N
            }
            mat2=FUN(additional.samples)
            inadmissible2=mat2[,'d']==0 | is.infinite(mat2[,'d'])
            mat=rbind(mat[!inadmissible,],mat2[!inadmissible2,])
        }
        mat
    }

    samp=if (!is.steady.state) sample.non.ss else if (sample.f0.in.ss) sample.ss.with.f0 else sample.ss.without.f0

    samp=resample(samp)

    N=nrow(samp)
    if (N==0) return(emptyres())

    sc=quantile(samp[,'s'],c(0.5-conf.int/2,0.5+conf.int/2))
    dc=quantile(samp[,'d'],c(0.5-conf.int/2,0.5+conf.int/2))

    re=list(
        s=mean(samp[,'s']),
        d=mean(samp[,'d']),
        s.conf.int=sc,
        d.conf.int=dc
    )
    if (return.samples)
        re$samples=data.frame(s=samp[,'s'],d=samp[,'d'])
    if (return.points)
        re$points=points

    re
}


#' Estimate parameters for a one-shot experiment.
#'
#' Under steady state conditions it is straight-forward to estimate s and d. Otherwise, the total levels at some other time point are needed.
#'
#' @param ntr the new to total RNA ratio (measured)
#' @param total the total level of RNA (measured)
#' @param t the labeling duration
#' @param t0 time before measurement at which f0 is total level (only necessary under non-steady-state conditions)
#' @param f0 total level at t0 (only necessary under non-steady-state conditions)
#'
#' @return a named vector for s and d
#'
#' @export
#'
TransformSnapshot=function(ntr,total,t,t0=NULL,f0=NULL) {
    if (is.null(f0)) {
        d=-1/t*log(1-ntr)
        s=total*d
    } else if (t0==t) {
        Fval=pmin(total*(1-ntr)/f0,1)
        d=ifelse(Fval>=1,0,-1/t*log(Fval))
        s=-1/t*total*ntr * ifelse(Fval>=1,-1,ifelse(is.infinite(Fval),0,log(Fval)/(1-Fval)))
    } else if (length(ntr)>1 || length(total)>1 || length(f0)>1) {
        m=cbind(ntr,total,f0)
        return(t(sapply(1:nrow(m),function(i) TransformSnapshot(m[i,1],m[i,2],t,t0,m[i,3]))))
    } else {
        new=total*ntr
        old=total-new
        f0new=f0-new
        inter=c(log(2)/48,log(2)/0.001)
        #print(c(total,ntr,f0))
        eq=function(d) total*exp(-t*d)-f0*exp(-t0*d)*exp(-t*d)+f0new*exp(-t0*d)-old
        #print((eq(inter[1])))
        #print((eq(inter[2])))
        if (eq(inter[1])<0) {
            d=0
            s=0
        } else if (eq(inter[2])>0) {
            d=Inf
            s=Inf
        } else {
            d=uniroot(eq,inter)$root
            s=new*d/(1-exp(-t*d))
        }
    }
    if (length(s)>1) cbind(s=unname(s),d=unname(d)) else c(s=unname(s),d=unname(d))
}


#' Generate plots about potential kinetic parameters and their changes for observed abundance and ntr estimates
#'
#' The standard mass action kinetics model of gene expression arises from the differential equation
#' \eqn{df/dt = s - d  f(t)}, with s being the constant synthesis rate, d the constant degradation rate and \eqn{f0=f(0)} (the abundance at time 0).
#' The RNA half-life is directly related to d via \eqn{HL=log(2)/d}.
#' This model dictates the time evolution of old and new RNA abundance after metabolic labeling starting at time t=0. However, the observed
#' abundance and NTR values alone do not allow to infer s,d and f0 simultaneously. Also changes from one condition to another of these parameters cannot be
#' estimated from the abundance and NTR value of both conditions. This function produces a plot showing the relation of these three parameters for given abundances and NTRs.
#'
#' @param total1 the total abundance of the RNA under condition 1 at time t
#' @param total2 the total abundance of the RNA under condition 2 at time t
#' @param ntr1 the new-to-total RNA ratio under condition 1 at time t
#' @param ntr2 the new-to-total RNA ratio under condition 2 at time t
#' @param f0.above.ss.factor the maximal value of f0 is computed according to total1 times this factor
#' @param t the time t
#' @param N the grid size to compute the heatmaps
#'
#' @return A patchworked set of ggplots
#'
#' @details f0 must be strictly greater than \eqn{total \cdot (1-ntr)} (as old RNA has to increase with d>0). If \eqn{f0>total}, the total RNA decreases.
#' Thus, the f0 range considered therefore goes from \eqn{total \cdot (1-ntr)} to \eqn{total \codt f0.above.ss.factor}.
#'
#' @export
#'
#' @examples
#' PlotSimpleKinetics(100,200,0.2,0.2)  # the total abundance increases with the ntr staying constant -> likely increase in synthesis
#' PlotSimpleKinetics(100,200,0.2,0.1)  # the total abundance increases but the new RNA stays constant -> likely increase in stability
#'
#'
PlotSimpleKinetics=function(total1,total2,ntr1,ntr2,f0.above.ss.factor=1,t=2,N=100) {
    old1=total1*(1-ntr1)
    old2=total2*(1-ntr2)
    new1=total1*ntr1
    new2=total2*ntr2

    k1=as.data.frame(TransformKineticParameters(old=old1,new=new1,t=t,f0=seq(old1+0.01,total1*f0.above.ss.factor,length.out=N),name.prefix = "1"))
    k2=as.data.frame(TransformKineticParameters(old=old2,new=new2,t=t,f0=seq(old2+0.01,total2*f0.above.ss.factor,length.out=N),name.prefix = "2"))
    df=data.frame(i_1=rep(1:N,N),i_2=rep(1:N,each=N))
    df=cbind(df,k1[df$i_1,])
    df=cbind(df,k2[df$i_2,])


    s_1=ggplot(k1,aes(f0_1,s_1))+geom_line()
    s_2=ggplot(k2,aes(f0_2,s_2))+geom_line()
    s_lfc=ggplot(df,aes(f0_1,f0_2,fill=log2(s_1/s_2)))+
        geom_tile()+
        scale_fill_viridis_c()

    HL_1=ggplot(k1,aes(f0_1,log(2)/d_1))+
        geom_line()+
        ylab("HL_1")+
        coord_cartesian(ylim=c(0,24))
    HL_2=ggplot(k2,aes(f0_2,log(2)/d_2))+
        geom_line()+
        ylab("HL_2")+
        coord_cartesian(ylim=c(0,24))
    HL_lfc=ggplot(df,aes(f0_1,f0_2,fill=log2(d_2/d_1)))+
        geom_tile()+
        scale_fill_viridis_c("log2(HL_1/HL_2)")

    t0_1=ggplot(k1,aes(f0_1,t0_1))+
        geom_line()+
        coord_cartesian(ylim=c(-24,0))
    t0_2=ggplot(k2,aes(f0_2,t0_2))+
        geom_line()+
        coord_cartesian(ylim=c(-24,0))
    t0_lfc=ggplot(df,aes(f0_1,f0_2,fill=log2(t0_1/t0_2)))+
        geom_tile()+
        scale_fill_viridis_c()

    (s_1 |  s_2 | s_lfc)  / (HL_1 |  HL_2 | HL_lfc) / (t0_1 |  t0_2 | t0_lfc)
}


#' Plot the abundance of new and old RNA and the fitted model over time for a single gene.
#'
#' For each \code{\link{Condition}} there will be one panel containing the values and the corresponding model fit.
#'
#' @param data a grandR object
#' @param gene the gene to be plotted
#' @param slot the data slot of the observed abundances
#' @param time the labeling duration column in the column annotation table
#' @param title the plot title
#' @param type how to fit the model (see link{FitKinetics})
#' @param bare.plot omit axis and main title?
#' @param exact.tics use axis labels directly corresponding to the labeling durations availanle?
#' @param return.tables also return the tables used for plotting
#' @param ... given to the fitting procedures
#'
#' @return either a ggplot object, or a list containing all tables used for plotting and the ggplot object.
#'
#' @seealso \link{FitKineticsGeneNtr}, \link{FitKineticsGeneLeastSquares}, \link{FitKineticsGeneLogSpaceLinear}
#'
#' @export
#'
#' @examples
#'
#' # straight-forward plot of a random gene
#' PlotGeneKinetics(sars,gene,type="ntr")
#'
#' # plot four ways of fitting the model: non-linear least squares, ntr fit, nlls fit without using the 3hpi time point for SARS, and nlls fit without assuming steady state for  SARS
#' (PlotGeneKinetics(sars,"SRSF6",bare.plot=T) |
#'   PlotGeneKinetics(sars,"SRSF6",type = "ntr",bare.plot=T))  /
#'    (PlotGeneKinetics(sars,"SRSF6",use.old=Coldata(sars)$Name!="SARS.no4sU.A",bare.plot=T) |
#'         PlotGeneKinetics(sars,"SRSF6",steady.state=list(Mock=TRUE,SARS=FALSE),bare.plot=T))
#'         REVIEWED BY Lygeri: changed nlls to nlls
PlotGeneKinetics=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU,title=Genes(data,genes=gene), type=c("nlls","ntr","lm"), bare.plot=FALSE,exact.tics=TRUE,show.CI=FALSE,return.tables=FALSE,...) {
    if (length(ToIndex(data,gene))==0) return(NULL)

    fit=switch(tolower(type[1]),
               ntr=FitKineticsGeneNtr(data,gene,slot=slot,time=time,...),
               nlls=FitKineticsGeneLeastSquares(data,gene,slot=slot,time=time,...),
               lm=FitKineticsGeneLogSpaceLinear(data,gene,slot=slot,time=time,...)
    )
    if (is.null(Coldata(data)$Condition)) fit=setNames(list(fit),gene)
    df=rbind(
        cbind(GetData(data,mode.slot=paste0("total.",slot),genes=gene),Type="Total"),
        cbind(GetData(data,mode.slot=paste0("new.",slot),genes=gene,ntr.na = FALSE),Type="New"),
        cbind(GetData(data,mode.slot=paste0("old.",slot),genes=gene,ntr.na = FALSE),Type="Old")
    )
    if (show.CI) {
        if (!all(c("lower","upper") %in% Slots(data))) stop("Compute lower and upper slots first! (ComputeNtrCI)")
        df=cbind(df,rbind(
            setNames(cbind(GetData(data,mode.slot=paste0("total.",slot),genes=gene,coldata=FALSE),
                     GetData(data,mode.slot=paste0("total.",slot),genes=gene,coldata=FALSE)),
                     c("lower","upper")),
            setNames(cbind(
                GetData(data,mode.slot=paste0("total.",slot),genes=gene,coldata=FALSE)*GetData(data,mode.slot="lower",genes=gene,coldata=FALSE),
                GetData(data,mode.slot=paste0("total.",slot),genes=gene,coldata=FALSE)*GetData(data,mode.slot="upper",genes=gene,coldata=FALSE)),
                c("lower","upper")),
            setNames(cbind(
                GetData(data,mode.slot=paste0("total.",slot),genes=gene,coldata=FALSE)*(1-GetData(data,mode.slot="upper",genes=gene,coldata=FALSE)),
                GetData(data,mode.slot=paste0("total.",slot),genes=gene,coldata=FALSE)*(1-GetData(data,mode.slot="lower",genes=gene,coldata=FALSE))),
                c("lower","upper"))
        ))
    }
    df$time=df[[time]]


    if (tolower(type[1])=="ntr") {
        #fac=unlist(lapply(as.character(df$Condition),function(n) fit[[n]]$f0))/df$Value[df$Type=="Total"]
        fac=unlist(lapply(1:nrow(df),function(i) {
            n=as.character(df$Condition)[i]
            tt=df$time[i]
            f.old.nonequi(tt,fit[[n]]$f0,fit[[n]]$Synthesis,fit[[n]]$Degradation)+f.new(tt,fit[[n]]$Synthesis,fit[[n]]$Degradation)
        }))/df$Value[df$Type=="Total"]
        df$Value=df$Value*fac
        if (show.CI) {
            df$lower=df$lower*fac
            df$upper=df$upper*fac
        }
    }

    df$Condition=if ("Condition" %in% names(df)) df$Condition else gene
    tt=seq(0,max(df$time),length.out=100)
    df.median=ddply(df,c("Condition","Type","duration.4sU","time"),function(s) data.frame(Value=median(s$Value)))
    fitted=ldply(fit,function(f) data.frame(time=c(tt,tt),Value=c(f.old.nonequi(tt,f$f0,f$Synthesis,f$Degradation),f.new(tt,f$Synthesis,f$Degradation)),Type=rep(c("Old","New"),each=length(tt))),.id="Condition")
    breaks=if (exact.tics) sort(unique(df[[time]])) else scales::breaks_extended(5)(df[[time]])

    g=ggplot(df,aes(time,Value,color=Type))
    if (show.CI) g=g+geom_errorbar(mapping=aes(ymin=lower,ymax=upper),width=0.1)
    g=g+geom_point()+
        geom_line(data=df.median[df.median$Type=="Total",])+
        scale_x_continuous(if (bare.plot) NULL else "4sU labeling",labels = scales::number_format(accuracy = max(0.01,my.precision(breaks)),suffix="h"),breaks=breaks)+
        scale_color_manual("RNA",values=c(Total="gray",New="red",Old="blue"),guide=if (bare.plot) "none" else "legend")+
        ylab(if (bare.plot) NULL else "Expression")+
        geom_line(data=fitted,aes(ymin=NULL,ymax=NULL),linetype=2)
    if (!is.null(Coldata(data)$Condition)) g=g+facet_wrap(~Condition,nrow=1)
    if (!bare.plot) g=g+ggtitle(title)
    if (return.tables) list(gg=g,df=df,df.median=df.median,fitted=fitted) else g
}




#' Simulate the kinetics of old and new RNA for given parameters.
#'
#' The standard mass action kinetics model of gene expression arises from the differential equation
#' \eqn{df/dt = s - d  f(t)}, with s being the constant synthesis rate, d the constant degradation rate and \eqn{f0=f(0)} (the abundance at time 0).
#' The RNA half-life is directly related to d via \eqn{HL=log(2)/d}.
#' This model dictates the time evolution of old and new RNA abundance after metabolic labeling starting at time t=0.
#' This function simulates data according to this model.
#'
#' @param s the synthesis rate
#' @param d the degradation rate
#' @param hl the RNA half-life
#' @param f0 the abundance at time t=0
#' @param min.time the start time to simulate
#' @param max.time the end time to simulate
#' @param N how many time points from min.time to max.time to simuate
#' @param name add a Name column to the resulting data frame
#' @param out which values to put into the data frame
#'
#' @return a data frame containing the simulated values
#' @export
#'
#' @seealso \link{PlotSimulation} for plotting the simulation
#'
#' @examples
#' SimulateKinetics(hl=2)   # simulate steady state kinetics for an RNA with half-life 2h
#'
#' @examples
#'
SimulateKinetics=function(s=100*d,d=log(2)/hl,hl=2,f0=s/d,min.time=-1,max.time=10,N = 1000,name=NULL,out=c("Old","New","Total","NTR")) {
    times=seq(min.time,max.time,length.out=N)
    ode.new=function(t,s,d) ifelse(t<0,0,s/d*(1-exp(-t*d)))
    ode.old=function(t,f0,s,d) ifelse(t<0,f0,f0*exp(-t*d))
    old=ode.old(times,f0,s,d)
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


#' Plot simulated data
#'
#' The input data is usually created by \code{\link{SimulateKinetics}}
#'
#' @param sim.df the input data frame
#' @param ntr show the ntr?
#' @param old show old RNA?
#' @param new show new RNA?
#' @param total show total RNA?
#'
#' @return a ggplot object
#'
#' @seealso \link{SimulateKinetics} for creating the input data frame
#' @export
#'
#' @examples
#' PlotSimulation(SimulateKinetics(hl=2))
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

#PlotCompareKinetics=function(type="NTR",...) {
#    dfs=list(...)
#    df=do.call("rbind",dfs)
#    ggplot(df[df$Type %in% type,],aes(Time,Value,color=Name))+
#        geom_line(size=1)+
#        facet_wrap(~Type)
#}

#sim.hl8.steady=SimulateKinetics(hl=8,name="Steady-state, 8h")
#sim.hl4.steady=SimulateKinetics(hl=4,name="Steady-state, 4h")
#sim.hl2.steady=SimulateKinetics(hl=2,name="Steady-state, 2h")
#sim.hl1.steady=SimulateKinetics(hl=1,name="Steady-state, 1h")
#sim.hl8.2x=SimulateKinetics(l0=100,hl=4,s=100*log(2)/8,name="2x down, 8h")
#sim.hl4.2x=SimulateKinetics(l0=100,hl=2,s=100*log(2)/4,name="2x down, 4h")
#sim.hl2.2x=SimulateKinetics(l0=100,hl=1,s=100*log(2)/2,name="2x down, 2h")
#sim.hl1.2x=SimulateKinetics(l0=100,hl=0.5,s=100*log(2)/1,name="2x down, 1h")
#sim.hl8.10x=SimulateKinetics(l0=100,hl=0.8,s=100*log(2)/8,name="10x down, 8h")
#sim.hl4.10x=SimulateKinetics(l0=100,hl=0.4,s=100*log(2)/4,name="10x down, 4h")
#sim.hl2.10x=SimulateKinetics(l0=100,hl=0.2,s=100*log(2)/2,name="10x down, 2h")
#sim.hl1.10x=SimulateKinetics(l0=100,hl=0.1,s=100*log(2)/1,name="10x down, 1h")



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
