

#' Fit kinetics using pulseR
#'
#' @param data A grandR object
#' @param name the user defined analysis name to store the results
#' @param time The column in the column annotation table representing the labeling duration
#'
#' @details This is adapted code from https://github.com/dieterich-lab/ComparisonOfMetabolicLabeling
#'
#' @return a new grandR object containing the pulseR analyses in a new analysis table
#' @export
#'
#' @concept kinetics
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

    deseq=function (x, loggeomeans)
    {
      finitePositive <- is.finite(loggeomeans) & x > 0
      if (any(finitePositive)) {
        res <- exp(median((log(x) - loggeomeans)[finitePositive],
                          na.rm = TRUE))
      }
      else {
        print(utils::head(x))
        stop("Can't normalise accross a condition.\n         Too many zero expressed genes. ")
      }
      res
    }
    findDeseqFactorsSingle=function (count_data)
    {
      loggeomeans <- rowMeans(log(count_data))
      deseqFactors <- apply(count_data, 2, deseq, loggeomeans = loggeomeans)
      deseqFactors
    }

    norms <- findDeseqFactorsSingle(totals)
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
            formulas = eval(substitute(alist(
              unlabelled = exp(mu1) + exp(mu2) + exp(mu3) * (1 + exp(-d * time)),
              labelled = exp(mu2) + exp(mu3) * (1 - exp(-d * time)))
            )),
#              pulseR::MeanFormulas(
#                unlabelled = exp(mu1) + exp(mu2) + exp(mu3) * (1 + exp(-d * time)),
#                labelled = exp(mu2) + exp(mu3) * (1 - exp(-d * time))),
            formulaIndexes = list(
                unlabelled = "unlabelled",
                labelled = "labelled"))

        pd <- do.call(getExportedValue("pulseR","PulseData"),list(
            counts,
            conditions[c("fraction", "time")], #  data.frame; the first column corresponds to the conditions given in formulas
            formulas$formulas,
            formulas$formulaIndexes,
            groups=~fraction+time
        ))
        pd$depthNormalisation <- norms1

        setOpts <- function(bounds, tolerance, cores = parallel::detectCores()-2, replicates = 5, normFactors = NULL) {
            opts <- do.call(getExportedValue("pulseR","setFittingOptions"),list(verbose = "verbose"))
            opts$cores <- cores
            opts$replicates <- replicates
            ## if rt-conversion data, we do not fit normalisation coefficients, because they
            ## are derived from the DESeq-like normalisation
            if (is.null(normFactors)) { opts$fixedNorms <- TRUE }
            opts <- do.call(getExportedValue("pulseR","setBoundaries"),list(bounds, normFactors = normFactors, options = opts))
            opts <- do.call(getExportedValue("pulseR","setTolerance"),list(
                params = tolerance$params,
                normFactors = tolerance$normFactors,
                logLik = tolerance$logLik,
                options = opts
            ))
            if (is.null(tolerance$normFactors)) {
                opts$tolerance$normFactors <- NULL
            }
            opts
        }

        opts <- setOpts(boundaries, tolerance)
        initf <- match.fun(init)
        initPars <- initf(pd$counts) # use pulseData here
        fit <- do.call(getExportedValue("pulseR","fitModel"),list(pd, initPars, opts))
        fit$d

 #       cis <- pulseR::ciGene("d",
 #                     par = re$fit,
 #                     geneIndexes = seq_along(re$fit$d),
 #                     pd = re$pd,
 #                     options = re$opts, confidence = CI.size)
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
#' @concept kinetics
f.old.equi=function(t,s,d) s/d*exp(-t*d)

#' @describeIn f.old.equi abundance of old RNA without assuming steady state
#' @export
f.old.nonequi=function(t,f0,s,d) f0*exp(-t*d)
#' @describeIn f.old.equi abundance of new RNA (steady state does not matter)
#' @export
f.new=function(t,s,d) s/d*(1-exp(-t*d))



#' Function to compute the abundance of new or old RNA at time t for non-constant rates.
#'
#' The standard mass action kinetics model of gene expression arises from the differential equation
#' \eqn{df/dt = s(t) - d(t)  f(t)}, with s(t) being the synthesis rate at time t, d(t) the degradation rate at time t and \eqn{f0=f(0)} (the abundance at time 0).
#' Here, both s and d have the following form \eqn{s(t)=so+sf \cdot t^{se}}.
#'
#'
#' @param t time in h (can be a vector)
#' @param so synthesis date offset
#' @param sf synthesis date factor
#' @param se synthesis date exponent
#' @param do degradation rate offset
#' @param df degradation rate factor
#' @param de degradation rate exponent
#' @param f0 the abundance at time t=0
#'
#' @return the RNA abundance at time t
#' @seealso \link{f.nonconst}
#'
#' @export
#' @concept kinetics
f.nonconst.linear=function(t,f0,so,sf,se,do,df,de) {
  if (length(t)>1) return(sapply(t,function(tt) f.nonconst.linear(tt,f0,so,sf,se,do,df,de)))

  ee=exp(-t*(df*t^de+do*de+do)/(de+1))
  ff=function(x) exp(x*((df*x^de)/(de+1)+do))*(so+sf*x^se)
  ii=integrate(ff,0,t)$value
  ee*(ii+f0)
}

#' Function to compute the abundance of new or old RNA at time t for non-constant rates.
#'
#' The standard mass action kinetics model of gene expression arises from the differential equation
#' \eqn{df/dt = s(t) - d(t)  f(t)}, with s(t) being the synthesis rate at time t, d(t) the degradation rate at time t and \eqn{f0=f(0)} (the abundance at time 0).
#' Here, both s and d have the following form \eqn{s(t)=so+sf \cdot t^{se}}.
#'
#'
#' @param t time in h (can be a vector)
#' @param f0 the abundance at time t=0
#' @param s the synthesis rate (see details)
#' @param d the degradation rate (see details)
#'
#' @details Both rates can be either (i) a single number (constant rate), (ii) a data frame with names "offset",
#' "factor" and "exponent" (for linear functions, see \link{ComputeNonConstantParam}; only one row allowed) or
#' (iii) a unary function time->rate. Functions
#'
#' @return the RNA abundance at time t
#' @seealso \link{f.nonconst.linear}
#'
#' @export
#' @concept kinetics
f.nonconst=function(t,f0,s,d) {

  if (is.numeric(s)) s = ComputeNonConstantParam(start=s)
  if (is.numeric(d)) d = ComputeNonConstantParam(start=d)
  if (is.data.frame(s) && nrow(s)!=1) stop("Only a single row data frame allowed!")
  if (is.data.frame(d) && nrow(d)!=1) stop("Only a single row data frame allowed!")

  if (is.function(s) || is.function(d)) {
    checkPackages("deSolve")

    sfun=if (is.data.frame(s)) function(t) s$offset+s$factor*t^s$exponent else s
    dfun=if (is.data.frame(d)) function(t) d$offset+d$factor*t^d$exponent else d

    re=deSolve::ode(y=f0,times=c(0,t),func=function(t,y,parms,...) list(sfun(t)-dfun(t)*y))
    return(re[-1,2])
  }

  so=s$offset
  sf=s$factor
  se=s$exponent
  do=d$offset
  df=d$factor
  de=d$exponent

  sapply(t,function(tt) {
  ee=exp(-tt*(df*tt^de+do*de+do)/(de+1))
  ff=function(x) exp(x*((df*x^de)/(de+1)+do))*(so+sf*x^se)
  ii=integrate(ff,0,tt)$value
  ee*(ii+f0)
  })
}

#' Compute and evaluate functions for non constant rates
#'
#' For simplicity, non constant rates here have the following form $o+f*t^e$.
#'
#'
#' @param start the value at t=0
#' @param end the value at t=end.time
#' @param exponent the exponent (e above)
#' @param end.time the end time
#' @param t vector of times
#' @param param output of \code{ComputeNonConstantParam()}, only a single row!
#'
#' @return data frame containing either the parameters o, f and e (ComputeNonConstantParam), or containing the value of $o+f*t^e$ for the given times (EvaluateNonConstantParam).
#'
#' @describeIn ComputeNonConstantParam compute a data frame containing the parameters for non constant rates
#' @export
#' @concept kinetics
ComputeNonConstantParam=function(start,end=start,exponent=1,end.time=2) data.frame(offset=start,factor=(end-start)/end.time^exponent,exponent=exponent)
#' @describeIn ComputeNonConstantParam compute a data frame containing the rates for the given parameter set (computed from \code{ComputeNonConstantParam})
#' @export
EvaluateNonConstantParam=function(t,param) data.frame(t=t,value=param$offset+param$factor*t^param$exponent)


#' Fit kinetic models to all genes.
#'
#' Fit the standard mass action kinetics model of gene expression by different methods. Some methods require steady state assumptions, for others
#' data must be properly normalized. The parameters are fit per \link{Condition}.
#'
#' @param data A grandR object
#' @param name.prefix the prefix of the analysis name to be stored in the grandR object
#' @param type Which method to use (either one of "full","ntr","lm", "chase")
#' @param slot The data slot to take expression values from
#' @param time The column in the column annotation table representing the labeling duration
#' @param CI.size A number between 0 and 1 representing the size of the confidence interval
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
#' @details Pulse-chase designs are fit using \link{FitKineticsGeneLeastSquares} while only considering the drop of labeled RNA. Note that in this case the notion "new" / "old" RNA is misleading,
#' since labeled RNA corresponds to pre-existing RNA!
#'
#' @details This function is flexible in what to put in the analysis table. You can specify the statistics using return.fields and return.extra (see \code{\link{kinetics2vector}})
#'
#' @seealso \link{FitKineticsGeneNtr}, \link{FitKineticsGeneLeastSquares}, \link{FitKineticsGeneLogSpaceLinear}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Cell",Design$dur.4sU,Design$Replicate))
#' sars <- FilterGenes(sars,use=1:10)
#' sars<-FitKinetics(sars,name="kinetics.ntr",type='ntr')
#' sars<-Normalize(sars)
#' sars<-FitKinetics(sars,name="kinetics.nlls",type='nlls')
#' sars<-FitKinetics(sars,name="kinetics.lm",type='lm')
#' head(GetAnalysisTable(sars,columns="Half-life"))
#'
#' @export
#'
#' @concept kinetics
FitKinetics=function(data,name.prefix="kinetics",type=c("nlls","ntr","lm","chase"),slot=DefaultSlot(data),time=Design$dur.4sU,CI.size=0.95,return.fields=c("Synthesis","Half-life"),return.extra=NULL,...) {

  fun=switch(tolower(type[1]),
             ntr=FitKineticsGeneNtr,
             nlls=FitKineticsGeneLeastSquares,
             chase=function(...) FitKineticsGeneLeastSquares(...,chase=TRUE),
             lm=FitKineticsGeneLogSpaceLinear)

  if (is.null(fun)) stop(sprintf("Type %s unknown!",type))
  result=plapply(Genes(data),
                 fun,
                 data=data,
                 slot=slot,time=time,CI.size=CI.size,
                 ...)

  adder=function(cond) {
    tttr=function(re) {
      if (is.matrix(re)) as.data.frame(t(re)) else setNames(data.frame(re),names(kinetics2vector(result[[1]],condition=cond,return.fields=return.fields,return.extra=return.extra)))
    }
    slam.param=tttr(sapply(result,kinetics2vector,condition=cond,return.fields=return.fields,return.extra=return.extra))
    rownames(slam.param)=Genes(data,use.symbols = FALSE)
    name=if (is.null(name.prefix)) cond else if(is.null(cond)) name.prefix else paste0(name.prefix,".",cond)
    AddAnalysis(data,name=name,slam.param)
  }

  if (!is.null(Condition(data))) {
    for (cond in levels(Condition(data))) data=adder(cond)
  } else {
    data=adder(NULL)
  }
  data
}



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
#' @param chase is this a pulse-chase experiment? (see details)
#' @param CI.size A number between 0 and 1 representing the size of the confidence interval
#' @param steady.state either a named list of logical values representing conditions in steady state or not, or a single logical value for all conditions
#' @param use.old a logical vector to exclude old RNA from specific time points
#' @param use.new a logical vector to exclude new RNA from specific time points
#' @param maxiter the maximal number of iterations for the Levenberg-Marquardt algorithm used to minimize the least squares
#' @param compute.residuals set this to TRUE to compute the residual matrix
#'
#' @return
#' A named list containing the model fit:
#' \itemize{
#'   \item{data: a data frame containing the observed value used for fitting}
#'   \item{residuals: the computed residuals if compute.residuals=TRUE, otherwise NA}
#'   \item{Synthesis: the synthesis rate (in U/h, where U is the unit of the slot)}
#'   \item{Degradation: the degradation rate (in 1/h)}
#'   \item{Half-life: the RNA half-life (in h, always equal to log(2)/degradation-rate}
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
#' @details For pulse-chase designs, only the drop of the labeled RNA is considered. Note that in this case the notion "new" / "old" RNA is misleading,
#' since labeled RNA corresponds to pre-existing RNA!
#'
#' @seealso \link{FitKinetics}, \link{FitKineticsGeneLogSpaceLinear}, \link{FitKineticsGeneNtr}
#'
#' @examples
#' sars <- ReadGRAND(system.file("extdata", "sars.tsv.gz", package = "grandR"),
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' sars <- Normalize(sars)
#' FitKineticsGeneLeastSquares(sars,"SRSF6",steady.state=list(Mock=TRUE,SARS=FALSE))
#'
#' @export
#'
#' @concept kinetics
FitKineticsGeneLeastSquares=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU,chase=FALSE,CI.size=0.95,steady.state=NULL,use.old=TRUE,use.new=TRUE, maxiter=250, compute.residuals=TRUE) {

    # residuals of the functions for usage with nls.lm
    res.fun.equi=function(par,old,new) {
        s=par[1]
        d=par[2]
        f=function(t) f.old.equi(t,s,d)
        g=function(t) f.new(t,s,d)
        c(old$Value-f(old$time),new$Value-g(new$time))
    }
    res.fun.chase=function(par,old,new) {
      s=par[1]
      d=par[2]
      f=function(t) f.old.equi(t,s,d)
      c(rep(0,nrow(old)),new$Value-f(new$time))
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

    lvl = median(GetData(data,mode.slot=slot,genes=gene,ntr.na = FALSE)$Value)

    newdf=GetData(data,mode.slot=if (chase) "ntr" else paste0("new.",slot),genes=gene,ntr.na = FALSE)
    newdf$use=1:nrow(newdf) %in% (1:nrow(newdf))[use.new]
    newdf$time=newdf[[time]]
    if (is.null(newdf$Condition)) {
        newdf$Condition="Data"
        if (length(steady.state)==1) names(steady.state)="Data"
    }
    if (chase) newdf=newdf[!newdf$no4sU,]
    newdf=plyr::dlply(newdf,"Condition",function(s) correct(s))

    correct=function(s) {
        if (max(s$Value)==0) s$Value=0.01
        s$Type="Old"
        s
    }
    olddf=GetData(data,mode.slot=paste0("old.",slot),genes=gene,ntr.na = FALSE)
    olddf$use=1:nrow(olddf) %in% (1:nrow(olddf))[use.old]
    olddf$time=olddf[[time]]
    if (is.null(olddf$Condition)) olddf$Condition="Data"
    if (chase) {
      olddf=olddf[!olddf$no4sU,]
      olddf$use=FALSE
    }
    olddf=plyr::dlply(olddf,"Condition",function(s) correct(s))


    fit.equi=function(ndf,odf) {

        res.fun = if (chase) res.fun.chase else res.fun.equi
        tinit=if (sum(ndf$Value>0)==0) min(ndf$time[ndf$time>0]) else min(ndf$time[ndf$Value>0])
        init.d=mean(-log(1-(0.1+ndf[ndf$time==tinit,"Value"])/(0.2+ndf[ndf$time==tinit,"Value"]+odf[odf$time==tinit,"Value"])))
        init.s=init.d*mean(odf[odf$time==tinit,"Value"])
        model.p=minpack.lm::nls.lm(c(init.s,init.d),lower=c(0,0.01),
                       fn=res.fun,
                       old=odf[odf$use,],
                       new=ndf[ndf$use,],
                       control=minpack.lm::nls.lm.control(maxiter = maxiter))

        if (model.p$niter==maxiter) return(list(data=NA,residuals=if (compute.residuals) data.frame(Name=c(as.character(odf$Name[odf$use]),as.character(ndf$Name[ndf$use])),Type=c(rep("old",nrow(ndf)),rep("new",nrow(odf))),Absolute=NA,Relative=NA) else NA,Synthesis=NA,Degradation=NA,`Half-life`=NA,conf.lower=c(NA,NA),conf.upper=c(NA,NA),f0=NA,logLik=NA,rmse=NA, rmse.new=NA, rmse.old=NA,total=NA,type="equi"))
        conf.p=try(confint.nls.lm(model.p,level=CI.size),silent = TRUE)
        if (!is.matrix(conf.p)) conf.p=matrix(c(NA,NA,NA,NA),nrow=2)
        par=setNames(model.p$par,c("s","d"))
        rmse=sqrt(model.p$deviance/(nrow(ndf)+nrow(odf)))
        fvec=res.fun(par,odf,ndf)
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

            resi=res.fun(model.p$par,odf,ndf)
            modval=c(odf[,"Value"],ndf[,"Value"])-resi
            residuals=data.frame(Name=c(as.character(odf$Name),as.character(ndf$Name)),Type=c(rep("old",nrow(ndf)),rep("new",nrow(odf))),Absolute=resi,Relative=resi/modval)
        }
        total=sum(ndf$Value)+sum(odf$Value)
        syn = if (chase) lvl*unname(par['d']) else unname(par['s']) # in chase designs, at time 0 new RNA might not be at steady state level!

        list(data=df,
             residuals=residuals,
             Synthesis=syn,
             Degradation=unname(par['d']),
             `Half-life`=log(2)/unname(par['d']),
             conf.lower=c(Synthesis=max(0,unname(conf.p[1,1])),Degradation=max(0,unname(conf.p[2,1])),`Half-life`=unname(log(2)/(conf.p[2,2]))),
             conf.upper=c(Synthesis=unname(conf.p[1,2]),Degradation=unname(conf.p[2,2]),`Half-life`=unname(log(2)/max(0,(conf.p[2,1])))),
             f0=unname(syn/par['d']),
             logLik=logLik.nls.lm(model.p),
             rmse=rmse, rmse.new=rmse.new, rmse.old=rmse.old,
             total=total,type="equi")
    }
    fit.nonequi=function(ndf,odf) {
        if (chase) stop("All conditions must be at steady state for pulse chase designs!")

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
        conf.m=try(confint.nls.lm(model.m,level=CI.size),silent = TRUE)
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
            resi=res.fun.nonequi(model.m$par,odf,ndf)
            modval=c(odf[,"Value"],ndf[,"Value"])-resi
            residuals=data.frame(Name=c(as.character(odf$Name),as.character(ndf$Name)),Type=c(rep("old",nrow(ndf)),rep("new",nrow(odf))),Absolute=resi,Relative=resi/modval)
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
#' @param CI.size A number between 0 and 1 representing the size of the confidence interval
#'
#' @return
#' A named list containing the model fit:
#' \itemize{
#'   \item{data: a data frame containing the observed value used for fitting}
#'   \item{Synthesis: the synthesis rate (in U/h, where U is the unit of the slot)}
#'   \item{Degradation: the degradation rate (in 1/h)}
#'   \item{Half-life: the RNA half-life (in h, always equal to log(2)/degradation-rate}
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
#' FitKineticsGeneLogSpaceLinear(sars,"SRSF6")   # fit per condition
#'
#' @export
#'
#' @concept kinetics
FitKineticsGeneLogSpaceLinear=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU,CI.size=0.95) {
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
    olddf=plyr::dlply(olddf,"Condition",function(s) correct(s))


    fit.lm=function(odf) {
        odf=odf[odf$Value>0,]
        fit=lm(log(Value)~time,data=odf)
        summ=summary(fit)

        par=setNames(c(exp(coef(fit)["(Intercept)"])*-coef(fit)["time"],-coef(fit)["time"]),c("s","d"))
        if (sum(abs(residuals(fit)))>0) {
            conf.p=confint(fit,level=CI.size)
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
             conf.lower=c(Synthesis=unname(conf.p[1,1]),Degradation=unname(conf.p[2,2]),`Half-life`=unname(log(2)/conf.p[2,1])),
             conf.upper=c(Synthesis=unname(conf.p[1,2]),Degradation=unname(conf.p[2,1]),`Half-life`=unname(log(2)/conf.p[2,2])),
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
#' @param CI.size A number between 0 and 1 representing the size of the credible interval
#' @param transformed.NTR.MAP Use the transformed NTR MAP estimator instead of the MAP of the transformed posterior
#' @param exact.ci compute exact credible intervals (see details)
#' @param total.fun use this function to summarize the expression values (only relevant for computing the synthesis rate s)
#'
#' @return
#' A named list containing the model fit:
#' \itemize{
#'   \item{data: a data frame containing the observed value used for fitting}
#'   \item{Synthesis: the synthesis rate (in U/h, where U is the unit of the slot)}
#'   \item{Degradation: the degradation rate (in 1/h)}
#'   \item{Half-life: the RNA half-life (in h, always equal to log(2)/degradation-rate}
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
#'                   design=c("Condition",Design$dur.4sU,Design$Replicate))
#' sars <- Normalize(sars)
#' sars <- subset(sars,columns=Condition=="Mock")
#' FitKineticsGeneNtr(sars,"SRSF6")
#'
#' @export
#'
#' @concept kinetics
FitKineticsGeneNtr=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU,CI.size=0.95,transformed.NTR.MAP=TRUE,exact.ci=FALSE,total.fun=median) {
    if (!all(c("alpha","beta") %in% Slots(data))) stop("Beta approximation data is not available in grandR object!")
    steady.state=TRUE

    crit <- qchisq(1-CI.size, df = 2, lower.tail = FALSE) / 2

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
    loglik=if (transformed.NTR.MAP) function(d,a,b,t) sum((a-1)*log1p(-exp(-t*d))-t*d*(b-1)) else function(d,a,b,t) sum((a-1)*log1p(-exp(-t*d))-t*d*b)

    t=a[,time]
    use=t>0
    names=a$Name[use]
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
        if (exact.ci) {
            plik=function(x) pmax(0,sapply(x,function(xx) exp(loglik(xx,a=a[ind],b=b[ind],t=t[ind])-max)))
            inte = function(u) integrate(plik,pbounds[1],u)$value
            lower1=uniroot.save(function(x) ploglik(x)-max+log(1E6),lower = pbounds[1],upper=d)
            upper1=uniroot.save(function(x) ploglik(x)-max+log(1E6),lower = d,upper=pbounds[2])

            total.inte=inte(upper1)
            lower=uniroot.save(function(x) inte(x)/total.inte-(1-CI.size)/2,lower = lower1,upper=d)
            upper=uniroot.save(function(x) inte(x)/total.inte-(1+CI.size)/2,lower = d,upper=upper1)
        } else {
            lower=uniroot.save(function(x) ploglik(x)-max+crit,lower = pbounds[1],upper=d)
            upper=uniroot.save(function(x) ploglik(x)-max+crit,lower = d,upper=pbounds[2])
        }

        rmse=sum(sqrt(((1-exp(-t[ind]*d)-ntr[ind]))^2))/sum(ind)
        if (!equi) {
            s=median(ntr/(1-exp(-t*d))*d*total)
            s.low=median(ntr/(1-exp(-t*lower))*lower*total)
            s.up=median(ntr/(1-exp(-t*upper))*upper*total)
        } else {
            s=f0*d
            s.low=f0*lower
            s.up=f0*upper
        }

        df=data.frame(Name=names,alpha=a,beta=b,t=t)
        df$Quantile = pbeta(1-exp(-df$t*d),df$alpha,df$beta)
        df$Expected.NTR=1-exp(-df$t*d)
        df$Observed.NTR=(df$alpha-1)/(df$alpha+df$beta-2)
        df$Absolute.Residual = df$Observed.NTR-df$Expected.NTR
        df$Relative.Residual = df$Absolute.Residual/df$Expected.NTR

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
#' @param CI.size The level for confidence intervals
#' @param compute.confidence should CIs be computed or not?
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
#' @export
#'
#' @concept recalibration
CalibrateEffectiveLabelingTimeKineticFit=function(data,slot=DefaultSlot(data),time=Design$dur.4sU,time.name="calibrated_time",time.conf.name="calibrated_time_conf",CI.size=0.95,compute.confidence=FALSE,n.estimate=1000, n.iter=10000, verbose=FALSE,...) {

    conds=Coldata(data)
    if (is.null(conds$Condition)) {
        conds$Condition=factor("Data")
    }
    re=matrix(NA,ncol=2,nrow=nrow(conds),dimnames = list(conds$Name,c(time.name,time.conf.name)))
    for (cond in levels(conds$Condition)) {
        if (verbose) cat(sprintf("Calibrating %s...\n",cond))
        sub=subset(data,columns=conds$Condition==cond)
        Condition(sub)=NULL
        sub=DropAnalysis(sub)

        # restrict to top n genes stratified by ntr
        totals=rowSums(GetTable(sub,type=slot))
        sub=FitKinetics(sub,type='nlls',slot=slot,time=time,return.fields="Half-life")
        HLs=GetAnalysisTable(sub,prefix.by.analysis = FALSE)$`Half-life`
        HL.cat=cut(HLs,c(seq(0,2*max(Coldata(data,time)),length.out=5),Inf),include.lowest = TRUE)
        fil=plyr::ddply(data.frame(Gene=Genes(sub),totals,HL.cat),plyr::.(HL.cat),function(s) {
           threshold=sort(s$totals,decreasing = TRUE)[min(length(s$totals),ceiling(n.estimate/length(unique(HL.cat))))]
           data.frame(Gene=s$Gene,use=s$totals>=threshold)
        })

        genes=as.character(fil$Gene[fil$use])
        sub=FilterGenes(sub,use=genes)
        sub=DropAnalysis(sub)

        n.eval=0
        opt.fun=function(times) {
            tt=init
            tt[use]=times
            Coldata(sub,Design$dur.4sU)=tt
            fit=FitKinetics(sub,return.fields='logLik',slot=slot,...)
        #    fit=FitKinetics(sub,return.fields='logLik',slot=slot,steady.state=steady.state[[cond]])
            re=sum(GetAnalysisTable(fit,gene.info=FALSE)[,1])
        #    fit=FitKinetics(sub,return.fields='Half-life',compute.residuals=TRUE,slot=slot,steady.state=steady.state[[cond]],return.extra = function(s) setNames(s$residuals$Relative,paste0("Residuals.",s$residuals$Name))[s$residuals$Type=="new"])
          #  fit=FitKinetics(sub,return.fields='Half-life',compute.residuals=TRUE,slot=slot,steady.state=steady.state[[cond]],return.extra = function(s) setNames(s$residuals$Relative,paste0("Residuals.",s$residuals$Name))[s$residuals$Type=="new"],...)
          #  hl=GetAnalysisTable(fit,gene.info = FALSE)[,1]
          #  res=GetAnalysisTable(fit,gene.info = FALSE)[,-1][,use]
          #  re=sum(sapply(res,function(rr) cor(hl,rr,method="spearman")^2))
            #cat(times)
            #cat(" ")
            #cat(re)
            #cat("\n")
            n.eval<<-n.eval+1
            if (verbose && n.eval%%10==0) cat(sprintf("Optimization round %d (current solution: %s, loglik: %.2f)...\n",n.eval,paste(sprintf("%.2f",times),collapse = ","),re))

            re
        }

        init=Coldata(sub)[[time]]
        init[init>0]=init[init>0]
        use=init>0 & init<max(init)

        opt.fun2=function(min) -opt.fun(init[use]-min)
        mini=optimize(opt.fun2,interval=c(0,min(init[use])))$minimum
        init[use]=init[use]-mini
        if (verbose) cat(sprintf("First step done.\n"))

        #fit=optim(init[use],fn=opt.fun,hessian=FALSE, control=list(fnscale=-1,maxit=n.iter),method="Nelder-Mead")
        ui=diag(length(init[use]))
        ci=rep(0,length(init[use]))
        fit=constrOptim(init[use],f=opt.fun,ui=ui,ci=ci,hessian=FALSE, control=list(fnscale=-1,maxit=n.iter),method="Nelder-Mead")
        if (fit$convergence!=0) warning(sprintf("Did not converge for %s. Maybe increase n.iter!",cond))

        tt=init
        tt[use]=fit$par
        re[as.character(Coldata(sub)$Name),time.name]=tt
        tt=rep(0,length(tt))

        if (compute.confidence) {
          fit$hessian=numDeriv::hessian(opt.fun,fit$par)
          conf=try(sd.from.hessian(fit$hessian)*qnorm(1-(1-CI.size)/2),silent=TRUE)
          tt[use]=conf
          re[as.character(Coldata(sub)$Name),time.conf.name]=tt
        }
    }
    data=Coldata(data,re)
    data
}

#' Calibrate the effective labeling time by matching half-lives to a .reference
#'
#' The NTRs of each sample might be systematically too small (or large). This function identifies such systematic
#' deviations and computes labeling durations without systematic deviations.
#'
#' @param data A grandR object
#' @param reference.halflives a vector of reference Half-lives named by genes
#' @param reference.columns the reference column description
#' @param slot The data slot to take expression values from
#' @param time.labeling the column in the column annotation table denoting the labeling duration or the labeling duration itself
#' @param time.experiment the column in the column annotation table denoting the experimental time point (can be NULL, see details)
#' @param time.name The name in the column annotation table to put the calibrated labeling durations
#' @param n.estimate the times are calibrated with the top n expressed genes
#' @param verbose verbose output
#'
#' @return
#' A new grandR object containing the calibrated durations in the column data annotation
#'
#' @seealso \link{FitKineticsGeneSnapshot}
#'
#'
#' @export
#'
#' @concept recalibration
CalibrateEffectiveLabelingTimeMatchHalflives=function(data,reference.halflives=NULL,reference.columns=NULL,slot=DefaultSlot(data),time.labeling=Design$dur.4sU,time.experiment=NULL,time.name="calibrated_time",n.estimate=1000,verbose=FALSE) {

  if (any(is.na(Genes(data,names(reference.halflives))))) stop("Not all names of reference.halflives are known gene names!")
  genes=names(reference.halflives)

  totals=rowSums(GetTable(data,type=slot,genes=genes))
  HL.cat=cut(reference.halflives,c(0,2,4,6,8,Inf),include.lowest = TRUE)
  fil=plyr::ddply(data.frame(Gene=genes,totals,HL.cat),plyr::.(HL.cat),function(s) {
    threshold=sort(s$totals,decreasing = TRUE)[min(length(s$totals),ceiling(n.estimate/length(unique(HL.cat))))]
    data.frame(Gene=s$Gene,use=s$totals>=threshold)
  })

  genes=as.character(fil$Gene[fil$use])
  reference.halflives=reference.halflives[genes]

  fit.column=function(column) {
    if (verbose) cat(sprintf("Starting %s...\n",column))
    refcol=reference.columns
    if (is.matrix(refcol)) refcol=apply(refcol[,Columns(data,column),drop=FALSE]==1,1,any)

    df=data.frame(
      ntr=GetTable(data,type ="ntr",columns = column,gene.info = FALSE,genes=genes)[,1],
      total=GetTable(data,type=slot,columns = column,gene.info = FALSE,genes=genes)[,1],
      ss=rowMeans(GetTable(data,type=slot,columns = refcol,gene.info = FALSE,genes=genes)),
      ref.HL=reference.halflives
    )

    if (!is.null(time.experiment) && length(unique(Coldata(data)[refcol,time.experiment]))!=1) stop("Steady state has to refer to a unique experimental time!")

    t=if (is.numeric(time.labeling)) time.labeling else Coldata(data)[column,time.labeling]
    if (t==0) return(0)

    t0=if (is.null(time.experiment)) t else Coldata(data)[column,time.experiment]-unique(Coldata(data)[refcol,time.experiment])
    is.steady.state=refcol[Columns(data,column)]
    if (is.steady.state) return(t)


    fun=function(t) {
      hl=log(2)/TransformSnapshot(df$ntr,df$total,t,t0,f0=df$ss)[,'d']
      re=median(log2(hl/df$ref.HL),na.rm=TRUE)
      #if (verbose) cat(sprintf("t=%.2f, mlfc=%.2f\n",t,re))
      re
    }
    upper=t
    while(fun(upper)<0) {upper=upper*2}
    re=uniroot(fun,interval=c(0,upper))$root
    if (verbose) cat(sprintf("Done  %s!\n",column))
    re
  }

  ctimes=psapply(Columns(data),fit.column)
  Coldata(data,time.name,ctimes)
}


#' Fits RNA kinetics from snapshot experiments
#'
#' Compute the posterior distributions of RNA synthesis and degradation from snapshot experiments for each condition
#'
#' @param data the grandR object
#' @param name.prefix the prefix for the new analysis name; a dot and the column names of the contrast matrix are appended; can be NULL (then only the contrast matrix names are used)
#' @param reference.columns a reference matrix usually generated by \link{FindReferences} to define reference samples for each sample (see details), can be NULL if all conditions are at steady state
#' @param slot the data slot to take f0 and totals from
#' @param conditions character vector of all condition names to estimate kinetics for; can be NULL (i.e. all conditions)
#' @param time.labeling the column in the column annotation table denoting the labeling duration or the labeling duration itself
#' @param time.experiment the column in the column annotation table denoting the experimental time point (can be NULL, see details)
#' @param sample.f0.in.ss whether or not to sample f0 under steady state conditions
#' @param N the sample size
#' @param N.max the maximal number of samples (necessary if old RNA > f0); if more are necessary, a warning is generated
#' @param CI.size A number between 0 and 1 representing the size of the credible interval
#' @param seed Seed for the random number generator
#' @param dispersion overdispersion parameter for each gene; if NULL this is estimated from data
#' @param sample.level Define how the NTR is sampled from the hierarchical Bayesian model (must be 0,1, or 2; see details)
#' @param correct.labeling Labeling times have to be unique; usually execution is aborted, if this is not the case; if this is set to true, the median labeling time is assumed
#' @param verbose Vebose output
#'
#' @details The kinetic parameters s and d are computed using \link{TransformSnapshot}. For that, the sample either must be in steady state
#' (this is the case if defined in the reference.columns matrix), or if the levels at an earlier time point are known from separate samples,
#' so called temporal reference samples. Thus, if s and d are estimated for a set of samples x_1,...,x_k (that must be from the same time point t),
#' we need to find (i) the corresponding temporal reference samples from time t0, and (ii) the time difference between t and t0.
#'
#' @details The temporal reference samples are identified by the reference.columns matrix. This is a square matrix of logicals, rows and columns correspond to all samples
#' and TRUE indicates that the row sample is a temporal reference of the columns sample. This time point is defined by \code{time.experiment}. If \code{time.experiment}
#' is NULL, then the labeling time of the A or B samples is used (e.g. useful if labeling was started concomitantly with the perturbation, and the steady state samples
#' are unperturbed samples).
#'
#' @details By default, the hierarchical Bayesian model is estimated. If sample.level = 0, the NTRs are sampled from a beta distribution
#' that approximates the mixture of betas from the replicate samples. If sample.level = 1, only the first level from the hierarchical model
#' is sampled (corresponding to the uncertainty of estimating the biological variability). If sample.level = 2, the first and second levels
#' are estimated (corresponding to the full hierarchical model).
#'
#' @details if N is set to 0, then no sampling from the posterior is performed, but the transformed MAP estimates are returned
#'
#' @return a new grandR object including new analysis tables (one per condition). The columns of the new analysis table are
#' \itemize{
#'  \item{"s"}{the posterior mean synthesis rate}
#'  \item{"HL"}{the posterior mean RNA half-life}
#'  \item{"s.cred.lower"}{the lower CI boundary of the synthesis rate}
#'  \item{"s.cred.upper"}{the upper CI boundary of the synthesis rate}
#'  \item{"HL.cred.lower"}{the lower CI boundary of the half-life}
#'  \item{"HL.cred.upper"}{the upper CI boundary of the half-life}
#' }
#'
#'
#' @export
#' @concept snapshot
FitKineticsSnapshot=function(data,name.prefix="Kinetics",reference.columns=NULL,
                             slot=DefaultSlot(data),conditions=NULL,
                             time.labeling=Design$dur.4sU,time.experiment=NULL,
                             sample.f0.in.ss=TRUE,N=10000,N.max=N*10,CI.size=0.95,seed=1337, dispersion=NULL,
                             sample.level=2, correct.labeling=FALSE, verbose=FALSE) {
  if (!check.slot(data,slot)) stop("Illegal slot definition!")
  if(!is.null(seed)) set.seed(seed)

  if (is.null(Condition(data))) {
    Condition(data)=""
  }

  if (is.null(conditions)) conditions=levels(Condition(data))

  for (n in conditions) {
    if (verbose) cat(sprintf("Computing snapshot kinetics for %s...\n",n))
    A=Condition(data)==n & !Coldata(data,"no4sU")

    if (is.null(reference.columns)) ss=A
    if (is.matrix(reference.columns)) ss=apply(reference.columns[,Columns(data,A),drop=FALSE]==1,1,any)
    if (length(Columns(data,ss))==0) stop("No reference columns found; check your reference.columns parameter!")
    dispersion = if (sum(ss)==1) rep(0.1,nrow(data)) else if (!is.null(dispersion)) rep(dispersion,length.out=nrow(data)) else estimate.dispersion(GetTable(data,type="count",columns = ss))
    if (verbose) {
      if (any(ss & A)) {
        cat(sprintf("Sampling from steady state for %s...\n",paste(colnames(data)[A],collapse = ",")))
      } else {
        cat(sprintf("Sampling from non-steady state for %s (reference: %s)...\n",paste(colnames(data)[A],collapse = ","),paste(colnames(data)[ss],collapse = ",")))
      }
    }

    # obtain prior from expression values
    get.beta.prior=function(columns,dispersion) {
      ex=rowMeans(GetTable(data,type = slot, columns = columns))
      ex=1/dispersion/(1/dispersion+ex)
      E.ex=mean(ex)
      V.ex=var(ex)
      c(
        shape1=(E.ex*(1-E.ex)/V.ex-1)*E.ex,
        shape2=(E.ex*(1-E.ex)/V.ex-1)*(1-E.ex)
      )
    }
    beta.prior=get.beta.prior(A,dispersion)
    if (verbose) cat(sprintf("Beta prior for %s: a=%.3f, b=%.3f\n",paste(colnames(data)[A],collapse = ","),beta.prior[1],beta.prior[2]))

    re=plapply(1:nrow(data),function(i) {
      #for (i in 1:nrow(data)) { print (i);
      fit.A=FitKineticsGeneSnapshot(data=data,gene=i,columns=A,dispersion=dispersion[i],reference.columns=reference.columns,slot=slot,time.labeling=time.labeling,time.experiment=time.experiment,sample.f0.in.ss=sample.f0.in.ss,sample.level=sample.level,beta.prior=beta.prior,return.samples=TRUE,N=N,N.max=N.max,CI.size=CI.size,correct.labeling=correct.labeling)
      samp.a=fit.A$samples

      N=nrow(samp.a)
      if (N==0 || is.null(fit.A$samples))
        return(c(
        s=unname(fit.A$s),
        HL=unname(log(2)/fit.A$d),
        s.cred.lower=-Inf,
        s.cred.upper=-Inf,
        HL.cred.lower=-Inf,
        HL.cred.upper=Inf
      ))

      samp.a=samp.a[1:N,,drop=FALSE]

      sc=quantile(samp.a[,'s'],c(0.5-CI.size/2,0.5+CI.size/2))
      hc=quantile(log(2)/samp.a[,'d'],c(0.5-CI.size/2,0.5+CI.size/2))

      return(
        c(
          s=mean(samp.a[,'s']),
          HL=log(2)/mean(samp.a[,'d']),
          s.cred.lower=unname(sc[1]),
          s.cred.upper=unname(sc[2]),
          HL.cred.lower=unname(hc[1]),
          HL.cred.upper=unname(hc[2])
        )
      )
    },seed=seed)

    re.df=as.data.frame(t(simplify2array(re)))
    rownames(re.df)=Genes(data)
    data=AddAnalysis(data,name = if (is.null(name.prefix)) n else if (n=="") name.prefix else paste0(name.prefix,".",n),table = re.df)
  }
  data
}




#' Compute the posterior distributions of RNA synthesis and degradation for a particular gene
#'
#' @param data the grandR object
#' @param gene a gene name or symbol or index
#' @param columns samples or cell representing the same experimental condition (must refer to a unique labeling duration)
#' @param reference.columns a reference matrix usually generated by \link{FindReferences} to define reference samples for each sample (see details)
#' @param dispersion dispersion parameter for the given columns (if NULL, this is estimated from the data, takes a lot of time!)
#' @param slot the data slot to take f0 and totals from
#' @param time.labeling the column in the column annotation table denoting the labeling duration or the labeling duration itself
#' @param time.experiment the column in the column annotation table denoting the experimental time point (can be NULL, see details)
#' @param sample.f0.in.ss whether or not to sample f0 under steady state conditions
#' @param sample.level Define how the NTR is sampled from the hierarchical Bayesian model (must be 0,1, or 2; see details)
#' @param beta.prior The beta prior for the negative binomial used to sample counts, if NULL, a beta distribution is fit to all expression values and given dispersions
#' @param return.samples return the posterior samples of the parameters?
#' @param return.points return the point estimates per replicate as well?
#' @param N the posterior sample size
#' @param N.max the maximal number of posterior samples (necessary if old RNA > f0); if more are necessary, a warning is generated
#' @param CI.size A number between 0 and 1 representing the size of the credible interval
#' @param correct.labeling whether to correct labeling times
#'
#' @details The kinetic parameters s and d are computed using \link{TransformSnapshot}. For that, the sample either must be in steady state
#' (this is the case if defined in the reference.columns matrix), or if the levels of reference samples from a specific prior time point are known. This time point is
#' defined by \code{time.experiment} (i.e. the difference between the reference samples and samples themselves). If
#' \code{time.experiment} is NULL, then the labeling time of the samples is used (e.g. useful if labeling was started concomitantly with
#' the perturbation, and the reference samples are unperturbed samples).
#'
#' @details By default, the hierarchical Bayesian model is estimated. If sample.level = 0, the NTRs are sampled from a beta distribution
#' that approximates the mixture of betas from the replicate samples. If sample.level = 1, only the first level from the hierarchical model
#' is sampled (corresponding to the uncertainty of estimating the biological variability). If sample.level = 2, the first and second levels
#' are estimated (corresponding to the full hierarchical model).
#'
#' @details Columns can be given as a logical, integer or character vector representing a selection of the columns (samples or cells).
#' The expression is evaluated in an environment having the \code{\link{Coldata}}, i.e. you can use names of \code{\link{Coldata}} as variables to
#' conveniently build a logical vector (e.g., columns=Condition=="x").
#'
#' @return a list containing the posterior mean of s and s, its credible intervals and,
#' if return.samples=TRUE a data frame containing all posterior samples
#'
#' @export
#'
#' @concept snapshot
FitKineticsGeneSnapshot=function(data,gene,columns=NULL,
                             reference.columns=NULL,
                             dispersion=NULL,
                             slot=DefaultSlot(data),time.labeling=Design$dur.4sU,time.experiment=NULL,
                             sample.f0.in.ss=TRUE,
                             sample.level=2,
                             beta.prior=NULL,
                             return.samples=FALSE,
                             return.points=FALSE,
                             N=10000,N.max=N*10,
                             CI.size=0.95,
                             correct.labeling=FALSE
                             ) {

  if (!sample.level %in% c(0,1,2)) stop("Sample level must be 0,1 or 2!")

  columns=substitute(columns)
  columns=if (is.null(columns)) colnames(data) else eval(columns,Coldata(data),parent.frame())
  columns=Columns(data,columns)
  columns=Columns(data) %in% columns

  if (is.null(reference.columns)) reference.columns=columns
  if (is.matrix(reference.columns)) reference.columns=apply(reference.columns[,Columns(data,columns),drop=FALSE]==1,1,any)

    ntr=GetData(data,mode.slot="ntr",genes=gene,columns = columns)
    if (correct.labeling) {
      if (!is.numeric(time.labeling) && length(unique(ntr[[time.labeling]]))!=1) ntr[[time.labeling]] = median(ntr[[time.labeling]])
      if (!is.null(time.experiment) && length(unique(ntr[[time.experiment]]))!=1) ntr[[time.experiment]] = median(ntr[[time.experiment]])
    } else {
      if (!is.numeric(time.labeling) && length(unique(ntr[[time.labeling]]))!=1) stop("Labeling duration has to be unique!")
      if (!is.null(time.experiment) && length(unique(ntr[[time.experiment]]))!=1) stop("Experimental time has to be unique!")
    }

    total=GetData(data,mode.slot=slot,genes=gene,columns = columns)
    ss=GetData(data,mode.slot=slot,genes=gene,columns = reference.columns)
    if (correct.labeling) {
      if (!is.null(time.experiment) && length(unique(ss[[time.experiment]]))!=1) ss[[time.experiment]] = median(ss[[time.experiment]])
    } else {
      if (!is.null(time.experiment) && length(unique(ss[[time.experiment]]))!=1) stop("Steady state has to refer to a unique experimental time!")
    }

    if (is.null(dispersion)) dispersion=estimate.dispersion(as.matrix(GetTable(data,type="count",columns=columns,gene.info = F)))[ToIndex(data,gene)]

    if (is.null(beta.prior)) {
      ex=rowMeans(GetTable(data,type = slot, columns = columns))
      ex=1/dispersion/(1/dispersion+ex)
      E.ex=mean(ex)
      V.ex=var(ex)
      beta.prior=c(
        shape1=(E.ex*(1-E.ex)/V.ex-1)*E.ex,
        shape2=(E.ex*(1-E.ex)/V.ex-1)*(1-E.ex)
      )
    }

    use=total$Value>0 & !is.na(ntr$Value)

    emptyres=function() {
        re=list(
            s=NA,
            d=NA,
            s.CI=c(-Inf,Inf),
            d.CI=c(-Inf,Inf)
        )
        if (return.samples)
            re$samples=data.frame(s=numeric(0),d=numeric(0),ntr=numeric(0),total=numeric(0),t=numeric(0),t0=numeric(0),f0=numeric(0))
        if (return.points)
            re$points=data.frame(s=numeric(0),d=numeric(0),ntr=numeric(0),total=numeric(0),t=numeric(0),t0=numeric(0),f0=numeric(0))
        re$N=0
        return(re)
    }



    ntr=ntr[use,]
    total=total[use,]
    ss=ss[ss$Value>0,]
    if ((nrow(ntr)==0)&(nrow(total)==0)) return(emptyres())

    t=if (is.numeric(time.labeling)) time.labeling else unique(ntr[[time.labeling]])
    t0=if (is.null(time.experiment)) t else unique(ntr[[time.experiment]])-unique(ss[[time.experiment]])
    is.steady.state=any(reference.columns & columns)

    if (return.points)
        points=as.data.frame(if (is.steady.state) TransformSnapshot(ntr=ntr$Value,total=total$Value,t=t) else TransformSnapshot(ntr=ntr$Value,total=total$Value,t=t,t0=t0,f0=mean(ss$Value)))

    if (N<1) {
        if (is.steady.state) {
            param=TransformSnapshot(ntr=mean(ntr$Value),total=mean(total$Value),t=t)
        } else {
            if (t0<=0) stop("Experimental time is not properly defined (the steady state sample must be prior to each of A and B)!")
            f0=mean(ss$Value)
            param=TransformSnapshot(ntr=mean(ntr$Value),total=mean(total$Value),t=t,t0=t0,f0=f0)
        }
        re=emptyres()
        re$s=param['s']
        re$d=param['d']
        if (return.points)
            re$points=points
        return(re)
    }

    if (sum(use)<2 || sum(ss$Value>0)==0) {
      warning("Less than 2 samples. Skipping...")
      return(emptyres())
    }
    alpha=GetData(data,mode.slot="alpha",genes=gene,columns = columns)
    beta=GetData(data,mode.slot="beta",genes=gene,columns = columns)
    alpha=alpha[use,]
    beta=beta[use,]

    if (sample.level==2) {
      mod=hierarchical.beta.posterior(alpha$Value,beta$Value,compute.marginal.likelihood = FALSE,compute.grid = TRUE,res=50)$sample
    } else if (sample.level==1) {
      mod=hierarchical.beta.posterior(alpha$Value,beta$Value,compute.marginal.likelihood = FALSE,compute.grid = TRUE,res=50)$sample.mu
    } else {
        fit=beta.approximate.mixture(alpha$Value,beta$Value)
        mod=function(N) rbeta(N,fit$a,fit$b)
    }
    #sample.counts=function(N,x) 0.1+rnbinom(N,size=1/dispersion,mu=mean(x))
    #sample.counts=function(N,x) rgamma(N,shape=1/dispersion,scale=mean(x)*dispersion)
    # the conjugate prior for the negative binomial p is the beta distribution, and the posterior after observing x_1,...,x_n, is Beta(a+n/dispersion,b+sum x_i)
    # pl
    sample.counts=function(N,x) {
      pp=rbeta(N,beta.prior[1]+length(x)*1/dispersion,beta.prior[2]+sum(x))
      (1-pp)/(pp*dispersion)
    }
    sample.ss.with.f0=function(N) {
        f0=sample.counts(N,total$Value)
        ntr=mod(N)
        total=sample.counts(N,total$Value)
        TransformSnapshot(ntr=ntr,total=total,t=t,f0=f0,t0=t,full.return=return.samples)
    }
    sample.ss.without.f0=function(N) {
        ntr=mod(N)
        total=sample.counts(N,total$Value)
        TransformSnapshot(ntr=ntr,total=total,t=t,full.return=return.samples)
    }
    sample.non.ss=function(N) {
        if (t0<=0) stop("Experimental time is not properly defined (the steady state sample must be prior to each of A and B)!")
        f0=sample.counts(N,ss$Value)
        total=sample.counts(N,total$Value)
        ntr=mod(N)
        TransformSnapshot(ntr=ntr,total=total,t=t,t0=t0,f0=f0,full.return=return.samples)
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

    sc=quantile(samp[,'s'],c(0.5-CI.size/2,0.5+CI.size/2))
    dc=quantile(samp[,'d'],c(0.5-CI.size/2,0.5+CI.size/2))

    re=list(
        s=mean(samp[,'s']),
        d=mean(samp[,'d']),
        s.CI=sc,
        d.CI=dc
    )
    if (return.samples)
        re$samples=as.data.frame(samp)
    if (return.points)
        re$points=points
    re$N=N
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
#' @param full.return also return the provided parameters
#'
#' @return a named vector for s and d
#'
#' @details t0 must be given as the total time in between the measurement of f0 and the given ntr and total values!
#'
#' @export
#'
#' @concept snapshot
TransformSnapshot=function(ntr,total,t,t0=NULL,f0=NULL,full.return=FALSE) {
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
        #eq=function(d) f0*exp(-t0*d) + new *(exp(-t*d)-exp(-t0*d))/(1-exp(-t*d))-old
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
    if (full.return) {
      if (is.null(t0)) t0=t
      if (is.null(f0)) f0=unname(s/d)
      if (length(s)>1) cbind(s=unname(s),d=unname(d),ntr=ntr,total=total,t=t,t0=t0,f0=f0) else c(s=unname(s),d=unname(d),ntr=ntr,total=total,t=t,t0=t0,f0=f0)
    } else {
      if (length(s)>1) cbind(s=unname(s),d=unname(d)) else c(s=unname(s),d=unname(d))
    }
}


#' Plot progressive labeling timecourses
#'
#' Plot the abundance of new and old RNA and the fitted model over time for a single gene.
#'
#' @param data a grandR object
#' @param gene the gene to be plotted
#' @param slot the data slot of the observed abundances
#' @param time the labeling duration column in the column annotation table
#' @param type how to fit the model (see link{FitKinetics})
#' @param exact.tics use axis labels directly corresponding to the available labeling durations?
#' @param show.CI show confidence intervals; one of TRUE/FALSE (default: FALSE)
#' @param return.tables also return the tables used for plotting
#' @param ... given to the fitting procedures
#'
#' @details For each \code{\link{Condition}} there will be one panel containing the values and the corresponding model fit.
#'
#' @return either a ggplot object, or a list containing all tables used for plotting and the ggplot object.
#'
#' @seealso \link{FitKineticsGeneNtr}, \link{FitKineticsGeneLeastSquares}, \link{FitKineticsGeneLogSpaceLinear}
#'
#' @export
#' @concept geneplot
PlotGeneProgressiveTimecourse=function(data,gene,slot=DefaultSlot(data),time=Design$dur.4sU, type=c("nlls","ntr","lm"),
                                       exact.tics=TRUE,show.CI=FALSE,return.tables=FALSE,...) {
  # R CMD check guard for non-standard evaluation
  Value <- Type <- lower <- upper <- NULL

  if (length(ToIndex(data,gene))==0) return(NULL)

    fit=switch(tolower(type[1]),
               ntr=FitKineticsGeneNtr(data,gene,slot=slot,time=time,...),
               nlls=FitKineticsGeneLeastSquares(data,gene,slot=slot,time=time,...),
               chase=FitKineticsGeneLeastSquares(data,gene,slot=slot,time=time,...,chase=TRUE),
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
            fit=if(is.null(Condition(data))) fit[[gene]] else fit[[as.character(df$Condition)[i]]]
            tt=df$time[i]
            f.old.nonequi(tt,fit$f0,fit$Synthesis,fit$Degradation)+f.new(tt,fit$Synthesis,fit$Degradation)
        }))/df$Value[df$Type=="Total"]
        df$Value=df$Value*fac
        if (show.CI) {
            df$lower=df$lower*fac
            df$upper=df$upper*fac
        }
    }

    if (tolower(type[1])=="chase") df=df[!df$no4sU,]

    df$Condition=if ("Condition" %in% names(df)) df$Condition else gene
    tt=seq(0,max(df$time),length.out=100)
    df.median=plyr::ddply(df,c("Condition","Type",time,"time"),function(s) data.frame(Value=median(s$Value)))

    if (tolower(type[1])=="chase") {
      fitted=plyr::ldply(fit,function(f) data.frame(time=c(tt,tt),Value=c(f$Synthesis/f$Degradation-f.old.equi(tt,f$Synthesis,f$Degradation),f.old.equi(tt,f$Synthesis,f$Degradation)),Type=rep(c("Old","New"),each=length(tt))),.id="Condition")
    } else {
      fitted=plyr::ldply(fit,function(f) data.frame(time=c(tt,tt),Value=c(f.old.nonequi(tt,f$f0,f$Synthesis,f$Degradation),f.new(tt,f$Synthesis,f$Degradation)),Type=rep(c("Old","New"),each=length(tt))),.id="Condition")
    }
    g=ggplot(df,aes(time,Value,color=Type))+cowplot::theme_cowplot()
    if (show.CI) g=g+geom_errorbar(mapping=aes(ymin=lower,ymax=upper),width=max(df$time)*0.02)
    g=g+geom_point(size=2)+
        geom_line(data=df.median[df.median$Type=="Total",],size=1)+
        scale_color_manual("RNA",values=c(Total="gray",New="#e34a33",Old="#2b8cbe"))+
        ylab("Expression")+
        geom_line(data=fitted,aes(ymin=NULL,ymax=NULL),linetype=2,size=1)

    if (exact.tics) {
      timeorig=paste0(time,".original")
      if (timeorig %in% names(df)) {
        brdf=unique(df[,c(time,timeorig)])
        brdf=brdf[order(brdf[[time]]),]
        brdf[[timeorig]]=gsub("_",".",brdf[[timeorig]])
        g=g+scale_x_continuous(NULL,labels = brdf[[timeorig]],breaks=brdf[[time]])+RotatateAxisLabels(45)
      } else {
        breaks=sort(unique(df$time))
        g=g+scale_x_continuous("4sU labeling [h]",labels = scales::number_format(accuracy = max(0.01,my.precision(breaks))),breaks=breaks)
      }
    } else {
      breaks=scales::breaks_extended(5)(df[[time]])
      g=g+scale_x_continuous("4sU labeling [h]",labels = scales::number_format(accuracy = max(0.01,my.precision(breaks))),breaks=breaks)
    }

    if (!is.null(Coldata(data)$Condition)) g=g+facet_wrap(~Condition,nrow=1)
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
#' @param s the synthesis rate (see details)
#' @param d the degradation rate (see details)
#' @param hl the RNA half-life
#' @param f0 the abundance at time t=0
#' @param min.time the start time to simulate
#' @param max.time the end time to simulate
#' @param N how many time points from min.time to max.time to simuate
#' @param name add a Name column to the resulting data frame
#' @param out which values to put into the data frame
#'
#' @details Both rates can be either (i) a single number (constant rate), (ii) a data frame with names "offset",
#' "factor" and "exponent" (for linear functions, see \link{ComputeNonConstantParam}) or (iii) a unary function time->rate. Functions
#'
#' @return a data frame containing the simulated values
#' @export
#'
#' @seealso \link{PlotSimulation} for plotting the simulation
#'
#' @examples
#' head(SimulateKinetics(hl=2))   # simulate steady state kinetics for an RNA with half-life 2h
#'
#' @concept kinetics
SimulateKinetics=function(s=100*d,d=log(2)/hl,hl=2,f0=NULL,min.time=-1,max.time=10,N = 1000,name=NULL,out=c("Old","New","Total","NTR")) {
    times=seq(min.time,max.time,length.out=N)
    if (is.numeric(s)) s = ComputeNonConstantParam(start=s)
    if (is.numeric(d)) d = ComputeNonConstantParam(start=d)
    if (is.null(f0)) f0=s$offset/d$offset
    #ode.new=function(t,s,d) ifelse(t<0,0,s/d*(1-exp(-t*d)))
    #ode.old=function(t,f0,s,d) ifelse(t<0,f0,f0*exp(-t*d))
    #old=ode.old(times,f0,s,d)
    #new=ode.new(times,s,d)
    #old=c(rep(f0,sum(times<0)),f.nonconst.linear(t = times[times>=0],f0 = f0,so = 0,sf=0,se=1,do = d$offset,df=d$factor,de=d$exponent))
    #new=c(rep(0,sum(times<0)),f.nonconst.linear(t = times[times>=0],f0 = 0,so = s$offset,sf=s$factor,se=s$exponent,do = d$offset,df=d$factor,de=d$exponent))
    old=c(rep(f0,sum(times<0)),f.nonconst(t = times[times>=0],f0 = f0,s=0,d=d))
    new=c(rep(0,sum(times<0)),f.nonconst(t = times[times>=0],f0 = 0,s=s,d=d))

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
#' @concept kinetics
PlotSimulation=function(sim.df,ntr=TRUE,old=TRUE,new=TRUE,total=TRUE) {
  # R CMD check guard for non-standard evaluation
  Time <- Value <- Type <- NULL

  if (!ntr) sim.df=sim.df[sim.df$Type!="NTR",]
    if (!old) sim.df=sim.df[sim.df$Type!="Old",]
    if (!new) sim.df=sim.df[sim.df$Type!="New",]
    if (!total) sim.df=sim.df[sim.df$Type!="Total",]
    sim.df$Type=droplevels(sim.df$Type)
    ggplot(sim.df,aes(Time,Value,color=Type))+
      cowplot::theme_cowplot()+
      geom_line(size=1)+
        scale_color_manual(NULL,values=c(Old="#54668d",New="#953f36",Total="#373737",NTR="#e4c534")[levels(sim.df$Type)])+
        facet_wrap(~ifelse(Type=="NTR","NTR","Timecourse"),scales="free_y",ncol=1)+
        ylab(NULL)+
        scale_x_continuous(breaks=scales::pretty_breaks())+
        theme(
	  strip.background = element_blank(),
	  strip.text.x = element_blank()
	)
}
