
sars <- ReadGRAND("https://zenodo.org/record/5834034/files/sars.tsv.gz",design=c(Design$Condition,Design$dur.4sU,Design$Replicate))
sars=subset(sars,columns=Condition(sars)=="Mock")
sars=FilterGenes(sars)
sars=Normalize(sars)
SetParallel()
sars=FitKinetics(sars,type="ntr")
seed=1337
noss=SimulateTimeCourse("noss",GeneInfo(sars),
                       s = GetAnalysisTable(sars)$Synthesis,
                       d = log(2)/GetAnalysisTable(sars)$`Half-life`,
                       f0 = rowMeans(GetTable(sars,type="norm")) *2^rnorm(nrow(sars),0,1),
                       dispersion=estimate.dispersion(GetTable(sars,type="count")),
                       s.variation = 1.05,
                       d.variation = 1.05,
                       timepoints = c(0,0,0,1,1,1,2,2,2,4,4,4,8,8,8),
                       beta.approx = TRUE,
                       conversion.reads = TRUE,
                       num.reads=1E7,
                       seed=seed)
ss=SimulateTimeCourse("ss",GeneInfo(sars),
                       s = GetAnalysisTable(sars)$Synthesis,
                       d = log(2)/GetAnalysisTable(sars)$`Half-life`,
                       dispersion=estimate.dispersion(GetTable(sars,type="count")),
                       s.variation = 1.05,
                       d.variation = 1.05,
                       timepoints = c(0,0,0,1,1,1,2,2,2,4,4,4,8,8,8),
                       beta.approx = TRUE,
                       conversion.reads = TRUE,
                       num.reads=1E7,
                       seed=seed)
m=merge(ss,noss)

extra=function(s) c(`Synthesis.lower`=unname(s$conf.lower["Synthesis"]),`Synthesis.upper`=unname(s$conf.upper["Synthesis"]),`Half-life.lower`=unname(s$conf.lower["Half-life"]),`Half-life.upper`=unname(s$conf.upper["Half-life"]))
m=FitKineticsPulseR(m,name="pulseR")
m=FitKinetics(m,name="nlls",type="nlls",steady.state=c(ss=TRUE,no.ss=FALSE),return.extra = extra)
m=FitKinetics(m,name="ntr",type="ntr",return.extra = extra,exact.ci=FALSE)
m=FitKinetics(m,name="entr",type="ntr",return.extra = extra,exact.ci=TRUE)
m=FitKinetics(m,name="lm",type="lm",return.extra = extra)

saveRDS(m,file="timecourse.rds")

m=readRDS("timecourse.rds")
