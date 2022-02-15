
library(grandR)
sars <- ReadGRAND("https://zenodo.org/record/5834034/files/sars.tsv.gz",design=c(Design$Condition,Design$dur.4sU,Design$Replicate))
sars=subset(sars,columns=Condition(sars)=="Mock")
sars=Normalize(sars)
SetParallel()
sars=FitKinetics(sars,type="ntr")
seed=1337
sim=SimulateTimeCourse("sim",GeneInfo(sars),
                       s = GetAnalysisTable(sars)$Synthesis,
                       d = log(2)/GetAnalysisTable(sars)$`Half-life`,
                       f0 = rowMeans(GetTable(sars,type="norm")) *2^rnorm(nrow(sars),0,1),
                       dispersion=estimate.dispersion(GetTable(sars,type="count")),
                       s.variation = 1.05,
                       d.variation = 1.05,
                       timepoints = c(0,0,0,1,1,1,2,2,2,4,4,4,8,8,8),
                       beta.approx = TRUE,
                       num.reads=1E7,
                       seed=seed)
saveRDS(sim,file="timecourse.rds")


