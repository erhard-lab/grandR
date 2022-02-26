

sars <- ReadGRAND("https://zenodo.org/record/5834034/files/sars.tsv.gz",design=c(Design$Condition,Design$dur.4sU,Design$Replicate))
sars=subset(sars,columns=Condition(sars)=="Mock")
sars=FilterGenes(sars)
sars=Normalize(sars)
SetParallel()
sars=FitKinetics(sars,type="ntr")
seed=1337
f0=GetAnalysisTable(sars)$Synthesis/log(2)*GetAnalysisTable(sars)$`Half-life`
control=SimulateTimeCourse("control",GeneInfo(sars),s = GetAnalysisTable(sars)$Synthesis,d = log(2)/GetAnalysisTable(sars)$`Half-life`,dispersion=estimate.dispersion(GetTable(sars,type="count")),f0 = f0,s.variation = 1.05,d.variation = 1.05,timepoints = c(2,2,2,2,2),beta.approx = TRUE,num.reads=5E7,seed=seed)
unperturbed=SimulateTimeCourse("unperturbed",GeneInfo(sars),s = GetAnalysisTable(sars)$Synthesis,d = log(2)/GetAnalysisTable(sars)$`Half-life`,dispersion=estimate.dispersion(GetTable(sars,type="count")),f0 = f0,s.variation = 1.05,d.variation = 1.05,timepoints = c(2,2,2,2,2),beta.approx = TRUE,num.reads=5E7,seed=seed*13)
perturbed.d=SimulateTimeCourse("HL",GeneInfo(sars),s = GetAnalysisTable(sars)$Synthesis,d = log(2)/GetAnalysisTable(sars)$`Half-life` *2^rnorm(nrow(sars),0,1) ,dispersion=estimate.dispersion(GetTable(sars,type="count")),f0 = f0,s.variation = 1.05,d.variation = 1.05,timepoints = c(2,2,2,2,2),beta.approx = TRUE,num.reads=5E7,seed=seed*13*13)
perturbed.s=SimulateTimeCourse("s",GeneInfo(sars),s = GetAnalysisTable(sars)$Synthesis *2^rnorm(nrow(sars),0,1) ,d = log(2)/GetAnalysisTable(sars)$`Half-life`,dispersion=estimate.dispersion(GetTable(sars,type="count")),f0 = f0,s.variation = 1.05,d.variation = 1.05,timepoints = c(2,2,2,2,2),beta.approx = TRUE,num.reads=5E7,seed=seed*13*13*13)

m=merge(control,unperturbed,perturbed.d,perturbed.s)


m=PairwiseDESeq2(m,"total",GetContrasts(m,contrast = c("Condition","control")),mode="total",verbose=T)
m=LFC(m,"total",GetContrasts(m,contrast = c("Condition","control")),mode="total")
m=PairwiseDESeq2(m,"new",GetContrasts(m,contrast = c("Condition","control")),mode="new",verbose=T)
m=LFC(m,"new",GetContrasts(m,contrast = c("Condition","control")),mode="new")
m=PairwiseDESeq2(m,"old",GetContrasts(m,contrast = c("Condition","control")),mode="old",verbose=T)
m=LFC(m,"old",GetContrasts(m,contrast = c("Condition","control")),mode="old")

m=Normalize(m)

SetParallel()
m=EstimateRegulation(m,"Regulation",contrasts = GetContrasts(m,contrast = c("Condition","control")),steady.state=FindReferences(m,Condition=="control",group=NULL),verbose=TRUE)
m=PairwiseNtrTest(m,"NTR",contrasts = GetContrasts(m,contrast = c("Condition","control")),verbose = T)


saveRDS(m,file="perturbed.rds")


m=readRDS("perturbed.rds")

power=function(cond) {
  c(
    total=sum(df[[sprintf("total.%s vs control.Q",cond)]]<0.05,na.rm=TRUE),
    new=sum(df[[sprintf("new.%s vs control.Q",cond)]]<0.05,na.rm=TRUE),
    old=sum(df[[sprintf("old.%s vs control.Q",cond)]]<0.05,na.rm=TRUE),
    ntr=sum(df[[sprintf("NTR.%s vs control.Q",cond)]]<0.05),
    combined=sum(p.adjust(pchisq(-2*(log(df[[sprintf("NTR.%s vs control.P",cond)]])+log(df[[sprintf("total.%s vs control.P",cond)]])),df=4,lower.tail = FALSE),method="BH")<0.05)
)
}



VulcanoPlot(m,analysis="total.s vs control",lfc.cutoff = 0.25) | VulcanoPlot(m,analysis="total.HL vs control",lfc.cutoff = 0.25)
VulcanoPlot(m,analysis="new.s vs control",lfc.cutoff = 0.25) | VulcanoPlot(m,analysis="new.HL vs control",lfc.cutoff = 0.25)
VulcanoPlot(m,analysis="old.s vs control",lfc.cutoff = 0.25) | VulcanoPlot(m,analysis="old.HL vs control",lfc.cutoff = 0.25)



df=cbind(GetAnalysisTable(m,"old|new|total",column="LFC|Q"),GetAnalysisTable(m,"Regulation",gene.info = FALSE),GetAnalysisTable(m,"NTR",gene.info = FALSE))
PlotScatter(df,x=-log2(true_d.perturbed.d/true_d),y=`old.HL vs control.LFC`,remove.outlier=F)+geom_abline()
PlotScatter(df,x=-log2(true_d.perturbed.d/true_d),y=`Regulation.HL vs control.HL.log2FC`,remove.outlier=F)+geom_abline()
PlotScatter(df,x=-log2(true_d.perturbed.d/true_d),y=-log10(`NTR.HL vs control.Q`),remove.outlier=F,ylim=c(-2,10))
PlotScatter(df,x=-log2(true_s.perturbed.s/true_s),y=-log10(`NTR.s vs control.Q`),remove.outlier=F,ylim=c(-2,10))


PlotScatter(df,x=rank(log(2)/true_d),y=`new.HL vs control.LFC`,ylim=c(-1,1))
PlotScatter(df,x=rank(log(2)/true_d),y=`Regulation.HL vs control.s.log2FC`,ylim=c(-1,1))







