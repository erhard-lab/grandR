

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

m=DropAnalysis(m)
m=subset(m,columns=Coldata(m)$Replicate %in% c("A","B","C"))

m=PairwiseDESeq2(m,"total",GetContrasts(m,contrast = c("Condition","control")),mode="total",verbose=T)
m=LFC(m,"total",GetContrasts(m,contrast = c("Condition","control")),mode="total")
m=PairwiseDESeq2(m,"new",GetContrasts(m,contrast = c("Condition","control")),mode="new",verbose=T)
m=LFC(m,"new",GetContrasts(m,contrast = c("Condition","control")),mode="new")
m=PairwiseDESeq2(m,"old",GetContrasts(m,contrast = c("Condition","control")),mode="old",verbose=T)
m=LFC(m,"old",GetContrasts(m,contrast = c("Condition","control")),mode="old")

m=Normalize(m)

SetParallel()
m=EstimateRegulation(m,"Regulation",contrasts = GetContrasts(m,contrast = c("Condition","control")),steady.state=FindReferences(m,Condition=="control",group=NULL),verbose=TRUE)
m=EstimateRegulation(m,"Regulation1",ROPE.max.log2FC = 1,contrasts = GetContrasts(m,contrast = c("Condition","control")),steady.state=FindReferences(m,Condition=="control",group=NULL),verbose=TRUE)
m=EstimateRegulation(m,"Regulationsampf0",sample.f0.in.ss = TRUE,contrasts = GetContrasts(m,contrast = c("Condition","control")),steady.state=FindReferences(m,Condition=="control",group=NULL),verbose=TRUE)
#m=PairwiseNtrTest(m,"NTR",contrasts = GetContrasts(m,contrast = c("Condition","control")),verbose = T)

saveRDS(m,file="perturbed_3rep.rds")
m=readRDS("perturbed_3rep.rds")

pp=function(min) ggroc(list(
  ROPE=roc(abs(log2(df$true_s/df$true_s.perturbed.s))>min,df$`Regulation.s vs control.s.ROPE`),
  ROPE1=roc(abs(log2(df$true_s/df$true_s.perturbed.s))>min,df$`Regulation1.s vs control.s.ROPE`),
  ROPEAs025=roc(abs(log2(df$true_s/df$true_s.perturbed.s))>min,df$`RegulationAs025.s vs control.s.ROPE`),
  ROPEAs=roc(abs(log2(df$true_s/df$true_s.perturbed.s))>min,df$`RegulationAs.s vs control.s.ROPE`),
  new=roc(abs(log2(df$true_s/df$true_s.perturbed.s))>min,df$`new.s vs control.Q`),
  total=roc(abs(log2(df$true_s/df$true_s.perturbed.s))>min,df$`total.s vs control.Q`),
  old=roc(abs(log2(df$true_s/df$true_s.perturbed.s))>min,df$`old.s vs control.Q`)
))

pp=function(min) ggroc(list(
  ROPE=roc(abs(log2(df$true_s/df$true_s.perturbed.s))>min,df$`Regulation.s vs control.s.ROPE`),
  new=roc(abs(log2(df$true_s/df$true_s.perturbed.s))>min,df$`new.s vs control.Q`),
  total=roc(abs(log2(df$true_s/df$true_s.perturbed.s))>min,df$`total.s vs control.Q`),
  old=roc(abs(log2(df$true_s/df$true_s.perturbed.s))>min,df$`old.s vs control.Q`)
))


pp=function(min) ggroc(list(
  ROPE=roc(abs(log2(df$true_d/df$true_d.perturbed.d))>min,df$`Regulation.HL vs control.HL.ROPE`),
  ROPE1=roc(abs(log2(df$true_d/df$true_d.perturbed.d))>min,df$`Regulation1.HL vs control.HL.ROPE`),
  ROPEAs025=roc(abs(log2(df$true_d/df$true_d.perturbed.d))>min,df$`RegulationAs025.HL vs control.HL.ROPE`),
  ROPEAs=roc(abs(log2(df$true_d/df$true_d.perturbed.d))>min,df$`RegulationAs.HL vs control.HL.ROPE`),
  new=roc(abs(log2(df$true_d/df$true_d.perturbed.d))>min,df$`new.HL vs control.Q`),
  total=roc(abs(log2(df$true_d/df$true_d.perturbed.d))>min,df$`total.HL vs control.Q`),
  old=roc(abs(log2(df$true_d/df$true_d.perturbed.d))>min,df$`old.HL vs control.Q`)
))




m=readRDS("perturbed.rds")

power=function(cond) {
  c(
    total=sum(df[[sprintf("total.%s vs control.Q",cond)]]<0.05,na.rm=TRUE),
    new=sum(df[[sprintf("new.%s vs control.Q",cond)]]<0.05,na.rm=TRUE),
    old=sum(df[[sprintf("old.%s vs control.Q",cond)]]<0.05,na.rm=TRUE),
    ntr=sum(df[[sprintf("NTR.%s vs control.Q",cond)]]<0.05),
    s.regu=sum(df[[sprintf("Regulation.%s vs control.s.regulated",cond)]]),
    HL.regu=sum(df[[sprintf("Regulation.%s vs control.HL.regulated",cond)]]),
    combined=sum(p.adjust(pchisq(-2*(log(df[[sprintf("NTR.%s vs control.P",cond)]])+log(df[[sprintf("total.%s vs control.P",cond)]])),df=4,lower.tail = FALSE),method="BH")<0.05)
)
}


contrasts=GetContrasts(m,contrast = c("Condition","control"))
n="s vs control"
A=contrasts[[n]]==1
B=contrasts[[n]]==-1



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







