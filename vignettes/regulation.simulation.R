

sars <- ReadGRAND("https://zenodo.org/record/5834034/files/sars.tsv.gz",design=c(Design$Condition,Design$dur.4sU,Design$Replicate))
sars=subset(sars,columns=Condition(sars)=="Mock")
sars=FilterGenes(sars)
sars=Normalize(sars)
SetParallel()
sars=FitKinetics(sars,type="ntr")
seed=1337
control=SimulateTimeCourse("control",GeneInfo(sars),s = GetAnalysisTable(sars)$Synthesis,d = log(2)/GetAnalysisTable(sars)$`Half-life`,dispersion=estimate.dispersion(GetTable(sars,type="count")),f0 = GetTable(sars,type="norm")$Mock.no4sU.A,s.variation=1.05,d.variation=1.05,timepoints = c(2,2,2,2,2),beta.approx = FALSE,num.reads=5E7,seed=seed)
perturbed.hl=SimulateTimeCourse("HL",GeneInfo(sars),s = GetAnalysisTable(sars)$Synthesis,d = log(2)/GetAnalysisTable(sars)$`Half-life` *2^rnorm(nrow(sars),0,1) ,dispersion=estimate.dispersion(GetTable(sars,type="count")),f0 = GetTable(sars,type="norm")$Mock.no4sU.A,s.variation=1.05,d.variation=1.05,timepoints = c(2,2,2,2,2),beta.approx = TRUE,num.reads=5E7,seed=seed)
perturbed.s=SimulateTimeCourse("s",GeneInfo(sars),s = GetAnalysisTable(sars)$Synthesis *2^rnorm(nrow(sars),0,1) ,d = log(2)/GetAnalysisTable(sars)$`Half-life`,dispersion=estimate.dispersion(GetTable(sars,type="count")),f0 = GetTable(sars,type="norm")$Mock.no4sU.A,s.variation=1.05,d.variation=1.05,timepoints = c(2,2,2,2,2),beta.approx = TRUE,num.reads=5E7,seed=seed)

m=merge(control,perturbed.hl,perturbed.s)
saveRDS(m,file="perturbed.rds")

m=readRDS("perturbed.rds")

m=PairwiseDESeq2(m,"total",GetContrasts(m,contrast = c("Condition","control")),mode="total",verbose=T)
m=LFC(m,"total",GetContrasts(m,contrast = c("Condition","control")),mode="total")
m=PairwiseDESeq2(m,"new",GetContrasts(m,contrast = c("Condition","control")),mode="new",verbose=T)
m=LFC(m,"new",GetContrasts(m,contrast = c("Condition","control")),mode="new")
m=PairwiseDESeq2(m,"old",GetContrasts(m,contrast = c("Condition","control")),mode="old",verbose=T)
m=LFC(m,"old",GetContrasts(m,contrast = c("Condition","control")),mode="old")

m=Normalize(m)

VulcanoPlot(m,name="total.s vs control",lfc.cutoff = 0.25) | VulcanoPlot(m,name="total.HL vs control",lfc.cutoff = 0.25)
VulcanoPlot(m,name="new.s vs control",lfc.cutoff = 0.25) | VulcanoPlot(m,name="new.HL vs control",lfc.cutoff = 0.25)
VulcanoPlot(m,name="old.s vs control",lfc.cutoff = 0.25) | VulcanoPlot(m,name="old.HL vs control",lfc.cutoff = 0.25)

m=PairwiseRegulation(m,"regu",contrasts = GetContrasts(m,contrast = c("Condition","control")),steady.state.columns = Coldata(m)$Condition=="control",verbose=TRUE)
PlotScatter(GetAnalysisTable(m,pattern="regu.HL vs control", regex=FALSE),"LFC.s","Prob.s",remove.outlier = F)
PlotScatter(GetAnalysisTable(m,pattern="regu.s vs control", regex=FALSE),"LFC.s","Prob.s",remove.outlier = F)
PlotScatter(GetAnalysisTable(m,pattern="regu.s vs control", regex=FALSE),"LFC.HL","Prob.HL",remove.outlier = F)
PlotScatter(GetAnalysisTable(m,pattern="regu.HL vs control", regex=FALSE),"LFC.HL","Prob.HL",remove.outlier = F)

df=cbind(GetAnalysisTable(m,"regu|simple"),GetAnalysisTable(m,"old|new",column="LFC",gene.info = F))
PlotScatter(df,x=log2(true_HL.perturbed.hl/true_HL),y=`old.HL vs control.LFC`,remove=F)+geom_abline()
PlotScatter(df,x=log2(true_HL.perturbed.hl/true_HL),y=`regu.HL vs control.LFC.HL`,remove=F)+geom_abline()
PlotScatter(df,x=true_HL,y=`new.HL vs control.LFC`,xlim=c(0,20),ylim=c(-1,1))
PlotScatter(df,x=true_HL,y=`regu.HL vs control.LFC.s`,xlim=c(0,20),ylim=c(-1,1))
