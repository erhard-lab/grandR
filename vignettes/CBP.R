
d=ReadGRAND("/mnt/hilbert/projects/erhard/grandslam/CBP_Diesch/CBP/grandslam_t15/CBP",design=c("Condition","Replicate"))
d=FilterGenes(d)
d=Normalize(d)
Coldata(d,"Experimental.time")=c(0,0,0,2,2,2)
Coldata(d,Design$dur.4sU)=1


contrasts=GetContrasts(d,contrast=c("Condition","DMSO"))
steady.state=FindReferences(d,Condition=="DMSO",group=NULL)

d=PairwiseDESeq2(d,"total",contrasts=contrasts,mode="total",logFC = TRUE)
d=PairwiseDESeq2(d,"new",contrasts=contrasts,mode="new",normalization="total",logFC = TRUE)
d=PairwiseDESeq2(d,"old",contrasts=contrasts,mode="old",normalization="total",logFC = TRUE)
d=LFC(d,"old",contrasts=contrasts,mode="old",normalization="total")
d=LFC(d,"new",contrasts=contrasts,mode="new",normalization="total")
d=LFC(d,"total",contrasts=contrasts,mode="total")

VulcanoPlot(d,"old.C646 vs DMSO")
VulcanoPlot(d,"new.C646 vs DMSO")
MAPlot(d,"new.C646 vs DMSO")
MAPlot(d,"total.C646 vs DMSO")

SetParallel()
d=EstimateRegulation(d,"Regulation",contrasts=contrasts,steady.state=steady.state,verbose=T,time.experiment = "Experimental.time",sample.f0.in.ss=FALSE)
d=EstimateRegulation(d,"RegulationSamp",contrasts=contrasts,steady.state=steady.state,verbose=T,time.experiment = "Experimental.time",sample.f0.in.ss=TRUE)

saveRDS(d,"CBP.rds")
d=readRDS("CBP.rds")




dispersion = estimate.dispersion(GetTable(d,type="count",columns = steady.state[,1]))
plotr=function(d,gene,time) {
  re=EstimateGeneRegulation(d,gene,A=contrasts[[1]]==1,B=contrasts[[1]]==-1,dispersion.A = dispersion[ToIndex(d,gene)],dispersion.B = dispersion[ToIndex(d,gene)],steady.state = steady.state, return.samples = TRUE,time.experiment = "Experimental.time",time.labeling = time,sample.f0.in.ss = TRUE)
  PlotScatter(re$samples)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)
}

ServeData(d,df.set=d,plot.single = list(DPlot(PlotGeneGroupsBars),DPlot(plotr,time=1),DPlot(plotr,time=0.5)))
