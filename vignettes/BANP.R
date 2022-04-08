
d=ReadGRAND("/mnt/hilbert/projects/erhard/grandslam/BANP_Schubeler/BANP/grandslam_t15/BANP",design=c("Condition","Experimental.time","Genotype",Design$dur.4sU,Design$has.4sU,"Replicate"))
d=FilterGenes(d)
d=Normalize(d)


contrasts=GetContrasts(d,contrast=c("Experimental.time.original","0h"),group = "Genotype")
reference.samples=FindReferences(d,reference.function=function(s) Experimental.time==max(c(0,Experimental.time[Experimental.time<s$Experimental.time])),group="Genotype")

d=PairwiseDESeq2(d,"total",contrasts=contrasts,mode="total")
d=PairwiseDESeq2(d,"new",contrasts=contrasts,mode="new",normalization="total")
d=PairwiseDESeq2(d,"old",contrasts=contrasts,mode="old",normalization="total")
d=LFC(d,"old",contrasts=contrasts,mode="old",normalization="total")
d=LFC(d,"new",contrasts=contrasts,mode="new",normalization="total")
d=LFC(d,"total",contrasts=contrasts,mode="total")


SetParallel()
d=EstimateRegulation(d,"Regulation",contrasts=contrasts[,1,drop=FALSE],reference.samples=reference.samples,verbose=T,time.experiment = "Experimental.time")

saveRDS(d,"BANP.rds")
d=readRDS("BANP.rds")




dispersion.0h = estimate.dispersion(GetTable(d,type="count",columns = Genotype=='dTag' & Experimental.time==0))
dispersion=lapply(c(1,2,4,6,20),function(time) estimate.dispersion(GetTable(d,type="count",columns = Genotype=='dTag' & Experimental.time==time)))
names(dispersion)=colnames(contrasts)

plotr=function(d,gene,i.contr) {
  re=EstimateGeneRegulation(d,gene,A=contrasts[[i.contr]]==1,B=contrasts[[i.contr]]==-1,dispersion.A = dispersion[[i.contr]][ToIndex(d,gene)],dispersion.B = dispersion.0h[ToIndex(d,gene)],reference.samples = reference.samples, return.samples = TRUE,time.experiment = "Experimental.time",sample.f0.in.ss = TRUE)
  PlotScatter(re$samples)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+ggtitle(i.contr)
}

ServeData(d,plot.single = list(DPlot(plotr,i.contr=colnames(contrasts)[1]),
                                        DPlot(plotr,i.contr=colnames(contrasts)[2]),
                                        DPlot(plotr,i.contr=colnames(contrasts)[3]),
                                        DPlot(plotr,i.contr=colnames(contrasts)[4]),
                                        DPlot(plotr,i.contr=colnames(contrasts)[5]),
                                        DPlot(PlotGeneGroupsBars,show.CI=TRUE)))
