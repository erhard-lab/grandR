d=ReadGRAND("/mnt/hilbert/projects/erhard/grandslam/muhar_science_2018/MYC/grandslam_t15/MYC",design=c("Condition","Treatment","Replicate"))
genes=FilterGenes(d,return.genes = TRUE)

d=ReadGRAND("/mnt/hilbert/projects/erhard/grandslam/muhar_science_2018/MYC/grandslam_t15_bias/MYC.high.tsv.gz",design=c("Condition","Treatment","Replicate"))
d=FilterGenes(d,use=genes)
Coldata(d,"Treatment2")=factor(gsub("IAA","Auxin",Coldata(d)$Treatment),levels=c("DMSO","Auxin"))
Coldata(d,Design$dur.4sU)=1
d=Normalize(d)
g.high=PlotGeneGroupsBars(d,"TWNK")
high=d

d=ReadGRAND("/mnt/hilbert/projects/erhard/grandslam/muhar_science_2018/MYC/grandslam_t15_bias/MYC.low.tsv.gz",design=c("Condition","Treatment","Replicate"))
d=FilterGenes(d,use=genes)
Coldata(d,"Treatment2")=factor(gsub("IAA","Auxin",Coldata(d)$Treatment),levels=c("DMSO","Auxin"))
Coldata(d,Design$dur.4sU)=1
d=Normalize(d)
g.low=PlotGeneGroupsBars(d,"TWNK")
low=d

d=ReadGRAND("/mnt/hilbert/projects/erhard/grandslam/muhar_science_2018/MYC/grandslam_t15_bias/MYC.tsv.gz",design=c("Condition","Treatment","Replicate"))
d=FilterGenes(d,minval=300)
Coldata(d,"Treatment2")=factor(gsub("IAA","Auxin",Coldata(d)$Treatment),levels=c("DMSO","Auxin"))
Coldata(d,Design$dur.4sU)=1
d=Normalize(d)
g.norm=PlotGeneGroupsBars(d,"TWNK")

wrap_plots(g.norm,g.high,g.low)

low=PairwiseDESeq2(low,"old",contrasts=contr,mode="old")
low=LFC(low,"old",contrasts=contr,mode="old")
VulcanoPlot(low,"old.Auxin vs DMSO.HCT116")
VulcanoPlot(low,"old.Auxin vs DMSO.K562",lfc.cutoff = 0.1)



d=ReadGRAND("/mnt/hilbert/projects/erhard/grandslam/muhar_science_2018/MYC/grandslam_t15/MYC",design=c("Condition","Treatment","Replicate"))
d=FilterGenes(d,minval=300)
Coldata(d,"Treatment2")=factor(gsub("IAA","Auxin",Coldata(d)$Treatment),levels=c("DMSO","Auxin"))
Coldata(d,Design$dur.4sU)=1
d=Normalize(d)
Coldata(d,"Experimental.time")=rep(c(0,1,0,1),each=3)



contr=GetContrasts(d,contrast=c("Treatment2","DMSO"),group = "Condition")
steady.state=FindReferences(d,Treatment2=="DMSO",group="Condition")

d=PairwiseDESeq2(d,"total",contrasts=contr,mode="total")
d=PairwiseDESeq2(d,"new",contrasts=contr,mode="new",separate = TRUE)
d=PairwiseDESeq2(d,"old",contrasts=contr,mode="old")
d=LFC(d,"old",contrasts=contr,mode="old")
d=LFC(d,"new",contrasts=contr,mode="new")
d=LFC(d,"total",contrasts=contr,mode="total")

d=LFC(d,"norm",contrasts=contr,slot = "norm",mode="total",LFC.fun = NormLFC, normalizeFun=function(a)a)


SetParallel()
d=EstimateRegulation(d,"Regulation1",contrasts=contr,steady.state=steady.state,verbose=T,time.experiment = "Experimental.time")
d=EstimateRegulation(d,"Regulation05",contrasts=contr,steady.state=steady.state,verbose=T,time.experiment = "Experimental.time", time.labeling = 0.5)
d=EstimateRegulation(d,"Regulationsamp1",sample.f0.in.ss = TRUE,contrasts=contr,steady.state=steady.state,verbose=T,time.experiment = "Experimental.time", time.labeling = 1)
d=EstimateRegulation(d,"Regulationsamp05",sample.f0.in.ss = TRUE,contrasts=contr,steady.state=steady.state,verbose=T,time.experiment = "Experimental.time", time.labeling = 0.5)
#d=PairwiseNtrTest(d,"NTR",contrasts = contr,verbose = T)
saveRDS(d,"MYC.rds")


df=cbind(GetAnalysisTable(d,"old|new|total",column="LFC|Q|P"),GetAnalysisTable(d,"Regulation",gene.info = FALSE),GetAnalysisTable(d,"NTR",gene.info = FALSE))



df=GetAnalysisTable(d)
PlotScatter(df,x=`new.Auxin vs DMSO.K562.LFC`,y=`old.Auxin vs DMSO.K562.LFC`,remove=F)

VulcanoPlot(d,"total.Auxin vs DMSO.HCT116")
VulcanoPlot(d,"total.Auxin vs DMSO.K562")
VulcanoPlot(d,"new.Auxin vs DMSO.HCT116")
VulcanoPlot(d,"new.Auxin vs DMSO.K562")
VulcanoPlot(d,"old.Auxin vs DMSO.HCT116")
VulcanoPlot(d,"old.Auxin vs DMSO.K562")


power=function(cond) {
  c(
    total=sum(df[[sprintf("total.Auxin vs DMSO.%s.Q",cond)]]<0.05,na.rm=TRUE),
    new=sum(df[[sprintf("new.Auxin vs DMSO.%s.Q",cond)]]<0.05,na.rm=TRUE),
    old=sum(df[[sprintf("old.Auxin vs DMSO.%s.Q",cond)]]<0.05,na.rm=TRUE),
    ntr=sum(df[[sprintf("NTR.Auxin vs DMSO.%s.Q",cond)]]<0.05,na.rm=TRUE),
    combined=sum(p.adjust(pchisq(-2*(log(df[[sprintf("NTR.Auxin vs DMSO.%s.P",cond)]])+log(df[[sprintf("total.Auxin vs DMSO.%s.P",cond)]])),df=4,lower.tail = FALSE),method="BH")<0.05,na.rm=TRUE)
  )
}
power("K562")
power("HCT116")



dispersion.K562 = estimate.dispersion(GetTable(d,type="count",columns = steady.state[,1]))
dispersion.HCT116 = estimate.dispersion(GetTable(d,type="count",columns = steady.state[,7]))
plot.HCT116=function(d,gene,time) {
  re=EstimateGeneRegulation(d,gene,A=contr$`Auxin vs DMSO.HCT116`==1,B=contr$`Auxin vs DMSO.HCT116`==-1,dispersion.A = dispersion.HCT116[ToIndex(d,gene)],dispersion.B = dispersion.HCT116[ToIndex(d,gene)],steady.state = steady.state, return.samples = TRUE,time.experiment = "Experimental.time",time.labeling = time,sample.f0.in.ss = TRUE)
  PlotScatter(re$samples)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+ggtitle("HCT116")
}
plot.K562=function(d,gene,time) {
  re=EstimateGeneRegulation(d,gene,A=contr$`Auxin vs DMSO.K562`==1,B=contr$`Auxin vs DMSO.K562`==-1,dispersion.A = dispersion.K562[ToIndex(d,gene)],dispersion.B = dispersion.K562[ToIndex(d,gene)],steady.state = steady.state, return.samples = TRUE,time.experiment = "Experimental.time",time.labeling = time,sample.f0.in.ss = TRUE)
  PlotScatter(re$samples)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+ggtitle("K562")
}

ServeData(d,df.set=d,plot.single = list(DPlot(plot.K562,time=1),DPlot(plot.K562,time=0.5),DPlot(plot.HCT116,time=1),DPlot(plot.HCT116,time=0.5),DPlot(PlotGeneGroupsBars)))
