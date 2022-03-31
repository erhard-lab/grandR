
d=ReadGRAND("/mnt/hilbert/projects/wolf/Spt6/SLAM-seq/Spt6_SLAM_fullrun_norescue/grandslam_t15/Spt6_SLAM_fullrun_norescue",design=c("Condition",Design$dur.4sU,"Replicate"))
d=subset(d,columns=Coldata(d)$duration.4sU<=1)
d=FilterGenes(d,minval = 200)
d=Normalize(d)
Coldata(d,"Experimental.time")=rep(c(6,7),each=6)



contrasts=GetContrasts(d,contrast=c("Condition","Spt6+"),columns = Coldata(d)$duration.4sU==1)
steady.state.mat=cbind(FindReferences(d,TRUE,group="duration.4sU", columns=Coldata(d)$Condition=="Spt6+"),FindReferences(d,duration.4sU==0,columns=Coldata(d)$Condition=="Spt6-"))

d=PairwiseDESeq2(d,"total",contrasts=contrasts,mode="total")
d=PairwiseDESeq2(d,"new",contrasts=contrasts,mode="new",separate = TRUE)
d=PairwiseDESeq2(d,"old",contrasts=contrasts,mode="old")
d=LFC(d,"old",contrasts=contrasts,mode="old")
d=LFC(d,"new",contrasts=contrasts,mode="new")
d=LFC(d,"total",contrasts=contrasts,mode="total")

SetParallel()
d=DropAnalysis(d,"Regu")
d=EstimateRegulation(d,"Regulation1",contrasts=contrasts,steady.state=steady.state.mat,verbose=T,time.experiment = "Experimental.time")
#d=EstimateRegulation(d,"Regulation05",contrasts=contr,steady.state=steady.state.mat,verbose=T,time.experiment = "Experimental.time", time.labeling = 0.5)
#d=PairwiseNtrTest(d,"NTR",contrasts = contr,verbose = T)

df=cbind(GetAnalysisTable(d,"old|new|total",column="LFC|Q|P"),GetAnalysisTable(d,"Regulation",gene.info = FALSE))


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



dispersion = estimate.dispersion(GetTable(d,type="count"))
plot.post=function(d,gene,time) {
  re=EstimateGeneRegulation(d,gene,A=contrasts$`Spt6- vs Spt6+`==1,B=contrasts$`Spt6- vs Spt6+`==-1,dispersion.A = dispersion,dispersion.B = dispersion,steady.state = steady.state.mat, return.samples = TRUE,time.experiment = "Experimental.time",time.labeling = time)
  PlotScatter(re$samples)+geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+ggtitle("HCT116")
}

ServeData(d,df=df,df.set=d,plot.single = list(DPlot(plot.post,time=1),DPlot(plot.post,time=0.5),DPlot(PlotGeneGroupsBars)))



ntrs=GetTable(d,type="ntr",columns=!no4sU & Condition=="Spt6+")
df=data.frame(mean=rowMeans(ntrs),apply(ntrs,1,var))
PlotScatter(df,x=mean,y=mean*(1-mean)/var-1,log=T,remove.outlier = F)

