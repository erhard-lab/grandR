
d=ReadGRAND("/mnt/hilbert/projects/erhard/grandslam/ameres/miRNA/grandslam_t15_notall/slam_ameres_mirna",design=c("Cell","Condition",Design$dur.4sU,"Replicate"), rename.sample = function(n) gsub(".pulse","",n))
d=subset(d,Coldata(d)$duration.4sU==3)
d=FilterGenes(d,minval = 500,mincol = 6)
d=Normalize(d)


contrasts=GetContrasts(d,contrast=c("Condition","wt"))
steady.state=FindReferences(d,TRUE,group="Condition")

d=PairwiseDESeq2(d,"total",contrasts=contrasts,mode="total")
d=PairwiseDESeq2(d,"new",contrasts=contrasts,mode="new",separate = TRUE)
d=PairwiseDESeq2(d,"old",contrasts=contrasts,mode="old")
d=LFC(d,"old",contrasts=contrasts,mode="old")
d=LFC(d,"new",contrasts=contrasts,mode="new")
d=LFC(d,"total",contrasts=contrasts,mode="total")

SetParallel()
d=EstimateRegulation(d,"Regulation",contrasts=contrasts,steady.state=steady.state,verbose=T)
df=GetAnalysisTable(d)

sm=read.delim("/mnt/hilbert/projects/erhard/grandslam/paper_plots/final/f7/seedmatches.tsv")
sm$Symbol=gsub("\\..*","",sm$Transcript)
sm$Seed=factor(sm$Seed,levels=c("No","6mer","7mer A1","7mer m8","8mer"))
sm=ddply(sm,.(Symbol),function(s) data.frame(Seed=s$Seed[order(s$Seed,decreasing=TRUE)][1]))
df=merge(df,sm,by="Symbol")

ggplot(df,aes(`Regulation.Xpo5KO vs wt.HL.log2FC`,color=Seed))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))
ggplot(df,aes(`Regulation.Xpo5KO vs wt.s.log2FC`,color=Seed))+stat_ecdf()+coord_cartesian(xlim=c(-2,2))

ggplot(df,aes(Seed,`Regulation.Xpo5KO vs wt.HL.log2FC`,color=Seed))+geom_boxplot()
ggplot(df,aes(Seed,`Regulation.Xpo5KO vs wt.s.log2FC`,color=Seed))+geom_boxplot()

ggplot(df[df$Seed=="No",],aes(x=`Regulation.Xpo5KO vs wt.s.A`,y=`Regulation.Xpo5KO vs wt.s.B`))+geom_point(color='gray',size=0.7)+geom_point(data=df[df$Seed=="8mer",],color='red',size=0.7)+coord_cartesian(xlim=c(1,5000),ylim=c(1,5000))+scale_x_log10()+scale_y_log10()

