
data=ReadGRAND("https://zenodo.org/record/5909576/files/zap_bac4.tsv.gz?download=1",design=c("Cell","hpi","Replicate"),Unknown = "HCMV")
data=NormalizeTPM(data)

cellular.genes=intersect(FilterGenes(data,type='tpm',minval=1,return.genes = TRUE),FilterGenes(filtered,use=filtered$gene.info$Type=="Cellular",return.genes = TRUE))
viral.genes=intersect(FilterGenes(data,use=data$gene.info$Type=="HCMV",return.genes = TRUE),FilterGenes(viral,type='tpm',minval=10,mincol=2,return.genes = TRUE))

data=FilterGenes(data,use=c(cellular.genes,viral.genes))


data=LFC(data,name = "total",GetContrasts(data,contrast=c("Cell","ko","wt"),group="hpi"),mode.slot="total.tpm",prior=c(1,1))
data=LFC(data,name = "new",GetContrasts(data,contrast=c("Cell","ko","wt"),group="hpi"),mode.slot="new.tpm",prior=c(1,1))

PlotScatter(GetAnalysisTable(data),
            "total.18hpi.Cell (ko vs wt).LFC","new.18hpi.Cell (ko vs wt).LFC",
            ylim=c(-1.3,4),xlim=c(-1.3,4),
            highlight = GeneInfo(data)$Type=="HCMV" & rowMeans(GetTable(data,type="tpm",columns=c("wt.18hpi.A","wt.18hpi.B")))>10,
            label=c("UL4","UL5"))+
  geom_abline()+
  geom_abline(intercept = c(-1,1),linetype=2)



PlotScatter(GetAnalysisTable(data),
            "total.72hpi.Cell (ko vs wt).LFC","new.72hpi.Cell (ko vs wt).LFC",
            ylim=c(-1.3,4),xlim=c(-1.3,4),
            highlight = GeneInfo(data)$Type=="HCMV" & rowMeans(GetTable(data,type="tpm",columns=c("wt.72hpi.A","wt.72hpi.B")))>10,
            label=c("UL4","UL5"))+
  geom_abline()+
  geom_abline(intercept = c(-1,1),linetype=2)



data=ReadGRAND("~/mcmv.tsv.gz")
data=FilterGenes(data,min.cond = 1)
data=FilterGenes(data,use=GeneInfo(data)$Type=="Cellular")
data=NormalizeTPM(data)
data=Normalize(data)
Coldata(data,Design$dur.4sU)=2

data=LFC(data,"new",GetContrasts(data,no4sU=FALSE),mode.slot="new.count")
data=LFC(data,"old",GetContrasts(data,no4sU=FALSE),mode.slot="old.count")

PlotScatter(GetAnalysisTable(data),
            "new.mcmv vs mock.LFC","old.mcmv vs mock.LFC",
            remove.outlier = FALSE,
            label=GetAnalysisTable(data)$`new.mcmv vs mock.LFC`>3|GetAnalysisTable(data)$`old.mcmv vs mock.LFC`>2)

data=DissectPairwiseRegulation(data,name="pair",contrasts=GetContrasts(data),steady.state.columns = Coldata(data)$Condition=="mock",seed=1337)

PlotScatter(GetAnalysisTable(data),
            "pair.mcmv vs mock.LFC.s","pair.mcmv vs mock.LFC.HL",
            remove.outlier = FALSE)+
  coord_fixed()

