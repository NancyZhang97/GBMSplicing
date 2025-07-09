####################################################
# Analysis of changes on cell-cell interaction during tumor recurrence for HighSF GBM patients using CellChat
# Include codes to plot Figure 4a-b, and supplementary Figure 7a-b.
# Necessary data can be found in ".../data/" directory
####################################################

datapath="/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/Manuscript/GBMSplicing_codes/data/"
# Load clinical information for gliomas in snRNA-seq dataset
meta_sn<-read.csv(paste0(datapath,"snRNA-seq.GBM.metadata.csv"))
# Load Seurat object of processed snRNA-seq data
sn_sum<-readRDS("/Volumes/Samsung_T5/plot_data/wlsc_data.rds")
# Load IDs of HighSF and LowSF NEU+ GBM patients
load(paste0(datapath,"HighSF_LowSF_snRNA_ID.RData"))

sn_sum$group<-"LowSF.NEU-"
sn_sum$group[which(sn_sum$pair%in%rec_HighSF_sn_patients)]<-"HighSF"
sn_sum$group[which(sn_sum$pair%in%rec_LowSF_NEU_sn_patients)]<-"LowSF.NEU+"
sn_sum$celltype_subtype<-factor(sn_sum$celltype_subtype,levels = c("Neuron","Tumor.NPC","Tumor.Cycle","Tumor.OPC","Tumor.AC","Tumor.MES","Tumor.Stem","Tumor.Fiber","Oligodendrocyte","Myeloid","Lymphocyte","Endothelial","Astrocyte"))

# Apply CellChat to estimate cell-cell interactions in primary and recurrent tumors of HighSF GBM patients
library(CellChat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
data.input<-subset(sn_sum,subset=condition=="Recurrent")
data.input<-subset(data.input,subset=group=="HighSF")
meta<-data.input@meta.data
data.input<-data.input@assays$RNA@data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype_subtype")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat<-netAnalysis_computeCentrality(cellchat)
df.net <- subsetCommunication(cellchat)
df.net.R<-df.net
cellchat.R<-cellchat

data.input<-subset(sn_sum,subset=condition=="Primary")
data.input<-subset(data.input,subset=group=="HighSF")
meta<-data.input@meta.data
data.input<-data.input@assays$RNA@data
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype_subtype")
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat<-netAnalysis_computeCentrality(cellchat)
df.net <- subsetCommunication(cellchat)
df.net.P<-df.net
cellchat.P<-cellchat

# Compare cell-cell interactions between primary and recurrent tumors of HighSF GBM patients
object.list <- list(Primary = cellchat.P, Recurrent = cellchat.R)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# Plot Fig. 4a, differential number of cell-cell interactions comparing HighSF recurrent tumors to their matched primary tumors
colors<-c("#00A087FF","#E7AC53","#C82B1E","#7A634B","#8E39C5","#4A885C","#AC9C88","#CD7692","#4DBBD5FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF")
netVisual_heatmap(cellchat,color.heatmap = c(blu = "royalblue", reddish = "tomato"),color.use = colors)

# Plot Fig. 4b, differential number of cell-cell interactions in neurotransmitter signaling pathways
data_P<-df.net.P[which(df.net.P$pathway_name=="Glutamate"),]
data<-as.data.frame(table(data_P$source,data_P$target))
colnames(data)<-c("Source","Target","Freq")
data_R<-df.net.R[which(df.net.R$pathway_name=="Glutamate"),]
data0<-as.data.frame(table(data_R$source,data_R$target))
colnames(data0)<-c("Source","Target","Freq")
diff_cellchat<-data
diff_cellchat$Freq<-data0$Freq-data$Freq
ggplot(diff_cellchat, aes(x=Target, y=Source, fill=Freq)) +
  geom_tile(color="white", size = 0.25) +
  scale_fill_gradient2(high = "red",low = "blue")+theme_minimal()

data_P<-df.net.P[which(df.net.P$pathway_name=="GABA-A"),]
data<-as.data.frame(table(data_P$source,data_P$target))
colnames(data)<-c("Source","Target","Freq")
data_R<-df.net.R[which(df.net.R$pathway_name=="GABA-A"),]
data0<-as.data.frame(table(data_R$source,data_R$target))
colnames(data0)<-c("Source","Target","Freq")
diff_cellchat<-data
diff_cellchat$Freq<-data0$Freq-data$Freq
ggplot(diff_cellchat, aes(x=Target, y=Source, fill=Freq)) +
  geom_tile(color="white", size = 0.25) +
  scale_fill_gradient2(high = "red",low = "blue")+theme_minimal()

# Plot supplementary Fig. 7a and 7b, comparison of relative information flow from neurons to Tumor. NPC cells and from Tumor. NPC cells to neurons
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = F,sources.use = "Neuron",targets.use = "Tumor.NPC",color.use = c("#008B8B","#7F17F5"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = F,sources.use = "Tumor.NPC",targets.use = "Neuron",color.use = c("#008B8B","#7F17F5"))
gg1 + gg2
