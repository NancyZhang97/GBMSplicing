####################################################
# Compare cell compositions for HighSF, LowSF NEU+ and LowSF NEU- GBM patients based on processed snRNA-seq data
# Include codes to plot Figure 3e-f, and supplementary Figure 6.
# Necessary data can be found in ".../data/" directory
####################################################

datapath="/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/Manuscript/GBMSplicing_codes/data/"
# Load clinical information for gliomas in snRNA-seq dataset
meta_sn<-read.csv(paste0(datapath,"snRNA-seq.GBM.metadata.csv"))
# Load Seurat object of processed snRNA-seq data
sn_sum<-readRDS("/Volumes/Samsung_T5/plot_data/wlsc_data.rds")
# Load IDs of HighSF and LowSF NEU+ GBM patients
load(paste0(datapath,"HighSF_LowSF_snRNA_ID.RData"))

# Plot Fig. 3e and supplementary Fig. 6a, changes of cell type proportion during tumor recurrence across three subgroups
data<-as.data.frame(proportions(table(sn_sum$celltype,sn_sum$orig.ident2),margin = 2))
data$condition<-substr(data$Var2,1,1)
data$pair<-gsub("I","",data$Var2)
data$pair<-gsub("R","",data$pair)
data$group<-"LowSF.NEU"
data$group[which(data$pair%in%rec_HighSF_sn_patients)]<-"HighSF"
data$group[which(data$pair%in%rec_LowSF_NEU_sn_patients)]<-"LowSF.NEU+"
data$group<-factor(data$group,levels = c("HighSF","LowSF.NEU+","LowSF.NEU-"))
data<-data[which(data$Var1%in%c("Myeloid","Oligodendrocyte","Neuron","Lymphocyte","Endothelial","Astrocyte")),]
library(ggpubr)
data<-data[order(data$pair),]
data<-data[order(data$condition),]
data0<-data[which(data$Var1=="Neuron"),]
ggpaired(data0,x="condition",y="Freq",fill = "condition",palette = c("#51B1B2","#8A00FF"),facet.by = "group")+
  stat_compare_means(paired = T, label = "p.format",method = "t.test")+
  ylab("Proportion of Neurons")+theme_classic()

# Plot supplementary Fig. 6b, changes of relative tumor substate abundance during tumor recurrence across three subgroups
data0<-sn_sum@meta.data
data0<-data0[which(data0$celltype=="Tumor"),]
data<-as.data.frame(proportions(table(data0$celltype_subtype,data0$orig.ident2),margin = 2))
data$condition<-substr(data$Var2,1,1)
data$pair<-gsub("I","",data$Var2)
data$pair<-gsub("R","",data$pair)
data$group<-"LowSF.NEU"
data$group[which(data$pair%in%rec_HighSF_sn_patients)]<-"HighSF"
data$group[which(data$pair%in%rec_LowSF_NEU_sn_patients)]<-"LowSF.NEU+"
data$group<-factor(data$group,levels = c("HighSF","LowSF.NEU+","LowSF.NEU-"))
data<-data[order(data$pair),]
data<-data[order(data$condition),]
data0<-data[which(data$Var1=="Tumor.NPC"),]
ggpaired(data0,x="condition",y="Freq",fill = "condition",palette = c("#51B1B2","#8A00FF"),facet.by = "group")+
  stat_compare_means(paired = T, label = "p.format")+
  ylab("Tumor abundance of Tumor.NPC")+theme_classic()

# Plot Fig. 3f, changes of relative tumor substate abundance during tumor recurrence across three subgroups
cell_meta<-sn_sum@meta.data
cell_meta$group<-"LowSF.NEU-"
cell_meta$group[which(cell_meta$pair%in%rec_HighSF_sn_patients)]<-"HighSF"
cell_meta$group[which(cell_meta$pair%in%rec_LowSF_NEU_sn_patients)]<-"LowSF.NEU+"
data0<-cell_meta[which(cell_meta$condition=="Primary"),]
data0<-data0[which(data0$celltype=="Tumor"),]
data<-as.data.frame(proportions(table(data0$group,data0$celltype_subtype),1))
data$condition<-"Primary"
data0<-cell_meta[which(cell_meta$condition=="Recurrent"),]
data0<-data0[which(data0$celltype=="Tumor"),]
data0<-as.data.frame(proportions(table(data0$group,data0$celltype_subtype),1))
data0$condition<-"Recurrent"
data<-rbind(data,data0)
data$condition2<-ifelse(data$condition=="Primary",0,1)
tumor_sub = c(Tumor.Cycle="#C82B1E",Tumor.OPC="#7A634B",Tumor.AC="#AC9C88",Tumor.NPC="#E7AC53",
             Tumor.MES="#4A885C",Tumor.Stem="#8E39C5",Tumor.Fiber="#CD7692")
data$Var2<-factor(data$Var2,levels = c("Tumor.NPC","Tumor.Cycle","Tumor.OPC","Tumor.AC","Tumor.MES","Tumor.Stem","Tumor.Fiber"))
data$Var1<-factor(data$Var1,levels = c("HighSF","LowSF.NEU+","LowSF.NEU-"))
ggplot(data, aes(x=condition2, y=Freq, fill=Var2)) + 
  geom_area(alpha=1 , size=0.5, colour="black")+scale_fill_manual(values = tumor_sub)+
  facet_wrap(~Var1)+theme_classic2()

# Plot Fig. 3g, violin plots of chr7 and chr10 CNV scores and dot plots of 9 DE SFs across different cell types
library(scCustomize)
sn_sum.sub<-sn_sum[,-which(is.na(sn_sum$chr7_cnv))]
Stacked_VlnPlot(seurat_object = sn_sum.sub, features = c("chr7_cnv","chr10_cnv"),x_lab_rotate = TRUE,
                group.by = "celltype_subtype",colors_use = c("#91D1C2FF","#8491B4FF","#F39B7FFF","#3C5488FF","#00A087FF","#4DBBD5FF","#AC9C88","#C82B1E","#CD7692","#4A885C","#E7AC53","#7A634B","#8E39C5"))
# Load 9 DE SFs' labels
sig_diff_sf<-readRDS(paste0(datapath,"DE_SFs.rds"))
library(viridis)
DotPlot(sn_sum,features = c(sig_diff_sf),group.by = "celltype_subtype") +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_viridis(option="magma") +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))+
  coord_flip()+rotate_x_text(angle = 45)

# Plot Fig. 3h, bar plot of GSEA results comparing Tumor. NPC with other tumor substates
sn_tumor<-subset(sn_sum,subset=celltype=="Tumor")
de_tumor_npc<-FindMarkers(sn_tumor,group.by = "celltype_subtype",ident.1 = "Tumor.NPC")
library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(tibble)
de_tumor_npc$gene<-row.names(de_tumor_npc)
gsea_genes<-de_tumor_npc %>%
  arrange(desc(p_val_adj), desc(avg_log2FC)) %>%
  dplyr::select(gene,avg_log2FC)
ranks <- deframe(gsea_genes)
m_df<- msigdbr(species = "Homo sapiens", category = "C5",subcategory = "BP")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks ,
                  minSize=5,
                  maxSize=500,
                  nPermSimple = 10000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>% arrange(desc(NES))
fgseaResTidy$pathway<-gsub("GOBP_","",fgseaResTidy$pathway)
data<-fgseaResTidy[which(fgseaResTidy$pathway%in%c("NEUROTRANSMITTER_SECRETION","CHEMICAL_SYNAPTIC_TRANSMISSION_POSTSYNAPTIC","REGULATION_OF_MEMBRANE_POTENTIAL",
                                                   "GLUTAMATE_RECEPTOR_SIGNALING_PATHWAY","DENDRITE_DEVELOPMENT","EXOCYTOSIS","ALTERNATIVE_MRNA_SPLICING_VIA_SPLICEOSOME",
                                                   "MRNA_SPLICE_SITE_SELECTION","CENTRAL_NERVOUS_SYSTEM_NEURON_DIFFERENTIATION","EXTERNAL_ENCAPSULATING_STRUCTURE_ORGANIZATION",              
                                                   "VASCULATURE_DEVELOPMENT","BLOOD_VESSEL_MORPHOGENESIS","INNATE_IMMUNE_RESPONSE","INFLAMMATORY_RESPONSE",                                      
                                                   "POSITIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION","POSITIVE_REGULATION_OF_CELL_CELL_ADHESION",                  
                                                   "CELL_SUBSTRATE_ADHESION","REGULATION_OF_MITOTIC_CELL_CYCLE")),]
ggplot(data, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill= NES < 0)) + scale_fill_manual(values = c("#B22222","#008080")) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GOBP pathways") +
  theme_minimal()
