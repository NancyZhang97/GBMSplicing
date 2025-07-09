####################################################
# Comparisons of clinical information and molecular features across HighSF, LowSF NEU+ and LowSF NEU- groups.
# Include codes to plot Figure 1e, Figure 2, supplementary Figure 2 and Figure 3a.
# Necessary data can be found in ".../data/" directory
####################################################

datapath="/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/Manuscript/GBMSplicing_codes/data/"
# Load patient's information and HighSF GBMs' IDs.
clinical<-read.csv(paste0(datapath,"CELLO2_patient_clinical_info.csv"))
IDHwt_tmz_sample<-intersect(clinical$ID_in_Cello2[which(clinical$IDH.mutation=="IDHwt")],clinical$ID_in_Cello2[which(clinical$TMZ.Treatment.Before.Recurrence==1)])
rec_HighSF_patients<-readRDS(paste0(datapath,"rec_HighSF_GBMs.rds"))

# Load TPM matrix of CELLO2 for IDHwt GBMs
tpm_sum<-readRDS(paste0(datapath,"CELLO2_TPM_matrix.rds"))
tpm_sum<-tpm_sum[,c(paste0(IDHwt_tmz_sample,"_I"),paste0(IDHwt_tmz_sample,"_R"))]

# Load marker genes of four GBM subtypes
genelist<-readRDS(paste0(datapath,"genelist_fourGBMsubtype.rds"))
# Classify IDHwt GBMs into four subtypes
library(GSVA)
gbm_es_tpm<-gsva(as.matrix(tpm_sum), genelist, mx.diff=FALSE,method = "ssgsea")
#Find biggest score
gbm_tpm_grp4<-apply(gbm_es_tpm, 2, function(x) which(x==max(x)))
grp4_ssgsea<-ifelse(gbm_tpm_grp4==1,"Classical",ifelse(gbm_tpm_grp4==2,"Mesenchymal",ifelse(gbm_tpm_grp4==3,"Neural","Proneural")))
grp4_ssgsea<-as.data.frame(grp4_ssgsea)
grp4_ssgsea$sampleID<-row.names(grp4_ssgsea)
grp4_ssgsea$Patient<-substr(grp4_ssgsea$sampleID,1,5)
grp4_ssgsea$condition<-substr(grp4_ssgsea$sampleID,7,7)
library(reshape2)
grp4_ssgsea_short<-dcast(grp4_ssgsea,Patient~condition,value.var = "grp4_ssgsea")

# Plot Fig. 1f, proportional barplot of four-subtype classification of recurrent GBMs between HighSF and LowSF groups
data<-grp4_ssgsea_short
data$group<-ifelse(data$Patient%in%rec_HighSF_patients,"HighSF","LowSF")
data0<-as.data.frame(table(data$group,data$R))
colnames(data0)<-c("group","subtype(R)","Count")
data0$`subtype(R)`<-factor(data0$`subtype(R)`,levels = c("Mesenchymal","Classical","Proneural","Neural"))
library(ggplot2)
ggplot(data0, aes(x = group, y = Count, fill = `subtype(R)`)) +
  geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#4682B4","#C71585","#9ACD32","#FF6347"))+theme_classic()

# Group LowSF patients into LowSF NEU+ and LowSF NEU-
rec_LowSF_NEU_patients<-grp4_ssgsea_short$Patient[which(grp4_ssgsea_short$R=="Neural")]
rec_LowSF_NEU_patients<-setdiff(rec_LowSF_NEU_patients,rec_HighSF_patients)

# Plot Fig. 2a, clinical information for three groups of GBM patients
IDHwt_tmz_clinical<-clinical[which(clinical$ID_in_Cello2%in%IDHwt_tmz_sample),]
IDHwt_tmz_clinical<-IDHwt_tmz_clinical[,c("ID_in_Cello2","Source","Grade.of.Initial.Tumor","Grade.of.Recurrent.Tumor","Gender","Race","Hypermutation.of.Recurrent.Tumor")]
row.names(IDHwt_tmz_clinical)<-IDHwt_tmz_clinical$ID_in_Cello2
IDHwt_tmz_clinical$group<-"LowSF_noNEU"
IDHwt_tmz_clinical$group[which(IDHwt_tmz_clinical$ID_in_Cello2%in%rec_LowSF_NEU_patients)]<-"LowSF_NEU"
IDHwt_tmz_clinical$group[which(IDHwt_tmz_clinical$ID_in_Cello2%in%rec_HighSF_patients)]<-"HighSF"
IDHwt_tmz_clinical$group<-factor(IDHwt_tmz_clinical$group,levels = c("HighSF","LowSF_NEU","LowSF_noNEU"))
data<-t(IDHwt_tmz_clinical[order(IDHwt_tmz_clinical$group,IDHwt_tmz_clinical$ID_in_Cello2),])
data<-data[2:7,]
library(ComplexHeatmap)
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
  CGGA = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#B22222", col = NA)),
  NG2016 = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#0066FFFF", col = NA)),
  NM2019 = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#FF7F50", col = NA)),
  SMCnew = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#6A5ACD", col = NA)),
  TCGA = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#FFD700", col = NA)),
  Asian = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#006400", col = NA)),
  `non-Asian` = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#DB7093", col = NA)),
  M = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#6495ED", col = NA)),
  `F` = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CD5C5C", col = NA)),
  Yes = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#BA55D3", col = NA)),
  No = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#00CED1", col = NA)),
  II = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#E6E6FA", col = NA)),
  III = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#8B008B", col = NA)),
  IV = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#9370DB", col = NA))
)
oncoPrint(data, alter_fun = alter_fun,row_order = 1:nrow(data),column_order = 1:ncol(data),
          col = c(CGGA = "#B22222", NG2016 = "#0066FFFF",NM2019="#FF7F50",
                  SMCnew ="#6A5ACD", TCGA ="#FFD700", Asian="#006400", `non-Asian`="#DB7093",
                  M="#6495ED", `F`="#CD5C5C", Yes="#BA55D3", 
                  No="#00CED1", II="#E6E6FA", III="#8B008B", IV="#9370DB"))

# Plot supplementary Fig. 2a, comparison of age at initial diagnosis across 3 GBM subgroups
IDHwt_tmz_clinical<-clinical[which(clinical$ID_in_Cello2%in%IDHwt_tmz_sample),]
IDHwt_tmz_clinical$Age.at.initial.diagnosis<-as.numeric(IDHwt_tmz_clinical$Age.at.initial.diagnosis)
IDHwt_tmz_clinical$group<-"LowSF_noNEU"
IDHwt_tmz_clinical$group[which(IDHwt_tmz_clinical$ID_in_Cello2%in%rec_LowSF_NEU_patients)]<-"LowSF_NEU"
IDHwt_tmz_clinical$group[which(IDHwt_tmz_clinical$ID_in_Cello2%in%rec_HighSF_patients)]<-"HighSF"
IDHwt_tmz_clinical$group<-factor(IDHwt_tmz_clinical$group,levels = c("LowSF_noNEU","LowSF_NEU","HighSF"))
library(ggpubr)
comparison<-list(c("LowSF_noNEU","LowSF_NEU"),c("LowSF_NEU","HighSF"),c("LowSF_noNEU","HighSF"))
ggboxplot(IDHwt_tmz_clinical,x="group",y="Age.at.initial.diagnosis",fill = "group",palette = c("#6C7FAE","#3069F6","#EE9238"))+
  stat_compare_means(comparisons = comparison,label = "p.format",method = "t.test")

# Plot supplementary Fig. 2b, flow chart of subtype transition from primary to recurrent GBMs
library(ggalluvial)
grp4_ssgsea_short$group2<-"LowSF_noNEU"
grp4_ssgsea_short$group2[which(grp4_ssgsea_short$R=="Neural")]<-"LowSF_NEU"
grp4_ssgsea_short$group2[which(grp4_ssgsea_short$Patient%in%rec_HighSF_patients)]<-"HighSF"
plot_table<-as.data.frame(table(grp4_ssgsea_short$I,grp4_ssgsea_short$R,grp4_ssgsea_short$group2))
colnames(plot_table)<-c("Subtype.I","Subtype.R","group","Count")
ggplot(plot_table,
       aes(axis1 = Subtype.I,
           axis2 = Subtype.R,
           y = Count)) +
  geom_alluvium(aes(fill = group)) +
  scale_fill_manual(values = c(HighSF = "#EE9238", LowSF_NEU = "#3069F6", LowSF_noNEU = "#6C7FAE"))+
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),size=3) +
  scale_x_discrete(limits = c("Primary", "Recurrent"),
                   expand = c(.1, .1)) +
  theme_minimal()

# Plot supplementary Fig. 2c, flow chart of subtype transition from primary to recurrent GBMs
data<-as.data.frame(table(grp4_ssgsea_short$I=="Mesenchymal",grp4_ssgsea_short$group2))
colnames(data)<-c("I","group","Count")
ggplot(data, aes(x = group, y = Count, fill = I)) +
  geom_col(colour = "black", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(`TRUE` = "#DC143C", `FALSE` = "#008080"))+theme_classic2()

# Plot Fig. 2d, changes on enrichment scores of mesenchymal and neural subtype
data<-grp4_ssgsea
data$group<-NA
data$group[which(data$Patient%in%rec_LowSF_NEU_patients)]<-"LowSF_NEU"
data$group[which(data$Patient%in%rec_HighSF_patients)]<-"HighSF"
data$group[which(is.na(data$group))]<-"LowSF_noNEU"
gbm_es_tpm<-as.data.frame(t(gbm_es_tpm))
gbm_es_tpm<-gbm_es_tpm[data$sampleID,]
data<-cbind(data,gbm_es_tpm)
ggpaired(data,x="condition",y="VERHAAK_GLIOBLASTOMA_MESENCHYMAL",facet.by = "group",fill = "condition",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T)+ theme_classic()+ ylab("Enrichment score of MES markers")
ggpaired(data,x="condition",y="VERHAAK_GLIOBLASTOMA_NEURAL",facet.by = "group",fill = "condition",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T)+ theme_classic()+ ylab("Enrichment score of NEU markers")

# Plot Fig. 2e and supplementary Fig. 2e, changes on expression level of TMZ-resistance markers
data$exp<-tpm_sum["MGMT",row.names(data)]
ggpaired(data,x="condition",y="exp",facet.by = "group",fill = "condition",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T)+ theme_classic()+ ylab("Normalized TPM of MGMT")

# Plot supplementary Fig. 2f, comparisions of expression level of TMZ-resistance markers for recurrent GBMs
data$exp<-tpm_sum["MGMT",row.names(data)]
data0<-data[which(data$condition=="R"),]
data0$group<-factor(data0$group,levels = c("LowSF_noNEU","LowSF_NEU","HighSF"))
comparison<-list(c("LowSF_noNEU","LowSF_NEU"),c("LowSF_NEU","HighSF"),c("LowSF_noNEU","HighSF"))
ggboxplot(data0,x="group",y="exp",fill = "group",palette = c("#6C7FAE","#3069F6","#EE9238"))+
  stat_compare_means(comparisons = comparison,label = "p.format",method = "t.test")

# Load RNA-seq count matrix for gliomas of CELLO2
count_sum<-readRDS(paste0(datapath,"CELLO2_count_matrix.rds"))
library(DESeq2)
# Differential gene expression for primary and recurrent tumors of IDHwt GBMs.
# HighSF patients: Recurrent vs Primary
data<-count_sum[,c(paste0(rec_HighSF_patients,"_I"),paste0(rec_HighSF_patients,"_R"))]
anno_sample<-data.frame(sampleID=colnames(data))
anno_sample$condition<-substr(anno_sample$sampleID,7,7)
anno_sample$patient<-substr(anno_sample$sampleID,1,5)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(data),
                              colData = anno_sample,
                              design = ~ condition+patient)
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]
dds <- DESeq(dds)
contrast_case_vs_ctrl <- c("condition", 'R', 'I')
DEresults <- results(dds, contrast = contrast_case_vs_ctrl)
DEresults <- DEresults[order(DEresults$padj),]
result_pair16<-as.data.frame(DEresults[!is.na(DEresults$padj),])
result_pair16$gene<-row.names(result_pair16)

# LowSF NEU+ patients: Recurrent vs Primary
data<-count_sum[,c(paste0(rec_LowSF_NEU_patients,"_I"),paste0(rec_LowSF_NEU_patients,"_R"))]
anno_sample<-data.frame(sampleID=colnames(data))
anno_sample$condition<-substr(anno_sample$sampleID,7,7)
anno_sample$patient<-substr(anno_sample$sampleID,1,5)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(data),
                              colData = anno_sample,
                              design = ~ condition+patient)
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]
dds <- DESeq(dds)
contrast_case_vs_ctrl <- c("condition", 'R', 'I')
DEresults <- results(dds, contrast = contrast_case_vs_ctrl)
DEresults <- DEresults[order(DEresults$padj),]
result_pair_LowSF_NEU<-as.data.frame(DEresults[!is.na(DEresults$padj),])
result_pair_LowSF_NEU$gene<-row.names(result_pair_LowSF_NEU)

# LowSF NEU- patients: Recurrent vs Primary
data<-count_sum[,c(paste0(setdiff(IDHwt_tmz_sample,union(rec_HighSF_patients,rec_LowSF_NEU_patients)),"_I"),paste0(setdiff(IDHwt_tmz_sample,union(rec_HighSF_patients,rec_LowSF_NEU_patients)),"_R"))]
anno_sample<-data.frame(sampleID=colnames(data))
anno_sample$condition<-substr(anno_sample$sampleID,7,7)
anno_sample$patient<-substr(anno_sample$sampleID,1,5)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(data),
                              colData = anno_sample,
                              design = ~ condition+patient)
keep <- rowSums(counts(dds) >= 10) >= 5
dds <- dds[keep,]
dds <- DESeq(dds)
contrast_case_vs_ctrl <- c("condition", 'R', 'I')
DEresults <- results(dds, contrast = contrast_case_vs_ctrl)
DEresults <- DEresults[order(DEresults$padj),]
result_pair_LowSF_noNEU<-as.data.frame(DEresults[!is.na(DEresults$padj),])
result_pair_LowSF_noNEU$gene<-row.names(result_pair_LowSF_noNEU)

# Plot Fig. 2c, GSEA results for primary-to-recurrent changes across three GBM subgroups
library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(tibble)
# HighSF patients: GSEA
result_pair16$gene<-row.names(result_pair16)
gsea_genes<-result_pair16 %>%
  arrange(desc(padj), desc(log2FoldChange)) %>%
  dplyr::select(gene,log2FoldChange)
ranks <- deframe(gsea_genes)
m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks ,
                  minSize=5,
                  maxSize=500,
                  nPermSimple = 10000)
fgseaResTidy_HighSF <- fgseaRes %>%
  as_tibble() %>% arrange(desc(NES))

# LowSF NEU+ patients: GSEA
result_pair_LowSF_NEU$gene<-row.names(result_pair_LowSF_NEU)
gsea_genes<-result_pair_LowSF_NEU %>%
  arrange(desc(padj), desc(log2FoldChange)) %>%
  dplyr::select(gene,log2FoldChange)
ranks <- deframe(gsea_genes)
m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks ,
                  minSize=5,
                  maxSize=500,
                  nPermSimple = 10000)
fgseaResTidy_LowSF_NEU <- fgseaRes %>%
  as_tibble() %>% arrange(desc(NES))

# LowSF NEU- patients: GSEA
result_pair_LowSF_noNEU$gene<-row.names(result_pair_LowSF_noNEU)
gsea_genes<-result_pair_LowSF_noNEU %>%
  arrange(desc(padj), desc(log2FoldChange)) %>%
  dplyr::select(gene,log2FoldChange)
ranks <- deframe(gsea_genes)
m_df<- msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name) 
fgseaRes <- fgsea(pathways = fgsea_sets,
                  stats = ranks ,
                  minSize=5,
                  maxSize=500,
                  nPermSimple = 10000)
fgseaResTidy_LowSF_noNEU <- fgseaRes %>%
  as_tibble() %>% arrange(desc(NES))

fgseaResTidy<-rbind(fgseaResTidy_HighSF,rbind(fgseaResTidy_LowSF_NEU,fgseaResTidy_LowSF_noNEU))
pathway_select<-c("KEGG_MISMATCH_REPAIR","HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","KEGG_CELL_CYCLE","KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
                  "KEGG_CALCIUM_SIGNALING_PATHWAY","KEGG_AXON_GUIDANCE","KEGG_GAP_JUNCTION","KEGG_FOCAL_ADHESION","KEGG_ECM_RECEPTOR_INTERACTION",
                  "KEGG_P53_SIGNALING_PATHWAY","KEGG_LONG_TERM_POTENTIATION","KEGG_PHOSPHATIDYLINOSITOL_SIGNALING_SYSTEM","KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION")
fgseaResTidy_HighSF<-fgseaResTidy_HighSF[order(fgseaResTidy_HighSF$NES),]
data<-fgseaResTidy[which(fgseaResTidy$pathway%in%pathway_select),]
data$log_pval<-(-log10(data$padj))
data$pathway<-factor(data$pathway,levels = fgseaResTidy_HighSF$pathway)
data<-data[order(data$pathway),]
library(forcats)
data %>%
  arrange(pathway) %>%
  ggplot(aes(NES, pathway)) +
  geom_line(aes(group = pathway)) +
  geom_point(aes(color = group, size = log_pval))+
  scale_color_manual(values = c(HighSF = "#EE9238", LowSF_NEU = "#3069F6", LowSF_noNEU = "#6C7FAE"))+theme_pubclean() + theme(legend.position = "right")

# Plot Fig. 3a, comparison of tumor purity across three subgroups of recurrent GBMs
# Load tumor purity information estimated by FACET
purity_facet<-read.csv(paste0(datapath,"CELLO2_tumor_purity_FACET.csv"))
data<-purity_facet[which(purity_facet$Patient%in%IDHwt_tmz_sample),]
data<-data[which(data$Timepoint=="R"),]
data$group<-"LowSF_noNEU"
data$group[which(data$Patient%in%rec_LowSF_NEU_patients)]<-"LowSF_NEU"
data$group[which(data$Patient%in%rec_HighSF_patients)]<-"HighSF"
data$group<-factor(data$group,levels = c("LowSF_noNEU","LowSF_NEU","HighSF"))
comparison<-list(c("LowSF_noNEU","LowSF_NEU"),c("LowSF_NEU","HighSF"),c("LowSF_noNEU","HighSF"))
ggbarplot(data, x = "group", y = "Purity",fill = "group", palette = c("#6C7FAE","#3069F6","#EE9238"),
          add = "mean_se")+stat_compare_means(comparisons = comparison,label = "p.format",label.y = c(0.65,0.7,0.75),method = "t.test")+
  ylab("Tumor Purity")+theme_classic2()+rotate_x_text(angle = 45)
