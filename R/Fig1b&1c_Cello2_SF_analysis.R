####################################################
# Differential expression analysis of SF genes on Cello2 primary-recurrent matched glioma samples
# Include codes to plot Figure 1b and 1c
# Necessary data can be found in ".../data/" directory
####################################################

datapath="/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/Manuscript/GBMSplicing_codes/data/"
# Load patient-level clinical information for CELLO2
clinical<-read.csv(paste0(datapath,"CELLO2_patient_clinical_info.csv"))
IDHwt_tmz_sample<-intersect(clinical$ID_in_Cello2[which(clinical$IDH.mutation=="IDHwt")],clinical$ID_in_Cello2[which(clinical$TMZ.Treatment.Before.Recurrence==1)])
IDHmut_tmz_sample<-intersect(clinical$ID_in_Cello2[which(clinical$IDH.mutation=="IDHmut")],clinical$ID_in_Cello2[which(clinical$TMZ.Treatment.Before.Recurrence==1)])

# Load RNA-seq count matrix for gliomas of CELLO2
count_sum<-readRDS(paste0(datapath,"CELLO2_count_matrix.rds"))
library(DESeq2)
# Differential gene expression for primary and recurrent tumors of IDHwt GBMs.
data<-count_sum[,c(paste0(IDHwt_tmz_sample,"_I"),paste0(IDHwt_tmz_sample,"_R"))]
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
DE_analysis_IDHwt<-as.data.frame(DEresults[!is.na(DEresults$padj),])
DE_analysis_IDHwt$gene<-row.names(DE_analysis_IDHwt)

# Differential gene expression for primary and recurrent tumors of IDHmut gliomas.
data<-count_sum[,c(paste0(IDHmut_tmz_sample,"_I"),paste0(IDHmut_tmz_sample,"_R"))]
anno_sample<-data.frame(sampleID=colnames(data))
anno_sample$condition<-substr(anno_sample$sampleID,7,7)
anno_sample$patient<-substr(anno_sample$sampleID,1,5)

dds<-DESeqDataSetFromMatrix(countData = as.matrix(data),
                            colData = anno_sample,
                            design = ~ condition+patient)
keep<-rowSums(counts(dds) >= 10) >= 5
dds<-dds[keep,]
dds<-DESeq(dds)
contrast_case_vs_ctrl <- c("condition", 'R', 'I')
DEresults<-results(dds, contrast=contrast_case_vs_ctrl)
DEresults<-DEresults[order(DEresults$padj),]
DE_analysis_IDHmut<-as.data.frame(DEresults[!is.na(DEresults$padj),])
DE_analysis_IDHmut$gene<-row.names(DE_analysis_IDHmut)

# Load list of splicing factor (SF) genes
data<-read.csv(paste0(datapath,"RBP.csv"))
sf_gene<-data$HGNC.symbol[which(data$type=="RBP/SF")]

# Plot Figure 1b: Volcano plots of differential expression level of 160 SFs after tumor recurrence for IDHwt and IDHmut gliomas.
library(EnhancedVolcano)
DE_analysis_IDHwt<-DE_analysis_IDHwt[which(DE_analysis_IDHwt$gene%in%sf_gene),]
DE_analysis_IDHmut<-DE_analysis_IDHmut[which(DE_analysis_IDHmut$gene%in%sf_gene),]
EnhancedVolcano(DE_analysis_IDHwt,lab = DE_analysis_IDHwt$gene,
                selectLab = NA,
                drawConnectors = T,
                labSize = 4.5,
                #labFace = 'bold',
                x="log2FoldChange",y="padj",
                pCutoff = 1e-3,
                FCcutoff = 1,
                pointSize = -log10(DE_analysis_IDHwt$padj),
                xlim = c(-2.5,2.5),
                ylim = c(0,10),
                colAlpha = 1,
                cutoffLineCol = '#C6C6C8',
                col = c("#C6C6C8","#C6C6C8","#C6C6C8","#FF4500"),
                gridlines.major = F,
                gridlines.minor = F,
                #border = 'full',
                borderWidth = 0.6,
                borderColour = 'black')

data<-DE_analysis_IDHwt[which(DE_analysis_IDHwt$gene%in%sf_gene),]
sig_diff_sf<-data$gene[intersect(which(data$log2FoldChange>1),which(data$padj<1e-3))]

# Calculate z-scores of 9 DE SFs of all IDHwt GBMs
library(preprocessCore)
new_count_sum<-count_sum[which(rowSums(count_sum)>500),]
new_count_sum<-log2(new_count_sum+1)
new_count_sum<-normalize.quantiles(as.matrix(new_count_sum),keep.names = T)
new_count_sum<-new_count_sum[,c(paste0(IDHwt_tmz_sample,"_I"),paste0(IDHwt_tmz_sample,"_R"))]
z_score_sum<-data.frame(PatientID=colnames(new_count_sum))
for (i in (1:length(sig_diff_sf))) {
  data<-as.data.frame(new_count_sum[sig_diff_sf[i],])
  colnames(data)<-"n_count"
  data$z_score<-(data$n_count-mean(data$n_count))/sd(data$n_count)
  colnames(data)<-c("n_count",sig_diff_sf[i])
  z_score_sum<-cbind(z_score_sum,data[,2])
}
colnames(z_score_sum)<-c("PatientID",sig_diff_sf)
z_score_sum$sum<-apply(z_score_sum[,c(2:10)],1,sum)
z_score_sum<-z_score_sum[order(z_score_sum$sum,decreasing = T),]
z_score_sum$condition<-substr(z_score_sum$PatientID,7,7)
library(mixtools)
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}

# Figure 1c: Fit Gaussian mixture curves of z-scores into 3 groups
set.seed(1)
mixmdl <- normalmixEM(z_score_sum$sum, k = 3) #clustering according to the distribution of "sum" value in "z_score_sum" dataframe, clustering into 3 groups
ggplot(z_score_sum, aes(x = sum,fill=condition)) + #here condition means "Primary" or "Recurrent"
  geom_histogram(bins=50, colour = "black",binwidth = 0.7) +
  scale_fill_manual(values = c("#CEF3EF","#D2B2FC")) + 
  mapply(function(mean, sd, lambda, n, binwidth) {
    stat_function(
      fun = function(x) {
        (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
      },colour = "#3069F6",lwd = 0.8
    )
  },
  mean = mixmdl[["mu"]], 
  sd = mixmdl[["sigma"]],
  lambda = mixmdl[["lambda"]],
  n = length(z_score_sum$sum),
  binwidth = 0.7
  )+
  mapply(function(mean, sd, lambda, n, binwidth) {
    stat_function(
      fun = function(x) {
        (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
      },colour = "#EE9238",lwd = 0.8
    )
  },
  mean = mixmdl[["mu"]][3], 
  sd = mixmdl[["sigma"]][3],
  lambda = mixmdl[["lambda"]][3],
  n = length(z_score_sum$sum),
  binwidth = 0.7
  )+
  ylab("Count") + xlab("Sum of z-scores for 9 DE SFs") + labs(title = "Fitted k=3 Guassian Mixture model for 204 samples", subtitle = " ")+
  theme_classic()

# Extract Patient IDs with HighSF recurrent GBMs.
rec_HighSF_patients<-z_score_sum$PatientID[intersect(which(z_score_sum$condition=="R"),which(z_score_sum$sum>12))]
rec_HighSF_patients<-gsub("_R","",rec_HighSF_patients)
# Save necessary information for the following analysis (DE SFs' gene name and HighSF GBMs' patient IDs)
saveRDS(sig_diff_sf,file = paste0(datapath,"DE_SFs.rds"))
saveRDS(rec_HighSF_patients,file = paste0(datapath,"rec_HighSF_GBMs.rds"))
