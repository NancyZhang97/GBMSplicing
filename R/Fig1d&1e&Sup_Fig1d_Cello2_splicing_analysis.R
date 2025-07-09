####################################################
# Differential splicing analysis of HighSF and LowSF recurrent GBMs from Cello2 cohort
# Include codes to plot Figure 1d, 1e and supplementary Figure 1d
# Necessary data can be found in ".../data/" directory
####################################################

datapath="/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/Manuscript/GBMSplicing_codes/data/"
# Load patient's information and HighSF GBMs' IDs.
clinical<-read.csv(paste0(datapath,"CELLO2_patient_clinical_info.csv"))
IDHwt_tmz_sample<-intersect(clinical$ID_in_Cello2[which(clinical$IDH.mutation=="IDHwt")],clinical$ID_in_Cello2[which(clinical$TMZ.Treatment.Before.Recurrence==1)])
rec_HighSF_patients<-readRDS(paste0(datapath,"rec_HighSF_GBMs.rds"))

# Load junction reads and splicing event lists for four types of events (SE, MXE, A5SS and A3SS)
event_list_sum<-readRDS(paste0(datapath,"CELLO2_splicing_events.rds"))
junc_count_sum<-readRDS(paste0(datapath,"CELLO2_junction_count.rds"))

# Differential splicing analysis for skipped exon (SE) events:
data<-junc_count_sum$SE[which(junc_count_sum$SE$Sample_ID%in%IDHwt_tmz_sample),]
data<-data[which(data$sum_count_R>=10),]
event_list_se_R<-event_list_sum$SE
diff_se_R<-data.frame(Event_ID=event_list_se_R$Event_ID,Diff_PSI=NA,Pval_ttest=NA)
for (i in (1:nrow(event_list_se_R))) {
  data0<-data[which(data$Event_ID==event_list_se_R$Event_ID[i]),]
  H<-data0$psi_R[which(data0$Sample_ID%in%rec_HighSF_patients)]
  L<-data0$psi_R[-which(data0$Sample_ID%in%rec_HighSF_patients)]
  tryCatch(expr = { 
    res<-t.test(H,L)
    diff_se_R$Diff_PSI[i]<-(mean(H)-mean(L))
    diff_se_R$Pval_ttest[i]<-res$p.value}, error = function(e){message(paste("An error occurred for item", i,":\n"), e)})
  print(i)
}
diff_se_R<-cbind(diff_se_R,event_list_se_R[2:12])
sig_diff_se_R<-diff_se_R[intersect(which(diff_se_R$Pval_ttest<1e-10),which(abs(diff_se_R$Diff_PSI)>0.2)),]
sig_psi_se<-matrix(nrow = nrow(sig_diff_se_R),ncol = 102)
colnames(sig_psi_se)<-IDHwt_tmz_sample
for (i in (1:nrow(sig_psi_se))) {
  data0<-data[which(data$Event_ID==sig_diff_se_R$Event_ID[i]),]
  sig_psi_se[i,data0$Sample_ID]<-data0$psi_R
}

# Differential splicing analysis for mutually exclusive exon (MXE) events:
data<-junc_count_sum$MXE[which(junc_count_sum$MXE$Sample_ID%in%IDHwt_tmz_sample),]
data<-data[which(data$sum_count_R>=15),]
event_list_mxe_R<-event_list_sum$MXE
diff_mxe_R<-data.frame(Event_ID=event_list_mxe_R$Event_ID,Diff_PSI=NA,Pval_ttest=NA)
for (i in (1:nrow(event_list_mxe_R))) {
  data0<-data[which(data$Event_ID==event_list_mxe_R$Event_ID[i]),]
  H<-data0$psi_R[which(data0$Sample_ID%in%rec_HighSF_patients)]
  L<-data0$psi_R[-which(data0$Sample_ID%in%rec_HighSF_patients)]
  tryCatch(expr = { 
    res<-t.test(H,L)
    diff_mxe_R$Diff_PSI[i]<-(mean(H)-mean(L))
    diff_mxe_R$Pval_ttest[i]<-res$p.value}, error = function(e){message(paste("An error occurred for item", i,":\n"), e)})
  print(i)
}
diff_mxe_R<-cbind(diff_mxe_R,event_list_mxe_R[2:13])
sig_diff_mxe_R<-diff_mxe_R[intersect(which(diff_mxe_R$Pval_ttest<1e-10),which(abs(diff_mxe_R$Diff_PSI)>0.2)),]
sig_psi_mxe<-matrix(nrow = nrow(sig_diff_mxe_R),ncol = 102)
colnames(sig_psi_mxe)<-IDHwt_tmz_sample
for (i in (1:nrow(sig_psi_mxe))) {
  data0<-data[which(data$Event_ID==sig_diff_mxe_R$Event_ID[i]),]
  sig_psi_mxe[i,data0$Sample_ID]<-data0$psi_R
}

# Differential splicing analysis for alternative 3’ splice site (A3SS) and alternative 5’ splice site (A5SS) events:
event_list_A53SS_R<-event_list_sum$A53SS
psi_A53SS_R<-junc_count_sum$A53SS
diff_A53SS_R<-data.frame(Event_ID=event_list_A53SS_R$Event_ID,Diff_PSI=NA,Pval_ttest=NA)
diff_A53SS_R$Pval_ttest<-apply(psi_A53SS_R, 1, function(t) t.test(as.numeric(t[rec_HighSF_patients]),as.numeric(t[unsig_sample]))$p.value)
diff_A53SS_R$Diff_PSI<-apply(psi_A53SS_R, 1, function(t) mean(as.numeric(t[rec_HighSF_patients]))-mean(as.numeric(t[unsig_sample])))
diff_A53SS_R<-cbind(diff_A53SS_R,event_list_A53SS_R)
sig_diff_A53SS_R<-diff_A53SS_R[intersect(which(diff_A53SS_R$Pval_ttest<1e-10),which(abs(diff_A53SS_R$Diff_PSI)>0.2)),]
sig_diff_A53SS_R$gene<-gsub(";","",sig_diff_A53SS_R$gene)
diff_A53SS_R$geneSymbol<-gsub(";","",diff_A53SS_R$gene)
sig_psi_A53SS_R<-psi_A53SS_R[sig_diff_A53SS_R$Event_ID,]

# Plot Fig. 1d, heatmap of PSI values of differential splicing events
col_annotation<-data.frame(sample=IDHwt_tmz_sample,condition="Recurrent")
row.names(col_annotation)<-col_annotation$sample
col_annotation$group<-ifelse(col_annotation$sample%in%rec_HighSF_patients,"HighSF","LowSF")
colors<-list(group=c(HighSF="#EE9235",LowSF="#4A68DA"),
             condition=c(Initial="#40E0D0",Recurrent="#8A00FF"))
pheatmap(rbind(sig_psi_se,sig_psi_mxe,sig_psi_A53SS_R),annotation_col = col_annotation[,2:3],show_colnames = F,show_rownames = F,
         cutree_cols = 2,color = colorRampPalette(c("navy", "white", "firebrick3"))(50), annotation_colors = colors)

data<-rbind(diff_se_R[,c("Event_ID","Diff_PSI","Pval_ttest","geneSymbol")],diff_mxe_R[,c("Event_ID","Diff_PSI","Pval_ttest","geneSymbol")])
data<-rbind(data,diff_A53SS_R[,c("Event_ID","Diff_PSI","Pval_ttest","geneSymbol")])
# Plot supplementary Fig. 1d, volcano plot of differential splicing results
library(EnhancedVolcano)
EnhancedVolcano(data,lab = data$geneSymbol,
                #drawConnectors = T,
                labSize = 3.5,
                #labFace = 'bold',
                x="Diff_PSI",y="Pval_ttest",
                pCutoff = 1e-10,
                FCcutoff = 0.2,
                pointSize = 1,
                xlim = c(-0.8,0.8),
                ylim = c(0,40),
                colAlpha = 1,
                cutoffLineCol = '#C6C6C8',
                col = c("#C6C6C8","#C6C6C8","#C6C6C8","#C795A7"),
                gridlines.major = F,
                gridlines.minor = F,
                border = 'full',
                borderWidth = 0.6,
                borderColour = 'black')

# Plot Fig. 1e, barplot of enrichment analysis results for differential spliced genes
gene<-unique(data$geneSymbol[intersect(which(data$Pval_ttest<1e-10),which(abs(data$Diff_PSI)>0.2))])
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
go_enrich_up<-enrichGO(gene = gene,
                       universe =unique(data$geneSymbol),#unique(nicheDE_result_sum$geneName)
                       OrgDb = org.Hs.eg.db, 
                       keyType = 'SYMBOL',
                       readable = F,
                       ont = "BP",
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.1)
categorys <- go_enrich_up@result$Description[c(1:3,5,6,8,12,13,19,26)]
barplot(go_enrich_up, showCategory=categorys, label_format = 60)+theme_classic()
