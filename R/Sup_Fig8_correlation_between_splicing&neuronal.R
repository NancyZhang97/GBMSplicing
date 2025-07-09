####################################################
# Analysis of correlation between specific splicing patterns and neuronal transitions in GBMs
# Include codes to plot supplementary Figure 8.
# Necessary data can be found in ".../data/" directory
####################################################

datapath="/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/Manuscript/GBMSplicing_codes/data/"
# Load patient's information and HighSF GBMs' IDs.
clinical<-read.csv(paste0(datapath,"CELLO2_patient_clinical_info.csv"))
IDHwt_tmz_sample<-intersect(clinical$ID_in_Cello2[which(clinical$IDH.mutation=="IDHwt")],clinical$ID_in_Cello2[which(clinical$TMZ.Treatment.Before.Recurrence==1)])
rec_HighSF_patients<-readRDS(paste0(datapath,"rec_HighSF_GBMs.rds"))
rec_LowSF_NEU_patients<-readRDS(paste0(datapath,"rec_LowSF_NEU_GBMs.rds"))

# Load estimated cell proportions of CELLO2 GBM samples by CIBERSORTx
prop_cello2<-read.csv(paste0(datapath,"CELLO2_CIBERSORTx_cell_proportions.csv"))

# Plot supplementary Fig. 8a, changes of cell proportions during tumor recurrence across three GBM subgroups
prop_cello2<-prop_cello2[which(prop_cello2$Patient%in%IDHwt_tmz_sample),]
prop_cello2$group<-"LowSF.NEU-"
prop_cello2$group[which(prop_cello2$Patient%in%rec_HighSF_patients)]<-"HighSF"
prop_cello2$group[which(prop_cello2$Patient%in%rec_LowSF_NEU_patients)]<-"LowSF.NEU+"
prop_cello2<-prop_cello2[order(prop_cello2$Patient),]
prop_cello2<-prop_cello2[order(prop_cello2$condition),]
prop_cello2$group<-factor(prop_cello2$group,levels = c("HighSF","LowSF.NEU+","LowSF.NEU-"))
library(ggpubr)
ggpaired(prop_cello2,x="condition",y="Myeloid",fill = "condition",palette = c("#51B1B2","#8A00FF"),facet.by = "group")+
  stat_compare_means(label = "p.format",method = "t.test")+theme_classic2()

# Plot supplementary Fig. 8b, changes of relative tumor abundance during tumor recurrence across three GBM subgroups
prop_cello2_tumor<-prop_cello2[,c("Tumor.Stem","Tumor.Fiber","Tumor.AC","Tumor.Cycle","Tumor.MES","Tumor.OPC","Tumor.NPC")]
prop_cello2_tumor<-t(apply(prop_cello2_tumor,1,function(t) t/sum(t)))
prop_cello2_tumor<-cbind(prop_cello2_tumor,prop_cello2[,18:20])
ggpaired(prop_cello2,x="condition",y="Tumor.MES",fill = "condition",palette = c("#51B1B2","#8A00FF"),facet.by = "group")+
  stat_compare_means(label = "p.format",method = "t.test")+theme_classic2()

# Load splicing events' junction reads and PSI values for CELLO2 GBM patients
junc_count_sum<-readRDS(paste0(datapath,"CELLO2_junction_count.rds"))

psi_to_matrix<-function(sample,psi_data){
  x<-rep(NA,length(sample))
  names(x)<-sample
  x[psi_data$Sample_ID]<-psi_data$psi_R
  return(x)
}

# Process SE events and do correlation test with neuronal proportions
# Read differential testing results of SE events
library(openxlsx)
diff_se_R<-read.xlsx(paste0(datapath,"diff_splicingEvents_HighSFvsLowSF.xlsx"),sheet = "SE")
data<-junc_count_sum$SE[which(junc_count_sum$SE$Sample_ID%in%IDHwt_tmz_sample),]
data<-data[which(data$sum_count_R>=10),]
psi_se_R<-matrix(nrow = nrow(diff_se_R),ncol = 102)
colnames(psi_se_R)<-IDHwt_tmz_sample
row.names(psi_se_R)<-diff_se_R$Event_ID

data0<-apply(diff_se_R, 1, function(t) psi_to_matrix(sample = IDHwt_tmz_sample, psi_data = data[which(data$Event_ID==t[1]),]))
psi_se_R<-t(data0)
colnames(psi_se_R)<-IDHwt_tmz_sample
cor_se_R<-data.frame(Event_ID=diff_se_R$Event_ID,cor_pval_neuron=NA,cor_pval_NPC=NA)
data0<-prop_cello2[which(prop_cello2$condition=="R"),]
row.names(data0)<-data0$Patient
data0<-data0[colnames(psi_se_R),]
cor_se_R$cor_pval_neuron<-apply(psi_se_R,1,function(t) cor.test(t,data0$Neuron,method = "spearman")$p.value)
data0<-prop_cello2_tumor[which(prop_cello2_tumor$condition=="R"),]
row.names(data0)<-data0$Patient
data0<-data0[colnames(psi_se_R),]
cor_se_R$cor_pval_NPC<-apply(psi_se_R,1,function(t) cor.test(t,data0$Tumor_Neuron.like,method = "spearman")$p.value)
cor_se_R<-cbind(cor_se_R,diff_se_R[,c(2:13)])
cor_se_R$log_p_NPC<-(-log10(cor_se_R$cor_pval_NPC))
cor_se_R$log_p_diff<-(-log10(cor_se_R$Pval_ttest))

# Process MXE events and do correlation test with neuronal proportions
# Read differential testing results of MXE events
diff_mxe_R<-read.xlsx(paste0(datapath,"diff_splicingEvents_HighSFvsLowSF.xlsx"),sheet = "MXE")
data<-junc_count_sum$MXE[which(junc_count_sum$MXE$Sample_ID%in%IDHwt_tmz_sample),]
data<-data[which(data$sum_count_R>=15),]
data0<-apply(diff_mxe_R, 1, function(t) psi_to_matrix(sample = IDHwt_tmz_sample, psi_data = data[which(data$Event_ID==t[1]),]))
psi_mxe_R<-t(data0)
cor_mxe_R<-data.frame(Event_ID=diff_mxe_R$Event_ID,cor_pval_neuron=NA,cor_pval_NPC=NA)
data0<-prop_cello2[which(prop_cello2$condition=="R"),]
row.names(data0)<-data0$Patient
data0<-data0[colnames(psi_se_R),]
cor_mxe_R$cor_pval_neuron<-apply(psi_mxe_R,1,function(t) cor.test(t,data0$Neuron,method = "spearman")$p.value)
data0<-prop_cello2_tumor[which(prop_cello2_tumor$condition=="R"),]
row.names(data0)<-data0$Patient
data0<-data0[colnames(psi_se_R),]
cor_mxe_R$cor_pval_NPC<-apply(psi_mxe_R,1,function(t) cor.test(t,data0$Tumor_Neuron.like,method = "spearman")$p.value)
cor_mxe_R<-cbind(cor_mxe_R,diff_mxe_R[,c(2:15)])

# Process A53SS events and do correlation test with neuronal proportions
# Read differential testing results of A53SS events
diff_A53SS_R<-read.xlsx(paste0(datapath,"diff_splicingEvents_HighSFvsLowSF.xlsx"),sheet = "A53SS")
cor_A53SS_R<-data.frame(Event_ID=diff_A53SS_R$Event_ID,cor_pval_neuron=NA,cor_pval_NPC=NA)
psi_A53SS_R<-junc_count_sum$A53SS
data0<-prop_cello2[which(prop_cello2$condition=="R"),]
row.names(data0)<-data0$Patient
data0<-data0[colnames(psi_se_R),]
cor_A53SS_R$cor_pval_NPC<-apply(psi_A53SS_R,1,function(t) cor.test(as.numeric(t),data0$Tumor_Neuron.like,method = "spearman")$p.value)
data0<-prop_cello2_tumor[which(prop_cello2_tumor$condition=="R"),]
row.names(data0)<-data0$Patient
data0<-data0[colnames(psi_se_R),]
cor_A53SS_R$cor_pval_neuron<-apply(psi_A53SS_R,1,function(t) cor.test(as.numeric(t),data0$Neuron,method = "spearman")$p.value)
cor_A53SS_R<-cbind(cor_A53SS_R,diff_A53SS_R[,c(2:14)])

sum_cor_diff_R<-cor_se_R[,c("Event_ID","Diff_PSI","Pval_ttest","cor_pval_NPC","cor_pval_neuron","geneSymbol")]
sum_cor_diff_R$Type<-"SE"
data<-cor_mxe_R[,c("Event_ID","Diff_PSI","Pval_ttest","cor_pval_NPC","cor_pval_neuron","geneSymbol")]
data$Type<-"MXE"
sum_cor_diff_R<-rbind(sum_cor_diff_R,data)
data<-cor_A53SS_R[,c("Event_ID","Diff_PSI","Pval_ttest","cor_pval_NPC","cor_pval_neuron","geneSymbol")]
data$Type<-"A53SS"
sum_cor_diff_R<-rbind(sum_cor_diff_R,data)
sum_cor_diff_R$sig_diff<-"none"
sum_cor_diff_R$sig_diff[intersect(which(sum_cor_diff_R$Pval_ttest<1e-3),which(abs(sum_cor_diff_R$Diff_PSI)>0.2))]<-"sig"
sum_cor_diff_R$sig_cor<-"none"
sum_cor_diff_R$sig_cor[intersect(which(sum_cor_diff_R$cor_pval_NPC<1e-3),which(sum_cor_diff_R$cor_pval_neuron<1e-5))]<-"sig"

# Plot supplementary Fig. 8c, splicing events' correlation coefficients with neuron and Tumor.NPC cells
sum_cor_diff_R$log_p_NPC<-(-log10(sum_cor_diff_R$cor_pval_NPC))
sum_cor_diff_R$log_p_neuron<-(-log10(sum_cor_diff_R$cor_pval_neuron))
ggscatter(sum_cor_diff_R,x="log_p_NPC",y="log_p_neuron",color = "sig_cor",label = "geneSymbol",label.select=c("CD47","ITGA6","NIN"),size = 0.1,palette = c("#C0C0C0","#DC143C")) 
