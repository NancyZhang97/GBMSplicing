####################################################
# Analysis of correlation between specific splicing patterns and patients' survival
# Include codes to plot Figure 5a and 5b.
# Necessary data can be found in ".../data/" directory
####################################################

datapath="/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/Manuscript/GBMSplicing_codes/data/"
# Load patient's information and HighSF GBMs' IDs.
clinical<-read.csv(paste0(datapath,"CELLO2_patient_clinical_info.csv"))
IDHwt_tmz_sample<-intersect(clinical$ID_in_Cello2[which(clinical$IDH.mutation=="IDHwt")],clinical$ID_in_Cello2[which(clinical$TMZ.Treatment.Before.Recurrence==1)])
clinical<-clinical[which(clinical$ID_in_Cello2%in%IDHwt_tmz_sample),]
clinical$time<-as.numeric(clinical$Overall.Survival..OS..months.)
clinical$status<-as.numeric(clinical$Censor.of.OS..0.alive.1.dead.)

# Load junction reads and splicing event lists for four types of events (SE, MXE, A5SS and A3SS)
event_list_sum<-readRDS(paste0(datapath,"CELLO2_splicing_events.rds"))
junc_count_sum<-readRDS(paste0(datapath,"CELLO2_junction_count.rds"))

library(survival)
library(survminer)
# Grouped CELLO2 GBM patients by PSI values of ITGA6 and compare patients' survival between two groups
data<-junc_count_sum$SE[which(junc_count_sum$SE$Sample_ID%in%IDHwt_tmz_sample),]
data<-data[which(data$Event_ID=="154929"),]
data<-merge(clinical,data,by.x="ID_in_Cello2",by.y="Sample_ID")
data$psi<-data$psi_R
data<-data[-which(data$sum_count_R<10),]
quantile(data$psi)
cut_num<-seq(0.18,0.94,by=0.001)
best_p=1
best_cut=1
for (i in (1:length(cut_num))) {
  data$cut_psi<-NA
  data$cut_psi[which(data$psi>cut_num[i])]<-"H"
  data$cut_psi[which(data$psi<=cut_num[i])]<-"L"
  fit1 <- survfit(Surv(time, status) ~ cut_psi,data = data)
  p_val<-surv_pvalue(fit1,data = data)$pval
  if (p_val<best_p){
    if (min(table(data$cut_psi))>10){
      best_p<-p_val
      best_cut<-cut_num[i]
    }
  }
}
# Plot the right part of Fig. 5d, survival curves of ITGA6 spliced and unspliced groups of CELLO2
data$cut_psi<-NA
data$cut_psi[which(data$psi>best_cut)]<-"H"
data$cut_psi[which(data$psi<=best_cut)]<-"L"
fit1 <- survfit(Surv(time, status) ~ cut_psi,data = data)
ggsurvplot(fit1, data = data, pval = TRUE, palette = "Set1",pval.coord=c(90,0.75),legend="right")

# Load CD47 PSI values and junction counts for exon9 and exon10 of CELLO2.
CD47_twoExons<-read.csv(paste0(datapath,"CELLO2_CD47_exon9&10_junc_count.csv"))
# Grouped CELLO2 GBM patients by PSI values of ITGA6 and compare patients' survival between two groups
data<-CD47_twoExons[which(CD47_twoExons$ID_in_Cello2%in%IDHwt_tmz_sample),]
data<-merge(clinical,data,by="ID_in_Cello2")
data$psi<-data$R_psi
quantile(data$psi)
cut_num<-seq(0.01,0.82,by=0.001)
best_p=1
best_cut=1
for (i in (1:length(cut_num))) {
  data$cut_psi<-NA
  data$cut_psi[which(data$psi>cut_num[i])]<-"H"
  data$cut_psi[which(data$psi<=cut_num[i])]<-"L"
  fit1 <- survfit(Surv(time, status) ~ cut_psi,data = data)
  p_val<-surv_pvalue(fit1,data = data)$pval
  if (p_val<best_p){
    if (min(table(data$cut_psi))>10){
      best_p<-p_val
      best_cut<-cut_num[i]
    }
  }
}
# Plot the right part of Fig. 5c, survival curves of CD47 spliced and unspliced groups of CELLO2
data$cut_psi<-NA
data$cut_psi[which(data$psi>best_cut)]<-"H"
data$cut_psi[which(data$psi<=best_cut)]<-"L"
fit1 <- survfit(Surv(time, status) ~ cut_psi,data = data)
ggsurvplot(fit1, data = data, pval = TRUE, palette = "Set1",pval.coord=c(90,0.75),legend="right")

# Get clinical data of TCGA-GBM cohort from: https://pmc.ncbi.nlm.nih.gov/articles/PMC3910500/
tcga_clinical<-read.csv(paste0(datapath,"TCGA_clinical_NIHMS530933.csv"))
tcga_clinical$TMZ_treat<-ifelse(tcga_clinical$Therapy.Class%in%c("Alkylating Chemo","Nonstandard Radiation, TMZ Chemo",
                                                                 "Standard Radiation, Alkylating Chemo","Standard Radiation, TMZ Chemo",
                                                                 "TMZ Chemoradiation, TMZ Chemo","TMZ Chemo"),TRUE,FALSE)

# Load CD47's PSI values and junction counts for exon9 and exon10 of TCGA-GBM.
tcga_CD47_twoExon<-read.csv(paste0(datapath,"TCGA_CD47_exon9&10_junc_count.csv"))
data<-merge(tcga_CD47_twoExon,tcga_clinical,by.x="ID",by.y="Case.ID")
data<-data[which(data$IDH1..status=="WT"),]
data<-data[which(data$TMZ_treat==TRUE),]
data$time<-as.numeric(data$OS..days.)
data$status<-ifelse(data$Vital.Status=="LIVING",0,1)

quantile(data$psi)
cut_num<-seq(0.01,0.64,by=0.001)
best_p2=1
best_cut2=1
for (i in (1:length(cut_num))) {
  data$cut_psi<-NA
  data$cut_psi[which(data$psi>cut_num[i])]<-"H"
  data$cut_psi[which(data$psi<=cut_num[i])]<-"L"
  fit1 <- survfit(Surv(time, status) ~ cut_psi,data = data)
  p_val<-surv_pvalue(fit1,data = data)$pval
  if (p_val<best_p2){
    if (min(table(data$cut_psi))>10){
      best_p2<-p_val
      best_cut2<-cut_num[i]
    }
  }
}
# Plot the right part of Fig. 5c, survival curves of CD47 spliced and unspliced groups of TCGA-GBM
data$cut_psi<-NA
data$cut_psi[which(data$psi>best_cut2)]<-"H"
data$cut_psi[which(data$psi<=best_cut2)]<-"L"
fit1 <- survfit(Surv(time, status) ~ cut_psi,data = data)
ggsurvplot(fit1, data = data, pval = TRUE, palette = "Set1",pval.coord=c(1200,0.75),legend="right")

# Load ITGA6's PSI values and junction counts for exon25 of TCGA-GBM.
tcga_ITGA6_exon25<-read.csv(paste0(datapath,"TCGA_ITGA6_exon25_junc_count.csv"))
data<-merge(tcga_ITGA6_exon25,tcga_clinical,by.x="ID",by.y="Case.ID")
data$sum_count<-data$count_StoExon+data$count_StoEnd+data$count_ExontoEnd
data<-data0[which(data0$sum_count>20),]
data<-data[which(data$IDH1..status=="WT"),]
data<-data[which(data$TMZ_treat==TRUE),]
data$time<-as.numeric(data$OS..days.)
data$status<-ifelse(data$Vital.Status=="LIVING",0,1)

quantile(data$psi)
cut_num<-seq(0.27,0.85,by=0.001)
best_p2=1
best_cut2=1
for (i in (1:length(cut_num))) {
  data$cut_psi<-NA
  data$cut_psi[which(data$psi>cut_num[i])]<-"H"
  data$cut_psi[which(data$psi<=cut_num[i])]<-"L"
  fit1 <- survfit(Surv(time, status) ~ cut_psi,data = data)
  p_val<-surv_pvalue(fit1,data = data)$pval
  if (p_val<best_p2){
    if (min(table(data$cut_psi))>10){
      best_p2<-p_val
      best_cut2<-cut_num[i]
    }
  }
}
# Plot the right part of Fig. 5c, survival curves of CD47 spliced and unspliced groups of TCGA-GBM
data$cut_psi<-NA
data$cut_psi[which(data$psi>best_cut2)]<-"H"
data$cut_psi[which(data$psi<=best_cut2)]<-"L"
fit1 <- survfit(Surv(time, status) ~ cut_psi,data = data)
ggsurvplot(fit1, data = data, pval = TRUE, palette = "Set1",pval.coord=c(1200,0.75),legend="right")
