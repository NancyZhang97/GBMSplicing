####################################################
# Analysis of pseudo bulk-level data gotten from snRNA-seq dataset
# Include codes to plot supplementary Figure 3.
# Necessary data can be found in ".../data/" directory
####################################################

datapath="/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/Manuscript/GBMSplicing_codes/data/"
# the raw sequencing matrices of snRNA-seq data for gliomas are gotten from the article: https://pmc.ncbi.nlm.nih.gov/articles/PMC9767870/
sn_pathway<-"pathway of directory that contains count matrices of gliomas"
# Load clinical information for gliomas in snRNA-seq dataset
meta_sn<-read.csv(paste0(datapath,"snRNA-seq.GBM.metadata.csv"))

pseudo_mtx<-matrix(data = NA,nrow = nrow(pbmc),ncol = nrow(meta_sn))
row.names(pseudo_mtx) <- row.names(pbmc)
colnames(pseudo_mtx) <- meta_sn$Prefix
for (i in (1:nrow(meta_sn))) {
  data_dir <- paste0(sn_pathway,'/GSE174554_RAW/',meta_sn$Prefix[i],"/")
  pbmc <- Read10X(data_dir,gene.column=1)
  pseudo_mtx[,i] <- as.numeric(apply(pbmc, 1, sum))
}

pseudo_mtx_trans<-pseudo_mtx[which(rowSums(pseudo_mtx)>10),]
pseudo_mtx_trans<-log2(pseudo_mtx_trans+1)
boxplot(pseudo_mtx_trans)
library(preprocessCore)
pseudo_mtx_trans<-normalize.quantiles(pseudo_mtx_trans,keep.names = T)

geneInfo<-"dataframe that contains gene names with length annotation"
# compute TPM
# find gene length normalized values 
rpk <- apply(pseudo_mtx_trans[intersect(row.names(pseudo_mtx_trans),row.names(geneInfo)),], 2, 
             function(x) x/(geneInfo[intersect(row.names(pseudo_mtx_trans),row.names(geneInfo)),"Length"]/1000))
# normalize by the sample size using rpk values
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
tpm_trans <- tpm[which(rowSums(tpm)>10),]
dim(tpm_trans)
tpm_trans<-log2(tpm_trans+1)
boxplot(tpm_trans)
tpm_trans<-normalize.quantiles(tpm_trans,keep.names = T)

# Directly load TPM matrix for pseudo bulk-level snRNA-seq data
tpm_trans<-readRDS(paste0(datapath,"snRNA_pseudo_bulk_tpm.rds"))
zscore_tpm_trans<-tpm_trans
zscore_tpm_trans<-t(apply(tpm_trans, 1, function(x) (x-mean(x))/sd(x)))

# Plot supplementary Fig. 3b, sum of z-scores for 9 DE SFs across IDHwt GBMs of snRNA-seq dataset
# Load 9 DE SFs' labels
sig_diff_sf<-readRDS(paste0(datapath,"DE_SFs.rds"))
zscore_sf_sum<-as.data.frame(t(zscore_tpm_trans[sig_diff_sf,]))
zscore_sf_sum$sum<-apply(zscore_sf_sum[,1:9],1,sum)
zscore_sf_sum$Prefix<-row.names(zscore_sf_sum)
zscore_sf_sum<-merge(zscore_sf_sum,meta_sn,by.x="Prefix",by.y = "Type")
library(mixtools)
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}
set.seed(5)
mixmdl <- normalmixEM(zscore_sf_sum$sum, k = 2)
ggplot(zscore_sf_sum, aes(x = sum,fill=Stage)) + #here condition means "Initial" or "Recurrent"
  geom_histogram(bins=50, colour = "black",binwidth = 0.7) +
  scale_fill_manual(values = c("#E0FFFF","#E6E6FA")) + #c("#F0FFFF","#E6E6FA")
  mapply(function(mean, sd, lambda, n, binwidth) {
    stat_function(
      fun = function(x) {
        (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
      },colour = "#3069F6",lwd = 0.6
    )
  },
  mean = mixmdl[["mu"]], 
  sd = mixmdl[["sigma"]],
  lambda = mixmdl[["lambda"]],
  n = length(zscore_sf_sum$sum),
  binwidth = 0.7
  )+
  mapply(function(mean, sd, lambda, n, binwidth) {
    stat_function(
      fun = function(x) {
        (dnorm(x, mean = mean, sd = sd)) * n * binwidth * lambda
      },colour = "#EE9238",lwd = 0.6
    )
  },
  mean = mixmdl[["mu"]][2], 
  sd = mixmdl[["sigma"]][2],
  lambda = mixmdl[["lambda"]][2],
  n = length(zscore_sf_sum$sum),
  binwidth = 0.7
  )+
  ylab("Count") + xlab("Sum of z-scores for 9 DE SF genes") + labs(title = "Fitted k=2 Guassian Mixture model for 58 samples (28 pairs)", subtitle = " ")+
  theme_bw()

# Extract patients with HighSF recurrent GBMs
rec_HighSF_sn_patients<-zscore_sf_sum$Pair.[which(zscore_sf_sum$sum>13)]

# Load marker genes of four GBM subtypes
genelist<-readRDS(paste0(datapath,"genelist_fourGBMsubtype.rds"))
# Classify IDHwt GBMs into four subtypes
library(GSVA)
gbm_4type<-gsva(zscore_tpm_trans, genelist, mx.diff=FALSE,method = "ssgsea")
gbm_zscore_grp4<-apply(gbm_4type, 2, function(x) which(x==max(x)))
grp4_ssgsea<-ifelse(gbm_zscore_grp4==1,"Classical",ifelse(gbm_zscore_grp4==2,"Mesenchymal",ifelse(gbm_zscore_grp4==3,"Neural","Proneural")))
grp4_ssgsea<-as.data.frame(grp4_ssgsea)
grp4_ssgsea$Type<-row.names(grp4_ssgsea)
grp4_ssgsea<-merge(meta_sn[,c("Pair.","Prefix","Stage","Type")],grp4_ssgsea,by="Type")
library(reshape2)
grp4_ssgsea_short<-dcast(grp4_ssgsea,Pair.~Stage,value.var = "grp4_ssgsea")

# Group LowSF patients into LowSF NEU+ and LowSF NEU-
rec_LowSF_NEU_sn_patients<-grp4_ssgsea_short$Pair.[which(grp4_ssgsea_short$Recurrent=="Neural")]
rec_LowSF_NEU_sn_patients<-setdiff(rec_LowSF_NEU_sn_patients,rec_HighSF_sn_patients)

# Plot supplementary Fig. 3c, changes of subtypes during tumor recurrence for GBMs in snRNA-seq dataset
library(ggalluvial)
grp4_ssgsea_short$group<-"LowSF_noNEU"
grp4_ssgsea_short$group[which(grp4_ssgsea_short$Pair.%in%rec_LowSF_NEU_sn_patients)]<-"LowSF_NEU"
grp4_ssgsea_short$group[which(grp4_ssgsea_short$Pair.%in%rec_HighSF_sn_patients)]<-"HighSF"
plot_table<-as.data.frame(table(grp4_ssgsea_short$Primary,grp4_ssgsea_short$Recurrent,grp4_ssgsea_short$group))
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

# Plot supplementary Fig. 3d, changes of enrichment scores of Neural and Mesenchymal subtypes across three subgroups of GBM patients
data<-grp4_ssgsea
data$group<-NA
data$group[which(data$Pair.%in%rec_LowSF_NEU_sn_patients)]<-"LowSF_NEU"
data$group[which(data$Pair.%in%rec_HighSF_sn_patients)]<-"HighSF"
data$group[which(is.na(data$group))]<-"LowSF_noNEU"
gbm_4type<-as.data.frame(t(gbm_4type))
gbm_4type<-gbm_4type[data$Type,]
data<-cbind(data,gbm_4type)
ggpaired(data,x="Stage",y="VERHAAK_GLIOBLASTOMA_MESENCHYMAL",facet.by = "group",fill = "Stage",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T)+ theme_classic()+ ylab("Enrichment score of MES markers")
ggpaired(data,x="Stage",y="VERHAAK_GLIOBLASTOMA_NEURAL",facet.by = "group",fill = "Stage",palette = c("#51B1B2","#8A00FF"))+
  stat_compare_means(label = "p.format",method = "t.test",paired = T)+ theme_classic()+ ylab("Enrichment score of NEU markers")

# Save IDs' for HighSF and LowSF NEU+ patients
save(rec_HighSF_sn_patients,rec_LowSF_NEU_sn_patients,file = paste0(datapath,"HighSF_LowSF_snRNA_ID.RData"))
