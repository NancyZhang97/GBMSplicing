####################################################
# Analysis of spatial transcriptomics data for IDHwt GBM samples
# Include codes to plot Figure 4c-f.
# Necessary data can be found in ".../data/" directory
####################################################

datapath="/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/Manuscript/GBMSplicing_codes/data/"
# Load Seurat object of processed snRNA-seq data
sn_sum<-readRDS("/Volumes/Samsung_T5/plot_data/wlsc_data.rds")

# Load ST count matrix of GBM samples and estimate cell proportions using snRNA-seq data as reference
ST_datapath<-"pathway of directory that contains ST data objects of GBM samples"
# Load clinical information of ST data
meta_ST<-read.csv(paste0(datapath,"CC2022_spatial_clinical.csv"))
meta_ST<-meta_ST[which(meta_ST$Tumor=="IDH-WT"),]
ST_id<-gsub("UKF","",meta_ST$ID)
ST_id<-gsub("#","",ST_id)
ST_id<-unique(ST_id)

library(SPATA2)
library(spacexr)
# Build up reference matrix for cell deconvolution by RCTD
data<-as.matrix(sn_sum@assays$RNA$counts)
sc_lable<-sn_sum$celltype_subtype
names(sc_lable) <- colnames(data)
sc_lable <- as.factor(sc_lable)
reference_sc <- Reference(data, sc_lable)
# Load ST count matrix of GBM samples and do cell deconvolution one by one
prop_estimate_ST<-list()
i=ST_id[1]
spata_obj <- loadSpataObject(paste0(ST_datapath,i,"_T.rds"))
count_data<-getCountMatrix(spata_obj) #extract count matrix from SPATA-object
coordinate<-as.data.frame(getCoordsDf(spata_obj)) #extract coordinate information from SPATA-object
row.names(coordinate)<-coordinate$barcodes
coordinate<-coordinate[,c("x","y")]
spata_pre <- SpatialRNA(coordinate, count_data)
myRCTD <- create.RCTD(spata_pre, reference_sc, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
results <- myRCTD@results
prop_estimate = as.data.frame(normalize_weights(results$weights))
prop_estimate<-prop_estimate[row.names(coordinate),]
prop_estimate<-prop_estimate[row.names(coordinate),]
prop_estimate$segmentation<-spata_obj@fdata[[paste0(i,"_T")]]$segmentation
prop_estimate_ST[[paste0("T",i)]]<-prop_estimate
prop_estimate$sample<-paste0("T",i)
prop_estimate_long<-prop_estimate
for (i in ST_id[2:17]) {
  spata_obj <- loadSpataObject(paste0(ST_datapath,i,"_T.rds"))
  count_data<-getCountMatrix(spata_obj) #extract count matrix from SPATA-object
  coordinate<-as.data.frame(getCoordsDf(spata_obj)) #extract coordinate information from SPATA-object
  row.names(coordinate)<-coordinate$barcodes
  coordinate<-coordinate[,c("x","y")]
  spata_pre <- SpatialRNA(coordinate, count_data)
  myRCTD <- create.RCTD(spata_pre, reference_sc, max_cores = 1)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  results <- myRCTD@results
  prop_estimate = as.data.frame(normalize_weights(results$weights))
  prop_estimate<-prop_estimate[row.names(coordinate),]
  prop_estimate<-prop_estimate[row.names(coordinate),]
  prop_estimate$segmentation<-spata_obj@fdata[[paste0(i,"_T")]]$segmentation
  prop_estimate_ST[[paste0("T",i)]]<-prop_estimate
  prop_estimate$sample<-paste0("T",i)
  prop_estimate_long<-rbind(prop_estimate_long,prop_estimate)
}

# Plot Fig. 4c, histological segmentation for two IDHwt GBMs
library(ggplot2)
data<-prop_estimate_ST[["T260"]]
data$segmentation<-factor(data$segmentation,levels = c("Cellular","Infiltartive","Necrotic_Edge","Vascular_Hyper","none"))
ggscatter(data,x="x",y="y",color = "segmentation",palette = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF"))
data<-prop_estimate_ST[["T304"]]
data$segmentation<-factor(data$segmentation,levels = c("Cellular","Infiltartive","Necrotic_Edge","Vascular_Hyper","none"))
ggscatter(data,x="x",y="y",color = "segmentation",palette = c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#8491B4FF"))

# Plot Fig .4e, estimated cell proportion of neurons and Tumor.NPC cells for two IDHwt GBMs
library(RColorBrewer)
data<-prop_estimate_ST[["T260"]]
ggscatter(data,x="x",y="y",color = "Neuron")&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggscatter(data,x="x",y="y",color = "Tumor.NPC")&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
data<-prop_estimate_ST[["T304"]]
ggscatter(data,x="x",y="y",color = "Neuron")&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
ggscatter(data,x="x",y="y",color = "Tumor.NPC")&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

# Plot Fig. 4d, average expression level of 9 DE SFs for two IDHwt GBMs
# Load 9 DE SFs' labels
sig_diff_sf<-readRDS(paste0(datapath,"DE_SFs.rds"))
spata_obj <- loadSpataObject(paste0(ST_datapath,"260_T.rds"))
data0<-spata_obj@data[["260_T"]]$scaled[sig_diff_sf,]
data<-prop_estimate_ST[["T260"]]
data$exp<-rowMeans(data0)
ggscatter(data,x="x",y="y",color = "exp")&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))
spata_obj <- loadSpataObject(paste0(ST_datapath,"304_T.rds"))
data0<-spata_obj@data[["304_T"]]$scaled[sig_diff_sf,]
data<-prop_estimate_ST[["T304"]]
data$exp<-rowMeans(data0)
ggscatter(data,x="x",y="y",color = "exp")&scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdYlBu")))

# Calculate correlation of cell proportions by each two cell-cell pair across all GBM samples
celltype<-unique(sn_sum$celltype_subtype)
data<-combn(celltype,2)
cell_comb<-c()
for (i in (1:ncol(data))) {
  cell_comb<-c(cell_comb,paste0(data[1,i], " & ", data[2,i]))
}
cor_ST_sum<-data.frame(cell.comb=cell_comb,cor.coef=NA,cor.pval=NA)
for (i in (1:ncol(data))) {
  res<-cor.test(prop_estimate_long[,data[1,i]],prop_estimate_long[,data[2,i]])
  cor_ST_sum$cor.coef[i]<-res$estimate
  cor_ST_sum$cor.pval[i]<-res$p.value
}
cor_ST_sum$avg_prop.x<-NA
cor_ST_sum$avg_prop.y<-NA
for (i in (1:ncol(data))) {
  cor_ST_sum$avg_prop.x[i]<-mean(prop_estimate_long[,data[1,i]],na.rm=T)
  cor_ST_sum$avg_prop.y[i]<-mean(prop_estimate_long[,data[2,i]],na.rm=T)
}
cor_ST_sum$avg_prop.product<-cor_ST_sum$avg_prop.x*cor_ST_sum$avg_prop.y

# Build up simulated ST count matrix to estimate general distribution of cell-cell correlations
Idents(sn_sum)<-sn_sum$celltype_subtype
cell.list <- WhichCells(sn_sum, downsample = 1000)
sn_sum.downsampled <- sn_sum[, cell.list]

spata_obj <- loadSpataObject(paste0(ST_datapath,"304_T.rds"))
prop_estimate<-prop_estimate_ST$T304
prop_estimate<-prop_estimate[-which(is.na(prop_estimate$Neuron)),]
coordinate<-prop_estimate[,c("x","y")]
prop_estimate<-prop_estimate[,celltype]
cell_pseudo<-apply(prop_estimate, 1, function(t) round(t*10))
cell_pseudo<-as.data.frame(t(cell_pseudo))
cell_count<-apply(cell_pseudo, 2, sum)
celltype<-colnames(prop_estimate)
cell_count<-cell_count[celltype]
cell_string<-c()
for (i in celltype) {
  x<-rep(i,cell_count[i])
  cell_string<-c(cell_string,x)
}

random_select_gene<-function(input,cell_labels,count_mtx){
  data<-as.data.frame(table(input))
  labels<-c()
  for (i in (1:nrow(data))) {
    label<-sample(names(cell_labels)[which(cell_labels==as.character(data$input[i]))],data$Freq[i])
    labels<-c(labels,label)
  }
  data0<-count_mtx[,labels]
  data0<-apply(data0, 1, sum)
  return(data0)
}

count<-as.matrix(sn_sum.downsampled@assays$RNA$counts)
common_gene<-intersect(row.names(sn_sum.downsampled),row.names(spata_obj@data[["304_T"]]$counts))
count<-count[common_gene,]
sc_lable<-as.character(sn_sum.downsampled$celltype_subtype)
names(sc_lable) <- colnames(count)
psd_cor_list<-list()
psd_prop_list<-list()
for (k in (1:1000)) {
  group<-split(sample(cell_string,size = length(cell_string)),1:nrow(prop_estimate))
  psd_count<-sapply(group,function(t) random_select_gene(input = t,cell_labels = sc_lable,count_mtx = count))
  colnames(psd_count)<-row.names(prop_estimate)
  spata_pre <- SpatialRNA(coordinate, psd_count)
  myRCTD <- create.RCTD(spata_pre, reference_sc, max_cores = 1)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  results <- myRCTD@results
  prop_estimate0 = as.data.frame(as.matrix(normalize_weights(results$weights)))
  prop_estimate0<-prop_estimate0[row.names(coordinate),]
  data<-combn(celltype,2)
  cell_comb<-c()
  for (i in (1:ncol(data))) {
    cell_comb<-c(cell_comb,paste0(data[1,i], " & ", data[2,i]))
  }
  cor_ST0<-data.frame(cell.comb=cell_comb,cor.coef=NA,cor.pval=NA)
  for (i in (1:ncol(data))) {
    res<-cor.test(prop_estimate0[,data[1,i]],prop_estimate0[,data[2,i]])
    cor_ST0$cor.coef[i]<-res$estimate
    cor_ST0$cor.pval[i]<-res$p.value
  }
  psd_cor_list[[k]]<-cor_ST0
  psd_prop_list[[k]]<-prop_estimate0
  print(k)
}

data<-combn(celltype,2)
for (k in (1:length(psd_cor_list))) {
  prop_estimate0<-psd_prop_list[[k]]
  cor_ST0<-data.frame(cell.comb=cell_comb,cor.coef=NA,cor.pval=NA)
  for (i in (1:ncol(data))) {
    res<-cor.test(prop_estimate0[,data[1,i]],prop_estimate0[,data[2,i]])
    cor_ST0$cor.coef[i]<-res$estimate
    cor_ST0$cor.pval[i]<-res$p.value
  }
  psd_cor_list[[k]]<-cor_ST0
}

comb_cor_mtx<-matrix(nrow=length(cell_comb),ncol = length(psd_cor_list))
row.names(comb_cor_mtx)<-cell_comb
for (i in (1:length(psd_cor_list))) {
  comb_cor_mtx[,i]<-psd_cor_list[[i]]$cor.coef
}

#generate spot matrix with similar cell distribution in the real sample
cell_pseudo_list<-list()
for (i in (1:nrow(cell_pseudo))) {
  x<-cell_pseudo[i,]
  y<-c()
  for (j in celltype) {
    y<-c(y,rep(j,x[j]))
  }
  cell_pseudo_list[[i]]<-y
}
psd_count<-sapply(cell_pseudo_list,function(t) random_select_gene(input = t,cell_labels = sc_lable,count_mtx = count))
colnames(psd_count)<-row.names(prop_estimate)
spata_pre <- SpatialRNA(coordinate, psd_count)
myRCTD <- create.RCTD(spata_pre, reference_sc, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
results <- myRCTD@results
prop_estimate0 = as.data.frame(as.matrix(normalize_weights(results$weights)))
prop_estimate0<-prop_estimate0[row.names(coordinate),]
data<-combn(celltype,2)
cell_comb<-c()
for (i in (1:ncol(data))) {
  cell_comb<-c(cell_comb,paste0(data[1,i], " & ", data[2,i]))
}
cor_ST0<-data.frame(cell.comb=cell_comb,cor.coef=NA,cor.pval=NA)
for (i in (1:ncol(data))) {
  res<-cor.test(prop_estimate0[,data[1,i]],prop_estimate0[,data[2,i]])
  cor_ST0$cor.coef[i]<-res$estimate
  cor_ST0$cor.pval[i]<-res$p.value
}
prop_estimate_real<-prop_estimate0
cor_ST_real<-cor_ST0

library(MASS)
cor_ST_real$psd.pval<-NA
cor_ST_real$log_psd.pval<-NA
for (i in (1:length(cell_comb))) {
  fit <- fitdistr(comb_cor_mtx[i,], "normal")
  para <- fit$estimate
  if (cor_ST_real$cor.coef[i]>0){
    cor_ST_real$psd.pval[i]<-pnorm(cor_ST_real$cor.coef[i],mean = para[1],sd = para[2],lower.tail = F)
    cor_ST_real$log_psd.pval[i]<-(-pnorm(cor_ST_real$cor.coef[i],mean = para[1],sd = para[2],log.p = T,lower.tail = F))
  } else {
    cor_ST_real$psd.pval[i]<-pnorm(cor_ST_real$cor.coef[i],mean = para[1],sd = para[2],lower.tail = T)
    cor_ST_real$log_psd.pval[i]<-(-pnorm(cor_ST_real$cor.coef[i],mean = para[1],sd = para[2],log.p = T,lower.tail = T))
  }
}
cor_ST_sum$psd.pval<-cor_ST_real$psd.pval
cor_ST_sum$log_psd.pval<-cor_ST_real$log_psd.pval
cor_ST_sum$log_cor.pval<-(-log10(cor_ST_sum$cor.pval))

# Plot Fig. 4f, results of cell-cell co-existing test using ST data. 
# Results of cell-cell co-existing analysis could be loaded from ".../data/".
cor_ST_sum<-read.csv(paste0(datapath,"ST_coexist_result.csv"))
x<-cor_ST_sum$cell.comb[intersect(which(abs(cor_ST_sum$cor.coef)>0.15),which(cor_ST_sum$log_psd.pval>20))]
ggscatter(cor_ST_sum,x="log_psd.pval",y="cor.coef",label = "cell.comb",label.rectangle = T,label.select = x,repel = T,color = "cor.coef",size = "log_cor.pval",legend="right",font.label = 10) + 
  gradient_color(c("blue", "light grey", "red"))

