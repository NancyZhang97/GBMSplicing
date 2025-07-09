####################################################
# Prepossess snRNA-seq dataset for matched primary and recurrent GBMs. Include integration, annotation, clustering and CNV inference.
# Include codes to plot Figure 3c, 3d, supplementary Figure 4, and supplementary Figure 5.
# Necessary data can be found in ".../data/" directory
####################################################

datapath="/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/Manuscript/GBMSplicing_codes/data/"
# the raw sequencing matrices of snRNA-seq data for gliomas are gotten from the article: https://pmc.ncbi.nlm.nih.gov/articles/PMC9767870/
sn_pathway<-"pathway of directory that contains count matrices of gliomas"
# Load clinical information for gliomas in snRNA-seq dataset
meta_sn<-read.csv(paste0(datapath,"snRNA-seq.GBM.metadata.csv"))

library(Seurat)
# Read count matrix of snRNA-seq data for all GBM samples into one Seurat object
data_dir <- paste0(sn_pathway,'/GSE174554_RAW/',meta_sn$Prefix[1],"/")
sn <- Read10X(data_dir,gene.column=1)
colnames(sn)<-paste0(meta_sn$Type[1],"_",colnames(sn))
sn <- CreateSeuratObject(sn, names.field = 1, names.delim = "_")
sn$orig.ident<-meta_sn$ID[1]
sn$Prefix<-meta_sn$Prefix[1]
sn$condition<-meta_sn$Stage[1]
sn$pair<-meta_sn$Pair.[1]
sn_sum<-sn
for (i in (2:nrow(meta_sn))) {
  data_dir <- paste0(sn_pathway,'/GSE174554_RAW/',meta_sn$Prefix[i],"/")
  sn <- Read10X(data_dir,gene.column=1)
  colnames(sn)<-paste0(meta_sn$Type[i],"_",colnames(sn))
  sn <- CreateSeuratObject(sn, names.field = 1, names.delim = "_")
  sn$orig.ident<-meta_sn$ID[i]
  sn$Prefix<-meta_sn$Prefix[i]
  sn$condition<-meta_sn$Stage[i]
  sn$pair<-meta_sn$Pair.[i]
  sn_sum<-merge(x=sn_sum,y=sn)
}
sn_sum$orig.ident2<-ifelse(sn_sum$condition=="Primary",paste0("I",sn_sum$pair),paste0("R",sn_sum$pair))

sn_sum <- subset(sn_sum, subset = nCount_RNA >= 200 & percent.mt < 2.5)
pct_expressed <- rowSums(sn_sum@assays$RNA@counts>0)
sn_sum <- subset(sn_sum, features = names(pct_expressed[pct_expressed>=50]))
dim(sn_sum)
Idents(sn_sum)<-sn_sum$orig.ident2
# batch correction (using standard parameters in Seurat)
data.list <- SplitObject(sn_sum, split.by = "orig.ident2")
# normalize and identify variable features for each dataset independently
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = data.list)
data.list <- lapply(X = data.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
# Use reciprocal PCA to do dimension reduction for very large dataset
anchors <- FindIntegrationAnchors(object.list = data.list, reduction = "rpca", 
                                  dims = 1:30)
sn_sum <- IntegrateData(anchorset = anchors, dims = 1:30)

# Analysis with batch corrected data
DefaultAssay(sn_sum) <- 'integrated'
sn_sum <- SetAssayData(sn_sum, slot = "counts", assay = 'integrated', new.data = sn_sum@assays$RNA@counts[rownames(sn_sum),])
sn_sum <- FindVariableFeatures(object = sn_sum, selection.method = "vst", nfeatures = 2000)
sn_sum <- ScaleData(sn_sum)
sn_sum <- RunPCA(object = sn_sum, features = VariableFeatures(object = sn_sum), npcs = 50)
sn_sum <- RunUMAP(object = sn_sum, dims = 1:20)

# Start to cluster the cells
sn_sum <- FindNeighbors(object = sn_sum, dim = 1:10)
DefaultAssay(sn_sum) <- 'integrated'
sn_sum <- FindClusters(object = sn_sum, resolution = 2)

# Seurat object "sn_sum" with processed count matrix could be loaded from ".../data/".
sn_sum<-readRDS("/Volumes/Samsung_T5/plot_data/wlsc_data.rds")
# Plot supplementary Fig. 4a, distribution of cells from different samples after integration
DimPlot(object = sn_sum, reduction = 'umap',group.by = "orig.ident")
# Plot Fig. 3c, distribution of cells from primary and recurrent GBMs after integration and dimension reduction
DimPlot(object = sn_sum, reduction = 'umap',group.by = "condition",cols = c("#73DDD0","#7F17F5"))
# Plot supplementary Fig. 4b, distribution of clusters with resolution = 2
library(scCustomize)
DimPlot(object = sn_sum, reduction = 'umap',group.by = "integrated_snn_res.2",label = T,cols = DiscretePalette_scCustomize(54, palette = "varibow", shuffle = FALSE))

# Cell annotation and CNV estimation
# Plot supplementary Fig. 4c, violin plot showing expression level of cell type markers across clusters
Stacked_VlnPlot(seurat_object = sn_sum, features = c("MOG","VWF","SNAP25","PTPRC","CD163","CD96","SLC1A2","PTPRZ1"),x_lab_rotate = TRUE,
                group.by = "integrated_snn_res.2",colors_use = DiscretePalette_scCustomize(54, palette = "varibow", shuffle = FALSE))
sn_sum$celltype<-NA
sn_sum$celltype[which(sn_sum$integrated_snn_res.2%in%c(6,9,20,2,23,18,25))]<-"Oligodendrocyte"
sn_sum$celltype[which(sn_sum$integrated_snn_res.2%in%c(38,27,31,41,44,39,46))]<-"Neuron"
sn_sum$celltype[which(sn_sum$integrated_snn_res.2%in%c(14,19,8,10,45,42))]<-"Myeloid"
sn_sum$celltype[which(sn_sum$integrated_snn_res.2%in%c(43))]<-"Endothelial"
sn_sum$celltype[which(sn_sum$integrated_snn_res.2%in%c(36))]<-"Astrocyte"
sn_sum$celltype[which(sn_sum$integrated_snn_res.2%in%c(35))]<-"Lymphocyte"
sn_sum$celltype[which(is.na(sn_sum$celltype))]<-"Tumor"
# Cluster tumor cells into different subgroups
sn_tumor<-subset(sn_sum,subset=celltype=="Tumor")
DefaultAssay(sn_tumor) <- 'integrated'
sn_tumor <- NormalizeData(sn_tumor, normalization.method = "LogNormalize", scale.factor = 10000)
sn_tumor <- SetAssayData(sn_tumor, slot = "counts", assay = 'integrated', new.data = sn_tumor@assays$RNA@counts[rownames(sn_tumor),])
sn_tumor <- FindVariableFeatures(object = sn_tumor, selection.method = "vst", nfeatures = 2000)
sn_tumor <- ScaleData(sn_tumor)
sn_tumor <- FindNeighbors(object = sn_tumor, dim = 1:20)
sn_tumor <- FindClusters(object = sn_tumor, resolution = 0.1)
sn_tumor$tumor_subtype<-NA
sn_tumor$tumor_subtype[which(sn_tumor$integrated_snn_res.0.1==0)]<-"Tumor.AC"
sn_tumor$tumor_subtype[which(sn_tumor$integrated_snn_res.0.1==1)]<-"Tumor.Stem"
sn_tumor$tumor_subtype[which(sn_tumor$integrated_snn_res.0.1==2)]<-"Tumor.MES"
sn_tumor$tumor_subtype[which(sn_tumor$integrated_snn_res.0.1==3)]<-"Tumor.OPC"
sn_tumor$tumor_subtype[which(sn_tumor$integrated_snn_res.0.1==4)]<-"Tumor.Cycle"
sn_tumor$tumor_subtype[which(sn_tumor$integrated_snn_res.0.1==5)]<-"Tumor.NPC"
sn_tumor$tumor_subtype[which(sn_tumor$integrated_snn_res.0.1==6)]<-"Tumor.Fiber"

# Find marker genes for each tumor subgroup
Idents(sn_tumor)<-sn_tumor$tumor_subtype
DefaultAssay(sn_tumor) <- "RNA"
find_marker_tumor<-FindAllMarkers(sn_tumor,only.pos = T)
marker_selected<-find_marker_tumor[which(find_marker_tumor$p_val_adj<1e-50),]
marker_selected<-marker_selected[which(marker_selected$avg_log2FC>0.5),]
# Further reduce marker genes by the list of tumor markers of Neftel, et al. and genes' correlation with tumor characterstics
# Load the final list of marker genes for 7 tumor subgroups
tumor_marker<-readRDS(paste0(datapath,"snRNA_tumor_markers.rds"))
# Plot supplementary Fig. 5b, marker genes' expression across different tumor subgroups
data<-sn_tumor@assays$RNA@data[unlist(tumor_marker),]
avg_tumor_marker<-matrix(nrow = length(unlist(tumor_marker)),ncol = length(tumor_marker))
colnames(avg_tumor_marker)<-names(tumor_marker)
row.names(avg_tumor_marker)<-row.names(data)
for (i in names(tumor_marker)) {
  avg_tumor_marker[,i]<-rowMeans(data[,which(sn_tumor$tumor_subtype==i)])
}
library(pheatmap)
ann_colors<-list(
  celltype = c(Tumor.Cycle="#C82B1E",Tumor.OPC="#7A634B",Tumor.AC="#AC9C88",Tumor.NPC="#E7AC53",
               Tumor.MES="#4A885C",Tumor.Stem="#8E39C5",Tumor.Fiber="#CD7692"))
anno<-data.frame(celltype=names(tumor_marker))
row.names(anno)<-names(tumor_marker)
pheatmap(avg_tumor_marker,annotation_col = anno,scale = "row",cluster_rows = F,cluster_cols = F,
         color=colorRampPalette(c("navy", "white", "red"))(50),annotation_colors = ann_colors)

# Plot supplementary Fig. 4e, distribution of cell types
DimPlot(sn_sum,reduction = 'umap',group.by = "celltype",cols = 
          c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF","#8491B4FF","#91D1C2FF"))
# Plot Fig. 3d and supplementary Fig. 5a, distribution of cell types and tumor substates
sn_sum$celltype_subtype<-as.character(sn_sum$celltype)
sn_sum$celltype_subtype[which(sn_sum$celltype=="Tumor")]<-sn_tumor$tumor_subtype
DimPlot(sn_sum,reduction = 'umap',group.by = "celltype_subtype",cols = 
          c("#91D1C2FF","#8491B4FF","#F39B7FFF","#3C5488FF","#00A087FF","#4DBBD5FF","#AC9C88","#C82B1E","#4A885C","#CD7692","#E7AC53","#7A634B","#8E39C5"))

# Use infercnv to infer copy number change for tumor cells 
# Downsample seurat object into 50,000 cells
cell.list <- WhichCells(wlsc, downsample = 50000)
sn_sum.downsampled <- sn_sum[, cell.list]
sn_sum.downsampled <- subset(sn_sum.downsampled,subset=celltype%in%c("Tumor","Oligodendrocyte","Neuron","Myeloid"))
# prepare input files
input_expression <- paste0(datapath,"InferCNV_input.expression.txt")
input_annotation <- paste0(datapath,"InferCNV_input.annotation.txt")
write.table(as.matrix(sn_sum.downsampled@assays$RNA@counts), file = input_expression, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(as.data.frame(sn_sum.downsampled$celltype_subtype, row.names = colnames(sn_sum.downsampled)), file = input_annotation, sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
# run inferCNV
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix=input_expression,
  annotations_file=input_annotation,
  delim="\t",
  gene_order_file="hg38.gene_ordering_file.nodup.txt",
  ref_group_names=c("Oligodendrocyte","Neuron","Myeloid"))
infercnv_obj = infercnv::run(
  infercnv_obj, cutoff=0.1,
  out_dir=paste0(datapath,"Result_inferCNV"), analysis_mode='subclusters', tumor_subcluster_pval=0.05, cluster_by_groups=TRUE, 
  plot_steps=FALSE, denoise = TRUE, HMM=FALSE, no_prelim_plot = TRUE, png_res = 300)
rm(infercnv_obj, infercnv_obj)
cnv_score<-as.matrix(infercnv_obj@expr.data)
gene_order<-infercnv_obj@gene_order
cnv_chr_avg<-as.data.frame(matrix(nrow = 23,ncol=ncol(cnv_score)))
row.names(cnv_chr_avg)<-as.character(unique(gene_order$chr))
colnames(cnv_chr_avg)<-colnames(cnv_score)
for (i in (1:ncol(cnv_score))) {
  df<-data.frame(cnv=cnv_score[,i],chr=gene_order$chr)
  df<-df %>%
    group_by(chr) %>%
    summarise_at(vars(cnv), list(name = mean))
  cnv_chr_avg[,i]<-df$name
}
data<-t(cnv_chr_avg)
cnv_chr_avg2<-data
cnv_chr_avg2<-as.data.frame(cnv_chr_avg2)
cnv_chr_avg2<-cnv_chr_avg2[colnames(pbmc),]
sn_sum.downsampled$chr7_cnv<-as.numeric(cnv_chr_avg2$chr7)
sn_sum.downsampled$chr10_cnv<-as.numeric(cnv_chr_avg2$chr10)

# Plot supplementary Fig. 4d, violin plot for estimated CNV scores across clusters
Stacked_VlnPlot(seurat_object = sn_sum.downsampled, features = c("chr7_cnv","chr10_cnv"),x_lab_rotate = TRUE,
                group.by = "integrated_snn_res.2",colors_use = DiscretePalette_scCustomize(54, palette = "varibow", shuffle = FALSE))

