####################################################
# Comparisons of genetic variants across HighSF, LowSF NEU+ and LowSF NEU- groups.
# Include codes to plot Figure 2b and supplementary Figure 2d.
# Necessary data can be found in ".../data/" directory
####################################################

datapath="/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/Manuscript/GBMSplicing_codes/data/"
# Load patient's information, HighSF GBMs' IDs and LowSF NEU+ GBMs' IDs
clinical<-read.csv(paste0(datapath,"CELLO2_patient_clinical_info.csv"))
IDHwt_tmz_sample<-intersect(clinical$ID_in_Cello2[which(clinical$IDH.mutation=="IDHwt")],clinical$ID_in_Cello2[which(clinical$TMZ.Treatment.Before.Recurrence==1)])
rec_HighSF_patients<-readRDS(paste0(datapath,"rec_HighSF_GBMs.rds"))
rec_LowSF_NEU_patients<-readRDS(paste0(datapath,"rec_LowSF_NEU_GBMs.rds"))

# Load genetic variants information for CELLO2 gliomas
mutation_info<-read.csv("/Users/xiaomeng/Documents/Lab_Material/GBM_splicing/data/mutation_info_183samples.csv",row.names = 1)
mutation_info<-mutation_info[which(mutation_info$ID_in_Cello2%in%IDHwt_tmz_sample),]
mutation_info$group<-"LowSF_noNEU"
mutation_info$group[which(mutation_info$ID_in_Cello2%in%rec_HighSF_patients)]<-"HighSF"
mutation_info$group[which(mutation_info$ID_in_Cello2%in%rec_LowSF_NEU_patients)]<-"LowSF_NEU"

# Compare genetic variants for primary tumors of GBM patients
data<-mutation_info[which(mutation_info$condition=="Initial"),]
x<-colnames(data)[23:64]
mut_test_I<-data.frame(event=x,pval=NA,OR=NA,pval_NEU=NA,OR_NEU=NA,pval_HighSF=NA,OR_HighSF=NA)
for (i in (1:length(x))) {
  tryCatch(expr = {
    y<-fisher.test(table(data$group,data[,x[i]]))
    mut_test_I$pval[i]<-y$p.value
    mut_test_I$OR[i]<-y$estimate
  }, error = function(e){message(paste("An error occurred for item", i,":\n"), e)})
  tryCatch(expr = {
    data0<-data[-which(data$group=="LowSF_noNEU"),]
    y<-fisher.test(table(data0$group,data0[,x[i]]))
    mut_test_I$pval_HighSF[i]<-y$p.value
    mut_test_I$OR_HighSF[i]<-y$estimate
  }, error = function(e){message(paste("An error occurred for item", i,":\n"), e)})
  tryCatch(expr = {
    data0<-data[-which(data$group=="LowSF_NEU"),]
    y<-fisher.test(table(data0$group,data0[,x[i]]))
    mut_test_I$pval_NEU[i]<-y$p.value
    mut_test_I$OR_NEU[i]<-y$estimate
  }, error = function(e){message(paste("An error occurred for item", i,":\n"), e)})
}

# Compare genetic variants for recurrent tumors of GBM patients
data<-mutation_info[which(mutation_info$condition=="Recurrent"),]
x<-colnames(data)[23:64]
mut_test_R<-data.frame(event=x,pval=NA,OR=NA,pval_NEU=NA,OR_NEU=NA,pval_HighSF=NA,OR_HighSF=NA)
for (i in (1:length(x))) {
  tryCatch(expr = {
    y<-fisher.test(table(data$group,data[,x[i]]))
    mut_test_R$pval[i]<-y$p.value
    mut_test_R$OR[i]<-y$estimate
  }, error = function(e){message(paste("An error occurred for item", i,":\n"), e)})
  tryCatch(expr = {
    data0<-data[-which(data$group=="LowSF_noNEU"),]
    y<-fisher.test(table(data0$group,data0[,x[i]]))
    mut_test_R$pval_HighSF[i]<-y$p.value
    mut_test_R$OR_HighSF[i]<-y$estimate
  }, error = function(e){message(paste("An error occurred for item", i,":\n"), e)})
  tryCatch(expr = {
    data0<-data[-which(data$group=="LowSF_NEU"),]
    y<-fisher.test(table(data0$group,data0[,x[i]]))
    mut_test_R$pval_NEU[i]<-y$p.value
    mut_test_R$OR_NEU[i]<-y$estimate
  }, error = function(e){message(paste("An error occurred for item", i,":\n"), e)})
}

# Plot Fig. 2b, pie charts for significantly different genetic variants across three groups
data<-mutation_info[which(mutation_info$condition=="Initial"),]
data0<-as.data.frame(proportions(table(data$PTENdel,data$group),margin = 2))
colnames(data0)<-c("Event","group","Freq")
library(ggplot2)
ggplot(data=data0, aes(x=" ", y=Freq, group=group, fill=Event)) + 
  geom_bar(width = 1, stat = "identity") + 
  #scale_color_manual(values = c(`TRUE` = "#AD5BCD", `FALSE` = "#3DCBCF")) +
  scale_fill_manual(values = c(`0` = "#D3D3D3", `1` = "#00CED1")) +
  coord_polar("y", start=0) +  
  facet_grid(.~ group) +theme_void() 

data<-mutation_info[which(mutation_info$condition=="Recurrent"),]
data0<-as.data.frame(proportions(table(data$CDKN2Adel,data$group),margin = 2))
colnames(data0)<-c("Event","group","Freq")
ggplot(data=data0, aes(x=" ", y=Freq, group=group, fill=Event)) + 
  geom_bar(width = 1, stat = "identity") + 
  #scale_color_manual(values = c(`TRUE` = "#AD5BCD", `FALSE` = "#3DCBCF")) +
  scale_fill_manual(values = c(`0` = "#D3D3D3", `1` = "#BA55D3")) +
  coord_polar("y", start=0) +  
  facet_grid(.~ group) +theme_void() 

# Plot supplementary Fig. 2d
x<-colnames(mutation_info)[25:63]
mut_patient<-unique(mutation_info$ID_in_Cello2)
mut_matrix<-matrix(nrow = length(x),ncol = length(mut_patient))
row.names(mut_matrix)<-x
colnames(mut_matrix)<-mut_patient
for (i in mut_patient) {
  data<-mutation_info[which(mutation_info$ID_in_Cello2==i),]
  data<-data[25:63]
  data<-apply(data, 1, as.numeric)
  data[,2]<-data[,2]*10
  z<-data[,1]+data[,2]
  y<-ifelse(z==0,"Both_neg",ifelse(z==11, "Both_pos",ifelse(z==10,"Recurrent_only","Initial_only")))
  mut_matrix[,i]<-y
}
x<-apply(mut_matrix, 2, function(t) length(which(is.na(t))))
mut_matrix<-mut_matrix[,-which(x>10)]
x<-ifelse(colnames(mut_matrix)%in%rec_HighSF_patients,"HighSF",ifelse(colnames(mut_matrix)%in%rec_LowSF_NEU_patients,"LowSF_NEU","LowSF_noNEU"))
x<-matrix(x)
colnames(x)<-"Group"
mut_matrix<-rbind(mut_matrix,t(x))
mut_matrix<-mut_matrix[,order(mut_matrix["Group",])]
library(ComplexHeatmap)
alter_fun = list(
  background = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#CCCCCC", col = NA)),
  # red rectangles
  Initial_only = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#B22222", col = NA)),
  # blue rectangles
  Recurrent_only = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#0066FFFF", col = NA)),
  # dots
  Both_pos = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#FF7F50", col = NA)),
  Both_neg = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#E6E6FA", col = NA)),
  HighSF = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#EE9238", col = NA)),
  LowSF_NEU = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#3069F6", col = NA)),
  LowSF_noNEU = function(x, y, w, h) 
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "#6C7FAE", col = NA))
)
oncoPrint(mut_matrix, alter_fun = alter_fun,row_order = 1:nrow(mut_matrix),column_order = 1:ncol(mut_matrix),
          col = c(Initial_only = "#B22222", Recurrent_only = "#0066FFFF",Both_pos="#FF7F50",
                  Both_neg ="#E6E6FA", HighSF ="#EE9238", LowSF_NEU="#3069F6",LowSF_noNEU="#6C7FAE"),show_column_names = TRUE)

