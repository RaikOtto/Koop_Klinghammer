library("ggplot2")
library("stringr")
library("grid")
library("umap")
library("dplyr")
source("~/Koop_Klinghammer/Scripts/Visualization_colors.R")

expr_raw = read.table(
  "~/Koop_Klinghammer/Data/S103.tsv",
  sep ="\t",
  stringsAsFactors = F,
  header = T,
  row.names = 1
)
colnames(expr_raw ) = str_replace_all(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
dim(expr_raw)

### Prep

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Sample_ID

matcher = match(colnames(expr_raw), meta_info$SampleID, nomatch = 0)
meta_data = meta_info[matcher,]
rownames(meta_data) = meta_data$SampleID

### Figure SM 1 UMAP with labels

## umap

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
custom.config$random_state = 281
custom.config$n_components= 2

cor_mat = cor(expr_raw);pcr = prcomp(t(cor_mat))
umap_result = umap::umap(
  cor_mat,
  colvec = meta_data$Subtype,
  preserve.seed = TRUE,
  config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point(size = .5, aes(  color = as.character(meta_data$Subtype) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$Subtype), level=.5, type ="t", size=1.5)
umap_p = umap_p + scale_color_manual( values = c("black","darkgreen","blue")) ##33ACFF ##FF4C33
umap_p = umap_p + theme(legend.position = "top") + xlab("") + ylab("")
umap_p = umap_p + geom_text(aes(label = meta_data$SampleID, color = meta_data$Subtype),hjust=0, vjust=0)

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_1.svg", width = 10, height = 10)
umap_p 
dev.off()

### Figure 2 PCA with labels

cor_mat = cor(expr_raw);pcr = prcomp(t(cor_mat))
col_vec = as.character(meta_data$Subtype)
col_vec[col_vec == "BA"] = "black"
col_vec[col_vec == "MS"] = "blue"
col_vec[col_vec == "CL"] = "green"

#svg("~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_2.svg", width = 10, height = 10)
ggbiplot::ggbiplot(
    pcr,
    groups = as.character(meta_data$Subtype),
    ellipse = TRUE,
    circle = TRUE,
    var.axes = F,
    labels = meta_data$SampleID
)  + scale_color_manual(name="Clusters", values=c("Black", "darkgreen", "blue"))
dev.off()

#### Figure 3 clinical metadata

#	Budding_10HPF		Zellnestgröße_zentral		Mitosen_10HPF		Kerngröße	Stroma_vitalerTumor		Nekrose		Entzündung  WHO-Grad L1 und Pn1
selectors = c("Stroma_vitalerTumor","Budding_10HPF","Mitosen_10HPF","Nekrose","Zellnestgröße_zentral","Entzündung","L1","Pn1")
selectors[!(selectors %in% colnames(meta_data))]

vis_mat = meta_data[!is.na(meta_data$Entzündung_ROC),selectors]
vis_mat = vis_mat[,selectors]
vis_mat = scale(t(vis_mat))

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_3.svg", width = 10, height = 10)
pheatmap::pheatmap(
  vis_mat,
  annotation_col = meta_data[c("Subtype")],
  annotation_colors = aka3,
  show_colnames = FALSE,
  treeheight_row = 0,
  cluster_rows = FALSE,
  fontsize_col = 7,
  clustering_method = "ward.D2"
)
dev.off()

### Figure 4 clinical metadata correlation matrix

#	Budding_10HPF		Zellnestgröße_zentral		Mitosen_10HPF		Kerngröße	Stroma_vitalerTumor		Nekrose		Entzündung  WHO-Grad L1 und Pn1
selectors = c("Stroma_vitalerTumor","Budding_10HPF","Mitosen_10HPF","Nekrose","Zellnestgröße_zentral","Entzündung","L1","Pn1")
selectors[!(selectors %in% colnames(meta_data))]

vis_mat = meta_data[!is.na(meta_data$Entzündung),selectors]
vis_mat = vis_mat[,selectors]
vis_mat = cor(t(vis_mat))

annotation_mat = meta_data
annotation_mat$OS_ab_ED = log(annotation_mat$OS_ab_ED+1)
annotation_mat$PFS_Monate_ab_Einschluss = log(annotation_mat$PFS_Monate_ab_Einschluss+1)

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_4.svg", width = 10, height = 10)
pheatmap::pheatmap(
  vis_mat,
  annotation_col = annotation_mat[c("Subtype","OS_ab_ED","PFS_Monate_ab_Einschluss")],
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = FALSE,
  treeheight_row = 0,
  fontsize_col = 7,
  clustering_method = "ward.D"
)
dev.off()

#### SM Figure 5, barplots of phenotypes

selectors = c("Subtype","Stroma_vitalerTumor","Budding_10HPF","Mitosen_10HPF","Nekrose","Zellnestgröße_zentral","Entzündung","Keratinisierung")
selectors[!(selectors %in% colnames(meta_data))]
vis_mat_pre_pre = meta_data[!is.na(meta_data$Entzündung),selectors]

vis_mat_pre = reshape2::melt(vis_mat_pre_pre)
colnames(vis_mat_pre) = c("Subtype","Characteristic","Value")
vis_mat_pre = vis_mat_pre %>% dplyr::filter (!is.na(Value))
vis_mat_pre$Value = log(vis_mat_pre$Value+1)

pheno_plot = ggplot(vis_mat_pre, aes( x = Characteristic, y = Value, fill = Subtype) )
pheno_plot = pheno_plot + geom_boxplot()
pheno_plot = pheno_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
pheno_plot = pheno_plot + scale_fill_manual(values = c("black","darkgreen","blue"))
pheno_plot = pheno_plot + ylab("Logarithm of measurement") + theme(legend.position = "top")
pheno_plot = pheno_plot + coord_flip()  + theme(axis.text.x = element_text(angle =360))

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_5.svg", width = 10, height = 10)
pheno_plot
dev.off()

#### SM Figure 6, barplots of categorial phenotypes

selectors = c("Subtype","Budding_10HPF_ROC","Zellnestgröße_ROC","Mitosen_HPF_ROC","Stroma_ROC","Nekrose_ROC","Entzündung_ROC")
selectors[!(selectors %in% colnames(meta_data))]
vis_mat_pre_pre = meta_data[!is.na(meta_data$Entzündung),selectors]
vis_mat_pre = reshape2::melt(vis_mat_pre_pre)
colnames(vis_mat_pre) = c("Subtype","Characteristic","Count")
vis_mat_pre_mean =  vis_mat_pre %>% dplyr::group_by(Subtype,Characteristic) %>% dplyr::summarise(Mean = mean(Count))
vis_mat_pre_sd = vis_mat_pre %>% dplyr::group_by(Subtype,Characteristic) %>% dplyr::summarise(SD = sd(Count))

pheno_plot = ggplot(vis_mat_pre_mean, aes( x = Characteristic, y = Mean, fill = Subtype) )
pheno_plot = pheno_plot + geom_bar(stat="identity", position=position_dodge(), color = "black")
pheno_plot = pheno_plot + geom_errorbar(aes(ymin = Mean, ymax = Mean + vis_mat_pre_sd$SD),  position = "dodge")
pheno_plot = pheno_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
pheno_plot = pheno_plot + scale_fill_manual(values = c("black","darkgreen","blue"))
pheno_plot = pheno_plot + ylab("Mean of counts") + theme(legend.position = "top")
pheno_plot = pheno_plot + theme(axis.text.x = element_text(angle = 45, vjust = 1))

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_6.Phenotype_categorial.svg", width = 10, height = 10)
pheno_plot
dev.off()

####### pathway plots

genes_of_interest_hgnc_t = read.table("~/Koop_Klinghammer/Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1

# EMT Budding 1 HPF threshold

i = 13
genes_of_interest_hgnc_t$V1[i]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]
sad_genes = sad_genes[which(sad_genes %in% rownames(expr_raw) )]

annotation_mat = meta_data[meta_data$Survivalstatistik_neuestes_Sample == 1,]
dim(annotation_mat)

#expr = expr_raw[sad_genes,colnames(expr_raw) %in% annotation_mat$SampleID]
expr = expr_raw[sad_genes,]
cor_mat = cor(expr)

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_6.Alpha-Keratin.svg", width = 10, height = 10)
pheatmap::pheatmap(
  cor_mat,
  #scale(expr),
  annotation_col = meta_data[,c("Subtype","Budding_1_hpf_threshold","Stroma_vitalerTumor_threshold","Entzündung_ROC","L1")],
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = FALSE,
  treeheight_row = 0,
  fontsize_col = 7,
  clustering_method = "complete"
)
dev.off()

# LM22 Entzündung_ROC

meta_data[is.na(meta_data$Entzündung_ROC),"Entzündung_ROC"] = "Unknown"

i = 11
genes_of_interest_hgnc_t$V1[i]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]
sad_genes = sad_genes[which(sad_genes %in% rownames(expr_raw) )]

annotation_mat = meta_data[meta_data$Survivalstatistik_neuestes_Sample == 1,]
dim(annotation_mat)

#expr = expr_raw[sad_genes,colnames(expr_raw) %in% annotation_mat$SampleID]
expr = expr_raw[sad_genes,]
cor_mat = cor(expr)

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_6.Alpha-Keratin.svg", width = 10, height = 10)
pheatmap::pheatmap(
  cor_mat,
  #scale(expr),
  annotation_col = meta_data[,c("Subtype","Budding_1_hpf_threshold","Stroma_vitalerTumor_threshold","Entzündung_ROC","L1")],
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = FALSE,
  treeheight_row = 0,
  fontsize_col = 7,
  clustering_method = "ward.D"
)
dev.off()

## umap

meta_data[is.na(meta_data$L1),"L1"] = "Unknown"
meta_data[is.na(meta_data$Entzündung_ROC),"Entzündung_ROC"] = "Unknown"

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
custom.config$random_state = 281
custom.config$n_components= 2

cor_mat = cor(expr);pcr = prcomp(t(cor_mat))
umap_result = umap::umap(
  cor_mat,
  colvec = meta_data$Budding_1_hpf_threshold,
  preserve.seed = TRUE,
  config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point(size = 4, aes(  color = as.character(meta_data$Subtype) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$Entzündung_ROC), type ="t", size=1.5, level = .5)
umap_p = umap_p + scale_color_manual( values = c("white","black","blue","black","darkgreen","gray")) ##33ACFF ##FF4C33
umap_p = umap_p + theme(legend.position = "top") + xlab("") + ylab("")
#umap_p = umap_p + geom_text(aes(label = meta_data$SampleID, color = meta_data$Subtype),hjust=0, vjust=0)

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_1.svg", width = 10, height = 10)
umap_p 
dev.off()


### Statistical tests

stat_mat = meta_data[!is.na(meta_data$Entzündung),]
stat_mat = stat_mat[stat_mat$Survivalstatistik_neuestes_Sample == "1",]
stat_mat = stat_mat[!is.na(stat_mat$Budding_10HPF_ROC),]

selectors = c("Budding_10HPF_ROC","Budding_1HPF_ROC","Grading_WHO","Zellnestgröße_ROC","Mitosen_HPF_ROC","Mitosen_10HPF_ROC","Stroma_ROC","BestRespnse","Entzündung_ROC")
#selectors = c("Budding_10HPF_ROC","Budding_1HPF_ROC","Zellnestgröße_ROC","Mitosen_HPF_ROC","Mitosen_10HPF_ROC","Stroma_ROC","Entzündung_ROC")

for (selector in selectors){
  
  print(selector)
  
  BA_vec = stat_mat[stat_mat$Subtype == "BA", selector]
  CL_vec = stat_mat[stat_mat$Subtype == "CL", selector]
  MS_vec = stat_mat[stat_mat$Subtype == "MS", selector]
  
  print(binom.test(x = sum(BA_vec> 1), n = length(BA_vec), p =  mean(CL_vec-1))$p.value)
  print(binom.test(x = sum(BA_vec> 1), n = length(BA_vec), p =  mean(MS_vec-1))$p.value)
  print(binom.test(x = sum(MS_vec> 1), n = length(MS_vec), p =  mean(CL_vec-1))$p.value)
  
  if(FALSE){
  chi_mat = table(
    stat_mat[stat_mat$Subtype %in% c("BA","CL"),selector],
    stat_mat[stat_mat$Subtype %in% c("BA","CL"),"Subtype"])
  print(chisq.test(chi_mat)$p.value)
  
  chi_mat = table(
    stat_mat[stat_mat$Subtype %in% c("BA","MS"),selector],
    stat_mat[stat_mat$Subtype %in% c("BA","MS"),"Subtype"])
  print(chisq.test(chi_mat)$p.value)
  
  chi_mat = table(
    stat_mat[stat_mat$Subtype %in% c("CL","MS"),selector],
    stat_mat[stat_mat$Subtype %in% c("CL","MS"),"Subtype"])
  print(chisq.test(chi_mat)$p.value)}
}

for (selector in selectors){
  
  print(selector)
  
}
BA_vec = vis_mat_pre$Keratinisierung[vis_mat_pre$Subtype == "BA"]
CL_vec = vis_mat_pre$Keratinisierung[vis_mat_pre$Subtype == "CL"]
MS_vec = vis_mat_pre$Keratinisierung[vis_mat_pre$Subtype == "MS"]

BA_vec = vis_mat_pre$Nekrose[vis_mat_pre$Subtype == "BA"]
CL_vec = vis_mat_pre$Nekrose[vis_mat_pre$Subtype == "CL"]
MS_vec = vis_mat_pre$Nekrose[vis_mat_pre$Subtype == "MS"]


t.test(BA_vec, CL_vec)
t.test(BA_vec, MS_vec)
anova_1_p_value = TukeyHSD(anova_1)$`as.factor(as.character(meta_grading$Grading))`

###

selectors = c("Budding_1HPF","Budding_10HPF","Zellnestgröße_zentral","Mitosen_HPF","Mitosen_10HPF","Nekrose","Entzündung","Stroma_vitalerTumor")
stat_mat = meta_data[,selectors]
stat_mat = meta_data[!is.na(meta_data$Entzündung),]
stat_mat = stat_mat[stat_mat$Survivalstatistik_neuestes_Sample == "1",]

for (selector in selectors){
  
  print(selector)
  
  BA_vec = stat_mat[stat_mat$Subtype == "BA", selector]
  CL_vec = stat_mat[stat_mat$Subtype == "CL", selector]
  MS_vec = stat_mat[stat_mat$Subtype == "MS", selector]
  
  print(t.test(BA_vec,MS_vec )$p.value)
  print(t.test(BA_vec, MS_vec)$p.value)
  print(t.test(MS_vec, CL_vec)$p.value)
  
}
