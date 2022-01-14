library("ggplot2")
library("stringr")
library("grid")
library("umap")
library("dplyr")
source("~/Koop_Klinghammer/Scripts/Visualization_colors.R")

expr_raw = read.table(
  "~/Koop_Klinghammer/Data/Data.S104.tsv",
  sep ="\t",
  stringsAsFactors = F,
  header = T
  #row.names = F
)
colnames(expr_raw ) = str_replace_all(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
dim(expr_raw)

### Prep

## Figure 1

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Sample_ID

meta_data = meta_info[colnames(expr_raw),]
#expr_raw = expr_raw[ ,meta_data$Included == "TRUE"  ]
#meta_data
#meta_data = meta_info[colnames(expr_raw),]
#dim(expr_raw)

###
i = 13
genes_of_interest_hgnc_t = read.table("~/Koop_Klinghammer/Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
genes_of_interest_hgnc_t$V1[i]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]

genes_of_interest_hgnc_t[i,1]

sad_genes[which(!(sad_genes %in% rownames(expr_raw)))]
table(sad_genes %in% rownames(expr_raw) )

expr = expr_raw
meta_data = meta_info[colnames(expr),]

### Figure 2 PCA

cor_mat = cor(expr);pcr = prcomp(t(cor_mat))
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
selectors = c("Subtype","OS_Monate","PFS_Monate","Budding_10HPF","Budding_1HPF","Zellnestgröße_zentral","Mitosen_HPF","Mitosen_10HPF","Nekrose","Entzündung","Stroma_vitalerTumor")
selection[!(selection %in% colnames(meta_data))]

vis_mat_pre = meta_data[!is.na(meta_data$Entzündung),selection]
vis_mat = reshape2::melt(vis_mat_pre)
colnames(vis_mat) = c("Subtype","Characteristic","Value")
vis_mat = vis_mat %>% dplyr::filter (!is.na(Value))
vis_mat = vis_mat %>% dplyr::group_by(Subtype, Characteristic) %>% dplyr::summarise( Value = log(Value+1) )

pheno_plot = ggplot(vis_mat, aes( x = Characteristic, y = Value, fill = Subtype) )
pheno_plot = pheno_plot + geom_boxplot(notch = TRUE,outlier.colour = "red", outlier.shape = 1)
pheno_plot = pheno_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
pheno_plot = pheno_plot + scale_fill_manual(values = c("black","darkgreen","blue"))
pheno_plot = pheno_plot + ylab("Log of measurement") + theme(legend.position = "top")
pheno_plot
#pheno_plot + ylim(c(0,4))

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


###

selection = c("Subtype","Grading_WHO","L1","Pn1","Keratinisierung","Kerngröße")
vis_mat_pre = meta_data[!is.na(meta_data$Entzündung),selection]
vis_mat = reshape2::melt(vis_mat_pre)
colnames(vis_mat) = c("Subtype","Characteristic","Value")
vis_mat = vis_mat %>% dplyr::filter (!is.na(Value))

vis_mat = vis_mat %>% dplyr::group_by(Subtype, Characteristic) %>% dplyr::summarise( Value = Value )
vis_mat_sd = vis_mat %>% dplyr::group_by(Subtype, Characteristic) %>% dplyr::summarise( SD = sd(Value) )

pheno_plot = ggplot(vis_mat, aes( x = Characteristic, y = Value, fill = Subtype) )
pheno_plot = pheno_plot + geom_bar(stat="identity", position=position_dodge())
pheno_plot = pheno_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
pheno_plot = pheno_plot + scale_fill_manual(values = c("black","darkgreen","blue"))
pheno_plot = pheno_plot + ylab("Measurement") + theme(legend.position = "top")
pheno_plot

#######

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

#####

selectors = c("Subtype","PFS_Monate","Budding_10HPF","Zellnestgröße_zentral","Mitosen_10HPF","Nekrose","Entzündung","Stroma_vitalerTumor","L1","Pn1","OS_Monate","PFS_Monate")
selectors %in% colnames(meta_data)

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

####

selectors = c("Budding_1HPF","Budding_10HPF","Zellnestgröße_zentral","Mitosen_HPF","Mitosen_10HPF","Nekrose","Entzündung","Stroma_vitalerTumor","L1","Pn1","OS_Monate","PFS_Monate")
vis_mat_pre = meta_data[(meta_data$Survivalstatistik_neuestes_Sample == "1"),selectors]
vis_mat_pre = vis_mat_pre %>% filter(!is.na(OS_Monate))
scale_data = scale(t(vis_mat_pre))

pheatmap::pheatmap(
  #cor(vis_mat_pre),
  scale_data,
  annotation_col = meta_data["Subtype"],
  annotation_colors = aka3,
  show_rownames = TRUE,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = FALSE,
  fontsize_col = 7,
  clustering_method = "average"#,
  #cellheight = 15
)
dev.off()

######

selectors = c("Budding_10HPF","Zellnestgröße_zentral","Mitosen_10HPF","Nekrose","Entzündung","Stroma_vitalerTumor","L1","Pn1","OS_Monate")
vis_mat_pre = meta_data[(meta_data$Survivalstatistik_neuestes_Sample == "1"),selectors]
vis_mat_pre = vis_mat_pre %>% filter(!is.na(OS_Monate))

annotation_mat = meta_data[,c("Subtype","OS_Monate")]
annotation_mat = annotation_mat[!is.na(annotation_mat$OS_Monate),]
annotation_mat$OS_Monate = log(annotation_mat$OS_Monate +1)
annotation_mat = annotation_mat[rownames(vis_mat_pre),]

pheatmap::pheatmap(
  cor(t(vis_mat_pre)),
  annotation_col = annotation_mat,
  annotation_colors = aka3,
  show_rownames = FALSE,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = FALSE,
  fontsize_col = 7,
  clustering_method = "average"#,
  #cellheight = 15
)

## umap
library("umap")

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
#custom.config$random_state = 281
custom.config$n_components= 2

umap_result = umap::umap(
  cor(t(vis_mat_pre)),
  colvec = annotation_mat$Subtype,
  preserve.seed = TRUE,
  config=custom.config
)

survival_vec = survival_vec_ori = annotation_mat$OS_Monate
percentile_thresh_high = quantile(survival_vec_ori,  seq(0,1,by=.01))[81]
survival_vec[survival_vec_ori < percentile_thresh_high] = "Low"
survival_vec[survival_vec_ori >= percentile_thresh_high] = "High"
annotation_mat$survival_vec = survival_vec

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point(size = 4, aes(  color = as.character(annotation_mat$survival_vec) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = annotation_mat$survival_vec), level=.5, type ="t", size=1.5)
#umap_p = umap_p + scale_color_manual( values = c("black","darkgreen","blue")) ##33ACFF ##FF4C33
umap_p = umap_p + theme(legend.position = "none") + xlab("") + ylab("")

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_1.svg", width = 10, height = 10)
umap_p# +geom_text(aes(label=vis_mat$SampleID, color = vis_mat$Subtype),hjust=0, vjust=0)
dev.off()


