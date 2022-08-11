library("ggplot2")
library("stringr")
library("grid")
library("umap")

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

## Figure 1
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Sample_ID

meta_data = meta_info[colnames(expr_raw),]
dim(meta_data)
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

expr = expr_raw#[ rownames(expr_raw) %in% sad_genes,]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

#svg("~/Koop_Klinghammer/Results/23_10_2019/Heatmap_94.svg")
selection = c("Subtype","Grading_WHO","Keratinisierung","Budding_10HPF","Zellnestgröße_zentral","Mitosen_10HPF","Nekrose","Entzündung")
selection = c("Subtype","Grading_WHO")
selection[!(selection %in% colnames(meta_data))]

vis_mat_cor_plot = meta_data[,selection]
#not_na_vec = !is.na(vis_mat_cor_plot[,"Budding_10HPF"])
#vis_mat_cor_plot[not_na_vec,"Mitosen_10HPF"] = log(vis_mat_cor_plot[not_na_vec,"Mitosen_10HPF"]+1)
#vis_mat_cor_plot[not_na_vec,"Budding_10HPF"] = log(vis_mat_cor_plot[not_na_vec,"Budding_10HPF"]+1)
#vis_mat_cor_plot[not_na_vec,"Zellnestgröße_zentral"] = log(vis_mat_cor_plot[not_na_vec,"Zellnestgröße_zentral"]+1)
#vis_mat_cor_plot[not_na_vec,"Nekrose"] = log(vis_mat_cor_plot[not_na_vec,"Nekrose"]+1)
#vis_mat_cor_plot[not_na_vec,"Entzündung"] = log(vis_mat_cor_plot[not_na_vec,"Entzündung"]+1)

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_1.Grading.svg", width = 10, height = 10)
pheatmap::pheatmap(
  cor_mat,
  #vis_mat_cor_plot,
  annotation_col = vis_mat_cor_plot,
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = TRUE,
  fontsize_col = 7,
  clustering_method = "ward.D"
)
dev.off()

aka3 = list(
  Group = c(refractive = "red", sensitive = "darkgreen", intermediate = "orange"),
  Subtype = c(BA = "black", CL = "darkgreen", MS = "blue"),
  OS = c(high = "red", medium = "orange", low = "green"),
  Budding_absentpresent = c("2" = "black", "1" = "white"),
  Zellnestgröße_zentral = c(low = "white", high = "black"),
  Mitosen_10HPF = c(low = "white", high = "black"),
  Nekrose = c(low = "white", high = "black"),
  Entzündung = c(low = "white", high = "black"),
  Budding_10HPF = c(low = "white", high = "black"),
  Keratinisierung = c(low = "white", high = "black"),
  Grading_WHO = c( "1" = "#CBDB34","2" = "#FFC000", "3" = "#EE5124")
)
#dev.off()

####

#Budding_1HPF		Budding_10HPF		Zellnestgröße_zentral		Mitosen_HPF		Mitosen_10HPF		Kerngröße	Stroma_vitalerTumor		Nekrose		Entzündung
selection_vis = c("Subtype","PFS_Monate","Keratinisierung","Budding_10HPF","Zellnestgröße_zentral","Mitosen_HPF","Mitosen_10HPF","Nekrose","Entzündung","Kerngröße","Stroma_vitalerTumor")
vis_mat_pre = meta_data[!is.na(meta_data$Entzündung),selection_vis]
vis_mat = reshape2::melt(vis_mat_pre)
colnames(vis_mat) = c("Subtype","Characteristic","Value")
vis_mat = vis_mat %>% dplyr::group_by(Subtype, Characteristic) %>% dplyr::summarise( Value = log(Value+1) )

pheno_plot = ggplot(vis_mat, aes( x = Characteristic, y = Value, fill = Subtype) )
pheno_plot = pheno_plot + geom_boxplot(notch = TRUE,outlier.colour = "red", outlier.shape = 1)
pheno_plot = pheno_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
pheno_plot = pheno_plot + scale_fill_manual(values = c("black","darkgreen","blue"))
pheno_plot = pheno_plot + ylab("Log of measurement") + theme(legend.position = "top")
pheno_plot

#####

selection_vis = c("Subtype","Grading_WHO")
vis_mat_pre = meta_data[!is.na(meta_data$Entzündung),selection_vis]
vis_mat = reshape2::melt(vis_mat_pre)
colnames(vis_mat) = c("Subtype","Characteristic","Value")
vis_mat = vis_mat %>% dplyr::group_by(Subtype, Characteristic) %>% dplyr::summarise( Value = log(Value+1) )

pheno_plot = ggplot(vis_mat, aes( x = Characteristic, y = Value, fill = Subtype) )
pheno_plot = pheno_plot + geom_boxplot(notch = TRUE,outlier.colour = "red", outlier.shape = 1)
pheno_plot = pheno_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
pheno_plot = pheno_plot + scale_fill_manual(values = c("black","darkgreen","blue"))
pheno_plot = pheno_plot + ylab("Log of measurement") + theme(legend.position = "top")
pheno_plot

####

scale_data = scale(t(vis_mat_pre[,colnames(vis_mat_pre) != "Subtype"]))
#data = apply(data, MARGIN = 1, FUN =scale)
pheatmap::pheatmap(
  scale_data,
  annotation_col = meta_data["Subtype"],
  annotation_colors = aka3,
  show_rownames = TRUE,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = TRUE,
  fontsize_col = 7,
  clustering_method = "ward.D2"
)

#
