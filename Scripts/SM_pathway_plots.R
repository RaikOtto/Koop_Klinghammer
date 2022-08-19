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

###
i = 1
genes_of_interest_hgnc_t = read.table("~/Koop_Klinghammer/Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
genes_of_interest_hgnc_t$V1[i]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]

genes_of_interest_hgnc_t[i,1]

sad_genes[which(!(sad_genes %in% rownames(expr_raw)))]
table(sad_genes %in% rownames(expr_raw) )

expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

#svg("~/Koop_Klinghammer/Results/23_10_2019/Heatmap_94.svg")
selectors = c("Subtype","Grading","Keratinization","Tumor_cell_budding","Cell_nest_size","Mitotic_Count","Nuclear_Size","Necrosis","Inflammatory_Infiltrate","Lymphangiosis","Perineural_Invasion","Overall_Survival_from_diagnosis")
selection = c("Subtype","Grading","Tumor_cell_budding_ROC","Nuclear_Size_ROC","Necrosis_ROC")
selection[!(selection %in% colnames(meta_data))]

vis_mat_cor_plot = meta_data[,selection]
#not_na_vec = !is.na(vis_mat_cor_plot[,"Budding_10HPF"])
#vis_mat_cor_plot[not_na_vec,"Mitosen_10HPF"] = log(vis_mat_cor_plot[not_na_vec,"Mitosen_10HPF"]+1)
#vis_mat_cor_plot[not_na_vec,"Budding_10HPF"] = log(vis_mat_cor_plot[not_na_vec,"Budding_10HPF"]+1)
#vis_mat_cor_plot[not_na_vec,"Zellnestgröße_zentral"] = log(vis_mat_cor_plot[not_na_vec,"Zellnestgröße_zentral"]+1)
#vis_mat_cor_plot[not_na_vec,"Nekrose"] = log(vis_mat_cor_plot[not_na_vec,"Nekrose"]+1)
#vis_mat_cor_plot[not_na_vec,"Entzündung"] = log(vis_mat_cor_plot[not_na_vec,"Entzündung"]+1)

source("~/Koop_Klinghammer/Scripts/Visualization_colors.R")
dim(cor_mat)

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
#dev.off()
genes_of_interest_hgnc_t$V1[i]
