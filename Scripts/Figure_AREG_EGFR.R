library("ggplot2")
library("stringr")
library("grid")
library("umap")

source("~/Koop_Klinghammer/Scripts/Visualization_colors.R")

expr_raw = read.table(
  "~/Koop_Klinghammer/Data/S103.tsv",
  sep ="\t",
  stringsAsFactors = F,
  header = T
  #row.names = F
)
colnames(expr_raw ) = str_replace_all(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
dim(expr_raw)

## Figure 1
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Sample_ID

meta_data = meta_info[colnames(expr_raw),]
expr = expr_raw
meta_data = meta_info[colnames(expr),]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

selection = c("Subtype","Grading")
selection[!(selection %in% colnames(meta_data))]
vis_mat_cor_plot = meta_data[,selection]
vis_mat_cor_plot$Grading = as.character(vis_mat_cor_plot$Grading)
vis_mat_cor_plot$Grading[is.na(vis_mat_cor_plot$Grading)] = "Unknown"

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_2.Grading.svg", width = 10, height = 10)
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

table(meta_data$Subtype,meta_data$CefcidNr)
