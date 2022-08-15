library("ggplot2")
library("stringr")
library("grid")
library("umap")

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

## Figure 1
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Sample_ID

cor_mat = cor(expr_raw);pcr = prcomp(t(cor_mat))
matcher = match(as.character(colnames(cor_mat)), as.character(meta_info$Sample_ID), nomatch = 0)
meta_data = meta_info[matcher,]
rownames(meta_data) = meta_data$Sample_ID

selection = c("Subtype","Inflammatory_Infiltrate_ROC","Necrosis_ROC","Lymphangiosis","Grading")
selection %in% colnames(meta_data)

p = pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[,selection],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = TRUE,
  treeheight_row = 0,
  legend = FALSE,
  fontsize_col = 7,
  clustering_method = "ward.D2"
)
dev.off()

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_1.svg", width = 10, height = 10)
print(p)
dev.off()

