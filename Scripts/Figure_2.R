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
matcher = match(as.character(colnames(cor_mat)), as.character(meta_info$SampleID), nomatch = 0)
meta_data = meta_info[matcher,]
rownames(meta_data) = meta_data$SampleID

selection = c("Subtype","Grading")
meta_data[is.na(meta_data$Grading),"Grading"] = "Unknown"
meta_data$Grading = factor(meta_data$Grading, levels = c("3","2","1","Unknown"))

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_2.Grading.svg", width = 10, height = 10)
pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[,selection],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = TRUE,
  fontsize_col = 7,
  clustering_method = "ward.D2"
)
dev.off()

#table(meta_data$Subtype,meta_data$CefcidNr)
