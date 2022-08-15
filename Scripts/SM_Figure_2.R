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

matcher = match(as.character(colnames(cor_mat)), as.character(meta_info$Sample_ID), nomatch = 0)
meta_data = meta_info[matcher,]
rownames(meta_data) = meta_data$Sample_ID

cor_mat = cor(expr_raw);pcr = prcomp(t(cor_mat))
col_vec = as.character(meta_data$Subtype)
col_vec[col_vec == "BA"] = "black"
col_vec[col_vec == "MS"] = "blue"
col_vec[col_vec == "CL"] = "green"

p = ggbiplot::ggbiplot(
  pcr,
  groups = as.character(meta_data$Subtype),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  labels = meta_data$SampleID
)
p = p + scale_color_manual(name="Clusters", values=c("Black", "darkgreen", "blue"))
p = p + theme(legend.position="top")

#svg("~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_2.svg", width = 10, height = 10)
p
dev.off()
