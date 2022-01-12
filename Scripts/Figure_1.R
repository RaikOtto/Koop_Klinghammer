library("ggplot2")
library("stringr")
library("grid")
library("umap")
library("dplyr")

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
expr = expr_raw
meta_data = meta_info[colnames(expr),]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

## umap

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
custom.config$random_state = 281
custom.config$n_components= 2

cor_mat = cor(expr)
vis_mat = meta_info[colnames(cor_mat),]

umap_result = umap::umap(
  cor_mat,
  colvec = vis_mat$Subtype,
  preserve.seed = TRUE,
  config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point(size = 4, aes(  color = as.character(vis_mat$Subtype) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = vis_mat$Subtype), level=.5, type ="t", size=1.5)
umap_p = umap_p + scale_color_manual( values = c("black","darkgreen","blue")) ##33ACFF ##FF4C33
umap_p = umap_p + theme(legend.position = "none") + xlab("") + ylab("")

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_1.svg", width = 10, height = 10)
umap_p# +geom_text(aes(label=vis_mat$SampleID, color = vis_mat$Subtype),hjust=0, vjust=0)
dev.off()