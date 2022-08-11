library("ggplot2")
library("stringr")
library("grid")
library("umap")
library("dplyr")

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

matcher = match(colnames(expr_raw), meta_info$Sample_ID,nomatch = 0)
colnames(expr_raw)[matcher == 0]

matcher_rev = match(meta_info$Sample_ID,colnames(expr_raw),nomatch = 0)
meta_info$Sample_ID[matcher_rev == 0]


meta_data = meta_info[matcher,]
dim(meta_data)
cor_mat = cor(expr_raw);pcr = prcomp(t(cor_mat))

## umap

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
custom.config$random_state = 314 # 314
custom.config$n_components= 2

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
  aes(x, y, label = meta_data$Sample_ID))
#umap_p = umap_p + geom_point()+geom_text(hjust=0, vjust=0)
umap_p = umap_p + geom_point(size = 4, aes(  color = as.character(meta_data$Subtype) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data$Subtype), level=.5, type ="t", size=1.5)
umap_p = umap_p + scale_color_manual( values = c("black","darkgreen","blue")) ##33ACFF ##FF4C33
umap_p = umap_p + theme(legend.position = "top") + xlab("") + ylab("")

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_1.svg", width = 10, height = 10)
umap_p# +geom_text(aes(label=vis_mat$SampleID, color = vis_mat$Subtype),hjust=0, vjust=0)
dev.off()
