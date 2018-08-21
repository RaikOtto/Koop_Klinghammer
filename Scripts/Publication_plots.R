library("stringr")
library("grid")

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

### Prep

aka3 = list(
  Subtype = c(BA = "blue", CL ="darkgreen", MS = "red", Not_sig = "gray" ),
  Location = c(Primary = "white", Metastasis = "black"),
  OS = c(high = "White", medium = "gray", low = "black"),
  Correlation = c(high = "Red", medium = "Yellow", low = "Green"),
  Included = c(Yes = "green", No = "red"))

###

cor_mat = cor(expr);pcr = prcomp(t(pure_data))

## Figure 1
meta_data_tmp = meta_data
meta_data_tmp$OS = log2(meta_data_tmp$OS + 1)
meta_data_tmp$Subtype[meta_data_tmp$Sig == "FALSE"] = "Not_sig"

pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data_tmp[c("Subtype","OS")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  legend = F,
  fontsize_col = 7,
  clustering_method = "ward.D2"
)

### Figure 2

ggbiplot::ggbiplot(
  pcr,
  groups = as.character(meta_data_tmp$Subtype),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  var.scale = meta_data_tmp$OS*10,
  labels = meta_data$Name
)  + geom_point( aes( size = as.double(meta_data$OS)**2, color = as.factor(meta_data_tmp$Subtype) ))

aggregate(meta_data$OS, FUN = mean, by = list(meta_data$Subtype))

### umap

umap_plot = umap::umap(t(pure_data))
vis_data = as.data.frame(umap_plot$layout)
colnames(vis_data) = c("x","y")
dist_mat = dist((vis_data))
p = ggplot2::qplot( x = vis_data$x, y = vis_data$y, color = meta_data$Subtype)
p

## km plots


fit <- survival::survfit( survival::Surv( OS) ~ Subtype,data = meta_data_tmp)

# Visualize with survminer

survminer::ggsurvplot(fit, data = meta_data_tmp, risk.table = T, pval = T)

## Same patient only

find_vec = meta_data$pID
match_vec = match(find_vec, unique(find_vec),nomatch = 0)
multi_match = which( table(match_vec) > 1  )
multi_match = multi_match[multi_match != 12]
true_match = meta_data$pID[which( match_vec %in% multi_match)]

multi_data = pure_data[,which( meta_data$pID %in% true_match )]
multi_cor = cor(multi_data)
meta_data$pID = as.character(meta_data$pID)

pheatmap::pheatmap(
  multi_cor,
  annotation_col = meta_data[c("Subtype","OS","pID")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  legend = F,
  fontsize_col = 7
)

## Selected patient only
# pID = c(2, 362, 365, 397)

multi_meta_data = as.data.frame(meta_data[which( meta_data$pID %in%  365),])
multi_data = pure_data[, multi_meta_data$Name]
multi_cor = cor(multi_data)
pcr = prcomp(t(multi_cor))

pheatmap::pheatmap(
  multi_cor,
  annotation_col = meta_data[c("Subtype","OS","pID")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  legend = F,
  fontsize_col = 7
)

ggbiplot::ggbiplot(
  pcr,
  groups = as.character(multi_meta_data$Subtype),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  labels = multi_meta_data$Name
) 
