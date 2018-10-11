library("stringr")
library("grid")

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

expr_raw = read.table(
  "/home/ottoraik/Generic_mRNA_Expression_Pipeline/Project_files/Output/ExpressionSet_Klinghammer_neun.tsv",
  sep ="\t",
  stringsAsFactors = F,
  header = T
  #row.names = F
)
colnames(expr_raw ) = str_replace_all(colnames(expr_raw), pattern = "^X", "")
# variance selection
hgnc_list = expr_raw[,1]
hgnc_list_uni = unique(expr_raw[,1])
expr_raw = expr_raw[,-1]
source("~/Koop_Klinghammer/Scripts/Variance_selection.R")

expr_raw[1:5,1:5]
dim(expr_raw)
expr = as.double(as.character(unlist(expr_raw)))
expr = matrix(expr, ncol = ncol(expr_raw),nrow = nrow(expr_raw))
colnames(expr) = colnames(expr_raw)
rownames(expr) = rownames(expr_raw)
### Prep

aka3 = list(
  Copan_Cetux = c(Resistant = "Black", Responder = "White"),
  Cetuximab = c(Resistant = "Black", Responder = "White"),
  Copanlisib = c(Resistant = "Black", Responder = "White"),
  Copanlisib_val = c(high = "Black", low = "White"),
  Cetuximab_val = c(high = "Black", low = "White"),
  Copan_Cetux_val = c(high = "Black", low = "White")
)

###

## Figure 1
meta_data = read.table("~/Generic_mRNA_Expression_Pipeline/Project_files/Klinghammer_neun/RTV.tsv",sep ="\t", stringsAsFactors = F)
colnames(meta_data) = meta_data[1,]
meta_data = meta_data[-1,]
rownames(meta_data) = meta_data$ID
meta_data$Copanlisib  = str_replace_all(meta_data$Copanlisib, "CTRL", "Responder")
meta_data$Copanlisib  = str_replace_all(meta_data$Copanlisib, "CASE", "Resistant")
meta_data$Cetuximab  = str_replace_all(meta_data$Cetuximab, "CTRL", "Responder")
meta_data$Cetuximab  = str_replace_all(meta_data$Cetuximab, "CASE", "Resistant")
meta_data$Copan_Cetux  = str_replace_all(meta_data$Copan_Cetux, "CTRL", "Responder")
meta_data$Copan_Cetux  = str_replace_all(meta_data$Copan_Cetux, "CASE", "Resistant")
meta_data$Copanlisib_val = log2(as.double(str_replace(meta_data$Copanlisib_val, pattern = ",", "."))+1)
meta_data$Cetuximab_val = log2(as.double(str_replace(meta_data$Cetuximab_val, pattern = ",", "."))+1)
meta_data$Copan_Cetux_val = log2(as.double(str_replace(meta_data$Copan_Cetux_val, pattern = ",", "."))+1)

sub_tab = read.table(
    #"/home/ottoraik/Generic_mRNA_Expression_Pipeline/Project_files/Output/Results_Klinghammer_neun/Dif_Copan_Cetux.tsv",
    #"/home/ottoraik/Generic_mRNA_Expression_Pipeline/Project_files/Output/Results_Klinghammer_neun/Dif_Cetuximab.tsv",
    "/home/ottoraik/Generic_mRNA_Expression_Pipeline/Project_files/Output/Results_Klinghammer_neun/Dif_Copanselib.tsv",
    sep = "\t",
    stringsAsFactors = F,
    header = T
)

expr_sub = expr[rownames(expr) %in% sub_tab$HGNC_symb,]
cor_mat = cor(expr_sub);pcr = prcomp(t(cor_mat))
#"Copan_Cetux","Cetuximab","Copanlisib"
pheatmap::pheatmap(
  cor_mat,
  #annotation_col = meta_data[c("Copan_Cetux")],
  annotation_col = meta_data[c("Copanlisib","Copanlisib_val")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  legend = F,
  fontsize_col = 7#,
  #clustering_method = "ward.D2"
)

### Figure 2

ggbiplot::ggbiplot(
  pcr,
  #groups = as.character(meta_data$Copanlisib),
  groups = as.character(meta_data$Copan_Cetux),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  #var.scale = ,
  labels = meta_data$ID
)#  + geom_point( aes( size = as.double(meta_data$OS)**2, color = as.factor(meta_data_tmp$Subtype) ))

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

copanlisib_cor = (apply(expr_raw, MARGIN = 1, FUN = function(vec){return(cor(meta_data$Copanlisib_val, vec))}))
cetuximab_cor = (apply(expr_raw, MARGIN = 1, FUN = function(vec){return(cor(meta_data$Cetuximab_val, vec))}))
copan_cetux_cor = (apply(expr_raw, MARGIN = 1, FUN = function(vec){return(cor(meta_data$Copan_Cetux_val, vec))}))

cor_t = data.frame(
  "Copanlisib" = copanlisib_cor,
  "Cetuximab" = cetuximab_cor,
  "copan_cetux" = copan_cetux_cor
)

pheatmap::pheatmap(
  cor_t[1:20,],
  cluster_rows = F,
  color = colorRampPalette((RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
)
vis_mat = cor_t[!is.na(cor_t)[,1],]
vis_mat = vis_mat[order(vis_mat[,2], decreasing = T),]
pheatmap::pheatmap(
  vis_mat[1:30,],
  cluster_rows = F,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
)
write.table(cor_t, "~/Koop_Klinghammer/Results/Klinghammer_Neun_11_10_2018/Drug_Exp_Correlations.tsv",sep="\t", quote =F, row.names = T)

expr_raw["TIMD4",]
