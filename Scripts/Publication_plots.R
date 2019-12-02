library("stringr")
library("grid")

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

expr_raw = read.table(
  "~/Koop_Klinghammer/Data/Pure_data.05_06_2018.tsv",
  sep ="\t",
  stringsAsFactors = F,
  header = T
  #row.names = F
)
colnames(expr_raw ) = str_replace_all(colnames(expr_raw), pattern = "^X", "")

# variance selection
hgnc_list = expr_raw[,1]
hgnc_list_uni = unique(expr_raw[,1])
#expr_raw = expr_raw[,-1]
#source("~/Koop_Klinghammer/Scripts/Variance_selection.R")

expr_raw[1:5,1:5]
dim(expr_raw)
#expr = as.double(as.character(unlist(expr_raw)))
#expr = matrix(expr, ncol = ncol(expr_raw),nrow = nrow(expr_raw))
#colnames(expr) = colnames(expr_raw)
#rownames(expr) = rownames(expr_raw)
### Prep

aka3 = list(
  Group = c(refractive = "red", sensitive = "darkgreen", intermediate = "orange"),
  Subtype = c(BA = "black", CL = "green", MS = "blue"),
  OS = c(high = "red", medium = "orange", low = "green")
)

###

## Figure 1
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_Information_Victoria.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Raw_Name
meta_data = meta_info[ colnames(expr_raw),]
meta_data = meta_data[ meta_data$Included,  ]

expr_raw = expr_raw[,as.character(meta_data$Sample_ID)]
expr_raw = expr_raw[,meta_data$Included]
dim(expr_raw)

###

subset_genes_t = read.table("~/Koop_Klinghammer/Misc/Stem_signatures.tsv", sep ="\t", header = F, stringsAsFactors = F)
subset_genes = as.character(subset_genes_t[7,3:ncol(subset_genes_t)])

expr_sub = expr_raw#[ rownames(expr_raw) %in% subset_genes,]
cor_mat = cor(expr_sub);pcr = prcomp(t(cor_mat))

pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[c("Best_response","Subtype")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = F,
  treeheight_col = 0,
  legend = F,
  fontsize_col = 7,
  clustering_method = "ward.D2"
)

### Figure 2

ggbiplot::ggbiplot(
    pcr,
    groups = as.character(meta_data$Subtype),
    ellipse = TRUE,
    circle = TRUE,
    var.axes = F,
    labels = meta_data$Sample_ID
)#  + geom_point( aes( size = as.double(meta_data$OS)**2, color = as.factor(meta_data_tmp$Subtype) ))

#fisher.test(meta_data$Group,meta_data$Subtype)

aggregate(meta_data$OS, FUN = mean, by = list(meta_data$Subtype))

### umap

umap_plot = umap::umap(t(expr_sub))
vis_data = as.data.frame(umap_plot$layout)
colnames(vis_data) = c("x","y")
dist_mat = dist((vis_data))
p = ggplot2::qplot( x = vis_data$x, y = vis_data$y, color = (meta_data[,"Subtype"]))
p

## km plots


fit <- survival::survfit( survival::Surv( OS) ~ Subtype,data = meta_data)

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
  groups = as.character(meta_data$Group),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  labels = meta_data$ID
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

## MEN1 transcripts

library(ggplot2)

meta_data$CDKN1B = rep("",nrow(meta_data))
meta_data$CDKN1B = as.double(expr_raw["CDKN1B", meta_data$ID])
  
vis_mat = reshape2::melt( meta_data[,c("ID","CDKN1B","Copanlisib","Copanlisib_val")], id = c("ID","Copanlisib","CDKN1B","Copanlisib_val"))
vis_mat = vis_mat[order(vis_mat$CDKN1B),]

# Basic barplot
label_vec = meta_data[rownames(vis_mat),"Copanlisib"]
label_vec = label_vec[order(meta_data$CDKN1B)]
label_vec[label_vec == "Resistant"] = "R"
label_vec[label_vec != "R"] = "S"
col_vec = label_vec
col_vec[col_vec == "S"] = "darkgreen"
col_vec[col_vec != "darkgreen"] = "darkred"

vis_mat$ID = factor(vis_mat$ID, levels = vis_mat$ID)
vis_mat$Copanlisib_val = log(vis_mat$Copanlisib_val)

p = ggplot( data = vis_mat)
p = p + geom_bar(aes( x = ID, y = CDKN1B, fill = Copanlisib_val ),stat="identity", colour="black")
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + scale_fill_gradientn(colours = c("white","yellow","red"), breaks = c(0.0,.5,1.0))
p = p + annotate("text", x=1:nrow(vis_mat),y = 12,parse=TRUE, label = label_vec, color = col_vec, size = 4.5 )
p = p + annotate("text", x=5,y = 10,parse=TRUE, label = "pearson-correlation: -0.66", color = "black", size = 4.5 )
p = p + xlab("") + ylab("CDKN1B expression versus Copalisib") + theme(legend.position = "top")
p
