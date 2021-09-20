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
dim(expr_raw)

### Prep

aka3 = list(
  Group = c(refractive = "red", sensitive = "darkgreen", intermediate = "orange"),
  Subtype = c(BA = "black", CL = "darkgreen", MS = "blue"),
  OS = c(high = "red", medium = "orange", low = "green")
)

###

## Figure 1
meta_info = read.table("~/Koop_Klinghammer/Misc/Datensatz_CEFCID.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Raw_Name
meta_data = meta_info[ colnames(expr_raw),]
meta_data = meta_data[ meta_data$Included %in% c("TRUE"),  ]

expr_raw = expr_raw[,as.character(meta_data$Sample_ID)]
meta_data = meta_info[colnames(expr_raw),]
dim(expr_raw)

###

expr_sub = expr_raw#[ rownames(expr_raw) %in% subset_genes,]
cor_mat = cor(expr_sub);pcr = prcomp(t(cor_mat))

#svg("~/Koop_Klinghammer/Results/23_10_2019/Heatmap_94.svg")
pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data[c("Subtype")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  legend = F,
  fontsize_col = 7,
  clustering_method = "ward.D2"
)
#dev.off()

### Figure 2

col_vec = as.character(meta_data$Subtype)
col_vec[col_vec == "BA"] = "black"
col_vec[col_vec == "MS"] = "blue"
col_vec[col_vec == "CL"] = "green"

#svg("~/Koop_Klinghammer/Results/23_10_2019/PCA_63.svg")
ggbiplot::ggbiplot(
    pcr,
    groups = as.character(meta_data$Subtype),
    ellipse = TRUE,
    circle = TRUE,
    var.axes = F,
    labels = meta_data$Sample_ID
)  + scale_color_manual(name="Clusters", values=c("Black", "darkgreen", "blue"))
#dev.off()

aggregate(meta_data$OS, FUN = mean, by = list(meta_data$Subtype))

### Gene_set_reductions

genes_of_interest_hgnc_t = read.table("~/Koop_Klinghammer/Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1

i = 10
genes_of_interest_hgnc_t[i,1]

sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) )
sad_genes = sad_genes[ sad_genes != ""]
length(sad_genes)

meta_data = meta_info[colnames(expr_raw),]
table(rownames(expr_raw) %in% sad_genes)
sad_genes[sad_genes %in% rownames(expr_raw)]
expr = expr_raw[which(rownames(expr_raw) %in% sad_genes),]
expr[1:nrow(expr),1:5]
dim(expr)

vis_mat = reshape2::melt(expr)
#vis_mat = cbind(vis_mat, "IL1B")
vis_mat = cbind(vis_mat, rownames(expr))
vis_mat = cbind(vis_mat, meta_data$Subtype)
colnames(vis_mat) = c("Sample","Expression","Gene","Cluster")
vis_mat$Expression = as.double(vis_mat$Expression)

### cor_mat

correlation_matrix = cor(expr)

# waterfall

vis_mat$Sample = factor(vis_mat$Sample, levels = vis_mat$Sample[order(vis_mat$Expression)])

p = ggplot( data = vis_mat,aes( x = Sample, y = Expression, fill = Cluster ))
p = p + geom_col(position = position_dodge2( preserve = "single")) 

p = p + geom_bar(stat="identity", position=position_dodge())
p = p+ theme_minimal() + xlab("") + ylab("") + theme(legend.position="top")
p = p + theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())# + ylim(0,6)
p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p = p + labs(fill = "Necros-associated expression of IL1B")
p = p + annotate("text", x=40,y = 9.75, label= "Wilcoxon-Smith test on rank MS vs. rank BA p-value = 0.0658", size = 4.5 )
p

cluster = vis_mat$Cluster[order(vis_mat$Expression)]
rank_MS = which(cluster == "MS")
rank_CL = which(cluster == "CL")
rank_BA = which(cluster == "BA")

rank_MS_CL = c(rank_MS, rank_CL)

wilcox.test(rank_MS,rank_BA)
wilcox.test(rank_MS,rank_CL)
wilcox.test(rank_BA,rank_CL)

wilcox.test(rank_MS_CL,rank_BA)

## km plots


fit <- survival::survfit( survival::Surv( OS) ~ Subtype,data = meta_data)

# Visualize with survminer

survminer::ggsurvplot(fit, data = meta_data, risk.table = T, pval = T)

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
  groups = as.character(meta_data$Subtype),
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  labels = meta_data$Sample_ID
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
#write.table(cor_t, "~/Koop_Klinghammer/Results/Klinghammer_Neun_11_10_2018/Drug_Exp_Correlations.tsv",sep="\t", quote =F, row.names = T)

expr_raw["TIMD4",]

## MEN1 transcripts

library(ggplot2)

meta_data$AREG = rep("",nrow(meta_data))
meta_data$EGFR = rep("",nrow(meta_data))
meta_data$AREG = as.double(expr_raw["AREG",])
meta_data$EGFR = as.double(expr_raw["EGFR",])
  
vis_mat = reshape2::melt(
  meta_data[,c("AREG","EGFR","Subtype")],
  id = c("Subtype")
)
colnames(vis_mat) = c("Subtype","Gene","Value")


# Basic barplot
#label_vec = meta_data[rownames(vis_mat),"Subtype"]
#label_vec = label_vec[order(meta_data$AREG)]
#label_vec[label_vec == "Resistant"] = "R"
#label_vec[label_vec != "R"] = "S"
col_vec = as.character(meta_data$Subtype)
col_vec[col_vec == "BA"] = "black"
col_vec[col_vec == "MS"] = "blue"
col_vec[col_vec == "CL"] = "green"

p = ggplot( data = vis_mat,aes( x = Gene, y = Value, fill = Subtype ))
p = p + geom_boxplot()
p = p + scale_fill_manual( values = c("Black","Darkgreen","Blue"))
p = p + theme(legend.position = "top")
p

#svg("~/Koop_Klinghammer/Results/23_10_2019/Bar_plot.svg")
p
#dev.off()
#p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))

p = p + annotate("text", x=1:nrow(vis_mat),y = 12,parse=TRUE, label = label_vec, color = col_vec, size = 4.5 )
p = p + annotate("text", x=5,y = 10,parse=TRUE, label = "pearson-correlation: -0.66", color = "black", size = 4.5 )
p = p + xlab("") + ylab("CDKN1B expression versus Copalisib") 
p
