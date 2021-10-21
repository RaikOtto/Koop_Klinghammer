library("ggplot2")
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

### Prep

###

## Figure 1
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_Information.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Sample_ID

meta_data = meta_info[colnames(expr_raw),]
expr_raw = expr_raw[ ,meta_data$Included == "TRUE"  ]
meta_data = meta_info[colnames(expr_raw),]
dim(expr_raw)

###
i =82
genes_of_interest_hgnc_t = read.table("~/Koop_Klinghammer/Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
#genes_of_interest_hgnc_t = read.table("~/Downloads/GOBP_CELL_CYCLE.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
genes_of_interest_hgnc_t$V1[i]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
#sad_genes = genes_of_interest_hgnc_t$V1[3:nrow(genes_of_interest_hgnc_t)]
sad_genes = sad_genes[sad_genes != ""]

genes_of_interest_hgnc_t[i,1]

sad_genes[which(!(sad_genes %in% rownames(expr_raw)))]
table(sad_genes %in% rownames(expr_raw) )

expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

#svg("~/Koop_Klinghammer/Results/23_10_2019/Heatmap_94.svg")
selection = c("Subtype","Grading_WHO","Keratinisierung","Budding_10HPF","Zellnestgröße_zentral","Mitosen_10HPF","Nekrose","Entzündung")
#selection = c("Subtype","Grading_WHO","Keratinisierung","Budding_10HPF","Zellnestgröße_ROC","Mitosen_10HPF","Nekrose","Entzündung")

vis_mat_cor_plot = meta_data[,selection]
vis_mat_cor_plot$Budding_10HPF = log(vis_mat_cor_plot$Budding_10HPF+1)
vis_mat_cor_plot$Zellnestgröße_zentral = log(vis_mat_cor_plot$Zellnestgröße_zentral+1)
vis_mat_cor_plot$Nekrose = log(vis_mat_cor_plot$Nekrose+1)
vis_mat_cor_plot$Entzündung = log(vis_mat_cor_plot$Entzündung+1)

pheatmap::pheatmap(
  cor_mat,
  #expr,
  annotation_col = vis_mat_cor_plot,
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = F,
  treeheight_row = 0,
  legend = F,
  fontsize_col = 7,
  clustering_method = "average"
)

aka3 = list(
  Group = c(refractive = "red", sensitive = "darkgreen", intermediate = "orange"),
  Subtype = c(BA = "black", CL = "darkgreen", MS = "blue"),
  OS = c(high = "red", medium = "orange", low = "green"),
  Budding_absentpresent = c("2" = "black", "1" = "white"),
  Zellnestgröße_zentral = c(low = "white", high = "black"),
  Mitosen_10HPF = c(low = "white", high = "black"),
  Nekrose = c(low = "white", high = "black"),
  Entzündung = c(low = "white", high = "black"),
  Budding_10HPF = c(low = "white", high = "black"),
  Keratinisierung = c(low = "white", high = "black"),
  Grading_WHO = c( "1" = "#CBDB34","2" = "#FFC000", "3" = "#EE5124")
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

dev.off()
#dev.off()

####

selection_vis = c("Subtype","Keratinisierung","Budding_10HPF","Zellnestgröße_zentral","Mitosen_10HPF","Nekrose","Entzündung","Kerngröße","L1","Pn1","Stroma_vitalerTumor")
vis_mat_pre = meta_data[!is.na(meta_data$Entzündung),selection_vis]
vis_mat_pre$Budding_10HPF = log(vis_mat_pre$Budding_10HPF+1)
vis_mat_pre$Zellnestgröße_zentral = log(vis_mat_pre$Zellnestgröße_zentral)
vis_mat_pre$Mitosen_10HPF = log(vis_mat_pre$Mitosen_10HPF)
vis_mat_pre$Nekrose  = log(vis_mat_pre$Nekrose+1)
vis_mat_pre$Entzündung = log(vis_mat_pre$Entzündung+1)
vis_mat_pre$Stroma_vitalerTumor = log(vis_mat_pre$Stroma_vitalerTumor+1)
vis_mat = reshape2::melt(vis_mat_pre)
colnames(vis_mat) = c("Subtype","Characteristic","Value")

p = ggplot( data = vis_mat,aes( x = Characteristic, y = Value, fill = Subtype ))
p = p +  geom_boxplot( ) + xlab("")
p

#

p_Keratinisierung= ggplot( data = vis_mat[vis_mat$Characteristic %in% "Keratinisierung",],aes( x = Characteristic, y = Value, fill = Subtype ))+  geom_boxplot( ) + xlab("") + scale_fill_manual(values=c("darkgreen","orange","blue"))
p_Budding_10HPF= ggplot( data = vis_mat[vis_mat$Characteristic %in% "Budding_10HPF",],aes( x = Characteristic, y = log(Value+1), fill = Subtype ))+  geom_boxplot( ) + xlab("") + scale_fill_manual(values=c("darkgreen","orange","blue"))
p_Zellnestgröße_zentral= ggplot( data = vis_mat[vis_mat$Characteristic %in% "Zellnestgröße_zentral",],aes( x = Characteristic, y = log(Value+1), fill = Subtype ))+  geom_boxplot( ) + xlab("") + scale_fill_manual(values=c("darkgreen","orange","blue"))
p_Mitosen_10HPF= ggplot( data = vis_mat[vis_mat$Characteristic %in% "Mitosen_10HPF",],aes( x = Characteristic, y = Value, fill = Subtype ))+  geom_boxplot( ) + xlab("") + scale_fill_manual(values=c("darkgreen","orange","blue"))
p_Nekrose= ggplot( data = vis_mat[vis_mat$Characteristic %in% "Nekrose",],aes( x = Characteristic, y = log(Value+1), fill = Subtype ))+  geom_boxplot( ) + xlab("") + scale_fill_manual(values=c("darkgreen","orange","blue"))
p_Entzündung= ggplot( data = vis_mat[vis_mat$Characteristic %in% "Entzündung",],aes( x = Characteristic, y = log(Value+1), fill = Subtype ))+  geom_boxplot( ) + xlab("") + scale_fill_manual(values=c("darkgreen","orange","blue"))
p_Kerngröße= ggplot( data = vis_mat[vis_mat$Characteristic %in% "Kerngröße",],aes( x = Characteristic, y = log(Value+1), fill = Subtype ))+  geom_boxplot( ) + xlab("") + scale_fill_manual(values=c("darkgreen","orange","blue"))
p_L1= ggplot( data = vis_mat[vis_mat$Characteristic %in% "L1",],aes( x = Characteristic, y = log(Value+1), fill = Subtype ))+  geom_boxplot( ) + xlab("") + scale_fill_manual(values=c("darkgreen","orange","blue"))
p_Pn1= ggplot( data = vis_mat[vis_mat$Characteristic %in% "Pn1",],aes( x = Characteristic, y = log(Value+1), fill = Subtype ))+  geom_boxplot( ) + xlab("") + scale_fill_manual(values=c("darkgreen","orange","blue"))
p_Stroma= ggplot( data = vis_mat[vis_mat$Characteristic %in% "Stroma_vitalerTumor",],aes( x = Characteristic, y = log(Value+1), fill = Subtype ))+  geom_boxplot( ) + xlab("") + scale_fill_manual(values=c("darkgreen","orange","blue"))

p = ggpubr::ggarrange(p_Keratinisierung, p_Budding_10HPF, p_Zellnestgröße_zentral, p_Mitosen_10HPF,p_Nekrose,p_Entzündung,
              labels = c("Keratinisierung", "Budding_10HPF", "Zellnestgröße","Mitosen_10HPF","Nekrose","Inflammation"),
              ncol = 3, nrow = 2,  common.legend = TRUE)
p
####

aggregate( vis_mat_pre$Budding_10HPF, FUN = mean, by = list(vis_mat_pre$Subtype))
aggregate( vis_mat_pre$Budding_10HPF, FUN = max, by = list(vis_mat_pre$Subtype))

#

Zellnestgröße_ROC = meta_data[ !is.na(meta_data$Budding_10HPF), "Zellnestgröße_zentral"]
aggregate( Zellnestgröße_ROC, FUN = mean, by = list(subtype))

#

Mitosen_10HPF = meta_data[ !is.na(meta_data$Mitosen_10HPF), "Mitosen_10HPF"]
aggregate( Mitosen_10HPF, FUN = mean, by = list(subtype))

#

Mitosen_1HPF = meta_data[ !is.na(meta_data$Mitosen_HPF), "Mitosen_HPF"]
aggregate( Mitosen_1HPF, FUN = mean, by = list(subtype))

#

Nekrose = meta_data[ !is.na(meta_data$Nekrose), "Nekrose"]
aggregate( Nekrose, FUN = mean, by = list(subtype))

#

Entzündung = meta_data[ !is.na(meta_data$Entzündung), "Entzündung"]
aggregate( Entzündung, FUN = mean, by = list(subtype))

# Grading

grading_t = as.data.frame(cbind(meta_data$Grading_WHO, meta_data$Subtype))
grading_t = grading_t[!is.na(grading_t$V1) ,]
table(grading_t)
chisq.test(table(grading_t))

# L1

L1_t = as.data.frame(cbind(meta_data$L1, meta_data$Subtype))
L1_t = L1_t[!is.na(L1_t$V1) ,]
table(L1_t)
chisq.test(table(L1_t))

# Pn1, Strom

Pn1_t = as.data.frame(cbind(meta_data$Pn1, meta_data$Subtype))
Pn1_t = Pn1_t[!is.na(Pn1_t$V1) ,]
table(Pn1_t)
chisq.test(table(Pn1_t))

# Stroma

Stroma_t = as.data.frame(cbind(meta_data$Stroma_vitalerTumor, meta_data$Subtype))
Stroma_t = Stroma_t[!is.na(Stroma_t$V1) ,]
t.test(as.double(Stroma_t$V1[Stroma_t$V2 == "CL"]),as.double(Stroma_t$V1[Stroma_t$V2 == "MS"]))

# Kerngröße

Kerngröße_t = as.data.frame(cbind(meta_data$Kerngröße, meta_data$Subtype))
Kerngröße_t = Kerngröße_t[!is.na(Kerngröße_t$V1) ,]
t.test(as.double(Kerngröße_t$V1[Kerngröße_t$V2 == "BA"]),as.double(Kerngröße_t$V1[Kerngröße_t$V2 == "CL"]))

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
