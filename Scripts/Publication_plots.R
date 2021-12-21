library("ggplot2")
library("stringr")
library("grid")

expr_raw = read.table(
  "~/Koop_Klinghammer/Data/Normalize_data.S128.tsv",
  sep ="\t",
  stringsAsFactors = F,
  header = T
  #row.names = F
)
colnames(expr_raw ) = str_replace_all(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
dim(expr_raw)

### Prep

###

## Figure 1
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_Information.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Sample_ID

meta_data = meta_info[colnames(expr_raw),]
#expr_raw = expr_raw[ ,meta_data$Included == "TRUE"  ]
#meta_data
#meta_data = meta_info[colnames(expr_raw),]
#dim(expr_raw)

###
i =82
genes_of_interest_hgnc_t = read.table("~/Koop_Klinghammer/Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
genes_of_interest_hgnc_t$V1[i]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]

genes_of_interest_hgnc_t[i,1]

sad_genes[which(!(sad_genes %in% rownames(expr_raw)))]
table(sad_genes %in% rownames(expr_raw) )

#expr = expr_raw[ rownames(expr_raw) %in% sad_genes,]
expr = expr_raw[ , meta_data$P_value < 5E-2]
exclusion_samples = c("56","54","23","151","71","127","105","98","66","90","83","62","75","130","4","110","21")
length(exclusion_samples)
expr = expr[,!(colnames(expr) %in% exclusion_samples)]
meta_data = meta_info[colnames(expr),]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

#svg("~/Koop_Klinghammer/Results/23_10_2019/Heatmap_94.svg")
selection = c("Subtype","Grading_WHO","Keratinisierung","Budding_10HPF","Zellnestgröße_zentral","Mitosen_10HPF","Nekrose","Entzündung")
#selection = c("Subtype","Grading_WHO","Keratinisierung","Budding_10HPF","Zellnestgröße_ROC","Mitosen_10HPF","Nekrose","Entzündung")

vis_mat_cor_plot = meta_data[,selection]
vis_mat_cor_plot$Budding_10HPF = log(vis_mat_cor_plot$Budding_10HPF+1)
vis_mat_cor_plot$Zellnestgröße_zentral = log(vis_mat_cor_plot$Zellnestgröße_zentral+1)
vis_mat_cor_plot$Nekrose = log(vis_mat_cor_plot$Nekrose+1)
vis_mat_cor_plot$Entzündung = log(vis_mat_cor_plot$Entzündung+1)

dim(expr)
pheatmap::pheatmap(
  cor_mat,
  #expr,
  #annotation_col = meta_data[,selection],
  annotation_col = meta_data["Subtype"],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = F,
  fontsize_col = 7,
  clustering_method = "ward.D"
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

### umap
library("umap")

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
#custom.config$random_state = 995
custom.config$n_components=2

cor_mat = cor(expr)
vis_mat = meta_info[colnames(cor_mat),]

umap_result = umap::umap(
  cor_mat,
  colvec = vis_mat$Subtype,
  #colvec = meta_data$Study,
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
umap_p# +geom_text(aes(label=meta_data$Sample_ID, color = meta_data$Subtype),hjust=0, vjust=0)
custom.config$random_state

#write.table(meta_data,"~/Downloads/Meta_info.S111.tsv",quote =FALSE,sep ="\t",row.names = FALSE)
