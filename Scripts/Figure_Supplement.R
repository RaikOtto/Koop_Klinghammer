library("ggplot2")
library("stringr")
library("grid")
library("umap")

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
#expr_raw = expr_raw[ ,meta_data$Included == "TRUE"  ]
#meta_data
#meta_data = meta_info[colnames(expr_raw),]
#dim(expr_raw)

###
i = 13
genes_of_interest_hgnc_t = read.table("~/Koop_Klinghammer/Misc/Stem_signatures.tsv",sep ="\t", stringsAsFactors = F, header = F)
genes_of_interest_hgnc_t$V1
genes_of_interest_hgnc_t$V1[i]
sad_genes = str_to_upper( as.character( genes_of_interest_hgnc_t[i,3:ncol(genes_of_interest_hgnc_t)]) ) # 13
sad_genes = sad_genes[sad_genes != ""]

genes_of_interest_hgnc_t[i,1]

sad_genes[which(!(sad_genes %in% rownames(expr_raw)))]
table(sad_genes %in% rownames(expr_raw) )

expr = expr_raw
meta_data = meta_info[colnames(expr),]

### Figure 2 PCA

cor_mat = cor(expr);pcr = prcomp(t(cor_mat))
col_vec = as.character(meta_data$Subtype)
col_vec[col_vec == "BA"] = "black"
col_vec[col_vec == "MS"] = "blue"
col_vec[col_vec == "CL"] = "green"

#svg("~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_2.svg", width = 10, height = 10)

ggbiplot::ggbiplot(
    pcr,
    groups = as.character(meta_data$Subtype),
    ellipse = TRUE,
    circle = TRUE,
    var.axes = F,
    labels = meta_data$SampleID
)  + scale_color_manual(name="Clusters", values=c("Black", "darkgreen", "blue"))

dev.off()

#### Figure 3 clinical metadata

#		Budding_10HPF		Zellnestgröße_zentral		Mitosen_10HPF		Kerngröße	Stroma_vitalerTumor		Nekrose		Entzündung  WHO-Grad L1 und Pn1
selection = c("Subtype","PFS_Monate","Budding_10HPF","Zellnestgröße_zentral","Mitosen_10HPF","Nekrose","Entzündung","Stroma_vitalerTumor")
selection[!(selection %in% colnames(meta_data))]

vis_mat_pre = meta_data[!is.na(meta_data$Entzündung),selection]
vis_mat = reshape2::melt(vis_mat_pre)
colnames(vis_mat) = c("Subtype","Characteristic","Value")
vis_mat = vis_mat %>% filter (!is.na(Value))

vis_mat = vis_mat %>% dplyr::group_by(Subtype, Characteristic) %>% dplyr::summarise( Value = log(Value+1) )

pheno_plot = ggplot(vis_mat, aes( x = Characteristic, y = Value, fill = Subtype) )
pheno_plot = pheno_plot + geom_boxplot(notch = TRUE,outlier.colour = "red", outlier.shape = 1)
pheno_plot = pheno_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
pheno_plot = pheno_plot + scale_fill_manual(values = c("black","darkgreen","blue"))
pheno_plot = pheno_plot + ylab("Log of measurement") + theme(legend.position = "top")
pheno_plot
#pheno_plot + ylim(c(0,4))

###

selection = c("Subtype","Grading_WHO","Grading_WHO","L1","Pn1","Keratinisierung","Kerngröße")
vis_mat_pre = meta_data[!is.na(meta_data$Entzündung),selection]
vis_mat = reshape2::melt(vis_mat_pre)
colnames(vis_mat) = c("Subtype","Characteristic","Value")
vis_mat = vis_mat %>% dplyr::filter (!is.na(Value))

vis_mat = vis_mat %>% dplyr::group_by(Subtype, Characteristic) %>% dplyr::summarise( Value = Value )

pheno_plot = ggplot(vis_mat, aes( x = Characteristic, y = Value, fill = Subtype) )
pheno_plot = pheno_plot + geom_boxplot(notch = TRUE,outlier.colour = "red", outlier.shape = 1)
pheno_plot = pheno_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
pheno_plot = pheno_plot + scale_fill_manual(values = c("black","darkgreen","blue"))
pheno_plot = pheno_plot + ylab("Measurement") + theme(legend.position = "top")
pheno_plot

###

chisq.test(table(vis_mat_pre$Subtype,vis_mat_pre$L1))
chisq.test(table(vis_mat_pre$Subtype,vis_mat_pre$Pn1))
chisq.test(table(vis_mat_pre$Subtype,vis_mat_pre$Grading_WHO))
tab = table(vis_mat_pre$Subtype,vis_mat_pre$Kerngröße)

BA_vec = vis_mat_pre$Keratinisierung[vis_mat_pre$Subtype == "BA"]
CL_vec = vis_mat_pre$Keratinisierung[vis_mat_pre$Subtype == "CL"]
MS_vec = vis_mat_pre$Keratinisierung[vis_mat_pre$Subtype == "MS"]

t.test(BA_vec, CL_vec)
t.test(BA_vec, MS_vec)
anova_1_p_value = TukeyHSD(anova_1)$`as.factor(as.character(meta_grading$Grading))`

###
scale_data = scale(t(vis_mat_pre[,colnames(vis_mat_pre) != "Subtype"]))
#data = apply(data, MARGIN = 1, FUN =scale)
pheatmap::pheatmap(
  scale_data,
  annotation_col = meta_data["Subtype"],
  annotation_colors = aka3,
  show_rownames = TRUE,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = TRUE,
  fontsize_col = 7,
  clustering_method = "ward.D2"
)

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

p = ggpubr::ggarrange(
  p_Keratinisierung, p_Budding_10HPF, p_Zellnestgröße_zentral, p_Mitosen_10HPF,p_Nekrose,p_Entzündung,
  labels = c("Keratinisierung", "Budding_10HPF", "Zellnestgröße","Mitosen_10HPF","Nekrose","Inflammation"),
  ncol = 3, nrow = 2,  common.legend = TRUE)
p
####
library(dplyr)

meta_data[,]

base_data_mean = vis_mat %>% 
  group_by(Subtype,Characteristic) %>% 
  summarize(mean = mean(Value))
base_data_SD = vis_mat %>% 
  group_by(Subtype,Characteristic) %>% 
  summarize(SD = sd(Value))

vis_mat = as.data.frame(vis_mat)
data_vec = vis_mat %>% filter(Characteristic == "Keratinisierung")#[,-2]
chisq.test()
#t.test(as.double(Kerngröße_t$V1[Kerngröße_t$V2 == "BA"]),as.double(Kerngröße_t$V1[Kerngröße_t$V2 == "CL"]))
