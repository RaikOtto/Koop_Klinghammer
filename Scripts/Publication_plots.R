library("ggplot2")
library("stringr")
library("grid")
library("umap")

expr_raw = read.table(
  "~/Koop_Klinghammer/Data/Normalized_data.S121.Significant.tsv",
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

expr = expr_raw#[ rownames(expr_raw) %in% sad_genes,]
exclusion_samples = c("36","33","66","89","61","82","74","21","2","55","22","53","70","125","104","64","94")
length(exclusion_samples)
expr = expr[,!(colnames(expr) %in% exclusion_samples)]
meta_data = meta_info[colnames(expr),]
cor_mat = cor(expr);pcr = prcomp(t(cor_mat))

#svg("~/Koop_Klinghammer/Results/23_10_2019/Heatmap_94.svg")
selection = c("Subtype","Grading_WHO","Keratinisierung","Budding_10HPF","Zellnestgröße_zentral","Mitosen_10HPF","Nekrose","Entzündung")
selection = c("Subtype","Grading_WHO")
selection[!(selection %in% colnames(meta_data))]

vis_mat_cor_plot = meta_data[,selection]
#not_na_vec = !is.na(vis_mat_cor_plot[,"Budding_10HPF"])
#vis_mat_cor_plot[not_na_vec,"Mitosen_10HPF"] = log(vis_mat_cor_plot[not_na_vec,"Mitosen_10HPF"]+1)
#vis_mat_cor_plot[not_na_vec,"Budding_10HPF"] = log(vis_mat_cor_plot[not_na_vec,"Budding_10HPF"]+1)
#vis_mat_cor_plot[not_na_vec,"Zellnestgröße_zentral"] = log(vis_mat_cor_plot[not_na_vec,"Zellnestgröße_zentral"]+1)
#vis_mat_cor_plot[not_na_vec,"Nekrose"] = log(vis_mat_cor_plot[not_na_vec,"Nekrose"]+1)
#vis_mat_cor_plot[not_na_vec,"Entzündung"] = log(vis_mat_cor_plot[not_na_vec,"Entzündung"]+1)

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_1.Grading.svg", width = 10, height = 10)
pheatmap::pheatmap(
  cor_mat,
  #vis_mat_cor_plot,
  annotation_col = vis_mat_cor_plot,
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = TRUE,
  fontsize_col = 7,
  clustering_method = "ward.D"
)
dev.off()

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

####

#Budding_1HPF		Budding_10HPF		Zellnestgröße_zentral		Mitosen_HPF		Mitosen_10HPF		Kerngröße	Stroma_vitalerTumor		Nekrose		Entzündung
selection_vis = c("Subtype","PFS_Monate","Keratinisierung","Budding_1HPF","Budding_10HPF","Zellnestgröße_zentral","Mitosen_HPF","Mitosen_10HPF","Nekrose","Entzündung","Kerngröße","Stroma_vitalerTumor")
vis_mat_pre = meta_data[!is.na(meta_data$Entzündung),selection_vis]
#scale_data = scale((vis_mat_pre[,colnames(vis_mat_pre) != "Subtype"]))

vis_mat = reshape2::melt(vis_mat_pre)
colnames(vis_mat) = c("Subtype","Characteristic","Value")
vis_mat = vis_mat %>% group_by(Subtype, Characteristic) %>% summarise( Proportion = log(Value+1) )

#vis_mat_pre$Budding_10HPF = log(vis_mat_pre$Budding_10HPF+1)
#vis_mat_pre$Zellnestgröße_zentral = log(vis_mat_pre$Zellnestgröße_zentral)
#vis_mat_pre$Mitosen_10HPF = log(vis_mat_pre$Mitosen_10HPF)
#vis_mat_pre$Nekrose  = log(vis_mat_pre$Nekrose+1)
#vis_mat_pre$Entzündung = log(vis_mat_pre$Entzündung+1)
#vis_mat_pre$Stroma_vitalerTumor = log(vis_mat_pre$Stroma_vitalerTumor+1)

#vis_mat = reshape2::melt(vis_mat_pre)
#colnames(vis_mat) = c("Subtype","Characteristic","Value")

#
pheno_plot = ggplot(vis_mat, aes( x = Characteristic, y = Proportion, fill = Subtype) )
#pheno_plot = pheno_plot + geom_boxplot(notch = TRUE,outlier.colour = "red", outlier.shape = 1)
pheno_plot = pheno_plot + geom_violin()
pheno_plot = pheno_plot + theme(axis.text=element_text(size=14)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + theme(legend.text = element_text(size=14))
pheno_plot = pheno_plot + scale_fill_manual(values = c("black","darkgreen","blue"))
pheno_plot = pheno_plot + ylab("Log of measurement") + theme(legend.position = "top")
pheno_plot
#pheno_plot + ylim(c(0,4))

vis_mat_pre = vis_mat_pre[!is.na(vis_mat_pre$Budding_10HPF),]

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

### umap

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
custom.config$random_state = 281
custom.config$n_components=2

#cor_mat = cor(t(vis_mat_pre[,colnames(vis_mat_pre) != "Subtype"]))
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

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_1.svg", width = 10, height = 10)
#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_2.svg", width = 10, height = 10)
umap_p# +geom_text(aes(label=vis_mat$SampleID, color = vis_mat$Subtype),hjust=0, vjust=0)
dev.off()

custom.config$random_state
#281

#write.table(meta_data,"~/Downloads/Meta_info.S111.tsv",quote =FALSE,sep ="\t",row.names = FALSE)
