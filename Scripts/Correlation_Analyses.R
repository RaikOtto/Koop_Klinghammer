library("ggplot2")
library("dplyr")
library("stringr")
library("survminer")
#library("umap")


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
rownames(meta_info) = meta_info$SampleID
matcher = match(colnames(expr_raw),meta_info$SampleID)

meta_data = meta_info[matcher,]
meta_data = meta_data %>% filter(Survivalstatistik_neuestes_Sample == 1)
table(meta_data$Survivalstatistik_neuestes_Sample)

meta_data_vis = meta_data # OS_Monate_ab_Einschluss BA CL PFS_Monate_ab_Einschluss BA CL

# OS_ab_ED # OS_Monate_ab_Einschluss # PFS_Monate_ab_Einschluss # PFS_ab_ED
#parameters = c("Tumorzellgehalt","Grading","Keratinisierung","Budding_1HPF","Budding_10HPF","Zellnestgröße_zentral","Mitosen_HPF","Mitosen_10HPF","Kerngröße","Stroma_vitalerTumor","Nekrose","Entzündung","L1","Pn1","Tumorstadium_codiert","Raucher","Alkohol")
parameters = c("OS_ab_ED","OS_Monate_ab_Einschluss","PFS_Monate_ab_Einschluss","PFS_ab_ED","Tumorzellgehalt","Grading","Keratinisierung","Budding_1HPF","Zellnestgröße_zentral","Mitosen_HPF","Kerngröße","Stroma_vitalerTumor","Nekrose","Entzündung","L1","Pn1","Tumorstadium_codiert","Raucher","Alkohol")
parameters[! parameters %in% colnames(meta_data)]
parameter_matrix = meta_data[,parameters]
parameter_matrix = parameter_matrix[!is.na(parameter_matrix[,1]),]

correlation_mat = cor(parameter_matrix  )
order_mat = as.dist((1-correlation_mat)/2)
hc <- hclust(order_mat)
correlation_mat = correlation_mat[hc$order, hc$order]

correlation_mat[upper.tri(correlation_mat)] = NA
diag(correlation_mat) = 0

melted_cormat <- reshape2::melt(correlation_mat, na.rm = TRUE)

text_size = 7
p = ggheatmap = ggplot(
  melted_cormat,
  aes( Var2, Var1, fill = value)
)
p = p + geom_tile(color = "white") + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1,1)) 
p = p + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = text_size),axis.text.y= element_text(size=text_size)) + coord_fixed()
p = p +   theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank()
    #legend.justification = c(1, 0),
    #legend.position = c(0.6,1),
    #legend.direction = "horizontal"
    )
p

#svg(filename = "~/Downloads/First_draft_correlation_matrix.svg", width = 10, height = 10)
print(p)
dev.off()


for (i in 1:ncol(parameter_matrix)){
  for (j in 1:ncol(parameter_matrix)){
    if (i == j)
      next()
    
    x = parameter_matrix[,i]
    y = parameter_matrix[,j]
    parameter_x = colnames(parameter_matrix)[i]
    parameter_y = colnames(parameter_matrix)[j]
    
    test_result = cor.test(x,y)$p.value
    
    if ( test_result < 0.05){
      print( c(parameter_x, parameter_y))
      print( test_result )
    }
  }
}

