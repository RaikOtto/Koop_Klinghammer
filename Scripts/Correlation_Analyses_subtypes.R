library("tibble")
library("ggplot2")
library("stringr")
library("grid")
library("umap")
library("dplyr")

### Prep

selectors = c("Subtype","Grading","Keratinization","Tumor_cell_budding","Cell_nest_size","Mitotic_Count","Nuclear_Size","Necrosis","Inflammatory_Infiltrate","Lymphangiosis","Perineural_Invasion","Overall_Survival_from_diagnosis")


## Figure 1

meta_data = meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
selectors[ ! (selectors %in% colnames(meta_data))]
# CATEGORICAL

parameters_categorical = c("Subtype","Histotyp","Grading","Tumor_cell_budding_ROC","Cell_nest_size_ROC","Mitotic_Count_ROC","Nuclear_Size_ROC","Stroma_Content_ROC","Necrosis_ROC","Inflammatory_Infiltrate_ROC","Lymphangiosis","Perineural_Invasion","Best_Response","Localization_primary_tumor")
parameters_double = c("Keratinization","Tumor_cell_budding","Cell_nest_size","Mitotic_Count","Nuclear_Size","Necrosis","Inflammatory_Infiltrate","Overall_Survival_from_diagnosis","Overall_Survival_from_Randomisation","Progression_free_Survival_from_Randomisation","Progression_free_Survival_from_Diagnosis","Tumorstadium_codiert")

meta_data_vis = meta_data
dim(meta_data_vis) # 113 37

###

result_mat <<- matrix(as.character(), ncol = 4)

for( selector_1 in selectors ){
  for( selector_2 in selectors ){
  
    ana_table = as_tibble(cbind(
      (meta_data[,selector_1]),
      (meta_data[,selector_2])
    ))
    ana_table = ana_table %>% filter( (!is.na(ana_table[,1]))  & ( ana_table[,1] != "Unknown" ) & ( ana_table[,1] != "" ))
    ana_table = ana_table %>% filter( (!is.na(ana_table[,2]))  & ( ana_table[,2] != "Unknown" ) & ( ana_table[,2] != "" ))
    ana_table = ana_table %>% as_tibble(ana_table)
    
    if ((selector_1 %in% parameters_double)      & (selector_2 %in% parameters_categorical)){
        p_value = as.double(unlist(summary(aov(as.double(ana_table$V1)~as.factor(ana_table$V2))))[9])
        #ltm::biserial.cor(x=as.double(ana_table$V1),y=ana_table$V2 )
        correlation = -1.0
    }
    
    if  ((selector_1 %in% parameters_categorical) & (selector_2 %in% parameters_double)){
      p_value = as.double(unlist(summary(aov(as.double(ana_table$V2)~as.factor(ana_table$V1))))[9])
      correlation = -1.0
    }      
    
    if ((selector_1 %in% parameters_categorical) & (selector_2 %in% parameters_categorical)) {
      p_value = chisq.test(ana_table$V1,ana_table$V2)$p.value
      correlation = -1.0 #as.double(rcompanion::cramerV(x=ana_table$V1,y=ana_table$V2))
      #DescTools::CorPolychor()
    }  
    
    if  ((selector_1 %in% parameters_double) & (selector_2 %in% parameters_double)){
      p_value = cor.test(ana_table$V1,ana_table$V2)$p.value
      correlation =  cor(ana_table$V1,ana_table$V2)
    } 
  
    if (selector_1 == selector_2){
      p_value = 0.0
      correlation = 1.0
    }
      
    result_vec = matrix(c(selector_1, selector_2, correlation,p_value), ncol = 4)
    result_mat = rbind(result_mat, result_vec)
  }
}
colnames(result_mat) = c("Feature_1","Feature_2","Correlation", "P-value")

#write.table(result_mat, "~/Koop_Klinghammer/Results/Phenotype_correlations.tsv", sep ="\t",row.names = FALSE)

### correlation matrix

vis_matrix = result_mat %>% as_tibble()
vis_matrix$Correlation = as.double(vis_matrix$Correlation)
vis_matrix$Correlation[vis_matrix$Correlation == -1] = 0.0
vis_matrix = as.data.frame(matrix( as.double(vis_matrix$Correlation),ncol = 12,nrow=12))
rownames(vis_matrix) = colnames(vis_matrix) = result_mat[1:12,2]

#vis_matrix[vis_matrix == 1] = 0

order_mat = as.dist((1-vis_matrix)/2)
hc <- hclust(order_mat)
vis_matrix = vis_matrix[hc$order, hc$order]

vis_matrix[upper.tri(vis_matrix)] = NA
diag(vis_matrix) = 0

text_size = 7
p = pheatmap::pheatmap(
  vis_matrix,
  show_rownames = TRUE,
  show_colnames = TRUE,
  treeheight_row = 0,
  treeheight_col = 0,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  fontsize_col = 7
)

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_2.svg", width = 10, height = 10)
print(p)
dev.off()


result_mat$Correlation = as.double(result_mat$`P-value`)
melted_cormat <- reshape2::melt(vis_matrix, na.rm = FALSE)
melted_cormat = cbind(colnames(vis_matrix),melted_cormat)
melted_cormat = melted_cormat[!is.na(melted_cormat$value),]
colnames(melted_cormat) = c("Feature_1","Feature_2","Correlation")

p = ggheatmap = ggplot(
  data = melted_cormat,
  aes( Feature_1, Feature_2, fill = Correlation)
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

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Supplement/SM_Figure_2.svg", width = 10, height = 10)
print(p)
dev.off()

