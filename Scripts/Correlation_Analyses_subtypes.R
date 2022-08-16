library("tibble")
library("ggplot2")
library("stringr")
library("grid")
library("umap")
library("dplyr")

### Prep

selectors = c("Subtype","Grading","Keratinization","Tumor_cell_budding","Cell_nest_size","Mitotic_Count","Nuclear_Size","Necrosis","Inflammatory_Infiltrate","Lymphangiosis","Perineural_Invasion","Overall_Survival_from_diagnosis","Progression_free_Survival_from_Diagnosis")


## Figure 1

meta_data = meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
selectors[ ! (selectors %in% colnames(meta_data))]
# CATEGORICAL

parameters_categorical = c("Subtype","Histotyp","Grading","Tumor_cell_budding_ROC","Cell_nest_size_ROC","Mitotic_Count_ROC","Nuclear_Size_ROC","Stroma_Content_ROC","Necrosis_ROC","Inflammatory_Infiltrate_ROC","Lymphangiosis","Perineural_Invasion","Best_Response","Localization_primary_tumor")
parameters_double = c("Keratinization","Tumor_cell_budding","Cell_nest_size","Mitotic_Count","Nuclear_Size","Necrosis","Inflammatory_Infiltrate","Overall_Survival_from_diagnosis","Overall_Survival_from_Randomisation","Progression_free_Survival_from_Randomisation","Progression_free_Survival_from_Diagnosis","Tumorstadium_codiert")

meta_data_vis = meta_data
dim(meta_data_vis) # 113 37

###

result_mat <<- matrix(as.character(), ncol = 3)

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
    }
    if  ((selector_1 %in% parameters_categorical) & (selector_2 %in% parameters_double))      p_value = as.double(unlist(summary(aov(as.double(ana_table$V2)~as.factor(ana_table$V1))))[9])
    if ((selector_1 %in% parameters_categorical) & (selector_2 %in% parameters_categorical))  p_value = chisq.test(ana_table$V1,ana_table$V2)$p.value
    if  ((selector_1 %in% parameters_categorical) & (selector_2 %in% parameters_categorical))  p_value = chisq.test(ana_table$V1,ana_table$V2)$p.value
  
    if (selector_1 == selector_2) p_value = 0
      
    result_vec = matrix(c(selector_1, selector_2, p_value), ncol = 3)
    result_mat = rbind(result_mat, result_vec)
  }
} 

write.table(result_mat, "~/Koop_Klinghammer/Results/Phenotype_correlations.tsv", sep ="\t",row.names = FALSE)

result_mat = result_mat %>% as_tibble()

### correlation matrix

table(result_mat$V1)

na_vec = apply(parameter_matrix, MARGIN = 1, FUN = function(vec){return(TRUE %in%is.na(vec))})
parameter_matrix = parameter_matrix[! na_vec,]

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

#### correlation heatmap samples

source("~/Koop_Klinghammer/Misc/Visualization_colors.R")
parameters_survival = c("Overall_survial","Progression_free_survial","Progression_free_survial_diagnosis","Overall_survival_diagnosis")
parameters = c(parameters_double,parameters_survival)
parameters[! parameters %in% colnames(meta_data)]
parameter_matrix = meta_data[,parameters]
na_vec = apply(parameter_matrix, MARGIN = 1, FUN = function(vec){return(TRUE %in%is.na(vec))})
parameter_matrix = parameter_matrix[! na_vec,]

correlation_mat = cor( t(parameter_matrix)  )
selection = c("EGFR","AREG","Subtype")
meta_data_vis = meta_data[colnames(correlation_mat),]
meta_data_vis$AREG = scale(as.double(expr_raw["AREG",colnames(correlation_mat)]))
meta_data_vis$EGFR = scale(as.double(expr_raw["EGFR",colnames(correlation_mat)]))

p = pheatmap::pheatmap(
  scale(t(parameter_matrix)),
  #correlation_mat,
  annotation_col = meta_data_vis[selection],
  annotation_colors = aka3,
  show_rownames = TRUE,
  show_colnames = FALSE,
  treeheight_row = 0,
  legend = FALSE,
  fontsize_col = 7,
  clustering_method = "average"
)


####

custom.config = umap.defaults
custom.config$random_state = sample(1:1000,size = 1)
custom.config$random_state = 281
custom.config$n_components= 2

umap_result = umap::umap(
  correlation_mat,
  colvec = meta_data_vis$Subtype,
  preserve.seed = TRUE,
  config=custom.config
)

umap_result$layout = as.data.frame(umap_result$layout)
colnames(umap_result$layout) = c("x","y")

umap_p = ggplot(
  umap_result$layout,
  aes(x, y))
umap_p = umap_p + geom_point(size = 4, aes(  color = as.character(meta_data_vis$Subtype) ))
umap_p = umap_p + stat_ellipse( linetype = 1, aes( color = meta_data_vis$Subtype), level=.5, type ="t", size=1.5)
umap_p = umap_p + scale_color_manual( values = c("black","darkgreen","blue")) ##33ACFF ##FF4C33
umap_p = umap_p + theme(legend.position = "top") + xlab("") + ylab("")

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_1.svg", width = 10, height = 10)
umap_p# +geom_text(aes(label=vis_mat$SampleID, color = vis_mat$Subtype),hjust=0, vjust=0)
