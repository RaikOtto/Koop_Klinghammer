library("ggplot2")
library("stringr")
library("grid")
library("umap")
library("dplyr")

### Prep

selectors = c("Subtype","Histotyp","Grading_WHO","Keratinization","Tumor_cell_budding","Tumor_cell_budding_2_tier","Cell_nest_size","Cell_nest_size_2_tier","Mitotic_Count","Mitotic_Count_2_tier","Nuclear_Size","Nuclear_Size_2_tier","Stroma_Content","Stroma_Content_2_tier","Necrosis","Necrosis_2_tier","Inflammatory_Infiltrate","Inflammatory_Infiltrate_2_tier","Lymphangiosis","Perineural_Invasion","DCR","Best_Response","Localization_primary_tumor","Tumorstadium_codiert","Overall_Survival_from_diagnosis","Overall_Survivall_from_Randomisation","Progression_free_Survivall_from_Randomisation","Progression_free_Survivall_from_Diagnosis")

## Figure 1

meta_data = meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
selectors[ ! (selectors %in% colnames(meta_data))]
# CATEGORICAL

parameters_categorical = c("Subtype","Histotyp","Grading_WHO","Tumor_cell_budding_2_tier","Cell_nest_size_2_tier","Mitotic_Count_2_tier","Nuclear_Size_2_tier","Stroma_Content_2_tier","Necrosis_2_tier","Inflammatory_Infiltrate_2_tier","Lymphangiosis","Perineural_Invasion","DCR","Best_Response","Localization_primary_tumor")
parameters_double = c("Keratinization","Tumor_cell_budding","Cell_nest_size","Mitotic_Count","Nuclear_Size","Stroma_Content","Necrosis","Inflammatory_Infiltrate","Overall_Survival_from_diagnosis","Overall_Survivall_from_Randomisation","Progression_free_Survivall_from_Randomisation","Progression_free_Survivall_from_Diagnosis","Tumorstadium_codiert")

meta_data_vis = meta_data#[meta_data$Subtype %in% c("BA","CL"),]
dim(meta_data_vis) # 113 37

result_mat_categorial <<- matrix(as.character(), ncol = 4)
for( i in seq(1,length(parameters_categorical)-1) ){
  for( j in seq(i+1,length(parameters_categorical)) ){
    
    parameter_1 = parameters_categorical[i]
    parameter_2 = parameters_categorical[j]
    
    ana_table = as.data.frame(cbind(
      (meta_data_vis[,parameter_1]),
      (meta_data_vis[,parameter_2])
    ))
    
    ana_table = ana_table[ (!is.na(ana_table[,1]) &  ( ana_table[,1] != "")), ]
    ana_table = ana_table[ (!is.na(ana_table[,2]) &  ( ana_table[,2] != "")), ]
    test_statistic = table(ana_table)
    
    p_value =  chisq.test(test_statistic)$p.value
    
    if ( p_value < 0.05){
      print( c(parameter_1, parameter_2, p_value) )
    }
    result_mat_categorial = rbind(result_mat_categorial, c("categorial",parameter_1, parameter_2, p_value))
  }
}
colnames(result_mat_categorial) = c("Type","Parameter_1","Parameter_2","P_value")

###### NUMERIC

result_mat_numeric <<- matrix(as.character(), ncol = 4)
for( i in seq(1,length(parameters_double)-1) ){
  for( j in seq(i+1,length(parameters_double)) ){
    
    parameter_1 = parameters_double[i]
    parameter_2 = parameters_double[j]
    
    ana_table = as.data.frame(cbind(
      (meta_data_vis[,parameter_1]),
      (meta_data_vis[,parameter_2])
    ))
    
    ana_table = ana_table[ (!is.na(ana_table[,1]) &  ( ana_table[,1] != "")), ]
    ana_table = ana_table[ (!is.na(ana_table[,2]) &  ( ana_table[,2] != "")), ]

    p_value =  cor.test(ana_table[,1],ana_table[,2])$p.value
    correlation =  cor(ana_table[,1],ana_table[,2])
    
    if ( p_value < 0.05){
      print( c(parameter_1, parameter_2, p_value) )
    }
    result_mat_numeric = rbind(result_mat_numeric, c("numeric",parameter_1, parameter_2, p_value))
  }
}
colnames(result_mat_numeric) = c("Type","Parameter_1","Parameter_2","P_value")

###### ANOVA

result_mat_anova <<- matrix(as.character(), ncol = 4)
for( i in seq(1,length(parameters_double)) ){
  for( j in seq(1,length(parameters_double)) ){
    
    parameter_1 = parameters_categorical[i]
    parameter_2 = parameters_double[j]
    
    ana_table = as.data.frame(cbind(
      (meta_data_vis[,parameter_1]),
      (meta_data_vis[,parameter_2])
    ))
    
    ana_table = ana_table[ (!is.na(ana_table[,1]) &  ( ana_table[,1] != "")), ]
    ana_table = ana_table[ (!is.na(ana_table[,2]) &  ( ana_table[,2] != "")), ]
    
    aov_analysis =  aov(ana_table[,2] ~ as.factor(ana_table[,1]))
    p_value_results = TukeyHSD(aov_analysis)
    
    # extract th class label identities
    res = p_value_results[names(p_value_results)]
    labels = rownames(as.data.frame(res))
    
    length_res = length(unlist(p_value_results))
    indices = length_res/4
    p_values = as.double(unlist(p_value_results)[seq(length_res-indices +1,length_res)])
    
    labels_para_2 = paste(parameter_2, as.vector(labels), sep ="_")
    anova_res_mat = matrix(c(rep("anova", length(p_values)),rep(parameter_1,length_res/4), labels_para_2, p_values), nrow = length(p_values))
    result_mat_anova = rbind(result_mat_anova, anova_res_mat)
  }
}
colnames(result_mat_anova) = c("Type","Parameter_1","Parameter_2","P_value")

### output everything

res_mat = rbind(result_mat_categorial,result_mat_numeric,result_mat_anova)
#write.table(res_mat, "~/Downloads/Phenotype_correlations.tsv", sep ="\t",row.names = FALSE)
### correlation matrix

parameters_survival = c("Overall_survial","Progression_free_survial","Progression_free_survial_diagnosis","Overall_survival_diagnosis")
parameters = c(parameters_double,parameters_survival)
parameters[! parameters %in% colnames(meta_data)]
parameter_matrix = meta_data[,parameters]
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

#### correlation heatmap sampels

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
