library("ggplot2")
library("stringr")
library("grid")
library("umap")
library("dplyr")


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

meta_info = read.table("~/Koop_Klinghammer/Misc/meta_info_table_klinghammer.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$SampleID
matcher = match(colnames(expr_raw),meta_info$SampleID, nomatch = 0)

meta_data = meta_info[matcher,]
#meta_data = meta_data %>% filter(Histotyp != "")
#meta_data = meta_data %>% filter(! is.na(Disease_control_rate ))

# CATEGORICAL

parameters_categorical = c("Grading","Histotyp","Budding_1HPC_ROC","Minimal_cell_nest_size_ROC","Mitotic_activity_1_HPF_ROC","Nuclear_size_ROC","Stroma_Tumor_Ratio_ROC","Necrosis_ROC","Inflammatory_infiltration_ROC","Lymphangiosis","Perineural_invasion","Disease_control_rate","Best_response","Localisation_primary_tumor")

meta_data_vis = meta_data[meta_data$Subtype %in% c("BA","CL"),]
dim(meta_data_vis) # BA 39

for( parameter in parameters_categorical ){
  
  print(parameter)
  
  ana_table = as.data.frame(cbind(
    (meta_data_vis[,"Subtype"]),
    (meta_data_vis[,parameter])
  ))
  ana_table = ana_table[!is.na(ana_table[,2]),]
  ana_table = ana_table[ ana_table[,2] != "" ,]
  test_statistic = table(ana_table)
  
  suppressWarnings(
    print( chisq.test(test_statistic)$p.value )
  )
}

###### NUMERIC

parameters_double = c("Keratinization","Budding_1HPF","Minimal_cell_nest_size","Mitotic_activity_1_HPF","Nuclear_size","Stroma_Tumor_Ratio","Necrosis","Inflammatory_infiltration","Age_at_randomisation")
Subtype_1 = c("MS","CL")
Subtype_2 = c("BA")
dim(meta_data_vis) # BA 39

for( parameter in parameters_double ){
  
  ana_table = as.data.frame(cbind(
    (meta_data[,"Subtype"]),
    (meta_data[,parameter])
  ))
  ana_table = ana_table[!is.na(ana_table[,2]),]
  ana_table = ana_table[ ana_table[,2] != "" ,]
  subtypte_data_1 = as.double(ana_table[ana_table[,1] == Subtype_1,2])
  subtypte_data_2 = as.double(ana_table[ana_table[,1] == Subtype_2,2])
  
  result = t.test(subtypte_data_1, subtypte_data_2)$p.value
  
  if ( result < 0.05){
      print(c(Subtype_1,Subtype_2,parameter, round(result,2)))
  }
}

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
