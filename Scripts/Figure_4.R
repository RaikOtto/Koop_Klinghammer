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
rownames(meta_info) = meta_info$Sample_ID
matcher = match(colnames(expr_raw),meta_info$Sample_ID, nomatch = 0)

meta_data = meta_info[matcher,]
meta_data = meta_data %>% filter(Survivalstatistik_neuestes_Sample == 1)
table(meta_data$Survivalstatistik_neuestes_Sample)

#meta_data = meta_data %>%  mutate( Subtype = case_when( Subtype %in% c("BA","MS") ~ "BA_MS", Subtype == "CL" ~ "CL"))
meta_data = meta_data %>%  filter( Subtype == "CL")

##########

survival_vectors =  c("Overall_Survival_from_diagnosis","Overall_Survival_from_Randomisation","Progression_free_Survival_from_Diagnosis","Progression_free_Survival_from_Randomisation")
#survival_categories = c("Subtype","Histotyp","Grading","Keratinization","Tumor_cell_budding_ROC","Cell_nest_size_ROC","Mitotic_Count_ROC","Nuclear_Size_ROC","Necrosis_ROC","Inflammatory_Infiltrate_ROC","Lymphangiosis","Perineural_Invasion","Localization_primary_tumor")
survival_categories = c("Histotyp","Grading","Keratinization","Tumor_cell_budding_ROC","Cell_nest_size_ROC","Mitotic_Count_ROC","Nuclear_Size_ROC","Necrosis_ROC","Inflammatory_Infiltrate_ROC","Lymphangiosis","Perineural_Invasion","Localization_primary_tumor")
overall_results <<- matrix(character(), ncol = 4 )
roc_types = c("Cutoff_Finder","Mean","Median")

for (survival_type in survival_vectors){
  
  if( str_detect(survival_type, "Overall") ) censor_type = "censor_OS" else censor_type = "censor_PFS"
  
  for ( survival_category in survival_categories){
    
    meta_data_vis = meta_data[,c(survival_type,survival_category,censor_type,str_replace(survival_category, "_ROC",""))]
    meta_data_vis = meta_data_vis %>% filter(meta_data_vis[,survival_type] != "")
    meta_data_vis = meta_data_vis %>% filter(meta_data_vis[,survival_category] != "")
    
    survival_time = as.double(as.character(meta_data_vis[,survival_type]))
    censor_vector = as.integer(meta_data_vis[,censor_type])
    
    ###
    
    for ( roc_type in roc_types){
      
      if (roc_type == "Cutoff_Finder") {
        feature_type = as.factor(as.character(meta_data_vis[,survival_category]))
      } else if (str_detect(survival_category, "_ROC")){
        
        feature_type = values_original = meta_data_vis[,str_replace(survival_category, "_ROC","")]
        if ( roc_type == "Mean") Threshold = mean(values_original) else Threshold = median(values_original)
        feature_type[ feature_type > Threshold] = "high"
        feature_type[ feature_type <= Threshold] = "low"
      }
      print(c(survival_category, survival_type, roc_type))
      fit_os = survival::survfit(
        survival::Surv( 
          survival_time,
          censor_vector
        ) ~ feature_type)
      p_value = as.double(survminer::surv_pvalue(fit_os, data = meta_data_vis)$pval)
      results = c(survival_type, as.character(survival_category), roc_type, as.character(p_value))
      
      if ( is.na(p_value) ) next()
      if (p_value < 0.05) print(results)
      
      overall_results = rbind(overall_results, results)
    }
  }
}
colnames(overall_results) = c("Survival_time_type","Feature","ROC_type","P_value")

openxlsx::write.xlsx(as.data.frame(overall_results), "~/Downloads/survival_cl.xlsx", asTable = TRUE, overwrite = TRUE)
