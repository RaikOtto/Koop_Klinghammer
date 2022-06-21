library("ggplot2")
library("dplyr")
library("stringr")
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

meta_info = read.table("~/Koop_Klinghammer/Misc/Clinical_metadata_survival.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$SampleID
sum(meta_info$Survivalstatistik_neuestes_Sample)

#meta_data = meta_info[colnames(expr_raw),]

# Overall_Survival_from_diagnosis
# Overall_Survivall_from_Randomisation
# Progression_free_Survivall_from_Diagnosis
# Progression_free_Survivall_from_Randomisation
survival_selector = "Overall_Survival_from_diagnosis"
dim(meta_info)

##

selector = "Keratinization"
selectors = c("Subtype","Histotyp","Grading_WHO","Keratinization","Tumor_cell_budding_HPF_2_tier","Tumor_cell_budding_10HPF_2_tier","Cell_nest_size_2_tier","Mitotic_Count_HPF_2_tier","Mitotic_Count_10HPF_2_tier","Nuclear_Size_2_tier","Stroma_Content_2_tier","Necrosis_2_tier","Inflammatory_Infiltrate_2_tier","Lymphangiosis","Perineural_Invasion","Best_Response")
#selectors[which(! selectors %in% colnames(meta_data_survival))]

res_vec <<- matrix(as.character(),ncol = 3)
for (selector in selectors){
  
  phenotype_vector = meta_info[,selector]
  survival_vector = as.double(meta_data_survival_stat[,survival_selector])

  if (survival_selector %in% c("Overall_Survival_from_diagnosis","Overall_Survivall_from_Randomisation") ) {
    censor_vector = meta_data_survival_stat$censor_OS
  } else {
    censor_vector = meta_data_survival_stat$censor_PFS  
  }
  
  length(phenotype_vector)
  length(survival_vector)
  length(censor_vector)
  
  fit = survival::survfit( 
    survival::Surv( 
      time = survival_vector,
      censor_vector
    ) ~ phenotype_vector
  )
  p_value = survminer::surv_pvalue(fit, data = meta_data_survival_stat)$pval
  
  if( p_value < 0.05){
    print(c(selector,round(p_value,4), survival_selector,nrow(meta_data_survival_stat)))
    res_vec = rbind(res_vec,c(selector,round(p_value,4), survival_selector))
  }
}

data = data.frame(
  "Sample_ID" = meta_data_survival_stat$Sample_ID,
  "Phenotype" = phenotype_vector,
  "Survival" = survival_vector,
  "Censor" = censor_vector
)
#write.table(res_vec,"~/Downloads/Progression_free_Survivall_from_Randomisation.tsv", sep = "\t", row.names = FALSE, quote = F)

print(survminer::ggsurvplot(fit, data = meta_data_survival_stat, risk.table = T, pval = T, censor.size = 10))

#write.table(meta_info,"~/Koop_Klinghammer/Misc/Meta_information.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
