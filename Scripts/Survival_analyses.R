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

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$SampleID
sum(meta_info$Survivalstatistik_neuestes_Sample)

#meta_data = meta_info[colnames(expr_raw),]

# Overall_Survival_from_diagnosis
# Overall_Survivall_from_Randomisation
# Progression_free_Survivall_from_Diagnosis
# Progression_free_Survivall_from_Randomisation
survival_selector = "Progression_free_Survivall_from_Randomisation"
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

meta_data = meta_info[colnames(expr_raw),]
meta_data = meta_data %>% filter(as.character(Survivalstatistik_neuestes_Sample) == "1")
length(meta_data$Survivalstatistik_neuestes_Sample)
dim(meta_data)

meta_data_vis = meta_data %>% filter(!is.na(meta_data$PFS_Monate_ab_Einschluss))
#meta_data_vis = meta_data %>% filter(!is.na(meta_data$OS_Monate_ab_Einschluss))

##

selector = "Budding_1HPF"
meta_data_vis_pfs = meta_data_vis[!is.na( meta_data_vis[,selector]),]
schwellert = meta_data_vis_pfs[,selector]
schwellert[schwellert %in% c(0,1)] = 0
schwellert[schwellert > 1 ] = 1

schwellert = meta_data_vis_pfs$Stroma_vitalerTumor
schwellert[schwellert < 31] = 0
schwellert[schwellert >= 31 ] = 1

selector = "Grading" # 0.034
selector = "Budding_1HPF_ROC"
selector = "Budding_10HPF_ROC"
selector = "Keratinisierung"
selector = "Subtype"
selector = "Zellnestgröße_zentral"
selector = "Zellnestgröße_ROC"
selector = "Mitosen_HPF_ROC"
selector = "Stroma_ROC"
selector = "Nekrose_ROC"
selector = "Entzündung_ROC"
selector = "L1"
selector = "Pn1"
meta_data_vis_pfs = meta_data_vis[!is.na( meta_data_vis_pfs[,selector]),]
schwellert = meta_data_vis_pfs[,selector]

fit = survival::survfit( survival::Surv( as.double(meta_data_vis_pfs$PFS_ab_ED), censor_PFS ) ~ schwellert, data = meta_data_vis_pfs)
survminer::surv_pvalue(fit, data = meta_data_vis_pfs)$pval
fit = survival::survfit( survival::Surv( as.double(meta_data_vis_pfs$PFS_Monate_ab_Einschluss), censor_PFS ) ~ schwellert, data = meta_data_vis_pfs)
survminer::surv_pvalue(fit, data = meta_data_vis_pfs)$pval
fit = survival::survfit( survival::Surv( as.double(meta_data_vis_pfs$OS_ab_ED), censor_PFS ) ~ schwellert, data = meta_data_vis_pfs)
survminer::surv_pvalue(fit, data = meta_data_vis_pfs)$pval
fit = survival::survfit( survival::Surv( as.double(meta_data_vis_pfs$OS_Monate_ab_Einschluss), censor_PFS ) ~ schwellert, data = meta_data_vis_pfs)
survminer::surv_pvalue(fit, data = meta_data_vis_pfs)$pval

print(survminer::ggsurvplot(fit_pfs, data = meta_data_vis_pfs, risk.table = F, pval = T, censor.size = 10))

#svg(filename = "~/Deko_Projekt/Results/Images/Figure_5_survival_ductal_three.svg", width = 10, height = 10)

#############

# OS_Monate_ab_Einschluss OS_ab_ED OS_Monate

meta_data_vis = meta_data

### OS_ab_ED

meta_data_vis_os = meta_data_vis %>% filter(!is.na(OS_ab_ED))

schwellert = meta_data_vis_os$Budding_1HPF
schwellert[schwellert %in% c(0,1)] = 0
schwellert[schwellert > 1 ] = 1

schwellert = meta_data_vis_os$Stroma_vitalerTumor
schwellert[schwellert < 31] = 0
schwellert[schwellert >= 31 ] = 1

####

selector = "Grading" # 0.034
#selector = "Budding_1HPF"
selector = "Budding_1HPF_ROC"
#selector = "Budding_10HPF"
selector = "Budding_10HPF_ROC"
#selector = "Keratinisierung"
selector = "Subtype"
#selector = "Zellnestgröße_zentral"
selector = "Zellnestgröße_ROC"
#selector = "Mitosen_HPF"
selector = "Mitosen_HPF_ROC"
selector = "Stroma_ROC"
#selector = "Nekrose"
selector = "Nekrose_ROC"
#selector = "Entzündung"
selector = "Entzündung_ROC"
selector = "L1"
selector = "Pn1"


meta_data_vis_os = meta_data_vis %>% filter(!is.na(selector))
meta_data_vis_os$schwellert = meta_data_vis_os[,selector]
selector

#fit_os = survival::survfit( survival::Surv( as.double(meta_data_vis_os$OS_ab_ED), meta_data_vis_os$censor_OS ) ~ meta_data_vis_os$schwellert)
#fit_os = survival::survfit( survival::Surv( as.double(meta_data_vis_os$OS_Monate_ab_Einschluss), meta_data_vis_os$censor_OS ) ~ meta_data_vis_os$schwellert)
#fit_os = survival::survfit( survival::Surv( as.double(meta_data_vis_os$PFS_ab_ED), meta_data_vis_os$censor_OS ) ~ meta_data_vis_os$schwellert)
#fit_os = survival::survfit( survival::Surv( as.double(meta_data_vis_os$PFS_Monate_ab_Einschluss), meta_data_vis_os$censor_OS ) ~ meta_data_vis_os$schwellert)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval

#pdf("~/Downloads/L1.OS_Monate_ab_Einschluss.pdf")
print(survminer::ggsurvplot(fit_os, data = meta_data_vis_os, risk.table = F, pval = T, censor.size = 10))
#dev.off()

## Schwellwert

schwellert = meta_data_vis_os$Budding_1HPF
schwellert[schwellert %in% c(0,1)] = 0
schwellert[schwellert > 1 ] = 1
meta_data_vis_os$schwellert = schwellert

schwellert = meta_data_vis_os$Stroma_vitalerTumor
schwellert[schwellert < 31] = 0
schwellert[schwellert >= 31 ] = 1
meta_data_vis_os$schwellert = schwellert

#write.table(meta_info,"~/Koop_Klinghammer/Misc/Meta_information.tsv", quote = FALSE, sep = "\t", row.names = FALSE)
