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

#meta_data$Subtype[meta_data$Subtype %in% c("MS","BA")] = "Non_MS"

meta_data_vis = meta_data # OS_Monate_ab_Einschluss BA CL PFS_Monate_ab_Einschluss BA CL
#meta_data_vis = meta_data %>% filter(Arm_codiert == 1)
#meta_data_vis = meta_data %>% filter(Arm_codiert == 2) # PFS_Monate_ab_Einschluss BA CL
#meta_data_vis = meta_data_vis %>% filter(Subtype %in% c("BA","CL"))

# OS_ab_ED

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Overall_Survival_from_diagnosis))
fit_os = survival::survfit( survival::Surv( as.double(Overall_Survival_from_diagnosis), censor_OS ) ~ Subtype, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("black","darkgreen","blue")
)

## OS_Monate_ab_Einschluss

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Overall_Survivall_from_Randomisation))
fit_os = survival::survfit( survival::Surv( as.double(Overall_Survivall_from_Randomisation), censor_OS ) ~ Subtype, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("black","darkgreen","blue")
)

# PFS_ab_ED

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Progression_free_Survival_from_Diagnosis))
fit_os = survival::survfit( survival::Surv( Progression_free_Survival_from_Diagnosis, censor_OS ) ~ Subtype, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("black","darkgreen","blue")
)

## PFS_Monate_ab_Einschluss

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Progression_free_Survival_from_Randomisation))
fit_os = survival::survfit( survival::Surv( Progression_free_Survival_from_Randomisation, censor_OS ) ~ Subtype, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("black","darkgreen","blue")
)

##################

selector = "L1" # 0.054


selection_vector = meta_data_vis_os[,selector]

# Budding_1HPF

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Tumor_cell_budding))
threshold = threshold_ori = meta_data_vis_os$Tumor_cell_budding
threshold[ threshold_ori < 10] = "low"
threshold[ threshold_ori >= 10] = "high"

fit_os = survival::survfit( survival::Surv( as.double(Overall_Survival_from_diagnosis), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Overall_Survivall_from_Randomisation), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Diagnosis), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Randomisation), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("red","darkgreen")
)

# Zellnestgröße_zentral

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Cell_nest_size))
threshold = threshold_ori = meta_data_vis_os$Cell_nest_size
threshold[ threshold_ori < 10] = "low" # 10
threshold[ threshold_ori >= 10] = "high" # 10

fit_os = survival::survfit( survival::Surv( as.double(Overall_Survival_from_diagnosis), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Overall_Survivall_from_Randomisation), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Diagnosis), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Randomisation), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval

# Mitosen_HPF

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Mitotic_Count))
threshold = threshold_ori = meta_data_vis_os$Mitotic_Count
threshold[ threshold_ori < median(threshold_ori)] = "low" # 2
threshold[ threshold_ori >= median(threshold_ori)] = "high" # 2

fit_os = survival::survfit( survival::Surv( as.double(Overall_Survival_from_diagnosis), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Overall_Survivall_from_Randomisation), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Diagnosis), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Randomisation), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval

# Kerngröße

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Nuclear_Size))
threshold = threshold_ori = meta_data_vis_os$Nuclear_Size
threshold[ threshold_ori < 2] = "small"
threshold[ threshold_ori >= 2] = "big"

fit_os = survival::survfit( survival::Surv( as.double(Overall_Survival_from_diagnosis), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Overall_Survivall_from_Randomisation), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("red","darkgreen")
)

fit_os = survival::survfit( survival::Surv( as.double(Overall_Survival_from_diagnosis), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Overall_Survivall_from_Randomisation), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Diagnosis), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Randomisation), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval

# Stroma_vitalerTumor

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Stroma_vitalerTumor))
threshold = threshold_ori = meta_data_vis_os$Stroma_vitalerTumor
threshold[ threshold_ori < 33] = "low"
threshold[ threshold_ori >= 33] = "high"

fit_os = survival::survfit( survival::Surv( as.double(Overall_Survival_from_diagnosis), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Overall_Survivall_from_Randomisation), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Diagnosis), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Randomisation), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("red","darkgreen")
)

# Nekrose

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Necrosis))
threshold = threshold_ori = meta_data_vis_os$Necrosis_ROC
threshold[ threshold_ori < 4] = "low"
threshold[ threshold_ori >= 4] = "high"

fit_os = survival::survfit( survival::Surv( as.double(Overall_Survival_from_diagnosis), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Overall_Survivall_from_Randomisation), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Diagnosis), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Randomisation), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval

# Entzündung

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Entzündung))
threshold = threshold_ori = meta_data_vis_os$Entzündung
threshold[ threshold_ori < 2] = "low"
threshold[ threshold_ori >= 2] = "high"

fit_os = survival::survfit( survival::Surv( as.double(Overall_Survival_from_diagnosis), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Overall_Survivall_from_Randomisation), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Diagnosis), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Randomisation), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("red","darkgreen")
)

# L1

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Lymphangiosis))
threshold = threshold_ori = meta_data_vis_os$Lymphangiosis

fit_os = survival::survfit( survival::Surv( as.double(Overall_Survival_from_diagnosis), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Overall_Survivall_from_Randomisation), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Diagnosis), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(Progression_free_Survival_from_Randomisation), censor_PFS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
