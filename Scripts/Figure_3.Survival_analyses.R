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
#meta_data_vis = meta_data %>% filter(Arm_codiert == 1)
#meta_data_vis = meta_data %>% filter(Arm_codiert == 2) # PFS_Monate_ab_Einschluss BA CL
#meta_data_vis = meta_data_vis %>% filter(Subtype %in% c("BA","CL"))

# OS_ab_ED

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$OS_ab_ED))
fit_os = survival::survfit( survival::Surv( as.double(OS_ab_ED), censor_OS ) ~ selection_vector, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("black","darkgreen","blue")
)

## OS_Monate_ab_Einschluss

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$OS_Monate_ab_Einschluss))
fit_os = survival::survfit( survival::Surv( as.double(OS_Monate_ab_Einschluss), censor_OS ) ~ selection_vector, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("black","darkgreen","blue")
)

# PFS_ab_ED

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$PFS_ab_ED))
fit_os = survival::survfit( survival::Surv( PFS_ab_ED, censor_OS ) ~ selection_vector, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("black","darkgreen","blue")
)

## OS_Monate_ab_Einschluss

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$PFS_Monate_ab_Einschluss))
fit_os = survival::survfit( survival::Surv( PFS_Monate_ab_Einschluss, censor_OS ) ~ selection_vector, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("black","darkgreen","blue")
)

##################

selector = "Entzündung" 

selector = "L1" # 0.054


selection_vector = meta_data_vis_os[,selector]

# Budding_1HPF

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Budding_1HPF))
threshold = threshold_ori = meta_data_vis_os$Budding_1HPF
threshold[ threshold_ori < 2] = "low"
threshold[ threshold_ori >= 2] = "high"

fit_os = survival::survfit( survival::Surv( as.double(OS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(OS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("red","darkgreen")
)

# Zellnestgröße_zentral

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Zellnestgröße_zentral))
threshold = threshold_ori = meta_data_vis_os$Zellnestgröße_zentral
threshold[ threshold_ori < 3] = "low"
threshold[ threshold_ori >= 3] = "high"

fit_os = survival::survfit( survival::Surv( as.double(OS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(OS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval

# Mitosen_HPF

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Mitosen_HPF))
threshold = threshold_ori = meta_data_vis_os$Mitosen_HPF
threshold[ threshold_ori < 4] = "low"
threshold[ threshold_ori >= 4] = "high"

fit_os = survival::survfit( survival::Surv( as.double(OS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(OS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval

# Kerngröße

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Kerngröße))
threshold = threshold_ori = meta_data_vis_os$Kerngröße
threshold[ threshold_ori < 2] = "small"
threshold[ threshold_ori >= 2] = "big"

fit_os = survival::survfit( survival::Surv( as.double(OS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(OS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("red","darkgreen")
)

fit_os = survival::survfit( survival::Surv( as.double(PFS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval

# Stroma_vitalerTumor

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Stroma_vitalerTumor))
threshold = threshold_ori = meta_data_vis_os$Stroma_vitalerTumor
threshold[ threshold_ori < 33] = "low"
threshold[ threshold_ori >= 33] = "high"

fit_os = survival::survfit( survival::Surv( as.double(OS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(OS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("red","darkgreen")
)

# Nekrose

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Nekrose))
threshold = threshold_ori = meta_data_vis_os$Nekrose
threshold[ threshold_ori < 4] = "low"
threshold[ threshold_ori >= 4] = "high"

fit_os = survival::survfit( survival::Surv( as.double(OS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(OS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval

# Entzündung

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$Entzündung))
threshold = threshold_ori = meta_data_vis_os$Entzündung
threshold[ threshold_ori < 2] = "low"
threshold[ threshold_ori >= 2] = "high"

fit_os = survival::survfit( survival::Surv( as.double(OS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(OS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("red","darkgreen")
)

# L1

meta_data_vis_os = meta_data_vis %>% filter(!is.na(meta_data_vis$L1))
threshold = threshold_ori = meta_data_vis_os$L1

fit_os = survival::survfit( survival::Surv( as.double(OS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(OS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("red","darkgreen")
)

fit_os = survival::survfit( survival::Surv( as.double(PFS_ab_ED), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
fit_os = survival::survfit( survival::Surv( as.double(PFS_Monate_ab_Einschluss), censor_OS ) ~ threshold, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T,
  palette = c("red","darkgreen")
)
