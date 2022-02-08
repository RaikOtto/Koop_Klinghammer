library("ggplot2")
library("dplyr")
library("stringr")
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
matcher = match(colnames(expr_raw),meta_info$SampleID)

meta_data = meta_info[matcher,]
meta_data = meta_data %>% filter(Survivalstatistik_neuestes_Sample == 1)
table(meta_data$Survivalstatistik_neuestes_Sample)

meta_data_vis = meta_data # OS_Monate_ab_Einschluss BA CL PFS_Monate_ab_Einschluss BA CL
meta_data_vis = meta_data %>% filter(Arm_codiert == 1)
meta_data_vis = meta_data %>% filter(Arm_codiert == 2) # PFS_Monate_ab_Einschluss BA CL
#meta_data_vis = meta_data_vis %>% filter(Subtype %in% c("BA","CL"))

# OS_Monate_ab_Einschluss OS_ab_ED OS_Monate

### OS_ab_ED

meta_data_vis_os = meta_data_vis %>% filter(!is.na(OS_Monate_ab_Einschluss))
selector = "Subtype" # 
selector = "Budding_10HPF_ROC" # 0.054
selector = "Zellnestgröße_ROC" #
selector = "Mitosen_HPF_ROC" # 0.055 
selector = "Mitosen_10HPF_ROC" # 0.07 0.052
selector = "Stroma_ROC" #
selector = "Nekrose_ROC" #
selector = "Entzündung_ROC" # 0.018 0.014

selection_vector = meta_data_vis_os[,selector]

fit_os = survival::survfit( survival::Surv( as.double(OS_ab_ED), censor_OS ) ~ selection_vector, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
#print(survminer::ggsurvplot(fit_os, data = meta_data_vis_os, risk.table = F, pval = T, censor.size = 10,  palette = c("black","darkgreen","blue")))

## OS_Monate_ab_Einschluss

selection_vector = meta_data_vis_os[,selector]

fit_os = survival::survfit( survival::Surv( as.double(OS_Monate_ab_Einschluss), censor_OS ) ~ selection_vector, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
surv_plot = survminer::ggsurvplot(
  fit_os,
  data = meta_data_vis_os,
  risk.table = F,
  pval = T#,
  #palette = c("black","darkgreen","blue")
)

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_3.Candidate.OS_Monate_ab_Einschluss.Arm_AB.BA_CL.svg", width = 10, height = 10)
#print(surv_plot)
#dev.off()

# PFS_ab_ED PFS_Monate_ab_Einschluss PFS_Monate

meta_data_vis_pfs = meta_data_vis %>% filter(!is.na(PFS_Monate_ab_Einschluss))
selection_vector = meta_data_vis_pfs[,selector]
  
## PFS_ab_ED

fit_pfs = survival::survfit( survival::Surv( as.double(PFS_ab_ED), censor_PFS ) ~ selection_vector, data = meta_data_vis_pfs)
survminer::surv_pvalue(fit_pfs, data = meta_data_vis_pfs)$pval
#print(survminer::ggsurvplot(fit_pfs, data = meta_data_vis_pfs, risk.table = F, pval = T, censor.size = 10,  palette = c("black","darkgreen","blue")))

## PFS_Monate_ab_Einschluss | Arm 2 Entzündung ROC

fit_pfs = survival::survfit( survival::Surv( as.double(PFS_Monate_ab_Einschluss), censor_PFS ) ~ selection_vector, data = meta_data_vis_pfs)
survminer::surv_pvalue(fit_pfs, data = meta_data_vis_pfs)$pval

#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_3.Candidate.PFS_Monate_ab_Einschluss.Arm_AB.BA_CL_MS.svg", width = 10, height = 10)
#svg(filename = "~/Koop_Klinghammer/Results/Figures/Figure_3.Candidate.PFS_Monate_ab_Einschluss.Arm_AB.BA_CL.svg", width = 10, height = 10)
#print(survminer::ggsurvplot(fit_pfs, data = meta_data_vis_pfs, risk.table = F, pval = T, censor.size = 10,  palette = c("black","darkgreen","blue")))
#dev.off()

