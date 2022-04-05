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

meta_info = read.table("~/Koop_Klinghammer/Misc/Masterliste.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$SampleID

meta_data = meta_info[colnames(expr_raw),]
meta_data = meta_data %>% filter(as.character(Survivalstatistik_neuestes_Sample) == "1")
length(meta_data$Survivalstatistik_neuestes_Sample)
dim(meta_data)

meta_data_vis = meta_data %>% filter(Arm_codiert == 1)
meta_data_vis = meta_data %>% filter(Arm_codiert == 2)
meta_data_vis = meta_data_vis %>% filter(Subtype %in% c("MS","CL"))

# OS_Monate_ab_Einschluss OS_ab_ED OS_Monate

### OS_ab_ED

meta_data_vis_os = meta_data_vis %>% filter(!is.na(OS_Monate))
meta_data_vis_os = meta_data_vis %>% filter(!is.na(PFS_Monate_ab_Einschluss))

fit_os = survival::survfit( survival::Surv( as.double(OS_ab_ED), censor_OS ) ~ Subtype, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
print(survminer::ggsurvplot(fit_os, data = meta_data_vis_os, risk.table = F, pval = T, censor.size = 10))

## OS_Monate_ab_Einschluss

fit_os = survival::survfit( survival::Surv( as.double(OS_Monate_ab_Einschluss), censor_OS ) ~ Subtype, data = meta_data_vis_os)
survminer::surv_pvalue(fit_os, data = meta_data_vis_os)$pval
print(survminer::ggsurvplot(fit_os, data = meta_data_vis_os, risk.table = F, pval = T, censor.size = 10))

# PFS_ab_ED PFS_Monate_ab_Einschluss PFS_Monate

meta_data_vis_pfs = meta_data_vis %>% filter(!is.na(meta_data_vis$PFS_Monate_ab_Einschluss))

## PFS_ab_ED

fit_pfs = survival::survfit( survival::Surv( as.double(PFS_ab_ED), censor_PFS ) ~ Subtype, data = meta_data_vis_pfs)
survminer::surv_pvalue(fit_pfs, data = meta_data_vis_pfs)$pval
print(survminer::ggsurvplot(fit_pfs, data = meta_data_vis_pfs, risk.table = F, pval = T, censor.size = 10))

## PFS_Monate_ab_Einschluss

schwellert = meta_data_vis_pfs$Budding_1HPF
schwellert[schwellert %in% c(0,1)] = 0
schwellert[schwellert > 1 ] = 1

schwellert = meta_data_vis_pfs$Stroma_vitalerTumor
schwellert[schwellert < 31] = 0
schwellert[schwellert >= 31 ] = 1

fit_pfs = survival::survfit( survival::Surv( as.double(meta_data_vis_pfs$PFS_Monate_ab_Einschluss), censor_PFS ) ~ schwellert, data = meta_data_vis_pfs)
survminer::surv_pvalue(fit_pfs, data = meta_data_vis_pfs)$pval
print(survminer::ggsurvplot(fit_pfs, data = meta_data_vis_pfs, risk.table = F, pval = T, censor.size = 10))

#svg(filename = "~/Deko_Projekt/Results/Images/Figure_5_survival_ductal_three.svg", width = 10, height = 10)
#print(survminer::ggsurvplot(fit, data = ratio_m, risk.table = F, pval = T, censor.size = 10))
#dev.off()