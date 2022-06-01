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

