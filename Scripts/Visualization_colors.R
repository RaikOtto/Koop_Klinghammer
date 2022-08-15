aka3 = list(
  Group = c(refractive = "red", sensitive = "darkgreen", intermediate = "orange"),
  Subtype = c(BA = "black", CL = "darkgreen", MS = "blue"),
  OS = c(high = "red", medium = "orange", low = "green"),
  Budding_absentpresent = c("2" = "black", "1" = "white"),
  Zellnestgröße_zentral = c(low = "white", high = "black"),
  Mitosen_10HPF = c(low = "white", high = "black"),
  Nekrose = c(low = "white", high = "black"),
  Entzündung = c(low = "white", high = "black"),
  Budding_10HPF = c(low = "white", high = "black"),
  Keratinisierung = c(low = "white", high = "black"),
  Grading = c( G1 = "#d1e6e6",G2 = "#78b8b3", G3 = "#f54b17", "Unknown" = "gray"),
  survival_vec = c("Low"="darkgreen","High" = "red"),
  PFS_Monate_ab_Einschluss = c(high="black",low = "white"),
  OS_ab_ED = c(high="black",low = "white"),
  #Entzündung_ROC = c("2" = "#efad16", "1" = "#0a7e8c","Unknown" = "#aaa9ad")
  Entzündung_ROC = c("2" = "black", "1" = "white","Unknown" = "#aaa9ad"),
  Budding_1_hpf_threshold= c("1" = "black", "0" = "white","Unknown" = "gray"),
  Lymphangiosis= c("1" = "black", "0" = "white","Unknown" = "gray"),
  Inflammatory_Infiltrate_ROC= c("high" = "black", "low" = "white","Unknown" = "gray"),
  Necrosis_ROC= c("high" = "black", "low" = "white","Unknown" = "gray")
)
