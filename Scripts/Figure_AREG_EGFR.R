library("ggplot2")
library("stringr")
library("grid")
library("umap")

source("~/Koop_Klinghammer/Scripts/Visualization_colors.R")

expr_raw = read.table(
  "~/Koop_Klinghammer/Data/S103.tsv",
  sep ="\t",
  stringsAsFactors = F,
  header = T,
  row.names = 1
)
colnames(expr_raw ) = str_replace_all(colnames(expr_raw), pattern = "^X", "")

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Sample_ID
matcher = match( colnames(expr_raw), meta_info$Sample_ID, nomatch = 0)
dim(expr_raw)

## Figure AREG EGRF
meta_data = meta_info[matcher,]
candidates = meta_data$Sample_ID[! is.na(meta_data$Subtype)]#[which(meta_data$Survivalstatistik_neuestes_Sample == 1)]
length(candidates)

expr = expr_raw[c("EGFR","AREG"),match(candidates, colnames(expr_raw))]
vis_mat = reshape2::melt(t(expr))
colnames(vis_mat) = c("Sample","Gene","Value")
Subtype = meta_data[ match(vis_mat$Sample, meta_data$Sample_ID),"Subtype"]
vis_mat = cbind(Subtype,vis_mat)
vis_mat$Gene = factor(vis_mat$Gene, levels = c("AREG","EGFR"))

p = ggplot(
  vis_mat,
  aes( x = Gene, y = Value, fill = Subtype )
)
p = p + geom_boxplot(notch = FALSE,outlier.colour = "red", outlier.shape = 1)
p = p + scale_fill_manual(values = c("black","darkgreen","blue"))
p = p + ylab("Expression") + theme(legend.position = "top")

#svg(filename = "~/Downloads/AREG_EGFR_plot.svg", width = 10, height = 10)
print(p)
dev.off()

###

areg_ms = vis_mat[(vis_mat$Subtype == "MS")&(vis_mat$Gene == "AREG"),"Value"]
areg_cl = vis_mat[(vis_mat$Subtype == "CL")&(vis_mat$Gene == "AREG"),"Value"]
areg_ba = vis_mat[(vis_mat$Subtype == "BA")&(vis_mat$Gene == "AREG"),"Value"]

t.test(areg_ms, areg_cl)$p.value # 0.03
t.test(areg_ms, areg_ba)$p.value # 1.8E-12
t.test(areg_cl, areg_ba)$p.value # 4E-4

#

egfr_ms = vis_mat[(vis_mat$Subtype == "MS")&(vis_mat$Gene == "EGFR"),"Value"]
egfr_cl = vis_mat[(vis_mat$Subtype == "CL")&(vis_mat$Gene == "EGFR"),"Value"]
egfr_ba = vis_mat[(vis_mat$Subtype == "BA")&(vis_mat$Gene == "EGFR"),"Value"]

t.test(egfr_ms, egfr_cl)$p.value # 0.03
t.test(egfr_ms, egfr_ba)$p.value # 1.8E-12
t.test(egfr_cl, egfr_ba)$p.value # 4E-4

#
