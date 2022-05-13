library("stringr")
library("limma")
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
expr_raw[1:5,1:5]
dim(expr_raw)

## Figure 1
meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Sample_ID

cor_mat = cor(expr_raw);pcr = prcomp(t(cor_mat))
matcher = match(as.character(colnames(cor_mat)), as.character(meta_info$SampleID), nomatch = 0)
meta_data = meta_info[matcher,]
rownames(meta_data) = meta_data$SampleID

design <- model.matrix(~0 + meta_data$Subtype)
#design <- model.matrix(~0 + as.factor(meta_data$Subtype))
colnames(design) = c("BA","CL","MS")

#vfit <- lmFit(t(r_mat[,c(-1,-2)]),design)
vfit <- lmFit(expr_raw,design)
#vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
#plotSA(efit)

#contr.matrix = makeContrasts( contrast = CL - MS,  levels = design )
contr.matrix = makeContrasts( contrast = BA - CL ,  levels = design )
vfit <- contrasts.fit( vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)

summary(decideTests(efit))

result_t = topTable( efit, coef = "contrast", number  = nrow(expr_raw), adjust  ="none", p.value = 1, lfc = 0)
result_t$hgnc_symbol = rownames(result_t)
colnames(result_t) = c("Log_FC","Average_Expr","t","P_value","adj_P_value","B","HGNC")

result_t = result_t[c("HGNC","Log_FC","Average_Expr","P_value","adj_P_value")]
result_t = result_t[order(result_t$P_value, decreasing = F),]
result_t$Log_FC = round(result_t$Log_FC, 1)
result_t$Average_Expr = round(result_t$Average_Expr, 1)
result_t = result_t[order(result_t$Log_FC,decreasing = T),]

#write.table("~/Koop_Klinghammer/Results/05_06_2018/CL_minus_MS.tsv", x = result_t, sep = "\t", quote = F, row.names = F)

###

