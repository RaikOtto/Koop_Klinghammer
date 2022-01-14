library("ggplot2")
library("stringr")
#library("umap")


expr_raw = read.table(
  "~/Koop_Klinghammer/Data/Data.S104.tsv",
  sep ="\t",
  stringsAsFactors = F,
  header = T
  #row.names = F
)
colnames(expr_raw ) = str_replace_all(colnames(expr_raw), pattern = "^X", "")
expr_raw[1:5,1:5]
dim(expr_raw)

### Prep

## Figure 1

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Sample_ID

meta_data = meta_info[colnames(expr_raw),]