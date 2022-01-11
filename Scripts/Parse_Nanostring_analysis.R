require('NanoStringNorm')
library("stringr")

meta_info = read.table("~/Koop_Klinghammer/Misc/Datensatz_CEFCID.tsv",sep="\t",header =T, stringsAsFactors = F)
meta_info$Raw_name = str_replace_all(meta_info$Raw_name, pattern = "-", "_")

rcc_files <- dir( "~/Koop_Klinghammer/Data/Raw_data_new/", full.names = TRUE, pattern = "*")
length(rcc_files)
raw_data = NanoStringNorm::read.markup.RCC("~/Koop_Klinghammer/Data/Raw_data_new/" , rcc.pattern =  "*.RCC")

#norm.comp.results.test = norm.comp(raw_data, verbose = T)
#sum_mean.2sd_none_vsn

normalized_data = NanoStringNorm::NanoStringNorm( 
  raw_data,
  CodeCount.methods = "sum",
  Background.methods = "mean.2sd",
  SampleContent.methods = "none",
  OtherNorm.methods = "vsn",
  take.log = T,
  round.values = T,
  return.matrix.of.endogenous.probes = F,
  verbose = TRUE
)
### normalization

source_mat = normalized_data$normalized.data
m    = matrix( as.character(unlist( source_mat)), nrow=  dim(source_mat)[1], ncol = dim(source_mat)[2])
info = m[,seq(3)]
data = matrix( as.double(m[,-seq(3)]), nrow=  dim(source_mat)[1], ncol = dim(source_mat)[2]-3)
data = round(data,1)
dim(data)

rownames(data) = rownames(source_mat)
col_labels = str_replace( colnames(source_mat)[-seq(3)], pattern = "^X", "") 
colnames(data) = col_labels

matcher = match(colnames(data), meta_info$Raw_name, nomatch = 0)
colnames(data) = meta_info$Sample_ID[matcher]
#write.table(data, "~/Koop_Klinghammer/Data/New_data.S130.tsv", quote= F, row.names = T, sep = "\t")

