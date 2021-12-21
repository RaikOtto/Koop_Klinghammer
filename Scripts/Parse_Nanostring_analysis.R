require('NanoStringNorm')
library("stringr")

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_Information.tsv",sep="\t",header =T, stringsAsFactors = F)
meta_info$Raw_name = str_replace_all(meta_info$Raw_name, pattern = "-", "_")

excluded_files = meta_info$Raw_name[meta_info$Included == FALSE]
#excluded_files = meta_info$Raw_name[as.character(meta_info$pID) != "365"]
excluded_files = excluded_files[!is.na(excluded_files)]
for ( file in excluded_files){
  
  ori_file = paste( c( "~/Koop_Klinghammer/Data/Raw_data/",file, ".RCC"), collapse = "" )
  print(c(file, file.exists(ori_file)))
  file.copy(ori_file, 
            paste( c( "~/Koop_Klinghammer/Data/Excluded/",file, ".RCC"), collapse = "" )
  )
  file.remove(paste( c( "~/Koop_Klinghammer/Data/Raw_data/",file, ".RCC"), collapse = "" ))
}
included_files = meta_info$Raw_name[meta_info$Included == TRUE]
for ( file in list.files("~/Koop_Klinghammer/Data/Raw_data/")){
  
  ori_file = paste( c( "~/Koop_Klinghammer/Data/Raw_data/",file, ".RCC"), collapse = "" )
  print(c(file, file.exists(ori_file)))
  if (!( file %in% included_files)){
    file.copy(ori_file, 
              paste( c( "~/Koop_Klinghammer/Data/Excluded/",file, ".RCC"), collapse = "" )
    )
    file.remove(paste( c( "~/Koop_Klinghammer/Data/Raw_data/",file, ".RCC"), collapse = "" ))
  }
}

rcc_files <- dir( "~/Koop_Klinghammer/Data/Raw_data_new/", full.names = TRUE, pattern = "*")
length(rcc_files)
raw_data = NanoStringNorm::read.markup.RCC("~/Koop_Klinghammer/Data/Raw_data_new/" , rcc.pattern =  "*.RCC")

norm.comp.results.test = norm.comp(raw_data, verbose = T)
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
  verbose = T
)

### tmp

expr_raw = normalized_data$normalized.data
dim(expr_raw)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X","")
colnames(expr_raw) = str_replace_all(colnames(expr_raw), pattern = "\\.","_")

meta_match = match( colnames(expr_raw), meta_info$Raw_name, nomatch = 0)
sample_names[meta_match == 0]
colnames(expr_raw)[-seq(3)] = meta_info$Sample_ID[meta_match]
colnames(expr_raw)[1:5]
dim(expr_raw)

#write.table(expr_raw, "~/Koop_Klinghammer/Data/New_data.S141.Kulbe.tsv",quote = FALSE, sep ="\t", row.names = TRUE)
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

matcher = match(colnames(expr_raw), meta_info$Raw_name, nomatch = 0)
colnames(expr_raw) = meta_info$Sample_ID[matcher]
#write.table(expr_raw, "~/Koop_Klinghammer/Data/New_data.S138.Kulbe.tsv", quote= F, row.names = T, sep = "\t")

