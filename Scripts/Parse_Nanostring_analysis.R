require('NanoStringNorm')
library("stringr")

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_Information.tsv",sep="\t",header =T, stringsAsFactors = F)
meta_info$Raw_name = str_replace_all(meta_info$Raw_name, pattern = "-", "_")
meta_info$OS = str_replace_all(meta_info$OS, pattern = ",", ".")
meta_info$OS = as.double(meta_info$OS)
meta_info$OS[is.na(meta_info$OS)] = ""

excluded_files = meta_info$Raw_name[meta_info$Included == FALSE]
#excluded_files = meta_info$Raw_name[as.character(meta_info$pID) != "365"]
excluded_files = excluded_files[!is.na(excluded_files)]
for ( file in excluded_files){
  
  ori_file = paste( c( "/home/ottoraik/Koop_Klinghammer/Data/Raw_data/",file, ".RCC"), collapse = "" )
  print(c(file, file.exists(ori_file)))
  file.copy(ori_file, 
            paste( c( "/home/ottoraik/Koop_Klinghammer/Data/Excluded/",file, ".RCC"), collapse = "" )
  )
  file.remove(paste( c( "/home/ottoraik/Koop_Klinghammer/Data/Raw_data/",file, ".RCC"), collapse = "" ))
}
included_files = meta_info$Raw_name[meta_info$Included == FALSE]
for ( file in list.files("~/Koop_Klinghammer/Data/Raw_data/")){
  
  ori_file = paste( c( "/home/ottoraik/Koop_Klinghammer/Data/Raw_data/",file, ".RCC"), collapse = "" )
  print(c(file, file.exists(ori_file)))
  if (!( file %in% included_files)){
      file.copy(ori_file, 
          paste( c( "/home/ottoraik/Koop_Klinghammer/Data/Excluded/",file, ".RCC"), collapse = "" )
      )
      file.remove(paste( c( "/home/ottoraik/Koop_Klinghammer/Data/Raw_data/",file, ".RCC"), collapse = "" ))
  }
}

raw_data  = read.markup.RCC( rcc.path = "~/Koop_Klinghammer//Data/Raw_data/", rcc.pattern = "*.RCC")

### tmp
sample_names = names(raw_data$header)
sample_names = str_replace(sample_names, pattern = "^X","")
sample_names = str_replace_all(sample_names, pattern = "\\.","_")

meta_match = match( sample_names, meta_info$Raw_name, nomatch = 0)
sample_names[meta_match == 0]
colnames(raw_data$x)[-seq(3)] = meta_info$Sample_ID[meta_match]

### normalization

#norm.comp.results.test = norm.comp(raw_data, verbose = T)
eset = NanoStringNorm::NanoStringNorm( 
  raw_data,
  CodeCount.methods = "sum",
  Background.methods = "mean.2sd",
  SampleContent.methods = "housekeeping.sum",
  OtherNorm.methods = "vsn",
  take.log = T,
  round.values = T,
  return.matrix.of.endogenous.probes = F,
  verbose = T
)

source_mat = eset$normalized.data
m    = matrix( as.character(unlist( source_mat)), nrow=  dim(source_mat)[1], ncol = dim(source_mat)[2])
info = m[,seq(3)]
data = matrix( as.double(m[,-seq(3)]), nrow=  dim(source_mat)[1], ncol = dim(source_mat)[2]-3)
data = round(data,1)

rownames(data) = rownames(source_mat)
col_labels = str_replace( colnames(source_mat)[-seq(3)], pattern = "^X", "") 

colnames(data) = col_labels
res = cbind( info,data )
res = cbind(rownames(source_mat), res)

pure_data = as.character(res)[-seq(dim(source_mat)[1]*4)]
pure_data = matrix( as.double( pure_data ), nrow = dim(source_mat)[1] )
rownames( pure_data ) = info[ ,2 ]
pure_data = pure_data[ info[,1] == "Endogenous"  ,]
colnames(pure_data) = colnames(data)

#pure_data = pure_data[,colnames(pure_data) != "25"]
s_match = match( as.character( colnames(pure_data)), as.character( meta_info$Sample_ID), nomatch = 0)
meta_data = meta_info[s_match,]
rownames(meta_data) = meta_data$Sample_ID
meta_data$OS = as.double(meta_data$OS)

pure_data = pure_data[,colnames(pure_data) %in% meta_info$Sample_ID[meta_info$Included]]
d = colMeans(pure_data)
boxplot(pure_data[,order( apply(pure_data, MARGIN = 2, FUN = function(vec){return(mean(vec))}) )])
#pure_data[,order( (colMeans(pure_data)) )]
#write.table(pure_data, "~/Koop_Klinghammer/Data/Pure_data.05_06_2018.tsv", quote= F, row.names = T, sep = "\t")

