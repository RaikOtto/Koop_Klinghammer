require('NanoStringNCTools')
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

rcc_files <- dir( "~/Koop_Klinghammer/Data/Raw_data/", full.names = TRUE, pattern = "*")
raw_data = readNanoStringRccSet(rcc_files, )
head(assayData(raw_data)[["exprs"]])
head(summary(raw_data, MARGIN = 1), 2)

demoData <- normalize(raw_data, type="nSolver", fromELT = "exprs", toELT = "exprs_norm")
assayDataElement(demoData, elt = "exprs_norm")[1:3, 1:2]

### tmp
sample_names = names(raw_data$raw)
sample_names = str_replace(sample_names, pattern = "^X","")
sample_names = str_replace_all(sample_names, pattern = "\\.","_")

meta_match = match( sample_names, meta_info$Raw_name, nomatch = 0)
sample_names[meta_match == 0]
names(raw_data$raw)[-seq(3)] = meta_info$Sample_ID[meta_match]
names(raw_data$raw)[1:3]
### normalization

#norm.comp.results.test = norm.comp(raw_data, verbose = T)
eset = HKnorm( 
  raw_data$raw)

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

d = colMeans(pure_data)
boxplot(pure_data[,order( apply(pure_data, MARGIN = 2, FUN = function(vec){return(mean(vec))}) )])
#pure_data[,order( (colMeans(pure_data)) )]
#write.table(pure_data, "~/Koop_Klinghammer/Data/Pure_data.05_06_2018.tsv", quote= F, row.names = T, sep = "\t")


####

rtv_t = read.table("~/Koop_Klinghammer/Misc/RTV.csv",sep="\t",header =T ,skip =0)
rownames(rtv_t) = rtv_t$Sample
rtv_t = rtv_t[,-1]
expr_raw = read.table("~/Koop_Klinghammer/Data/ExpressionSet_Klinghammer_neun.tsv",sep="\t",header =T ,skip =0)
colnames(expr_raw) = str_replace_all(colnames(expr_raw),"^X","")
rownames(expr_raw) = expr_raw$HGNC
expr_raw = expr_raw[,-1]
expr_raw[1:5,1:5]

match(colnames(expr_raw), rownames(rtv_t))
meta_data = rtv_t[,]
dim(meta_data)

meta_data = meta_data[meta_data$ENummer != "",]
enummer = str_replace_all(pattern = "^HN_",meta_data$ENummer,"")
enummer = str_replace_all(pattern = "^HNSCC_",enummer,"")
rownames(meta_data) = enummer
meta_data[rownames(rtv_t),]
