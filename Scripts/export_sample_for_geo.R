require('NanoStringNorm')
library("stringr")
library("dplyr")
library("stringr")


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

meta_info = read.table("~/Koop_Klinghammer/Misc/Archiv/Meta_information_save.tsv",sep ="\t", stringsAsFactors = F, header = T)
rownames(meta_info) = meta_info$Name

# raw_data

rcc_files <- dir( "~/Koop_Klinghammer/Data/Raw_data_new/", full.names = TRUE, pattern = "*")
length(rcc_files)
raw_data = NanoStringNorm::read.markup.RCC("~/Koop_Klinghammer/Data/Raw_data/" , rcc.pattern =  "*.RCC")$x
raw_data = raw_data %>% select(-CodeClass) 
rownames(raw_data) = raw_data$Name
raw_data = raw_data %>% select(-Name) 
colnames(raw_data) = str_replace_all(colnames(raw_data), "^X","")
raw_data[1:5,1:5]

#matcher = match(colnames(raw_data), meta_info$Raw_Name, nomatch = 0)
#sample_id = meta_info$Name[matcher]

matcher = match( colnames(expr_raw), meta_info$Name , nomatch = 0)
raw_id = meta_info$Raw_Name[matcher]
raw_id %in% colnames(raw_data)
raw_data_export = raw_data[,c("Accession",raw_id)]
colnames(raw_data_export) = colnames(expr_raw)

#write.table(raw_data_export,"~/Koop_Klinghammer/Data/Raw_data_dd.tsv",sep ="\t", row.names = TRUE, quote =F)

files = sapply(colnames(raw_data_export)[-1], FUN = function(name){return(paste0(c("~/Koop_Klinghammer/Data/Raw_data/",name,".RCC"), collapse = ""))})

file.exists(files)
file.copy(files, "~/Koop_Klinghammer/Data/Raw_data/Export/")