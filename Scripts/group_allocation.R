library(stringr)

expr_raw = read.table("~/Koop_Klinghammer/Data/Normalized_data.S114.20012022.tsv", sep ="\t", row.names = 1, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X","")
expr_raw[1:5,1:5]
dim(expr_raw)

balanced.centroid = read.table( "~/Koop_Klinghammer/Misc/balanced.centroid.txt", header=TRUE, row.names=1, sep="\t",stringsAsFactors = F)
balanced.centroid_importance = sort(rowSums(abs(balanced.centroid)), decreasing = T)
balanced.centroid = balanced.centroid[ match(names(balanced.centroid_importance),rownames(balanced.centroid)),]

### Preparation

rownames(expr_raw) = str_to_upper(rownames(expr_raw))
rownames(expr_raw) = str_replace_all(rownames(expr_raw), pattern = "\\_","" )
rownames(expr_raw) = str_replace_all(rownames(expr_raw), pattern = "-","" )
rownames(balanced.centroid) = str_to_upper(rownames(balanced.centroid))
rownames(balanced.centroid) = str_replace_all(rownames(balanced.centroid), pattern = "\\_","" )
rownames(balanced.centroid) = str_replace_all(rownames(balanced.centroid), pattern = "-","" )
colnames(expr_raw) = str_replace_all(colnames(expr_raw), pattern = "^X","" )

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_information.tsv",sep="\t",header =T, stringsAsFactors = F)
###

source("~/Koop_Klinghammer/Scripts/Classification_scripts.R")
centroid_genes = rownames(balanced.centroid)
genes_matching_centroid = rownames(expr_raw)[which(  (rownames(expr_raw) %in% rownames(balanced.centroid) ) ) ]
length(genes_matching_centroid)
genes_not_matching_centroid = rownames(expr_raw)[which(!(  (rownames(expr_raw) %in% rownames(balanced.centroid) ) )) ]
length(genes_not_matching_centroid)

### centroid classification

pub_cor <<- matrix( as.double(), ncol = length( colnames( balanced.centroid )  ) )

expr2bc = centroid2expr( balanced.centroid[,], expr_raw )
colnames(expr2bc$correlation) = c("Sample","Subtype","Correlation","P_value")
class_data = as.data.frame(expr2bc$correlation)
dim(class_data)

table(as.double(class_data$P_value) <= 0.05)
candidates = class_data$Sample[ which(as.double(class_data$P_value) <= 0.05)]
expr_raw = expr_raw[,match(candidates, colnames(expr_raw),nomatch = 0)]
dim(expr_raw)
subtype_vec = class_data[match(colnames(expr_raw),class_data$Sample),"Subtype"]
p_value_vec = class_data[match(colnames(expr_raw),class_data$Sample),"P_value"]

#write.table( class_data, "~/Koop_Klinghammer/Results/Data.S108.20012022.tsv",sep ="\t", quote =F , row.names = FALSE)

matcher = match(colnames(expr_raw),meta_info$SampleID,nomatch = 0)
meta_data = meta_info[matcher,]
dim(meta_data)

meta_data[matcher,"Subtype"] = subtype_vec
meta_data[matcher,"P_value"] = p_value_vec

#write.table(meta_data,"~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t",quote =F,row.names =F)
matcher = match(meta_data$SampleID, colnames(expr_raw),nomatch = 0)
expr_raw = expr_raw[,matcher]
#write.table(expr_raw,"~/Koop_Klinghammer/Data/S108.tsv",sep ="\t",quote =F,row.names =TRUE)
