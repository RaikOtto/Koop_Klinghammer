library(stringr)

expr_raw = read.table("~/Koop_Klinghammer/Data/New_data.S130.tsv", sep ="\t", row.names = 1, header = T)
colnames(expr_raw) = str_replace(colnames(expr_raw), pattern = "^X","")
expr_raw[1:5,1:5]
dim(expr_raw)

balanced.centroid = read.table( "~/Koop_Klinghammer/Misc/Archiv//balanced.centroid.txt", header=TRUE, row.names=1, sep="\t",stringsAsFactors = F)
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
meta_data = meta_info[match(colnames(expr_raw), as.character( meta_info$SampleID)),]
dim(meta_data)

###

source("~/Koop_Klinghammer/Scripts/Classification_scripts.R")
table( rownames(expr_raw) %in% rownames(balanced.centroid) )
table( rownames(balanced.centroid) %in% rownames(expr_raw) )

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

#write.table( class_data, "~/Koop_Klinghammer/Results/Normalized_classification_data.S130.tsv",sep ="\t", quote =F , row.names = FALSE)

meta_match = match( colnames(expr_raw), meta_info$SampleID, nomatch = 0 )
meta_info$Subtype = rep("",nrow(meta_info))
meta_info$Subtype[meta_match] = as.character( class_data$Subtype )
meta_info$P_value = rep("",nrow(meta_info))
meta_info$P_value[meta_match] = as.double( as.character( class_data$P_value ) )

meta_info[(as.double(meta_info$P_value[meta_match]) > 0.05),"Included"] = FALSE

write.table(meta_info,"~/Koop_Klinghammer/Misc/Meta_information.tsv",sep ="\t",quote =F,row.names =F)
rownames(meta_info) = meta_info$SampleID

meta_data = meta_info[colnames(expr_raw),]
meta_data$P_value = as.double(meta_data$P_value)
expr_raw = expr_raw[,meta_data$P_value < .05]
sum(meta_data$P_value < .05)
#SLC16A1, AKR1C1, HIF1A, AKR1C3 E2F2, TMUB2, PRPF38A, STAT1, EGFR, LAG3, SLAMF6, ITGB1, CXCL10, AMMECR1L, TNFRSF1A, AREG, E6 (HPV16), VEGF, IDO1, DDX50, RPA2

meta_info[exclusion_samples,"Included"] = "FALSE"
meta_info[colnames(expr),"Included"] = "TRUE"

#write.table(expr_raw,"~/Koop_Klinghammer/Data/Normalized_data.S121.Significant.tsv",quote =FALSE, sep ="\t", row.names = TRUE)
