library("stringr")
library("grid")

meta_info = read.table("~/Koop_Klinghammer/Misc/Meta_Information.tsv",sep="\t",header =T, stringsAsFactors = F)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 90, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

### Prep

aka3 = list(
  #Subtype = c(BA = "blue", CL ="darkgreen", MS = "red"),
  Subtype = c(BA = "blue", CL ="darkgreen", MS = "red", Not_sig = "gray"),
  Location = c(Primary = "white", Metastasis = "black"),
  OS = c(high = "White", medium = "gray", low = "black"),
  Correlation = c(high = "Red", medium = "Yellow", low = "Green"),
  Included = c(Yes = "green", No = "red"),
  OS_Log2_Monate = c(high = "White", medium = "gray", low = "black"),
  Best_response = c(PD = "Yellow", SD = "Darkred", PR = "Green")
)

###

## Figure 1

meta_info$OS = as.double(str_replace_all(meta_info$OS, pattern = ",", "."))
meta_info$Subtype[meta_info$Subtype == ""] = "Not_sig"

meta_data = meta_info[ meta_info$Sample_ID %in% colnames(pure_data),]
meta_data_tmp = meta_data[as.double(meta_data$P_value) < 0.001,]
pure_data_tmp = pure_data[,colnames(pure_data) %in% meta_data_tmp$Sample_ID ]
cor_mat = cor(pure_data_tmp);pcr = prcomp(t(cor_mat))

meta_data_tmp$OS_Log2_Monate = log2(as.double(meta_data_tmp$OS)+1)
meta_data_tmp$pID = as.factor(meta_data_tmp$pID)
pheatmap::pheatmap(
  cor_mat,
  annotation_col = meta_data_tmp[c("Subtype","Best_response","OS_Log2_Monate")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  legend = F,
  fontsize_col = 7,
  clustering_method = "average"
)

### Figure 2

ggbiplot::ggbiplot(
  pcr,
  groups = meta_data_tmp[,c("Subtype")],
  ellipse = TRUE,
  circle = TRUE,
  var.axes = F,
  var.scale = meta_data_tmp$OS**2#,
 # labels = meta_data_tmp$Sample_ID
)  + geom_point( aes( size = as.double(meta_data_tmp$OS)**2, color = as.factor(meta_data_tmp$Subtype) ))

aggregate(as.double(meta_data$OS), FUN = mean, by = list(meta_data$Subtype))

### umap

#umap_plot = umap::umap(t(pure_data))
#vis_data = as.data.frame(umap_plot$layout)
#colnames(vis_data) = c("x","y")
#dist_mat = dist((vis_data))
#p = ggplot2::qplot( x = vis_data$x, y = vis_data$y, color = meta_data$Subtype)
#p

## km plots

meta_data_SP = meta_data_tmp
meta_data_SP = meta_data_SP[match(unique(meta_data_SP$pID),meta_data_SP$pID),]

meta_data_SP$OS_Monate = as.double(str_replace_all(meta_data_SP$OS_Monate, pattern = ",", "."))
meta_data_SP$PFS_Monate = as.double(str_replace_all(meta_data_SP$PFS_Monate, pattern = ",", "."))

fit <- survival::survfit( survival::Surv( OS_Monate) ~ Subtype,data = meta_data_SP)
survminer::ggsurvplot(fit, data = meta_data_SP, risk.table = T, pval = T)

fisher.test( meta_data_SP$Subtype[meta_data_SP$Subtype %in% c("BA","CL")], meta_data_SP$OS_Monate[meta_data_SP$Subtype %in% c("BA","CL")])
t.test( meta_data_SP$OS_Monate[meta_data_SP$Subtype %in% c("BA")], meta_data_SP$OS_Monate[meta_data_SP$Subtype %in% c("CL")])

## Same patient only

find_vec = meta_data$pID
match_vec = match(find_vec, unique(find_vec),nomatch = 0)
multi_match = which( table(match_vec) > 1  )
multi_match = multi_match[multi_match != 12]
true_match = meta_data$pID[which( match_vec %in% multi_match)]

multi_data = pure_data[,which( meta_data$pID %in% true_match )]
multi_cor = cor(multi_data)
meta_data$pID = as.character(meta_data$pID)

pheatmap::pheatmap(
  multi_cor,
  annotation_col = meta_data[c("Subtype","OS","pID")],
  annotation_colors = aka3,
  show_rownames = F,
  show_colnames = T,
  #treeheight_col = 0,
  legend = F,
  fontsize_col = 7
)

## GOI expression

genes_of_interest = c("AKR1C1","AKR1C3","ALDH3A1","E2F2","MCM2","CDC7","CDKN2A")
genes_of_interest = c("CD74","LAG3","IL17A","IDO1","CXCL9","STAT1","IFNG","HLADRA","VIM")
genes_of_interest = c("VEGF","HIF1A","EGFR","AREG","KRT17","COL17A1","ITGB1","KRT19")
genes_of_interest = c("SNAI2","TGFBI","EMT","SPRR3","S100A4")
genes_of_interest = c("AKR1C1","AKR1C3","ALDH3A1","E2F2","MCM2","CDC7","CDKN2A","CD74","LAG3","IL17A","IDO1","CXCL9","STAT1","IFNG","HLADRA","VIM","VEGF","HIF1A","EGFR","AREG","KRT17","COL17A1","ITGB1","KRT19","SNAI2","TGFBI","SPRR3","S100A4")

genes_of_interest[!( genes_of_interest %in% rownames(pure_data) )]

vis_mat = t(pure_data)
vis_mat = reshape2::melt(vis_mat )
vis_mat = cbind(vis_mat, meta_data$Subtype[match(vis_mat[,1], meta_data$Sample_ID)])
colnames(vis_mat) = c("Sample","Gene","Expression","Subtype")
vis_mat$Expression = as.double(vis_mat$Expression)
vis_mat$Gene = as.character(vis_mat$Gene)
vis_mat = subset(vis_mat, Gene %in% genes_of_interest)
vis_mat = subset(vis_mat, Subtype != "Not_sig")

Subtype_col = vis_mat$Subtyp 
Subtype_col[Subtype_col == "BA"] = "red"
Subtype_col[Subtype_col == "CL"] = "Darkgreen"
Subtype_col[Subtype_col == "MS"] = "Blue"
Subtype_col[Subtype_col == "Control"] = "Yellow"
men1_plot = ggplot( data = vis_mat, aes ( x = Gene,  y = Expression))
men1_plot = men1_plot + geom_boxplot( aes(fill = Subtype_col))
men1_plot = men1_plot + theme(axis.text.x = element_text(angle = 45, hjust = 1))
men1_plot = men1_plot + scale_fill_manual(values=c("Red", "Darkgreen", "Blue","BLACK"), name="Subtype",labels=c("BA", "CL", "MS","Control"))
men1_plot

### Exp per Sample

## MEN1 transcripts

library(ggplot2)

# Basic barplot

for(gene in genes_of_interest){
    
    vis_mat = pure_data[rownames(pure_data) == gene,]
    vis_mat = cbind(vis_mat, colnames(pure_data))
    vis_mat = cbind(vis_mat, meta_data$Subtype)
    colnames(vis_mat) = c("Expression","Sample_ID","Subtype")
    vis_mat = as.data.frame(vis_mat)
    vis_mat$Sample_ID = factor(vis_mat$Sample_ID, levels = vis_mat$Sample_ID[order(vis_mat$Expression)] )
    vis_mat$Expression = as.double(vis_mat$Expression)  
  
    os_vec = meta_data$OS_Monate
    os_vec = os_vec[match(meta_data$Sample_ID,vis_mat$Sample_ID)]
    os_vec = as.double( str_replace_all( meta_data$OS_Monate, pattern = ",", "."))
    os_mean = mean(os_vec[!is.na(os_vec)])
    os_vec[is.na(os_vec)] = os_mean
    os_vec[os_vec > os_mean] = "A"
    os_vec[os_vec != "A"] = "B" 
    os_vec[levels(vis_mat$Sample_ID[order(vis_mat$Expression)]) == 12] = "C"
    
    col_vec = os_vec
    col_vec[col_vec == "A"] = "darkgreen"
    col_vec[col_vec != "darkgreen"] = "darkred"
    col_vec[levels(vis_mat$Sample_ID[order(vis_mat$Expression)]) == 12] = "Black"
    
    p = ggplot( data = vis_mat,aes( x = Sample_ID, y = Expression, fill = Subtype ))
    p = p + geom_bar(stat="identity", position=position_dodge())
    p = p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    p = p + scale_fill_manual(values=c("Red", "Darkgreen","BLACK", "Blue"), name="Subtype",labels=c("BA", "CL","Control", "MS"))
    p = p + annotate("text", x=1:ncol(pure_data),y = max(vis_mat$Expression) -1 ,parse=TRUE, label = os_vec,color = col_vec, size = 4.5 )
    p = p + xlab("") + ylab( paste0(gene, " expression in nCounts") ) + theme(legend.position = "top")
    
    file_name = paste0( collapse = "" ,c("~/Koop_Klinghammer/Results/Plots_Darbeit_Victoria/6_", gene, ".png"))
    png(file_name, width = 1024, height = 680)
        plot(p)
    dev.off()
}
