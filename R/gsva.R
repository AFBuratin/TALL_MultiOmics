library(GSVA)
library(GSEABase)
library(msigdbr)
# library(pd.clariom.s.human)
require(affy)
require(limma)
require(hgu95a.db)
require(annotate)

hallmark_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens", # Can change this to what species you need
  category = "H"#, # H for Only hallmark gene sets
  # subcategory = "BP"
)

hallmark_gene_sets <- msigdbr::msigdbr(
  species = "Homo sapiens", # Can change this to what species you need
  category = "C2" # H for Only hallmark gene sets
  # subcategory = "CP"
)

hallmarks_list <- split(
  hallmark_gene_sets$gene_symbol, # The genes we want split into pathways
  hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
)


library(data.table)
Cluster.data = rbind(mRNA_Sel_FC_CS1vsCS2[match(input.mRNA, mRNA_Sel_FC_CS1vsCS2$SYMBOL), c(1,2:3)],
  circRNA_Sel_FC_CS1vsCS2[match(input.circ, circRNA_Sel_FC_CS1vsCS2$SYMBOL), c(1, 2:3)])
Cluster.matrix = Cluster.data[,-1]
rownames(Cluster.matrix) = Cluster.data$SYMBOL

gsva_results_conc_sil <- gsva(
  as.matrix(Cluster.matrix),#[rownames(pat.matrix),],
  hallmarks_list, #$DANG_BOUND_BY_MYC,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Gaussian",
  # Minimum gene set size
  min.sz = 5,
  # Maximum gene set size
  max.sz = 1000,
  # Compute Gaussian-distributed scores
  mx.diff = FALSE,
  # Don't print out the progress bar
  verbose = FALSE
)


gsva_results_conc_sil[rownames(gsva_results_conc_sil)%in%"GSE5463_CTRL_VS_DEXAMETHASONE_TREATED_THYMOCYTE_DN",]
gsva_results_conc_sil[rownames(gsva_results_conc_sil)%in%"GSE5463_CTRL_VS_DEXAMETHASONE_TREATED_THYMOCYTE_UP",]

annot_df_cong = data.frame("Cluster" = c("CS1", "CS2"))
rownames(annot_df_cong) = colnames(Cluster.matrix)

ann_colors <- list(
  Cluster = c(`CS1` = "#619CFF", `CS2` = "#F8766D")
)

rownames(gsva_results_conc_sil) = substr(gsub("_", " ", gsub("^HALLMARK_", "", rownames(gsva_results_conc_sil))), 1, 35)
merge.dt = merge(gsva_results_pat_M[rownames(gsva_results_pat_M)%like%"NOTCH",],
      gsva_results_pat_S[rownames(gsva_results_pat_S)%like%"NOTCH",], by = "row.names")
plot.dt = merge.dt[,-1]
rownames(plot.dt) = merge.dt$Row.names
pathway_heatmap <- pheatmap::pheatmap(gsva_results_conc_sil, 
                                      clustering_distance_cols = "correlation", 
                                      #clustering_distance_rows = "correlation", 
                                      clustering_method = "ward.D2",
                                      # annotation_col = annot_df_cong, # Add metadata labels!
                                      # annotation_colors = ann_colors,
                                      show_colnames = T, # Don't show sample labels
                                      fontsize_row = 18 # Shrink the pathway labels a tad
)
pathway_heatmap <- pheatmap::pheatmap(gsva_results_conc_sil, 
                                      clustering_distance_cols = "correlation", 
                                      #clustering_distance_rows = "correlation", 
                                      clustering_method = "ward.D2",
                                      annotation_col = annot_df_cong, # Add metadata labels!
                                      annotation_colors = ann_colors,
                                      show_colnames = T, # Don't show sample labels
                                      fontsize_row = 20 # Shrink the pathway labels a tad
)
pathway_heatmap

library(matrixTests)
colnames(gsva_results_conc_sil) = c("sirNEG", "sir-cFBXW7")
row_t_welch(gsva_results_conc_sil[,c("sirNEG")], gsva_results_conc_sil[,c("sir-cFBXW7")])
Ttest.pvalue = data.frame(pvalue = unlist(sapply(1:nrow(gsva_results_conc_sil), 
                                                 function(i) t.test(as.numeric(as.character(unlist(gsva_results_conc_sil[i,c(1,3)]))), 
                                                                    as.numeric(as.character(unlist(gsva_results_conc_sil[i,c(2,4)]))),
                                                                    paired = T, alternative = "less")[c("p.value")])))

Ttest.pvalue$padj = p.adjust(Ttest.pvalue$pvalue, method = "BH")
Ttest.pvalue$HM = rownames(gsva_results_conc_sil)
MYC_ptw_sign = Ttest.pvalue[Ttest.pvalue$pvalue<=0.1&Ttest.pvalue$HM%like%"MYC",]

diff = rowMeans(gsva_results_conc_sil[rownames(gsva_results_conc_sil)%in%MYC_ptw_sign$HM,c(2,4)])+0.001 - rowMeans(gsva_results_conc_sil[rownames(gsva_results_conc_sil)%in%MYC_ptw_sign$HM,c(1,3)])+0.001
diff = gsva_results_conc_sil[,2]+0.001 - gsva_results_conc_sil[,1]+0.001


# Set stylings for row names and make our selected rows unique
row_idx <- which(rownames(gsva_results_conc_sil)%like% "NOTCH|MYC")
fontsizes <- rep(10, nrow(gsva_results_conc_sil))
fontsizes[row_idx] <- 10
fontcolors <- rep('black', nrow(gsva_results_conc_sil))
fontcolors[row_idx] <- 'black'
fontfaces <- rep('plain',nrow(gsva_results_conc_sil))
fontfaces[row_idx] <- 'bold'
rowAnno <- rowAnnotation(rows = anno_text(names(diff), #rot = -45,
                                          gp = gpar(fontsize = 6, fontface = fontfaces, col = fontcolors)))
library(ComplexHeatmap)
library(circlize)
ha = HeatmapAnnotation(df = data.frame(Treatment = c("sirNEG", "sir-cFBXW7")),
                       col = list(Treatment = c(`sirNEG` = "#C0C0C0", `sir-cFBXW7` = "#FF9933")),
                       show_annotation_name = F,
                       annotation_legend_param = list(Treatment = list(direction = "horizontal")))
m <- which(rownames(gsva_results_conc_sil) %like% 'MYC|NOTCH')
m = which(abs(diff)>=2)
ha.row = rowAnnotation(foo = anno_mark(at = m, labels = rownames(gsva_results_conc_sil)[m], labels_gp = gpar(fontsize = 8)))
library("RColorBrewer")
png(filename = "/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/ht_Boldhallmark_FC.png",
    width = 12, height = 15, units = "cm", res = 200)
draw(ComplexHeatmap::Heatmap(gsva_results_conc_sil,#[rownames(gsva_results_conc_sil)%in%MYC_ptw_sign$HM,], 
                             name = "GSVA score", 
                             # km = 2,
                             # column_km = 2,
                             # column_title = "circFBXW7-sil/RPMI",
                             # col = inferno(255),
                             col = rev(brewer.pal(9,"PRGn")),#colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                             top_annotation = ha,
                             # top_annotation_height = unit(2, "mm"),
                             clustering_distance_columns = "spearman",
                             cluster_columns = F,
                             clustering_method_columns = "ward.D2",
                             clustering_method_rows = "ward.D2",
                             # cluster_columns = dend,
                             clustering_distance_rows = "euclidean",
                             cluster_rows = T,
                             column_order =  c("sirNEG", "sir-cFBXW7"), #c("sirNEG", "sirNEG", "sir-cFBXW7", "sir-cFBXW7"),
                             row_dend_side = "left",
                             row_names_side = "right",
                             show_row_names = F, 
                             show_column_names = F, 
                             width = unit(1.5, "cm"),
                             show_row_dend = T,
                             # right_annotation = ha.row,
                             # show_column_dend = T) + 
                             # row_names_gp = gpar(fontsize = 8),
                             heatmap_legend_param = list(direction = "horizontal")) +
       Heatmap(diff, show_column_names = F, width = unit(2, "mm"),
               colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255),
               show_row_names = F,
               right_annotation = rowAnno, 
               heatmap_legend_param = list(title = "GSVA score diff.", direction = "horizontal")), #+
       # rowAnnotation(link = anno_mark(at = which(diff >= 2), 
       #                                labels = rownames(gsva_results_conc_sil)[diff >= 2], 
       #                                labels_gp = gpar(fontsize = 10), padding = unit(1, "mm"))),
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
     
dev.off()
# Print out heatmap here
png(filename = "/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/gsva_NOTCHpath_RPMI.png",
    width = 50, height = 31, units = "cm", res = 150)
pathway_heatmap
dev.off()

# Print out heatmap here
png(filename = "/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/gsva_NOTCH_SoulierMullighan.png",
    width = 50, height = 31, units = "cm", res = 150)
pathway_heatmap
dev.off()

# BiocManager::install("GeneAnswers")
library(GeneAnswers)
length(getSymbols(as.character(hallmarks_list$WP_NUCLEAR_RECEPTORS), 'org.Hs.eg.db'))

DEpwys_es <- gsva_results_pat
colorLegend <- c("chartreuse4", "burlywood3")
names(colorLegend) <- c("Control", "Silencing")
sample.color.map <- colorLegend[pData(rawData)[, "group"]]
names(sample.color.map) <- colnames(DEpwys_es)
sampleClustering <- hclust(as.dist(1-cor(DEpwys_es, method="spearman")),
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(DEpwys_es), method="spearman")),
                            method="complete")
heatmap(DEpwys_es, ColSideColors=sample.color.map, xlab="samples",
        ylab="Pathways", margins=c(2, 20),
        # labRow=substr(gsub("_", " ", gsub("^KEGG_|^REACTOME_|^BIOCARTA_", "",
        #                                   rownames(DEpwys_es))), 1, 35),
        # labRow=substr(gsub("_", " ", gsub("^HALLMARK_", "",
        #                                   rownames(DEpwys_es))), 1, 35),
        labRow=substr(gsub("_", " ", gsub("^GOBP_|^GOCC_|^GOMF_", "",
                                          rownames(DEpwys_es))), 1, 35),
        labCol="", scale="row", Colv=as.dendrogram(sampleClustering),
        Rowv=as.dendrogram(geneSetClustering))
legend("topright", names(colorLegend), fill=colorLegend, bg="white")
