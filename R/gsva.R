library(GSVA)
library(GSEABase)
library(GSVAdata)
library(pd.clariom.s.human)
require(affy)
require(limma)
require(hgu95a.db)
require(annotate)
# vst_df <- exprs(data_final) 
vst_df <- exp %>%
  as.data.frame() %>% # Make into a data frame
  tibble::rownames_to_column("PROBEID") # Make Gene IDs into their own column

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
genes_UPbyMYC = read.table("/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/genesetMYC_UP.txt")
genes_UPbyMYC = read.table("/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/harmonizome_c-Myc_MotifMap+Predicted+Transcription+Factor+Targets.txt")
genes_UPbyMYC_entrez_ids <- data.frame(mapIds(org.Hs.eg.db, genes_UPbyMYC$V1,
                     keytype = "SYMBOL",
                     column = "ENTREZID"))

hallmarks_list <- split(
  hallmark_gene_sets$gene_symbol, # The genes we want split into pathways
  hallmark_gene_sets$gs_name # The pathways made as the higher levels of the list
)
hallmarks_list[["GSE5463_CTRL_VS_DEXAMETHASONE_TREATED_THYMOCYTE_UP"]]
hallmarks_list[["YAGUE_PRETUMOR_DRUG_RESISTANCE_DN"]]
hallmarks_list$MYC_UP = genes_UPbyMYC$V1
library(data.table)
names(hallmarks_list)[names(hallmarks_list) %like% 'NOTCH']

mapped_df <- data.frame(
  "symbol_id" = mapIds(
    # Replace with annotation package for the organism relevant to your data
    clariomshumantranscriptcluster.db,
    keys = vst_df$PROBEID,
    # Replace with the type of gene identifiers in your data
    keytype = "PROBEID",
    # Replace with the type of gene identifiers you would like to map to
    column = "SYMBOL",
    # This will keep only the first mapped value for each Ensembl ID
    multiVals = "first"
  )
) %>%
  # If an Ensembl gene identifier doesn't map to a Entrez gene identifier,
  # drop that from the data frame
  dplyr::filter(!is.na(symbol_id)) %>%
  # Make an `Ensembl` column to store the row names
  tibble::rownames_to_column("PROBEID") %>%
  # Now let's join the rest of the expression data
  dplyr::inner_join(vst_df, by = c("PROBEID" = "PROBEID"))
head(mapped_df)
# First let's determine the gene means
gene_means <- rowMeans(mapped_df %>% dplyr::select(-PROBEID, -symbol_id))

# Let's add this as a column in our `mapped_df`.
mapped_df <- mapped_df %>%
  # Add gene_means as a column called gene_means
  dplyr::mutate(gene_means) %>%
  # Reorder the columns so `gene_means` column is upfront
  dplyr::select(PROBEID, symbol_id, gene_means, dplyr::everything())
filtered_mapped_df <- mapped_df %>%
  # Sort so that the highest mean expression values are at the top
  dplyr::arrange(dplyr::desc(gene_means)) %>%
  # Filter out the duplicated rows using `dplyr::distinct()`
  dplyr::distinct(symbol_id, .keep_all = TRUE)
filtered_mapped_matrix <- filtered_mapped_df %>%
  # GSVA can't the Ensembl IDs so we should drop this column as well as the means
  dplyr::select(-PROBEID, -gene_means) %>%
  # We need to store our gene identifiers as row names
  tibble::column_to_rownames("symbol_id") %>%
  # Now we can convert our object into a matrix
  as.matrix()

cong.matrix = data.frame(Control = rowMeans(filtered_mapped_matrix[,c(1,3)]),
                         Silencing = rowMeans(filtered_mapped_matrix[,c(2,4)]))
colnames(cong.matrix) = c("sirNEG", "sir-cFBXW7")
cong.matrix = filtered_mapped_matrix
cong.matrix["NR3C1",]
Sul_patients = norm.counts.m %>% melt() %>% merge(meta_tall, by.x = "sample_id", by.y = "sample_id") %>% 
  dplyr::group_by(SYMBOL, circKD) %>%
  dplyr::summarize(mean_count = mean(value, na.rm = TRUE)) %>% dcast(SYMBOL~circKD)
colnames(Sul_patients) = c("SYMBOL", paste0(colnames(Sul_patients)[-1],"_Solulier"))
Mul_patients = norm.counts.m %>% merge(meta_validation_complete, by.x = "sample_id", by.y = "sample_id") %>% 
  dplyr::group_by(SYMBOL, cFBXW7_group) %>%
  dplyr::summarize(mean_count = mean(norm.count, na.rm = TRUE)) %>% reshape2::dcast(SYMBOL~cFBXW7_group)
colnames(Mul_patients) = c("SYMBOL", paste0(colnames(Mul_patients)[-1],"_Mullighan"))
conc.patsil <- merge_pat[merge_pat$adj.P.Val<=0.1,]
conc.patsil <- merge_pat[merge_pat$adj.P.Val<=0.1 & abs(merge_pat$logFC_patients)>=0.2,]
patients = merge(Sul_patients, Mul_patients, by = "SYMBOL", all.x = T, all.y = T) %>% filter(SYMBOL%in%rownames(cong.matrix))
patients_Mul = Mul_patients %>% filter(SYMBOL%in%rownames(cong.matrix))
patients_Sul = Sul_patients %>% filter(SYMBOL%in%rownames(cong.matrix))
pat.matrix = as.matrix(patients[,-1])
pat.matrix["NR3C1",]
pat.matrix_M = as.matrix(patients_Mul[,-1])
pat.matrix_S = as.matrix(patients_Sul[,-1])
pat.entrez = data.frame("entrez" = mapIds(org.Hs.eg.db, as.character(patients$SYMBOL),
                              keytype = "SYMBOL",
                              column = "ENTREZID"))
entrex_to_symbol <- data.frame(symbol = mapIds(org.Hs.eg.db, rownames(cong.matrix),
                                              keytype = "ENTREZID",
                                              column = "SYMBOL"))
rownames(cong.matrix) = entrex_to_symbol$symbol
colnames(cong.matrix) = c("sirNEG", "sir-cFBXW7", "sirNEG", "sir-cFBXW7")

rownames(pat.matrix) = patients$SYMBOL #pat.entrez$entrez
pat.matrix = pat.matrix[!is.na(rownames(pat.matrix)),]
rownames(pat.matrix_S) = patients_Sul$SYMBOL #pat.entrez$entrez
pat.matrix_S = pat.matrix_S[!is.na(rownames(pat.matrix_S)),]
rownames(pat.matrix_M) = patients_Mul$SYMBOL #pat.entrez$entrez
pat.matrix_M = pat.matrix_M[!is.na(rownames(pat.matrix_M)),]

gsva_results_conc_sil <- gsva(
  as.matrix(cong.matrix),#[rownames(pat.matrix),],
  hallmarks_list, #$DANG_BOUND_BY_MYC,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Poisson",
  # Minimum gene set size
  min.sz = 5,
  # Maximum gene set size
  max.sz = 1000,
  # Compute Gaussian-distributed scores
  mx.diff = FALSE,
  # Don't print out the progress bar
  verbose = FALSE
)
save(gsva_results_conc_sil, file = "/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/gsva_NOTCH1_sil_RPMI.RData")
save(gsva_results_conc_sil, file = "/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/gsva_HM_sil_RPMI.RData")
load("/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/gsva_HM_sil_RPMI.RData")
gsva_results_pat_S <- gsva(
  pat.matrix_S,
  hallmarks_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Poisson",
  # Minimum gene set size
  min.sz = 2,
  # Maximum gene set size
  max.sz = 1000,
  # Compute Gaussian-distributed scores
  mx.diff = FALSE,
  # Don't print out the progress bar
  verbose = FALSE
)
gsva_results_pat_M <- gsva(
  pat.matrix_M,
  hallmarks_list,
  method = "gsva",
  # Appropriate for our vst transformed data
  kcdf = "Poisson",
  # Minimum gene set size
  min.sz = 2,
  # Maximum gene set size
  max.sz = 1000,
  # Compute Gaussian-distributed scores
  mx.diff = FALSE,
  # Don't print out the progress bar
  verbose = FALSE
)
hallmarks_list$"DANG_BOUND_BY_MYC"
hallmarks_list$DANG_MYC_TARGETS_UP

gsva_results_conc_sil[rownames(gsva_results_conc_sil)%in%"GSE5463_CTRL_VS_DEXAMETHASONE_TREATED_THYMOCYTE_DN",]
gsva_results_conc_sil[rownames(gsva_results_conc_sil)%in%"GSE5463_CTRL_VS_DEXAMETHASONE_TREATED_THYMOCYTE_UP",]

annot_df <- pData(rawData) %>%
  # We need the sample IDs and the main column that contains the metadata info
  dplyr::select(
    group,
    experiment) %>% dplyr::rename(Treatment = group,
                           Experiment = experiment)
annot_df_cong = data.frame("Treatment" = c("sirNEG", "sir-cFBXW7"))
annot_df_cong = data.frame("Treatment" = rep(c("sirNEG", "sir-cFBXW7"), 2))
rownames(annot_df_cong) = colnames(cong.matrix)

annot_df_cong = data.frame("CircFBXW7" = c("High", "Low"))
rownames(annot_df_cong) = colnames(pat.matrix)
  # ) %>%
  # # Create our `time_point` variable based on `refinebio_title`
  # dplyr::mutate(
  #   time_point = dplyr::case_when(
  #     # Create our new variable based whether the refinebio_title column
  #     # contains _AV_ or _CV_
  #     stringr::str_detect(refinebio_title, "_AV_") ~ "acute illness",
  #     stringr::str_detect(refinebio_title, "_CV_") ~ "recovering"
  #   )
  # ) %>%
  # # We don't need the older version of the variable anymore
  # dplyr::select(-refinebio_title)

ann_colors <- list(
  Treatment = c(`sirNEG` = "#C0C0C0", `sir-cFBXW7` = "#FF9933")
  # Experiment = c("1" = "blue4", "2" = "cadetblue2")
)
ann_colors <- list(
  CircFBXW7 = c(High = "chartreuse4", Low = "burlywood3")
  # Experiment = c("1" = "blue4", "2" = "cadetblue2")
)
rownames(gsva_results_conc_sil) = substr(gsub("_", " ", gsub("^HALLMARK_", "", rownames(gsva_results_conc_sil))), 1, 35)
merge.dt = merge(gsva_results_pat_M[rownames(gsva_results_pat_M)%like%"NOTCH",],
      gsva_results_pat_S[rownames(gsva_results_pat_S)%like%"NOTCH",], by = "row.names")
plot.dt = merge.dt[,-1]
rownames(plot.dt) = merge.dt$Row.names
pathway_heatmap <- pheatmap::pheatmap(gsva_results_conc_sil[rownames(gsva_results_conc_sil)%like%"MYC",], 
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
