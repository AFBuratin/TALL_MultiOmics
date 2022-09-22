table_Sil = read.csv("/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/Table_Sil_RPMI.csv")
table_Sil = tab_sil_KOPTK1
## sign genes
DE_genes <- subset(table_Sil, adj.P.Val <= 0.1)
geneList <- pat_enrich$logFC_mean
geneList = DE_genes$logFC
names(geneList) <- DE_genes$SYMBOL#circ_to_genes$gene_names[match(rownames(res[which(padj<0.05),]), circ_to_genes$circ_id)]
names(geneList) = pat_enrich$SYMBOL
geneList = sort(geneList, decreasing = TRUE)
gene <- names(geneList)[!is.na(names(geneList))]#[abs(geneList)>=0.5]

library(clusterProfiler)
library(enrichplot)
library(DOSE)
# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
library(org.Hs.eg.db)
AnnotationDbi::species(org.Hs.eg.db)

ego <- clusterProfiler::enrichGO(gene          = gene,
                                 universe      = merge_pat_conc2$SYMBOL[!is.na(merge_pat_conc2$SYMBOL)],
                                 OrgDb         = "org.Hs.eg.db",
                                 keyType       = "SYMBOL",
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE)
head(ego)
gse <- clusterProfiler::gseGO(geneList=geneList, 
                              ont ="All", 
                              keyType = "SYMBOL", 
                              nPerm = 10000,
                              minGSSize = 3,
                              maxGSSize = 800,
                              pvalueCutoff = 0.1,
                              verbose = TRUE, 
                              OrgDb = org.Hs.eg.db, 
                              pAdjustMethod = "none")
# save(gse, file = "/media/Data/TALL/FBXW7/Array_sicFBXW7_KOPTK1/gse_all.rda")
# write.csv(gse@result, "/media/Data/TALL/FBXW7/GSEA_meanlogFC_concPats.csv")
gse_old_res = read.csv("/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/gse.csv")
gse_old_res
gse@result[order(gse@result$enrichmentScore),][order(gse@result$NES),]
gse@result[gse@result$Description=="regulation of cell cycle G1/S phase transition",]

dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
library(data.table)
library(ggh4x)
activated = gse@result[gse@result$Description %like% 'G1/S|cell cycle|centrosome|microtubule',] %>% as.data.table() %>% 
  arrange(enrichmentScore) %>% top_n(20) %>% select(Description)

png(filename = "/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/gse_Silencing_ACT.png",
    width = 20, height = 22, units = "cm", res = 200)
dotplot(gse, showCategory=activated$Description, split=".sign") + facet_wrap(.~.sign, nrow = 1) +
  force_panelsizes(cols = c(0.1, 0.1)) +
  theme_bw() + # Just to show you can still add things
  theme(legend.position = "right",
        text = element_text(size=15),
        axis.title=element_blank(),
        # axis.ticks.length.y = unit(3, "cm"),
        axis.ticks=element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text.x=element_text(size=15)) 
dev.off()

suppressed = gse@result %>% filter(enrichmentScore<0) %>% as.data.table() %>% 
  arrange(enrichmentScore) %>% top_n(20, -p.adjust) %>% select(Description)

png(filename = "/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/gse_Silencing_SUPP.png",
    width = 20, height = 22, units = "cm", res = 200)
dotplot(gse, showCategory=suppressed$Description, split=".sign") + facet_wrap(.~.sign, nrow = 1) +
  force_panelsizes(cols = c(0.1, 0.1)) +
  theme_bw() + # Just to show you can still add things
  theme(legend.position = "right",
        text = element_text(size=15),
        axis.title=element_blank(),
        # axis.ticks.length.y = unit(3, "cm"),
        axis.ticks=element_blank(),
        strip.text.x = element_text(size = 20),
        axis.text.x=element_text(size=15)) 
dev.off()

load("/media/Data/TALL/FBXW7/Array_sicFBXW7_KOPTK1/data_final.rda")
back_genes_idx <- genefilter::genefinder(data_final,  
                                         as.character(DE_genes$PROBEID), 
                                         method = "manhattan", scale = "none")
back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)
back_genes <- featureNames(data_final)[back_genes_idx]
# back_genes <- setdiff(back_genes, DE_genes$PROBEID)
intersect(back_genes, DE_genes$PROBEID)

gene_IDs <- table_Sil$PROBEID
in_universe <- gene_IDs %in% c(DE_genes$PROBEID, back_genes)
in_selection <- gene_IDs %in% DE_genes$PROBEID 

all_genes <- in_selection[in_universe]
all_genes <- factor(as.integer(in_selection[in_universe]))
names(all_genes) <- gene_IDs[in_universe] 
sum(all_genes==1)

# BiocManager::install("topGO")
library(topGO)
top_GO_data <- new("topGOdata", description = "GO annotaion DEGs of circFBXW7 silencing in RPMI cell lines", 
                   ontology = "BP", allGenes = all_genes,
                   nodeSize = 10, annot = annFUN.db, affyLib = "clariomshumantranscriptcluster.db")
# save(top_GO_data, file = "/media/Data/TALL/FBXW7/Array_sicFBXW7_KOPTK1/top_go_data.rda")
result_top_GO_elim <- 
  runTest(top_GO_data, algorithm = "elim", statistic = "Fisher")
result_top_GO_classic <- 
  runTest(top_GO_data, algorithm = "classic", statistic = "Fisher")
res_top_GO <- GenTable(top_GO_data, Fisher.elim = result_top_GO_elim,
                       Fisher.classic = result_top_GO_classic,
                       orderBy = "Fisher.elim" , topNodes = 20)

genes_top_GO <- printGenes(top_GO_data, whichTerms = res_top_GO$GO.ID,
                           chip = "clariomshumantranscriptcluster.db", geneCutOff = 1000)
library(stringr)
res_top_GO$sig_genes <- sapply(genes_top_GO, function(x){
  str_c(paste0(x[x$'raw p-value' == 2, "Symbol.id"],";"), 
        collapse = "")
})

head(res_top_GO, 20)

allRes <- GenTable(top_GO_data, classic = result_top_GO_classic, topNodes = 10)
allRes

library(Rgraphviz)
png(filename = "/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/Graph_BP_DEGs.png",
    width = 20, height = 25, units = "cm", res = 300)
showSigOfNodes(top_GO_data, score(result_top_GO_classic), firstSigNodes = 5, useInfo = 'all')
dev.off()
score(result_top_GO_classic)[1:10]
top_GO_data@graph@nodeData
printGraph(top_GO_data, result_top_GO_classic, firstSigNodes = 5, fn.prefix = "sampleFile", useInfo = "all", pdfSW = TRUE)


entrez_ids <- mapIds(clariomshumantranscriptcluster.db, 
                     keys = DE_genes$PROBEID, 
                     keytype = "PROBEID",
                     column = "ENTREZID")
gene_enrezID <- mapIds(clariomshumantranscriptcluster.db, 
                       keys = table_Sil$PROBEID, 
                       keytype = "PROBEID",
                       column = "ENTREZID")
DE_genes$ENTREZID = entrez_ids[match(DE_genes$PROBEID, names(entrez_ids))]
table_Sil$ENTREZID = gene_enrezID

library(ReactomePA)
library(DOSE)

de = DE_genes$ENTREZID#[match(names(target_geneList), DE_genes$SYMBOL)]
de = de[!is.na(de)]

geneList_entrez = DE_genes$logFC[!is.na(DE_genes$ENTREZID)]
names(geneList_entrez) = de

geneList_entrez = sort(geneList_entrez, decreasing = TRUE)

edo <- enrichDGN(de)
library(enrichplot)
library(ggplot2)

barplot(edo, showCategory=20) 

edo2 <- gseDO(geneList_entrez, pAdjustMethod = "none")
edo2@result$Description
dotplot(gse, showCategory=30) + facet_grid(.~.sign) + ggtitle("dotplot for GSEA")

## convert gene ID to Symbol
library(ggnewscale)
# edox <- setReadable(x = gse, 'org.Hs.eg.db', 'SYMBOL')
categorys <- c("mitotic cell cycle",
               "cell cycle",
               "regulation of cell cycle process",
               "regulation of cell cycle G1/S phase transition")
png(filename = "/media/Data/TALL/FBXW7/Array_sicFBXW7_RPMI/CNETplot_BP_DEGs.png",
    width = 50, height = 45, units = "cm", res = 300)
cnetplot(gse, showCategory = categorys, foldChange = geneList)
dev.off()

p1 <- cnetplot(gse, foldChange=geneList_entrez)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(gse, categorySize="pvalue", foldChange=geneList)
gse@result[1:5,]

png(filename = "/media/Data/TALL/FBXW7/Array_sicFBXW7_KOPTK1/CNETplot_BP_DEGs.png",
    width = 28, height = 25, units = "cm", res = 300)
cnetplot(gse, showCategory = categorys, 
         foldChange=geneList, circular = TRUE, 
         layout = "kk", cex_category = 0.5, colorEdge= FALSE)

dev.off()

edox2 <- pairwise_termsim(gse)
treeplot(edox2, color = geneList)
treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

heatplot(gse, showCategory=5)
heatplot(gse, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

kk <- enrichKEGG(gene         = de,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1)
head(kk)
kk2 <- clusterProfiler::gseKEGG(geneList     = geneList_entrez,
                                organism     = 'hsa',
                                minGSSize    = 3,
                                pvalueCutoff = 0.1,
                                verbose      = FALSE,
                                pAdjustMethod = "none")
head(kk2)
library(enrichplot)
upsetplot(kk)
dotplot(kk)

reactome_enrich <- enrichPathway(gene = de,#[DE_genes$PROBEID], 
                                 universe = table_Sil$ENTREZID,
                                 organism = "human",
                                 pvalueCutoff = 0.1,
                                 # qvalueCutoff = 0.9,
                                 pAdjustMethod = "none",
                                 readable = TRUE)

reactome_enrich@result$Description <- paste0(str_sub(
  reactome_enrich@result$Description, 1, 20),
  "...")

as.data.frame(reactome_enrich)[1:20,1:7]

barplot(reactome_enrich)
