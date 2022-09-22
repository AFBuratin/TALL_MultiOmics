# BiocManager::install("ViSEAGO")
library(ViSEAGO)
# load genes background
background <- back_genes[!is.na(back_genes)]
# load gene selection
selection <- target_genes$SYMBOL[!is.na(target_genes$SYMBOL)]
table <- table_Sil[abs(table_Sil$logFC)>=0.1,c("SYMBOL","logFC")]
selection <- table_Sil$SYMBOL[table_Sil$logFC>=0.3]
# rank gene identifiers according statistical value
data.table::setorder(table,logFC)

## Not run: 
# connect to Bioconductor
Bioconductor<-ViSEAGO::Bioconductor2GO()

# load GO annotations from Bioconductor
myGENE2GO<-ViSEAGO::annotate(
  "org.Hs.eg.db",
  Bioconductor
)

entrez_ids <- mapIds(clariomshumantranscriptcluster.db, 
                     keys = table_Sil$PROBEID, 
                     keytype = "PROBEID",
                     column = "ENTREZID")
table_Sil$ENTREZID = entrez_ids[table_Sil$PROBEID]
# create topGOdata for BP
BP<-ViSEAGO::create_topGOdata(
  # geneSel=target_genes$ENTREZID[!is.na(target_genes$ENTREZID)],
  geneSel=table_Sil$ENTREZID[table_Sil$logFC>=0.3][!is.na(table_Sil$ENTREZID[table_Sil$logFC>=0.3])],
  allGenes=table_Sil$ENTREZID[!is.na(table_Sil$ENTREZID)],
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5
)
classic<-topGO::runTest(
  BP,
  algorithm = "classic",
  statistic = "fisher"
)
# merge results from topGO
BP_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("BP","classic")
  )
)


# initialyse 
myGOs<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults
)

# compute all available Semantic Similarity (SS) measures
myGOs<-ViSEAGO::compute_SS_distances(
  myGOs,
  distance="Wang"
)


# GOterms heatmap with the default parameters
Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=FALSE,
  showGOlabels=TRUE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic = list(
        pamStage = TRUE,
        pamRespectsDendro = TRUE,
        deepSplit = 2,
        minClusterSize = 5
      )
    )
  ),
  samples.tree=NULL
)

# Display the clusters-heatmap
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms"
)

# GOclusters heatmap
# Create GOterms heatmap
Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=FALSE,
  showGOlabels =TRUE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =2
      )
    )
  ),
  samples.tree=NULL
)

Wang_clusters_wardD2<-ViSEAGO::compute_SS_distances(
  Wang_clusters_wardD2,
  distance="BMA"
)

Wang_clusters_wardD2<-ViSEAGO::GOclusters_heatmap(
  Wang_clusters_wardD2,
  tree=list(
    distance="BMA",
    aggreg.method="ward.D2"
  )
)

ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOclusters"
)

# display the heatmap
ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOclusters"
)

# # print the clusters-heatmap
# ViSEAGO::show_heatmap(
#   Wang_clusters_wardD2,
#   "GOterms",
#   "cluster_heatmap_Wang_wardD2.png"
# )

# Display the clusters-heatmap table
ViSEAGO::show_table(Wang_clusters_wardD2)

# # Print the clusters-heatmap table
# ViSEAGO::show_table(
#   Wang_clusters_wardD2,
#   "cluster_heatmap_Wang_wardD2.xls"
# )

# display colored MDSplot
ViSEAGO::MDSplot(
  Wang_clusters_wardD2,
  "GOclusters"
)

# # print colored MDSplot
# ViSEAGO::MDSplot(
#   Wang_clusters_wardD2,
#   "GOterms",
#   file="mdsplot2.png"
# )


