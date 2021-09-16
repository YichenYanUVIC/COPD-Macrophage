library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
#data loading
{
  #FC score list with gene names in all clusters
  data_enrichment_featured_ALL <- readxl::read_xlsx("./data/feature_genes_log2FC_1.5.xlsx")
  data_enrichment_featured_ALL <- as.data.frame(data_enrichment_featured_ALL)
  
  #change the "cluster == n" for certain cluster (n = 0, 1, 2, ..., 15)
  FCList <- data_enrichment_featured_ALL[which(data_enrichment_featured_ALL$cluster == 15), "avg_log2FC"]
  names(FCList) <- data_enrichment_featured_ALL[which(data_enrichment_featured_ALL$cluster == 15), "gene"]
  FCList <- sort(FCList, decreasing = T)
  FCList
  
  #feature genes of cluster n (n = 0, 1, 2, ... ,15) (gene names)
  gene_enrichment <- data_enrichment_featured_ALL[which(data_enrichment_featured_ALL$cluster == 15), "gene"]
}


#KEGG
KEGG_data <- bitr(gene_enrichment, 
                  fromType = "SYMBOL", 
                  toType = "ENTREZID", 
                  OrgDb = "org.Hs.eg.db")

KEGG_result <- enrichKEGG(gene = KEGG_data$ENTREZID,
                          organism = "hsa",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05)

KEGG_result <- setReadable(KEGG_result,
                           OrgDb = "org.Hs.eg.db",
                           keyType = "ENTREZID")
head(KEGG_result)
#change n with the certain cluster you chose
write.csv(KEGG_result, file = "KEGG_pathway_n.csv",
          row.names = F)

library(enrichplot)
barplot(KEGG_result)
cnetplot(KEGG_result)
emapplot(pairwise_termsim(KEGG_result),
         showCategory = 10)

browseKEGG(KEGG_result, 
           pathID = "hsa05133") #set pathID according to the results


#GO enrichment for feature genes of cluster n
{
  #enrichGO (ALL)
  enrichGO <- enrichGO(gene = gene_enrichment,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "ALL",
                       pAdjustMethod = "BH", #default
                       pvalueCutoff = 0.05, #default
                       qvalueCutoff = 0.2) #default
  head(enrichGO)
  dim(enrichGO)
  
  #Results divided by category
  #data.frame, cannot use for analysis in functions
  {
    enrichGO_BP_result <- enrichGO[which(enrichGO$ONTOLOGY == "BP"), -c(1,2,6,10)]
    enrichGO_CC_result <- enrichGO[which(enrichGO$ONTOLOGY == "CC"), -c(1,2,6,10)]
    enrichGO_MF_result <- enrichGO[which(enrichGO$ONTOLOGY == "MF"), -c(1,2,6,10)]
  }
  
  #change n with the certain cluster you chose
  {
    write.csv(enrichGO_BP_result, file = "enrichGO_n_BP.csv") #Biological process
    write.csv(enrichGO_CC_result, file = "enrichGO_n_CC.csv") #Cellular component
    write.csv(enrichGO_MF_result, file = "enrichGO_n_MF.csv") #Molecular function
  }

  #Simplify GO enrichment results (MF category)
  enrichGO_MF <- enrichGO(gene = gene_enrichment,
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",
                          ont = "MF",
                          pAdjustMethod = "BH", #default
                          pvalueCutoff = 0.05, #default
                          qvalueCutoff = 0.2) #default

  enrichGO_MF_simp <- simplify(enrichGO_MF,
                               cutoff = 0.7, #default
                               by = "p.adjust",
                               select_fun = min,
                               measure = "Wang")
  head(enrichGO_MF_simp)
}


#GO Gene Set Enrichment Analysis (GSEA)
{
  GSEA_GO <- gseGO(FCList,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   keyType = "SYMBOL",
                   minGSSize = 10, #minimal size of each geneSet for analyzing
                   maxGSSize = 500, #maximal size of genes annotated for testing
                   pvalueCutoff = 0.05,
                   exponent = 1,
                   eps = 0, #warning recommendation
                   pAdjustMethod = "BH",
                   verbose = TRUE, #print message or not
                   #nPerm = 1000, # by DOSE only
                   by = "fgsea")

  head(GSEA_GO, 3)
  head(GSEA_GO[,2], 10)
  #change n with the certain cluster you chose
  write.csv(GSEA_GO, file = "GSEA_GO_n.csv")
}


#Visualization
{
  library(enrichplot)
  
  library(topGO)
  library(Rgraphviz)
  plotGOgraph(enrichGO_MF)
  
  #barplot
  barplot(enrichGO_MF_simp,
          showCategory = 10)
  
  #dotplot
  library(ggplot2)
  dotplot(enrichGO, #ALL
          showCategory = 10,
          title = "Top10 GO terms of each sub-class (BP, CC, MF)",
          split = "ONTOLOGY") +
    facet_grid(ONTOLOGY~.,
               scale = "free")
  #Enrichment map
  emapplot(pairwise_termsim(enrichGO_MF_simp),
           showCategory = 10)
  
  #Gene-Concept Network
  cnetplot(enrichGO_MF_simp,
           showCategory = 5,
           node_label = "all") #"category" / "gene" / "all" / "none"
  cnetplot(enrichGO, #ALL
           showCategory = 5,
           foldChange = FCList,
           circular = TRUE,
           colorEdge = TRUE)
  
  #UpSet Plot
  library(ggupset)
  upsetplot(enrichGO_MF_simp) #着重于不同基因集间基因的重叠情况
  upsetplot(GSEA_GO, n = 10) #对于GSEA结果 将绘制不同类别的fold change分布
  
  #Heatmap-like functional classification
  heatplot(enrichGO_MF_simp)
  heatplot(enrichGO,
           foldChange = FCList)
  
  #GO DAG graph
  goplot(enrichGO_MF_simp,
         showCategory = 10)
  
  #Ridgeline plot for expression distribution of GSEA result
  ridgeplot(GSEA_GO)
}



