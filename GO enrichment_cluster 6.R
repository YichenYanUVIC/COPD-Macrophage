library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
#data loading
{
  #FC score list with gene names (cluster 6)
  data_enrichment_featured_ALL <- readxl::read_xlsx("/Users/richard.yan/Desktop/COPD latest/feature_genes_log2FC_1.5.xlsx")
  data_enrichment_featured_ALL <- as.data.frame(data_enrichment_featured_ALL)
  
  FCList_6 <- data_enrichment_featured_ALL[which(data_enrichment_featured_ALL$cluster == 15), "avg_log2FC"]
  names(FCList_6) <- data_enrichment_featured_ALL[which(data_enrichment_featured_ALL$cluster == 15), "gene"]
  FCList_6 <- sort(FCList_6, decreasing = T)
  FCList_6
  
  #feature genes of cluster 6 (gene names)
  gene_enrichment_6 <- data_enrichment_featured_ALL[which(data_enrichment_featured_ALL$cluster == 15), "gene"]
}


#KEGG
KEGG_data <- bitr(gene_enrichment_6, 
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
setwd("/Users/richard.yan/Desktop/COPD\ latest/GO\ enrichment\ clusters/enrichment\ result\ data/KEGG\ pathway")
write.csv(KEGG_result, file = "KEGG_pathway_15.csv",
          row.names = F)

library(enrichplot)
barplot(KEGG_result)
cnetplot(KEGG_result)
emapplot(pairwise_termsim(KEGG_result),
         showCategory = 10)

browseKEGG(KEGG_result, 
           pathID = "hsa05133") #set pathID according to the results



#GO enrichment for feature genes of cluster 6 
{
  #enrichGO (ALL)
  enrichGO_6 <- enrichGO(gene = gene_enrichment_6,
                         OrgDb = org.Hs.eg.db,
                         keyType = "SYMBOL",
                         ont = "ALL",
                         pAdjustMethod = "BH", #default
                         pvalueCutoff = 0.05, #default
                         qvalueCutoff = 0.2) #default
  head(enrichGO_6)
  dim(enrichGO_6)
  
  #Results divided by category
  #data.frame, cannot use for analysis in functions
  {
    enrichGO_6_BP_result <- enrichGO_6[which(enrichGO_6$ONTOLOGY == "BP"), -c(1,2,6,10)]
    enrichGO_6_CC_result <- enrichGO_6[which(enrichGO_6$ONTOLOGY == "CC"), -c(1,2,6,10)]
    enrichGO_6_MF_result <- enrichGO_6[which(enrichGO_6$ONTOLOGY == "MF"), -c(1,2,6,10)]
  }
  
  setwd("/Users/richard.yan/Desktop/COPD\ latest/GO\ enrichment\ clusters/enrichment\ result\ data")
  {
    write.csv(enrichGO_6_BP_result, file = "enrichGO_15_BP.csv") #Biological process基因产物参与的生物路径或机制
    write.csv(enrichGO_6_CC_result, file = "enrichGO_15_CC.csv") #Cellular component基因产物在细胞内外的位置
    write.csv(enrichGO_6_MF_result, file = "enrichGO_15_MF.csv") #Molecular function基因产物分子层次的功能
  }

  #Simplify GO enrichment results (MF category)
  enrichGO_6_MF <- enrichGO(gene = gene_enrichment_6,
                            OrgDb = org.Hs.eg.db,
                            keyType = "SYMBOL",
                            ont = "MF",
                            pAdjustMethod = "BH", #default
                            pvalueCutoff = 0.05, #default
                            qvalueCutoff = 0.2) #default

  enrichGO_6_MF_simp <- simplify(enrichGO_6_MF,
                                 cutoff = 0.7, #default
                                 by = "p.adjust",
                                 select_fun = min,
                                 measure = "Wang")
  head(enrichGO_6_MF_simp)
  dim(enrichGO_6_MF_simp) # 31减少到19
}


#GO Gene Set Enrichment Analysis (GSEA)
{
  GSEA_GO_6 <- gseGO(FCList_6,
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
                     #nPerm = 1000, # by DOSE
                     by = "fgsea")

  head(GSEA_GO_6, 3)
  head(GSEA_GO_6[,2], 10)
  write.csv(GSEA_GO_6, file = "GSEA_GO_15.csv")
}


#Visualization
{
  library(enrichplot)
  
  library(topGO)
  library(Rgraphviz)
  plotGOgraph(enrichGO_6_MF)
  
  #barplot
  barplot(enrichGO_6_MF_simp,
          showCategory = 10)
  
  #dotplot
  library(ggplot2)
  dotplot(enrichGO_6, #ALL
          showCategory = 10,
          title = "Top10 GO terms of each sub-class (BP, CC, MF)",
          split = "ONTOLOGY") +
    facet_grid(ONTOLOGY~.,
               scale = "free")
  #Enrichment map
  emapplot(pairwise_termsim(enrichGO_6_MF_simp),
           showCategory = 10)
  
  #Gene-Concept Network
  cnetplot(enrichGO_6_MF_simp,
           showCategory = 5,
           node_label = "all") #"category" / "gene" / "all" / "none"
  cnetplot(enrichGO_6, #ALL
           showCategory = 5,
           foldChange = FCList_6,
           circular = TRUE,
           colorEdge = TRUE)
  
  #UpSet Plot
  library(ggupset)
  upsetplot(enrichGO_6_MF_simp) #着重于不同基因集间基因的重叠情况
  upsetplot(GSEA_GO_6, n = 10) #对于GSEA结果 将绘制不同类别的fold change分布
  
  #Heatmap-like functional classification
  heatplot(enrichGO_6_MF_simp)
  heatplot(enrichGO_6,
           foldChange = FCList_6)
  
  #GO DAG graph
  goplot(enrichGO_6_MF_simp,
         showCategory = 10)
  
  #Ridgeline plot for expression distribution of GSEA result
  ridgeplot(GSEA_GO_6)
  
  #Running score and pre-ranked list of GSEA result
  p001 <- gseaplot(GSEA_GO_6, geneSetID = 1, by = "runningScore", title = GSEA_GO_6$Description[1])
  p002 <- gseaplot(GSEA_GO_6, geneSetID = 1, by = "preranked", title = GSEA_GO_6$Description[1])
  p003 <- gseaplot(GSEA_GO_6, geneSetID = 1, title = GSEA_GO_6$Description[1]) #combined p001 & p002
  cowplot::plot_grid(p001, p002, p003, ncol = 2, labels = LETTERS[1:3])
  
  p004 <- gseaplot2(GSEA_GO_6, geneSetID = 3, title = GSEA_GO_6$Description[3])
  p005 <- gseaplot2(GSEA_GO_6, geneSetID = 2:4, subplots = 1)
  p006 <- gseaplot2(GSEA_GO_6, geneSetID = 2:4, pvalue_table = TRUE,
                    color = c("#E495A5", "#86B875", "#7DB0DD"), 
                    ES_geom = "dot")
  cowplot::plot_grid(p004, p005, p006, ncol = 1, labels = LETTERS[1:3])
}



