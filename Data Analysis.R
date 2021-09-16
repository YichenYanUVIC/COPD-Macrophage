library(Seurat)
library(monocle)

#Data loading
{
  raw <- readRDS("./data/Clustering_Data_136831.Rds") 
  Disease_Barcode <- readRDS("./data/Clustering_Labels_136831.Rds")
}

#Seurat Object Building
{
  Seurat_Validation <- CreateSeuratObject(raw, 
                                          min.cells = 5,
                                          min.features = 2000) #35603 genes 56318 cells
  cell_type <- as.factor(Disease_Barcode[which(Disease_Barcode[,1] %in% Seurat_Validation@assays$RNA@counts@Dimnames[[2]]),2])
  disease_type <- as.factor(Disease_Barcode[which(Disease_Barcode[,1] %in% Seurat_Validation@assays$RNA@counts@Dimnames[[2]]),3])
  Seurat_Validation@meta.data <- cbind(Seurat_Validation@meta.data, disease_type, cell_type)
}

#QC: remove high ^MT-  
{
  Seurat_Validation[["percent.mt"]] <- PercentageFeatureSet(Seurat_Validation, 
                                                            pattern = "^MT-")
  Seurat_Validation <- subset(Seurat_Validation, 
                              subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  dim(Seurat_Validation) #35603 genes 7929 cells
}

#Normalization and HVG selection
{
  #Normalization
  Seurat_Validation <- NormalizeData(Seurat_Validation, 
                                     normalization.method = "LogNormalize", 
                                     scale.factor = 10000)
  #HVG
  Seurat_Validation <- FindVariableFeatures(Seurat_Validation, 
                                            selection.method = "vst", 
                                            nfeatures = 1000)
}

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(Seurat_Validation), 20)
# plot variable features with labels
plot1 <- VariableFeaturePlot(Seurat_Validation)
plotHVG <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plotHVG

#Scaling the data
all.genes <- rownames(Seurat_Validation)
Seurat_Validation <- ScaleData(Seurat_Validation, 
                               features = all.genes)

#Perform linear dimensional reduction (PCA)
Seurat_Validation <- RunPCA(Seurat_Validation, 
                            npcs = 30,
                            features = VariableFeatures(object = Seurat_Validation))

#Examine and visualize PCA results in a few different ways
print(Seurat_Validation[["pca"]], 
      dims = 1:5, 
      nfeatures = 5)
VizDimLoadings(Seurat_Validation, 
               dims = 1:2, 
               reduction = "pca")
DimPlot(Seurat_Validation, 
        reduction = "pca",
        group.by = 'disease_type')
DimPlot(Seurat_Validation, 
        reduction = "pca",
        group.by = 'cell_type')
DimHeatmap(Seurat_Validation, 
           dims = 1, 
           cells = 500, 
           balanced = TRUE)
ElbowPlot(Seurat_Validation, 
          ndims = 30, 
          reduction = "pca")

#Cluster parameters setting
Seurat_Validation <- FindNeighbors(Seurat_Validation, dims = 1:20)
Seurat_Validation <- FindClusters(Seurat_Validation, resolution = 0.8)  # default 0.8
#Look at cluster IDs of the first 5 cells
head(Idents(Seurat_Validation), 5)

#UMAP
Seurat_Validation <- RunUMAP(Seurat_Validation, dims = 1:20, label = T)

#UMAP Visualization
DimPlot(Seurat_Validation, reduction = "umap", 
        shuffle = T, seed = 1,
        label = T, label.size = 10, repel = F,
        group.by = 'ident', pt.size = 0.5) + NoAxes() + NoLegend()

#group by disease type
DimPlot(Seurat_Validation, reduction = "umap", shuffle = T,
        group.by = 'disease_type') + NoAxes()

DimPlot(Seurat_Validation, reduction = "umap",
        split.by = 'disease_type',
        label = T, label.size = 7) + NoAxes()

DimPlot(Seurat_Validation, reduction = "umap",
        group.by = 'disease_type', split.by = 'ident',
        ncol = 4)

#group by cell type
DimPlot(Seurat_Validation, reduction = "umap", shuffle = T,
        group.by = 'cell_type', pt.size = 0.7,
        label = T, label.size = 7) + NoAxes() + NoLegend()

DimPlot(Seurat_Validation, reduction = "umap",
        split.by = 'cell_type',
        label = T, label.size = 6) + NoAxes() + NoLegend()

DimPlot(Seurat_Validation, reduction = "umap",
        group.by = 'cell_type', split.by = 'ident',
        ncol = 4)

#DE
{
  #Doing DE and getting results of all genes (both marker genes and other genes)
  Validation.markers_ALL_136831 <- FindAllMarkers(Seurat_Validation, 
                                                  logfc.threshold = 0,
                                                  only.pos = FALSE, 
                                                  min.pct = -Inf,
                                                  min.diff.pct = -Inf,
                                                  max.cells.per.ident = Inf,
                                                  return.thresh = 1,
                                                  test.use = 'MAST')
  write.csv(Validation.markers_ALL_136831, 
            file = "./data/DE_MAST_ALL_136831.csv")
}

#Top marker genes
library(dplyr)
#load DE file saved
Validation.markers_ALL_136831 <- read.csv("./data/DE_MAST_ALL_136831.csv")
#Select top 5 marker genes of each cluster according to log2FC
top5 <- Validation.markers_ALL_136831 %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

#DE Heatmap with top marker genes
library(DoMultiBarHeatmap)

#group by disease type
DoMultiBarHeatmap(Seurat_Validation, 
                  features = top5$gene,
                  group.by = 'seurat_clusters',
                  additional.group.by = 'disease_type',
                  group.bar = TRUE,
                  draw.lines = T)# + NoLegend()

#group by cell type
DoMultiBarHeatmap(Seurat_Validation, 
                  features = top5$gene,
                  group.by = 'seurat_clusters',
                  additional.group.by = 'cell_type',
                  group.bar = TRUE,
                  draw.lines = T)# + NoLegend()

#Crosstable
{
  DiseaseType <- Seurat_Validation@meta.data$disease_type
  CellType <- Seurat_Validation@meta.data$cell_type
  ClusterResults <- Seurat_Validation@meta.data$seurat_clusters 
}

library(gmodels)
CrossTable(x = DiseaseType, 
           y = ClusterResults, 
           prop.chisq = FALSE)
CrossTable(x = CellType, 
           y = ClusterResults, 
           prop.chisq = FALSE)

#Bar plots of percentage of diseases or cell types in each cluster
DTC <- table(DiseaseType, ClusterResults)
CTC <- table(CellType, ClusterResults)[c(1,4,2,3),]
barplot(DTC,
        col = colors()[c(12,555)], 
        font.axis = 2, 
        beside = TRUE,
        legend = rownames(DTC), 
        args.legend = list(x = "topleft"),
        xlab = "Cluster", 
        ylab = "Cell Number",
        ylim = c(0,1500),
        font.lab = 2)
barplot(CTC,
        col = colors()[c(14,400,99,465)], 
        font.axis = 2, 
        beside = TRUE,
        legend = rownames(CTC), 
        args.legend = list(x = "topleft"),
        xlab = "Cluster", 
        ylab = "Cell Number",
        ylim = c(0,1500),
        font.lab = 2)

#Pseudo-time Analysis with monocle
#Starting from UMAP finished

library(monocle)
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(Seurat_Validation@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = Seurat_Validation@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
pseudo_monocle <- newCellDataSet(data,
                                 phenoData = pd,
                                 featureData = fd,
                                 lowerDetectionLimit = 0.5,
                                 expressionFamily = negbinomial.size())
#Compute SizeFactor & Dispersions
pseudo_monocle <- estimateSizeFactors(pseudo_monocle)
pseudo_monocle <- estimateDispersions(pseudo_monocle)
head(pData(pseudo_monocle))
#saveRDS(pseudo_monocle, "pseudo_monocle.Rds") #save RDS data for time saving

diff_test_res <- differentialGeneTest(pseudo_monocle, 
                                      fullModelFormulaStr = "~seurat_clusters")
#saveRDS(diff_test_res, "diff_test_res.Rds") #save RDS data for time saving

#Apply saved RDS for downstream analysis
diff_test_res <- readRDS("./data/diff_test_res.Rds")
ordering_genes <- row.names(subset(diff_test_res, 
                                   qval < 0.01))

pseudo_monocle <- setOrderingFilter(pseudo_monocle, ordering_genes)
plot_ordering_genes(pseudo_monocle)
pseudo_monocle <- reduceDimension(pseudo_monocle, 
                                  max_components = 2, 
                                  reduction_method = "DDRTree")
#saveRDS(pseudo_monocle, "pseudo_monocle_reducedDimension.Rds") #save RDS data for time saving

#Apply saved RDS for downstream analysis
pseudo_monocle <- readRDS("./data/pseudo_monocle_reducedDimension.Rds")
pseudo_monocle <- orderCells(pseudo_monocle)

#GM_state function: return a fixed state as the root state of pseudo time process
#Here, find a state with the most Monocyte cells as the root state
#Since cluster 8 has the purest Monocyte cells (>99%), 
#Calculating how the cells in cluster 8 distributed in each state,
#Then choose the state with a highest number as the initial root state
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$seurat_clusters)[,"8"] # monocyte majority in cluster 8, set it as root
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
pseudo_monocle <- orderCells(pseudo_monocle, 
                             root_state = GM_state(pseudo_monocle))
plot_cell_trajectory(pseudo_monocle, color_by = "Pseudotime")

#Color in pseudotime plot with various standard
plot_cell_trajectory(pseudo_monocle, color_by = "seurat_clusters") + NoAxes()
plot_cell_trajectory(pseudo_monocle, color_by = "disease_type") + NoAxes()
plot_cell_trajectory(pseudo_monocle, color_by = "cell_type") + NoAxes()
plot_cell_trajectory(pseudo_monocle, color_by = "State") + NoAxes()

#group by State
plot_cell_trajectory(pseudo_monocle, color_by = "State") + NoAxes() + NoLegend() +
  facet_wrap(~State, nrow = 3, strip.position = "bottom") 

#cluster to diseases / cell types
plot_cell_trajectory(pseudo_monocle, color_by = "seurat_clusters") + NoAxes() +
  facet_wrap(~seurat_clusters, nrow = 4, strip.position = "bottom")
plot_cell_trajectory(pseudo_monocle, color_by = "disease_type") + NoAxes() +
  facet_wrap(~seurat_clusters, nrow = 4, strip.position = "bottom")
plot_cell_trajectory(pseudo_monocle, color_by = "cell_type") + NoAxes() +
  facet_wrap(~seurat_clusters, nrow = 4, strip.position = "bottom")

#cell types to diseases
plot_cell_trajectory(pseudo_monocle, color_by = "disease_type") + 
  facet_wrap(~cell_type, nrow = 1)
#disease to cell types
plot_cell_trajectory(pseudo_monocle, color_by = "cell_type") + 
  facet_wrap(~disease_type, nrow = 1)

#state to diseases / cell types
plot_cell_trajectory(pseudo_monocle, color_by = "disease_type") + NoAxes() +
  facet_wrap(~State, nrow = 3, strip.position = "bottom")
plot_cell_trajectory(pseudo_monocle, color_by = "cell_type") + NoAxes() +
  facet_wrap(~State, nrow = 3, strip.position = "bottom")

#tables: how diseases / cell types number vary with the states in pseudotime process
t(table(pData(pseudo_monocle)$State, pData(pseudo_monocle)$disease_type)[c(7,5,6,4,8,2,3,1,9), ])
t(table(pData(pseudo_monocle)$State, pData(pseudo_monocle)$cell_type)[c(7,5,6,4,8,2,3,1,9),
                                                                    c(1,4,2,3)])
CrossTable(y = pData(pseudo_monocle)$State,
           x = pData(pseudo_monocle)$disease_type,
           chisq = FALSE,
           prop.chisq = FALSE,
           prop.c = T,
           prop.r = F,
           prop.t = F)

CrossTable(y = pData(pseudo_monocle)$State,
           x = pData(pseudo_monocle)$cell_type,
           chisq = FALSE,
           prop.chisq = FALSE,
           prop.c = T,
           prop.r = F,
           prop.t = F)

#barplots of the tables data above
#par(mfrow = c(1,2))
barplot_data_celltype <- table(pData(pseudo_monocle)$cell_type, pData(pseudo_monocle)$State)[c(1,4,2,3), c(7,5,6,4,8,2,3,1,9)]
barplot(barplot_data_celltype,
        col = colors()[c(12,89,500,555)], 
        font.axis = 2, 
        beside = TRUE,
        legend = rownames(barplot_data_celltype), 
        args.legend = list(x = "topleft"),
        xlab = "States", 
        ylab = "Cell Number",
        ylim = c(0,1500),
        font.lab = 2)

barplot_data_diseasetype <- table(pData(pseudo_monocle)$disease_type, pData(pseudo_monocle)$State)[,c(7,5,6,4,8,2,3,1,9)]
barplot(barplot_data_diseasetype,
        col = colors()[c(12,555)], 
        font.axis = 2, 
        beside = TRUE,
        legend = rownames(barplot_data_diseasetype), 
        args.legend = list(x = "topleft"),
        xlab = "States", 
        ylab = "Cell Number",
        ylim = c(0,2000),
        font.lab = 2)

#calculating various cells number in each 0.5 pseudo-second gap and plot
#test data (not used) (don't delete!)
{
  head(table(pData(pseudo_monocle)$Pseudotime, pData(pseudo_monocle)$cell_type))
  pseudotime_values <- pData(pseudo_monocle)$Pseudotime # in used
  celltype_values <- as.character(pseudo_monocle$cell_type)
  celltype_pseudotime_plotdata <- data.frame(cbind(pseudotime_values, celltype_values))
  celltype_pseudotime_plotdata$pseudotime_values <- as.numeric(celltype_pseudotime_plotdata$pseudotime_values)
  celltype_pseudotime_plotdata <- celltype_pseudotime_plotdata[order(celltype_pseudotime_plotdata$pseudotime_values),]
  
  max(celltype_pseudotime_plotdata$pseudotime_values) #38.00622
  max(pData(pseudo_monocle)$Pseudotime) #38.00622
  
  #head(table(pData(pseudo_monocle)$Pseudotime, pData(pseudo_monocle)$cell_type))
}

#data in used
XXX <- cbind(pData(pseudo_monocle)$Pseudotime,
             table(pData(pseudo_monocle)$Pseudotime, pData(pseudo_monocle)$cell_type))
rownames(XXX) <- NULL
colnames(XXX) <- c("pseudotimes", "cMonocyte", "Macrophage", "Macrophage_Alveolar", "ncMonocyte")
XXX <- data.frame(XXX)
XXX <- XXX[order(XXX$pseudotimes),]
head(XXX)

#0.5s gap plot with legend
{
  #data for the plot
  XXX_result <- matrix(nrow = 78, ncol = 4) #39s / 0.5 == 78
  for (i in 0:77) {
    XXX_result[i+1,] <- colSums(XXX[which(XXX$pseudotimes >= 0.5*i & XXX$pseudotimes < 0.5*(i+1) ),
                                    c(2,3,4,5)]) 
    
  }
  XXX_result <- data.frame(XXX_result)
  rownames(XXX_result) <- seq(from = 1, to = 78, by = 1)
  colnames(XXX_result) <- c("cMonocyte", "Macrophage", "Macrophage_Alveolar", "ncMonocyte")
  XXX_result <- cbind(seq(from = 0.5, to = 39, by = 0.5), XXX_result)
  colnames(XXX_result) <- c("pseudotimes_0.5_gap", "cMonocyte", "Macrophage", "Macrophage_Alveolar", "ncMonocyte")
  head(XXX_result)
}
{
  #plot
  par(mfrow=c(2,2))
  plot(XXX_result$pseudotimes_0.5_gap, XXX_result$cMonocyte,
       type = "o",
       col = "darkred",
       cex = 0.5,
       pch = 19,
       ylim = c(0,150),
       xlab = "pseudotime (second)",
       ylab = "#cMonocyte per 0.5s")
  legend(x = 0,
         y = 150,
         fill = "darkred",
         legend = "cMonocyte",
         bty = "n",
         cex = 1.5,
         adj = c(0.35,0.5))
  
  plot(XXX_result$pseudotimes_0.5_gap, XXX_result$Macrophage_Alveolar,
       type = "o",
       col = "darkgreen",
       cex = 0.5,
       pch = 19,
       ylim = c(0,250),
       xlab = "pseudotime (second)",
       ylab = "#Macrophage Alveolar per 0.5s")
  legend(x = 0,
         y = 250,
         fill = "darkgreen",
         legend = "Macrophage_Alveolar",
         bty = "n",
         cex = 1.5,
         adj = c(0.2,0.5))
  
  plot(XXX_result$pseudotimes_0.5_gap, XXX_result$ncMonocyte,
       type = "o",
       col = "purple",
       cex = 0.5,
       pch = 19,
       ylim = c(0,100),
       xlab = "pseudotime (second)",
       ylab = "#ncMonocyte per 0.5s")
  legend(x = 0,
         y = 100,
         fill = "purple",
         legend = "ncMonocyte",
         bty = "n",
         cex = 1.5,
         adj = c(0.3,0.5))
  
  plot(XXX_result$pseudotimes_0.5_gap, XXX_result$Macrophage,
       type = "o",
       col = "blue",
       cex = 0.5,
       pch = 19,
       ylim = c(0,600),
       xlab = "pseudotime (second)",
       ylab = "#Macrophage per 0.5s")
  legend(x = 0,
         y = 600,
         fill = "blue",
         legend = "Macrophage",
         bty = "n",
         cex = 1.5,
         adj = c(0.3,0.5))
}

#4 curves in one (0.5s)
{
  plot(XXX_result$pseudotimes_0.5_gap, XXX_result$cMonocyte,
       type = "o",
       col = "darkred",
       cex = 0.5,
       pch = 19,
       ylim = c(0,600),
       xlab = "pseudotime (second)",
       ylab = "#cells per 0.5s")
  lines(XXX_result$pseudotimes_0.5_gap, XXX_result$Macrophage_Alveolar,
        type = "o",
        col = "darkgreen",
        cex = 0.5,
        pch = 19,
        ylim = c(0,600),
        xlab = "pseudotime (second)")#,
  #ylab = "#Macrophage Alveolar per 0.5s")
  lines(XXX_result$pseudotimes_0.5_gap, XXX_result$ncMonocyte,
        type = "o",
        col = "purple",
        cex = 0.5,
        pch = 19,
        ylim = c(0,600),
        xlab = "pseudotime (second)")#,
  #ylab = "#ncMonocyte per 0.5s")
  lines(XXX_result$pseudotimes_0.5_gap, XXX_result$Macrophage,
        type = "o",
        col = "blue",
        cex = 0.5,
        pch = 19,
        ylim = c(0,600),
        xlab = "pseudotime (second)")#,
  #ylab = "#Macrophage per 0.5s")
  legend(x = 0,
         y = 600,
         fill = c( "blue", "darkgreen", "darkred", "purple"),
         legend = c("Macrophage","Macrophage_Alveolar","cMonocyte","ncMonocyte"),
         bty = "n",
         cex = 1.2,
         adj = c(0,0.5))
}

#calculating disease number in each 0.5 pseudo-second gap and plot
#data in used
ZZZZZ <- cbind(pData(pseudo_monocle)$Pseudotime,
               table(pData(pseudo_monocle)$Pseudotime, pData(pseudo_monocle)$disease_type))
rownames(ZZZZZ) <- NULL
colnames(ZZZZZ) <- c("pseudotimes", "Control", "COPD")
ZZZZZ <- data.frame(ZZZZZ)
ZZZZZ <- ZZZZZ[order(ZZZZZ$pseudotimes),]
head(ZZZZZ)

#0.5s gap plot with legend
{
  #data for plot
  ZZZZZ_result <- matrix(nrow = 78, ncol = 2) #39s / 0.5 == 78
  for (i in 0:77) {
    ZZZZZ_result[i+1,] <- colSums(ZZZZZ[which(ZZZZZ$pseudotimes >= 0.5*i & ZZZZZ$pseudotimes < 0.5*(i+1) ),
                                        c(2,3)]) 
    
  }
  ZZZZZ_result <- data.frame(ZZZZZ_result)
  rownames(ZZZZZ_result) <- seq(from = 1, to = 78, by = 1)
  colnames(ZZZZZ_result) <- c("Control", "COPD")
  ZZZZZ_result <- cbind(seq(from = 0.5, to = 39, by = 0.5), ZZZZZ_result)
  colnames(ZZZZZ_result) <- c("pseudotimes_0.5_gap", "Control", "COPD")
  head(ZZZZZ_result)
}
{
  #plot
  par(mfrow=c(2,1))
  plot(ZZZZZ_result$pseudotimes_0.5_gap, ZZZZZ_result$Control,
       type = "o",
       col = "darkgreen",
       cex = 0.5,
       pch = 19,
       ylim = c(0,650),
       xlab = "pseudotime (second)",
       ylab = "#Control per 0.5s")
  legend(x = 0,
         y = 650,
         fill = "darkgreen",
         legend = "Control",
         bty = "n",
         cex = 1.5,
         adj = c(0.2,0.5))
  
  plot(ZZZZZ_result$pseudotimes_0.5_gap, ZZZZZ_result$COPD,
       type = "o",
       col = "darkred",
       cex = 0.5,
       pch = 19,
       ylim = c(0,350),
       xlab = "pseudotime (second)",
       ylab = "#COPD per 0.5s")
  legend(x = 0,
         y = 350,
         fill = "darkred",
         legend = "COPD",
         bty = "n",
         cex = 1.5,
         adj = c(0.2,0.5))
}

#2 curves in one plot (0.5s)
{
  plot(ZZZZZ_result$pseudotimes_0.5_gap, ZZZZZ_result$Control,
       type = "o",
       col = "darkgreen",
       cex = 0.5,
       pch = 19,
       ylim = c(0,650),
       xlab = "pseudotime (second)",
       ylab = "#cells per 0.5s")
  lines(ZZZZZ_result$pseudotimes_0.5_gap, ZZZZZ_result$COPD,
        type = "o",
        col = "darkred",
        cex = 0.5,
        pch = 19,
        ylim = c(0,650),
        xlab = "pseudotime (second)")
  legend(x = 0,
         y = 650,
         fill = c("darkgreen", "darkred"),
         legend = c("Control","COPD"),
         bty = "n",
         cex = 1.2,
         adj = c(0,0.5))
}


#pseudotime intervals and cell total numbers in each interval
table(pData(pseudo_monocle)$Pseudotime, pData(pseudo_monocle)$State)
state_values <- as.numeric(pData(pseudo_monocle)$State)
#pseudotime_values
YYY <- data.frame(cbind(pseudotime_values, state_values))
YYY <- YYY[order(YYY$pseudotime_values),]
head(YYY)

nrow(YYY[which(YYY$state_values == 7),]) #1343
min(YYY$pseudotime_values[which(YYY$state_values == 7)]) #0
max(YYY$pseudotime_values[which(YYY$state_values == 7)]) #11.99622

nrow(YYY[which(YYY$state_values == 5),]) #661
min(YYY$pseudotime_values[which(YYY$state_values == 5)]) #12.32711
max(YYY$pseudotime_values[which(YYY$state_values == 5)]) #19.78476

nrow(YYY[which(YYY$state_values == 6),]) #568
min(YYY$pseudotime_values[which(YYY$state_values == 6)]) #12.55471
max(YYY$pseudotime_values[which(YYY$state_values == 6)]) #14.37097

nrow(YYY[which(YYY$state_values == 4),]) #592
min(YYY$pseudotime_values[which(YYY$state_values == 4)]) #20.23017
max(YYY$pseudotime_values[which(YYY$state_values == 4)]) #23.84594

nrow(YYY[which(YYY$state_values == 8),]) #1234
min(YYY$pseudotime_values[which(YYY$state_values == 8)]) #21.31814
max(YYY$pseudotime_values[which(YYY$state_values == 8)]) #24.65926

nrow(YYY[which(YYY$state_values == 2),]) #66
min(YYY$pseudotime_values[which(YYY$state_values == 2)]) #23.84866
max(YYY$pseudotime_values[which(YYY$state_values == 2)]) #25.07432

nrow(YYY[which(YYY$state_values == 3),]) #196
min(YYY$pseudotime_values[which(YYY$state_values == 3)]) #24.01473
max(YYY$pseudotime_values[which(YYY$state_values == 3)]) #25.22413

nrow(YYY[which(YYY$state_values == 1),]) #2653
min(YYY$pseudotime_values[which(YYY$state_values == 1)]) #25.09848
max(YYY$pseudotime_values[which(YYY$state_values == 1)]) #38.00622

nrow(YYY[which(YYY$state_values == 9),]) #616
min(YYY$pseudotime_values[which(YYY$state_values == 9)]) #25.26338
max(YYY$pseudotime_values[which(YYY$state_values == 9)]) #27.64269

