library(Seurat)
library(monocle)

#Data loading: GSE 136831 (COPD Control)
{
  Validation_raw <- readRDS("/Users/richard.yan/Desktop/COPD latest/Clustering_Data_136831.Rds") 
  Validation_Disease_Barcode_ALL <- readRDS("/Users/richard.yan/Desktop/COPD latest/Clustering_Labels_136831.Rds")
}

#Seurat Object Building
{
  Seurat_Validation <- CreateSeuratObject(Validation_raw, 
                                          min.cells = 5,
                                          min.features = 2000) #35603 genes 56318 cells
  cell_type <- as.factor(Validation_Disease_Barcode_ALL[which(Validation_Disease_Barcode_ALL[,1] %in% Seurat_Validation@assays$RNA@counts@Dimnames[[2]]),2])
  disease_type <- as.factor(Validation_Disease_Barcode_ALL[which(Validation_Disease_Barcode_ALL[,1] %in% Seurat_Validation@assays$RNA@counts@Dimnames[[2]]),3])
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
                                            nfeatures = 1000) #don't skip
}

# Identify the 20 most highly variable genes
top20 <- head(VariableFeatures(Seurat_Validation), 20)
# plot variable features with labels
plot1 <- VariableFeaturePlot(Seurat_Validation)
plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
plot2

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
#head(Seurat_Validation@reductions$umap@cell.embeddings) # 提取UMAP坐标值。
#note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters

#Cluster result
write.csv(Seurat_Validation@active.ident, 
          file = "/Users/richard.yan/Desktop/COPD latest/clusters_results_136831.csv", 
          quote = F)

#Visualization
#png(file="/Users/richard.yan/Desktop/COPD latest/UMAP_COPD&Control.png")#图片输出为PNG文件

DimPlot(Seurat_Validation, reduction = "umap", 
        shuffle = T, seed = 1,
        label = T, label.size = 10, repel = F,
        group.by = 'ident', pt.size = 0.5) + NoAxes() + NoLegend()

DimPlot(Seurat_Validation, reduction = "umap", 
        shuffle = T, seed = 1,
        label = T, label.size = 8, repel = F,
        group.by = 'ident', pt.size = 0.5,
        cols.highlight = c("orange2","blue","green","yellow"),
        cells.highlight = c(names(Seurat_Validation$seurat_clusters[which(Seurat_Validation$seurat_clusters == 0)]),
                            names(Seurat_Validation$seurat_clusters[which(Seurat_Validation$seurat_clusters == 6)]),
                            names(Seurat_Validation$seurat_clusters[which(Seurat_Validation$seurat_clusters == 8)]),
                            names(Seurat_Validation$seurat_clusters[which(Seurat_Validation$seurat_clusters == 15)]))) + NoAxes() + NoLegend()

#group by disease type
DimPlot(Seurat_Validation, reduction = "umap", shuffle = T,
        group.by = 'disease_type') + NoAxes()

DimPlot(Seurat_Validation, reduction = "umap", #shuffle = T,
        split.by = 'disease_type',
        label = T, label.size = 7) + NoAxes()
DimPlot(Seurat_Validation, reduction = "umap", #shuffle = T,
        split.by = 'disease_type', label = T, label.size = 6,
        cols.highlight = "orange2",
        cells.highlight = names(Seurat_Validation$seurat_clusters[which(Seurat_Validation$seurat_clusters == 15)])) + NoAxes() + NoLegend()

DimPlot(Seurat_Validation, reduction = "umap", #shuffle = T,
        group.by = 'disease_type', split.by = 'ident',
        ncol = 4)

#group by cell type
DimPlot(Seurat_Validation, reduction = "umap", shuffle = T,
        group.by = 'cell_type', pt.size = 0.7,
        label = T, label.size = 7) + NoAxes() + NoLegend()

DimPlot(Seurat_Validation, reduction = "umap", #shuffle = T,
        split.by = 'cell_type',
        label = T, label.size = 6) + NoAxes() + NoLegend()

DimPlot(Seurat_Validation, reduction = "umap", #shuffle = T,
        group.by = 'cell_type', split.by = 'ident',
        ncol = 4)
dev.off()

#find all markers of cluster 1
#cluster1.markers <- FindMarkers(Seurat_Validation, 
#ident.1 = 1, 
#min.pct = 0.25,
#test.use = 'MAST')
#head(cluster1.markers, n = 30)

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
            file = "/Users/richard.yan/Desktop/COPD latest/DE_MAST_ALL_136831.csv")
}

#Top marker genes
library(dplyr)
#load DE file saved
Validation.markers_ALL_136831 <- read.csv("/Users/richard.yan/Desktop/COPD latest/DE_MAST_ALL_136831.csv")
#Select top 5 marker genes of each cluster according to log2FC
top5 <- Validation.markers_ALL_136831 %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)

#DE Heatmap with top marker genes
library(DoMultiBarHeatmap)
#png(file="/Users/richard.yan/Desktop/DE_COPD&Control.png")#图片输出为PNG文件

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
dev.off()

#crosstable
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

#各个cluster中的疾病类型和细胞类型比例图
KKK01 <- table(DiseaseType, ClusterResults)
KKK02 <- table(CellType, ClusterResults)[c(1,4,2,3),]
#par(mfrow = c(1,2))
barplot(KKK01,
        col = colors()[c(12,555)], 
        font.axis = 2, 
        beside = TRUE,
        legend = rownames(KKK01), 
        args.legend = list(x = "topleft"),
        xlab = "Cluster", 
        ylab = "Cell Number",
        ylim = c(0,1500),
        font.lab = 2)
barplot(KKK02,
        col = colors()[c(14,400,99,465)], 
        font.axis = 2, 
        beside = TRUE,
        legend = rownames(KKK02), 
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
#setwd("/Users/richard.yan/Desktop/COPD latest")
#saveRDS(pseudo_monocle, "pseudo_monocle.Rds") #raw data

diff_test_res <- differentialGeneTest(pseudo_monocle, 
                                      fullModelFormulaStr = "~seurat_clusters")
#setwd("/Users/richard.yan/Desktop/COPD latest")
#saveRDS(diff_test_res, "diff_test_res.Rds") #raw data

#使用本地保存好的diff_test_res文件
diff_test_res <- readRDS("/Users/richard.yan/Desktop/COPD latest/diff_test_res.Rds")
ordering_genes <- row.names(subset(diff_test_res, 
                                   qval < 0.01))

pseudo_monocle <- setOrderingFilter(pseudo_monocle, ordering_genes)
plot_ordering_genes(pseudo_monocle)
pseudo_monocle <- reduceDimension(pseudo_monocle, 
                                  max_components = 2, 
                                  reduction_method = "DDRTree")
#setwd("/Users/richard.yan/Desktop/COPD latest")
#saveRDS(pseudo_monocle, "pseudo_monocle_reducedDimension.Rds") #raw data

#读入本地保存的已降维的文件
pseudo_monocle <- readRDS("/Users/richard.yan/Desktop/COPD latest/pseudo_monocle_reducedDimension.Rds")
pseudo_monocle <- orderCells(pseudo_monocle)

#GM_state函数 返回一个指定的状态作为pseudo time的起始状态
#此处规定：找到包含最多monocyte cell的State作为pseudotime的起始类别
#计算cluster 8中的所有细胞分别有多少落在各个State中(cluster 8中的ncMonocyte + cMonocyte最纯：99%)
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

#单一着色
plot_cell_trajectory(pseudo_monocle, color_by = "seurat_clusters") + NoAxes()
plot_cell_trajectory(pseudo_monocle, color_by = "disease_type") + NoAxes()
plot_cell_trajectory(pseudo_monocle, color_by = "cell_type") + NoAxes()
plot_cell_trajectory(pseudo_monocle, color_by = "State") + NoAxes()#只是单纯定义了segments

#状态对状态(解释说明State的定义)
plot_cell_trajectory(pseudo_monocle, color_by = "State") + NoAxes() + NoLegend() +
  facet_wrap(~State, nrow = 3, strip.position = "bottom") 

#聚类结果对疾病类型/细胞类型
plot_cell_trajectory(pseudo_monocle, color_by = "seurat_clusters") + NoAxes() +
  facet_wrap(~seurat_clusters, nrow = 4, strip.position = "bottom")
plot_cell_trajectory(pseudo_monocle, color_by = "disease_type") + NoAxes() +
  facet_wrap(~seurat_clusters, nrow = 4, strip.position = "bottom")
plot_cell_trajectory(pseudo_monocle, color_by = "cell_type") + NoAxes() +
  facet_wrap(~seurat_clusters, nrow = 4, strip.position = "bottom")

#细胞类型对疾病类型
plot_cell_trajectory(pseudo_monocle, color_by = "disease_type") + 
  facet_wrap(~cell_type, nrow = 1)
#疾病类型对细胞类型(important)
plot_cell_trajectory(pseudo_monocle, color_by = "cell_type") + 
  facet_wrap(~disease_type, nrow = 1)

#状态对疾病类型/细胞类型(important)
plot_cell_trajectory(pseudo_monocle, color_by = "disease_type") + NoAxes() +
  facet_wrap(~State, nrow = 3, strip.position = "bottom")
plot_cell_trajectory(pseudo_monocle, color_by = "cell_type") + NoAxes() +
  facet_wrap(~State, nrow = 3, strip.position = "bottom")

#表格: 状态对疾病类型/ Monocyte/ Macrophage 按照时间顺序排列States
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

#各类型细胞在各状态中的数量(状态按伪时间开始节点从左到右排序)
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

#两种疾病类型细胞在各状态中的数量(状态按伪时间开始节点从左到右排序)
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

#计算四种免疫细胞每0.1s / 0.5s 伪时间间隔内的数量 作图

#测试数据构建(not used) (don't delete!)
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

#以0.5s为间隔(有图例)
{
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

#同一张图显示4条曲线(0.5s)
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

#计算两种疾病状态每0.5s伪时间间隔内的数量 作图
#data in used
ZZZZZ <- cbind(pData(pseudo_monocle)$Pseudotime,
               table(pData(pseudo_monocle)$Pseudotime, pData(pseudo_monocle)$disease_type))
rownames(ZZZZZ) <- NULL
colnames(ZZZZZ) <- c("pseudotimes", "Control", "COPD")
ZZZZZ <- data.frame(ZZZZZ)
ZZZZZ <- ZZZZZ[order(ZZZZZ$pseudotimes),]
head(ZZZZZ)

#以0.5s为间隔(有图例)
{
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
  par(mfrow=c(2,1))
  plot(ZZZZZ_result$pseudotimes_0.5_gap, ZZZZZ_result$Control,
       type = "o",
       col = "darkgreen",
       cex = 0.5,
       pch = 19,
       ylim = c(0,650),
       xlab = "pseudotime (second)",
       ylab = "#Control per 0.5s")
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 7)]),
           y0 = 100,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 7)]),
           y1 = 100,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 5)]),
           y0 = 180,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 5)]),
           y1 = 180,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 6)]),
           y0 = 250,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 6)]),
           y1 = 250,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 4)]),
           y0 = 220,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 4)]),
           y1 = 220,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 8)]),
           y0 = 600,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 8)]),
           y1 = 600,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 2)]),
           y0 = 500,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 2)]),
           y1 = 500,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 3)]),
           y0 = 430,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 3)]),
           y1 = 430,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 1)]),
           y0 = 270,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 1)]),
           y1 = 270,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 9)]),
           y0 = 350,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 9)]),
           y1 = 350,
           lwd = 2)
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
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 7)]),
           y0 = 80,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 7)]),
           y1 = 80,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 5)]),
           y0 = 130,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 5)]),
           y1 = 130,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 6)]),
           y0 = 160,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 6)]),
           y1 = 160,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 4)]),
           y0 = 260,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 4)]),
           y1 = 260,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 8)]),
           y0 = 230,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 8)]),
           y1 = 230,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 2)]),
           y0 = 200,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 2)]),
           y1 = 200,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 3)]),
           y0 = 170,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 3)]),
           y1 = 170,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 1)]),
           y0 = 120,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 1)]),
           y1 = 120,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 9)]),
           y0 = 150,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 9)]),
           y1 = 150,
           lwd = 2)
  legend(x = 0,
         y = 350,
         fill = "darkred",
         legend = "COPD",
         bty = "n",
         cex = 1.5,
         adj = c(0.2,0.5))
}

#同一张图显示2条曲线(0.5s)
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
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 7)]),
           y0 = 100,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 7)]),
           y1 = 100,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 5)]),
           y0 = 180,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 5)]),
           y1 = 180,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 6)]),
           y0 = 250,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 6)]),
           y1 = 250,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 4)]),
           y0 = 220,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 4)]),
           y1 = 220,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 8)]),
           y0 = 600,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 8)]),
           y1 = 600,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 2)]),
           y0 = 500,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 2)]),
           y1 = 500,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 3)]),
           y0 = 430,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 3)]),
           y1 = 430,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 1)]),
           y0 = 270,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 1)]),
           y1 = 270,
           lwd = 2)
  segments(x0 = min(YYY$pseudotime_values[which(YYY$state_values == 9)]),
           y0 = 350,
           x1 = max(YYY$pseudotime_values[which(YYY$state_values == 9)]),
           y1 = 350,
           lwd = 2)
  legend(x = 0,
         y = 650,
         fill = c("darkgreen", "darkred"),
         legend = c("Control","COPD"),
         bty = "n",
         cex = 1.2,
         adj = c(0,0.5))
}

#计算两种疾病状态的总量变化 每0.5s区间 作图
#data in used
ZZZZZ <- cbind(pData(pseudo_monocle)$Pseudotime,
               table(pData(pseudo_monocle)$Pseudotime, pData(pseudo_monocle)$disease_type))
rownames(ZZZZZ) <- NULL
colnames(ZZZZZ) <- c("pseudotimes", "Control", "COPD")
ZZZZZ <- data.frame(ZZZZZ)
ZZZZZ <- ZZZZZ[order(ZZZZZ$pseudotimes),]
head(ZZZZZ)

#以0.5s为间隔(有图例)
{
  ZZZZZ_result_total <- matrix(nrow = 78, ncol = 2) #39s / 0.5 == 78
  for (i in 0:77) {
    ZZZZZ_result_total[i+1,] <- colSums(ZZZZZ[which(ZZZZZ$pseudotimes >= 0 & ZZZZZ$pseudotimes < 0.5*(i+1) ),
                                              c(2,3)]) 
    
  }
  ZZZZZ_result_total <- data.frame(ZZZZZ_result_total)
  rownames(ZZZZZ_result_total) <- seq(from = 1, to = 78, by = 1)
  colnames(ZZZZZ_result_total) <- c("Control", "COPD")
  ZZZZZ_result_total <- cbind(seq(from = 0.5, to = 39, by = 0.5), ZZZZZ_result_total)
  colnames(ZZZZZ_result_total) <- c("pseudotimes_0.5_gap", "Control", "COPD")
  head(ZZZZZ_result_total)
}
{
  par(mfrow=c(2,1))
  plot(ZZZZZ_result_total$pseudotimes_0.5_gap, ZZZZZ_result_total$Control,
       type = "o",
       col = "darkgreen",
       cex = 0.5,
       pch = 19,
       ylim = c(0,6000),
       xlab = "pseudotime (second)",
       ylab = "#Control in total per 0.5s")
  legend(x = 0,
         y = 6000,
         fill = "darkgreen",
         legend = "Control",
         bty = "n",
         cex = 1.5,
         adj = c(0.2,0.5))
  
  plot(ZZZZZ_result_total$pseudotimes_0.5_gap, ZZZZZ_result_total$COPD,
       type = "o",
       col = "darkred",
       cex = 0.5,
       pch = 19,
       ylim = c(0,3000),
       xlab = "pseudotime (second)",
       ylab = "#COPD in total per 0.5s")
  legend(x = 0,
         y = 3000,
         fill = "darkred",
         legend = "COPD",
         bty = "n",
         cex = 1.5,
         adj = c(0.2,0.5))
}

#同一张图显示2条曲线(0.5s)
{
  plot(ZZZZZ_result_total$pseudotimes_0.5_gap, ZZZZZ_result_total$Control,
       type = "o",
       col = "darkgreen",
       cex = 0.5,
       pch = 19,
       ylim = c(0,6000),
       xlab = "pseudotime (second)",
       ylab = "#cells in total per 0.5s")
  lines(ZZZZZ_result_total$pseudotimes_0.5_gap, ZZZZZ_result_total$COPD,
        type = "o",
        col = "darkred",
        cex = 0.5,
        pch = 19,
        ylim = c(0,6000),
        xlab = "pseudotime (second)")
  legend(x = 0,
         y = 6000,
         fill = c("darkgreen", "darkred"),
         legend = c("Control","COPD"),
         bty = "n",
         cex = 1.2,
         adj = c(0,0.5))
}

#计算四种免疫细胞每0.5s 伪时间间隔内的总量变化 作图

#data in used
XXX <- cbind(pData(pseudo_monocle)$Pseudotime,
             table(pData(pseudo_monocle)$Pseudotime, pData(pseudo_monocle)$cell_type))
rownames(XXX) <- NULL
colnames(XXX) <- c("pseudotimes", "cMonocyte", "Macrophage", "Macrophage_Alveolar", "ncMonocyte")
XXX <- data.frame(XXX)
XXX <- XXX[order(XXX$pseudotimes),]
head(XXX)

#以0.5s为间隔(有图例)
{
  XXX_result_total <- matrix(nrow = 78, ncol = 4) #39s / 0.5 == 78
  for (i in 0:77) {
    XXX_result_total[i+1,] <- colSums(XXX[which(XXX$pseudotimes >= 0 & XXX$pseudotimes < 0.5*(i+1) ),
                                          c(2,3,4,5)]) 
    
  }
  XXX_result_total <- data.frame(XXX_result_total)
  rownames(XXX_result_total) <- seq(from = 1, to = 78, by = 1)
  colnames(XXX_result_total) <- c("cMonocyte", "Macrophage", "Macrophage_Alveolar", "ncMonocyte")
  XXX_result_total <- cbind(seq(from = 0.5, to = 39, by = 0.5), XXX_result_total)
  colnames(XXX_result_total) <- c("pseudotimes_0.5_gap", "cMonocyte", "Macrophage", "Macrophage_Alveolar", "ncMonocyte")
  head(XXX_result_total)
}
{
  par(mfrow=c(2,2))
  plot(XXX_result_total$pseudotimes_0.5_gap, XXX_result_total$cMonocyte,
       type = "o",
       col = "darkred",
       cex = 0.5,
       pch = 19,
       ylim = c(0,1500),
       xlab = "pseudotime (second)",
       ylab = "#cMonocyte in total per 0.5s")
  legend(x = 0,
         y = 1500,
         fill = "darkred",
         legend = "cMonocyte",
         bty = "n",
         cex = 1.5,
         adj = c(0.35,0.5))
  
  plot(XXX_result_total$pseudotimes_0.5_gap, XXX_result_total$Macrophage_Alveolar,
       type = "o",
       col = "darkgreen",
       cex = 0.5,
       pch = 19,
       ylim = c(0,1500),
       xlab = "pseudotime (second)",
       ylab = "#Macrophage Alveolar in total per 0.5s")
  legend(x = 0,
         y = 1500,
         fill = "darkgreen",
         legend = "Macrophage_Alveolar",
         bty = "n",
         cex = 1.5,
         adj = c(0.2,0.5))
  
  plot(XXX_result_total$pseudotimes_0.5_gap, XXX_result_total$ncMonocyte,
       type = "o",
       col = "purple",
       cex = 0.5,
       pch = 19,
       ylim = c(0,1000),
       xlab = "pseudotime (second)",
       ylab = "#ncMonocyte in total per 0.5s")
  legend(x = 0,
         y = 1000,
         fill = "purple",
         legend = "ncMonocyte",
         bty = "n",
         cex = 1.5,
         adj = c(0.3,0.5))
  
  plot(XXX_result_total$pseudotimes_0.5_gap, XXX_result_total$Macrophage,
       type = "o",
       col = "blue",
       cex = 0.5,
       pch = 19,
       ylim = c(0,5000),
       xlab = "pseudotime (second)",
       ylab = "#Macrophage in total per 0.5s")
  legend(x = 0,
         y = 5000,
         fill = "blue",
         legend = "Macrophage",
         bty = "n",
         cex = 1.5,
         adj = c(0.3,0.5))
}

#同一张图显示4条曲线(0.5s)
{
  plot(XXX_result_total$pseudotimes_0.5_gap, XXX_result_total$cMonocyte,
       type = "o",
       col = "darkred",
       cex = 0.5,
       pch = 19,
       ylim = c(0,5000),
       xlab = "pseudotime (second)",
       ylab = "#cells in total per 0.5s")
  lines(XXX_result_total$pseudotimes_0.5_gap, XXX_result_total$Macrophage_Alveolar,
        type = "o",
        col = "darkgreen",
        cex = 0.5,
        pch = 19,
        ylim = c(0,5000),
        xlab = "pseudotime (second)")#,
  #ylab = "#Macrophage Alveolar per 0.5s")
  lines(XXX_result_total$pseudotimes_0.5_gap, XXX_result_total$ncMonocyte,
        type = "o",
        col = "purple",
        cex = 0.5,
        pch = 19,
        ylim = c(0,5000),
        xlab = "pseudotime (second)")#,
  #ylab = "#ncMonocyte per 0.5s")
  lines(XXX_result_total$pseudotimes_0.5_gap, XXX_result_total$Macrophage,
        type = "o",
        col = "blue",
        cex = 0.5,
        pch = 19,
        ylim = c(0,5000),
        xlab = "pseudotime (second)")#,
  #ylab = "#Macrophage per 0.5s")
  legend(x = 0,
         y = 5000,
         fill = c("blue", "darkgreen", "darkred", "purple"),
         legend = c("Macrophage","Macrophage_Alveolar","cMonocyte","ncMonocyte"),
         bty = "n",
         cex = 1.2,
         adj = c(0,0.5))
}






#各状态伪时间区间和细胞总数
table(pData(pseudo_monocle)$Pseudotime, pData(pseudo_monocle)$State)
state_values <- as.numeric(pData(pseudo_monocle)$State)
#pseudotime_values
YYY <- data.frame(cbind(pseudotime_values, state_values))
YYY <- YYY[order(YYY$pseudotime_values),]
head(YYY)

nrow(YYY[which(YYY$state_values == 7),]) #1343
min(YYY$pseudotime_values[which(YYY$state_values == 7)]) #0
max(YYY$pseudotime_values[which(YYY$state_values == 7)]) #11.99622
#检查是否有遗漏细胞
YYY[which(YYY$pseudotime_values > 11.95 & YYY$pseudotime_values < 12.35), ]


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






