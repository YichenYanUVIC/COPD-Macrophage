#COPD & Control
#Macrophage & Monocyte

{
  library(Matrix)
  library(Seurat)
}
{
  matrix_dir = "/Users/richard.yan/Desktop/COPD latest/GSE 136831/"
  features.path <- paste0(matrix_dir, "GeneIDs.txt")
  barcode.path <- paste0(matrix_dir, "Barcodes.txt")
  matrix.path <- paste0(matrix_dir, "RawCounts_Sparse.mtx")
  celltype.path <- paste0(matrix_dir, "CellType.MetadataTable.txt")
}

feature.names = read.delim(features.path,
                           header = TRUE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
celltype.names = read.delim(celltype.path,
                            header = TRUE,
                            stringsAsFactors = FALSE)

mat <- readMM(file = matrix.path)

colnames(mat) = barcode.names$V1
rownames(mat) = feature.names$HGNC_EnsemblAlt_GeneID

Clustering_metadata <- celltype.names[c(which(celltype.names$Manuscript_Identity == "Macrophage"),
                                        which(celltype.names$Manuscript_Identity == "Macrophage_Alveolar"),
                                        which(celltype.names$Manuscript_Identity == "cMonocyte"),
                                        which(celltype.names$Manuscript_Identity == "ncMonocyte")), 
                                      c(1,5,7)]
Clustering_COPD_Control_barcode <- Clustering_metadata[which(Clustering_metadata$Disease_Identity != "IPF" ), c(1,3)]
levels(as.factor(Clustering_COPD_Control_barcode[,2]))

#Data in use
Clustering_Data <- mat[ , Clustering_COPD_Control_barcode[,1]]

setwd("/Users/richard.yan/Desktop/COPD latest")
#Disease Types: COPD + Control
saveRDS(Clustering_Data, "Clustering_Data_136831.Rds") #raw data
saveRDS(Clustering_metadata, "Clustering_Labels_136831.Rds") #raw disease and cell types labels



  