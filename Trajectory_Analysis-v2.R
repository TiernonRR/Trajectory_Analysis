library(SeuratDisk)
library(Seurat)
library(reticulate)
library(readr)
set.seed(42)

data_prefix <- "CMC"
seurat_obj_query <- LoadH5Seurat(paste0("/home/tr455/project/subset/", data_prefix, "_seurat_sub.h5seurat"))

#common_genes <- read_lines("/home/tr455/project/common_genes.txt")
#print(seurat_obj_query)
#seurat_obj_query <- seurat_obj_query[common_genes,]
#print(seurat_obj_query)
#### Output UMAP to file
# png(file=paste0(data_prefix,"_Test_UMAP_orig.png"))
# DimPlot(seurat_obj_query, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# dev.off()

#### BEGIN PHATE
#library(phateR)
# Since the dataset is so high dimensional, we will first run PHATE on PCA
#pca <- Embeddings(object=seurat_obj_query, reduction='harmony')

umap <- Embeddings(object=seurat_obj_query, reduction='umap')
#phate_object_PCA <- phate(as.matrix(pca))

#png(file=paste0("/home/tr455/project/plots/", data_prefix,"_PHATE_coding.png"))
#plot(phate_object_PCA, col=Idents(seurat_obj_query))
#dev.off()

# Run slingshot
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("slingshot")

library(slingshot)
library(viridis)

seurat_sling <- slingshot(as.matrix(umap), clusterLabels = Idents(seurat_obj_query), start.clus='L2/3 IT')
saveRDS(seurat_sling, paste0("/home/tr455/project/subset/", data_prefix, "_umap_slingResult_sub.rds"))

#seurat_sling <- readRDS(paste0("/home/tr455/project/", data_prefix, "_slingResultv3.rds"))
plotcol <- viridis::turbo(length(levels(seurat_obj_query)))



png(file=paste0("/home/tr455/project/plots/", data_prefix,"_umap_sub.png"))
plot(umap[,1], umap[,2], col=plotcol[Idents(seurat_obj_query)])
legend("bottomright",legend=levels(seurat_obj_query), col=plotcol, pch=16, bty='n', cex=3/4)
lines(SlingshotDataSet(seurat_sling), lwd=2, col='black')
dev.off()

