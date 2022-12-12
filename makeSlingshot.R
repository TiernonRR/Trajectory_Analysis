library(SeuratDisk)
library(Seurat)
library(reticulate)
library(readr)
library(slingshot)
library(viridis)
set.seed(42)

data_prefix <- "Ma_Sestan"

seurat_obj_query <-
  LoadH5Seurat(paste0(
    "/home/tr455/project/subset/",
    data_prefix,
    "_seurat_temp.h5seurat"
  ))


#### Run Slingshot on the UMAP embedding

umap <- Embeddings(object = seurat_obj_query, reduction = 'umap')
seurat_sling <- slingshot(
  as.matrix(umap),
  clusterLabels = Idents(seurat_obj_query),
  start.clus = 'L2/3 IT'
)
saveRDS(
  seurat_sling,
  paste0(
    "/home/tr455/project/subset/",
    data_prefix,
    "_umap_slingResult_temp.rds"
  )
)


plotcol <- viridis::turbo(length(levels(seurat_obj_query)))

umap <- Embeddings(object = seurat_obj_query, reduction = "umap")

png(file = paste0(
  "/home/tr455/project/subset/plots/",
  data_prefix,
  "_umap_subtmp2.png"
))
plot(umap[, 1], umap[, 2], col = plotcol[Idents(seurat_obj_query)])
legend(
  "bottomright",
  legend = levels(seurat_obj_query),
  col = plotcol,
  pch = 16,
  bty = 'n',
  cex = 3 / 4
)
lines(SlingshotDataSet(seurat_sling),
      lwd = 2,
      col = 'black')
dev.off()
