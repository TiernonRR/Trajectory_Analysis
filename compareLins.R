library(reticulate)
library(slingshot)
library(tradeSeq)
library(SeuratDisk)
library(Seurat)
library(reticulate)
library(SingleCellExperiment)
library(pheatmap)
library(UpSetR)
library(stringr)
library(cowplot)

#### This file may be used to compare the computed lineages from our pipelin


#### Read in the computed lineages

l1_k <- readRDS("/home/tr455/project/final/subset/Kriegstein_l1_cond")
l1_d <- readRDS("/home/tr455/project/final/subset/DevBrain_l1_cond")
l1_i <- readRDS("/home/tr455/project/final/subset/IsoHuB_l1")
l2_i <- readRDS("/home/tr455/project/final/subset/IsoHuB_l2")
l1_m <- readRDS("/home/tr455/project/final/subset/Ma_Sestan_l1")
l2_m <- readRDS("/home/tr455/project/final/subset/Ma_Sestan_l2")

u_k <- unique(l1_k)
u_d <- unique(l1_d)
u_i <- unique(union(l1_i, l2_i))
u_m <- unique(union(l1_m, l2_m))

#### Find the common genes
common_genes <- intersect(u_k, intersect(u_d, intersect(u_i, u_m)))

## Upset plot for lineages
png(file=paste0("/home/tr455/project/final/subset/plots/comparison_it_upset.png"), res=100)
UpSetR::upset(fromList(list(krie=u_k, dev=u_d, iso=u_i, ma_s=u_m)))
dev.off()

gam_K <- readRDS(paste0("/home/tr455/project/final/subset/", "Kriegstein", "_gam_final.rds"))
gam_I <- readRDS(paste0("/home/tr455/project/subset/", "IsoHuB", "_gam_final.rds"))
gam_M <- readRDS(paste0("/home/tr455/project/final/subset/", "Ma_Sestan", "_gam_final.rds"))
gam_D <- readRDS(paste0("/home/tr455/project/final/subset/", "DevBrain", "_gam_final.rds"))

### Generate heatmaps of the common genes for each dataset's GAM
yhatSmooth <- predictSmooth(gam_K, gene=common_genes, nPoints=50, tidy=F)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols=F,
                       show_rownames=F,show_colnames=F, filename = "/home/tr455/project/final/subset/plots/common_genes_it_Kriegstein.png")

yhatSmooth <- predictSmooth(gam_I, gene=common_genes, nPoints=50, tidy=F)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols=F,
                       show_rownames=F,show_colnames=F, filename = "/home/tr455/project/final/subset/plots/common_genes_it_IsoHuB.png")

yhatSmooth <- predictSmooth(gam_M, gene=common_genes, nPoints=50, tidy=F)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols=F,
                       show_rownames=F,show_colnames=F, filename = "/home/tr455/project/final/subset/plots/common_genes_it_Ma.png")

yhatSmooth <- predictSmooth(gam_D, gene=common_genes, nPoints=50, tidy=F)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols=F,
                       show_rownames=F,show_colnames=F, filename = "/home/tr455/project/final/subset/plots/common_genes_it_Dev.png")

fileConn <- file("/home/tr455/project/final/subset/common_genes_it.txt")
writeLines(common_genes, fileConn)
close(fileConn)

#### Find all of the genes across GAMs

gene_union <- c(unique(rownames(gam_D)),unique(rownames(gam_M)),
                unique(rownames(gam_K)), unique(rownames(gam_I)))

fileConn <- file("/home/tr455/project/final/subset/union_genes_it.txt")
writeLines(gene_union, fileConn)
close(fileConn)
