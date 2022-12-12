library(reticulate)
library(slingshot)
library(tradeSeq)
library(SingleCellExperiment)
library(pheatmap)
library(UpSetR)
library(stringr)
library(cowplot)

data_prefix <- "DevBrain"
seurat_sling <- readRDS(paste0("/home/tr455/project/coding/", data_prefix, "_umap_slingResult_coding.rds"))

gam <- readRDS(paste0("/home/tr455/project/coding/", data_prefix, "_gam_it.rds"))

to_remove <- list()

levels <- c("L2_3_IT", "L4_IT", "L5_IT", "L6_IT")

for (clust in levels){
  if(data_prefix == "Kriegstein"){
    t <- read.table(file=
        paste0("/home/tr455/project/DEGs_trajectory_Analysis/Kriegstein_", 
               clust, "_DEGs.tsv"), sep="\t", header=T)
  } else {
    t <- read.table(file=paste0("/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",
                                data_prefix, "_", clust, "_DEGs.csv"), 
                                sep="\t", header=T)
  }
  
  
  indices <- row.names(t[which(t['p_val_adj'] < 0.05),])
  to_remove <- append(to_remove, indices)
}

to_remove <- unique(to_remove)
sub <- rownames(gam) %in% to_remove
gam <- gam[which(sub==F),]

heatMap <- function(lineage, title, rows=F){
  yhatSmooth <- predictSmooth(gam, gene=lineage, nPoints=50, tidy=F)
  
  heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols=F, 
                         show_rownames=rows,show_colnames=F, 
                         filename = paste0("/home/tr455/project/coding/plots/",
                                           data_prefix, "_" ,title, ".png"))
 
}

### Association Test

rowData(gam)$assocRes <- associationTest(gam, l2fc = log2(2), lineages=T)

assocRes <- rowData(gam)$assocRes


l1 <- rownames(assocRes)[
  which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.05)
] 

# If the dataset has conditions
l1_cond <- rownames(assocRes)[
  which(p.adjust(assocRes$`pvalue_lineage1_conditionNot Applicable`, "fdr") <= 0.05)
]

### Finding if issues arise

yhatSmooth <- predictSmooth(gam, gene=l1_cond, nPoints=50, tidy=F)
scale_result <- t(scale(t(yhatSmooth[, 1:50])))
infs <- rep("", length(rownames(scale_result)))
NAs <- rep("", length(rownames(scale_result)))
zeros <- rep("", length(rownames(scale_result)))
names <- rownames(scale_result)
for(i in 1:length(names)){
  gene <- names[i]
  
  if(all(yhatSmooth[gene, ] == 0)){
    zeros[i] <- gene
  }
}

for(i in 1:length(names)){
  gene <- names[i]
  
  if(any(is.na(scale_result[gene,]) )){
    NAs <- gene
  } else if(any(scale_result[gene,] == Inf)){
    infs[i] <- gene
  }
}

unique(infs)
unique(NAs)
unique(zeros)

# If to_remove is non-empty, will have to rerun code above after remove
# from the gam
to_remove <- union(unique(infs), union(unique(NAs), unique(zeros)))

heatMap(l1, "l1")
heatMap(l1_cond, "l1_cond")

saveRDS(l1, paste0("/home/tr455/project/coding/", data_prefix, "_l1"))
saveRDS(l1_cond, paste0("/home/tr455/project/coding/", data_prefix, "_l1_cond"))

### Common
l1_k <- readRDS("/home/tr455/project/coding/Kriegstein_l1_cond")
l1_d <- readRDS("/home/tr455/project/coding/DevBrain_l1_cond")
l1_i <- readRDS("/home/tr455/project/coding/IsoHuB_l1")
l1_m <- readRDS("/home/tr455/project/coding/Ma_Sestan_l1")

u_k <- unique(l1_k)
u_d <- unique(l1_d)
u_i <- unique(l1_i)
u_m <- unique(l1_m)

common_genes <- intersect(u_k, intersect(u_d, intersect(u_i, u_m)))

png(file=paste0("/home/tr455/project/coding/plots/comparison_it_upset.png"), res=100)
UpSetR::upset(fromList(list(krie=u_k, dev=u_d, iso=u_i, ma_s=u_m)))
dev.off()

gam_K <- readRDS(paste0("/home/tr455/project/coding/", "Kriegstein", "_gam_it.rds"))
gam_I <- readRDS(paste0("/home/tr455/project/coding/", "IsoHuB", "_gam_it.rds"))
gam_M <- readRDS(paste0("/home/tr455/project/coding/", "Ma_Sestan", "_gam_it.rds"))
gam_D <- readRDS(paste0("/home/tr455/project/coding/", "DevBrain", "_gam_it.rds"))

### First with Kriegstein gam
yhatSmooth <- predictSmooth(gam_K, gene=common_genes, nPoints=50, tidy=F)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols=F, 
                       show_rownames=T,show_colnames=F, filename = "/home/tr455/project/coding/plots/common_genes_it_Kriegstein.png")
### IsoHuB
yhatSmooth <- predictSmooth(gam_I, gene=common_genes, nPoints=50, tidy=F)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols=F, 
                       show_rownames=T,show_colnames=F, filename = "/home/tr455/project/coding/plots/common_genes_it_IsoHuB.png")
yhatSmooth <- predictSmooth(gam_M, gene=common_genes, nPoints=50, tidy=F)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols=F, 
                       show_rownames=T,show_colnames=F, filename = "/home/tr455/project/coding/plots/common_genes_it_Ma.png")
yhatSmooth <- predictSmooth(gam_D, gene=common_genes, nPoints=50, tidy=F)
heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols=F, 
                       show_rownames=T,show_colnames=F, filename = "/home/tr455/project/coding/plots/common_genes_it_Dev.png")

fileConn <- file("/home/tr455/project/coding/common_genes_it.txt")
writeLines(common_genes, fileConn)
close(fileConn)

gene_union <- c(unique(rownames(gam_D)),unique(rownames(gam_M)), 
  unique(rownames(gam_K)), unique(rownames(gam_I)))

fileConn <- file("/home/tr455/project/coding/union_genes_it.txt")
writeLines(gene_union, fileConn)
close(fileConn)
