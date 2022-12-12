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

data_prefix <- "DevBrain"

# Flag True when using a dataset with conditions
conds <- T 

seurat_sling <- readRDS(paste0("/home/tr455/project/final/subset/", data_prefix, "_umap_slingResult_temp.rds"))

gam <- readRDS(paste0("/home/tr455/project/final/subset/", data_prefix, "_gam_final.rds"))

to_remove <- list()

levels <- c("L2_3_IT", "L4_IT", "L5_IT", "L6_IT")

## Remove DE genes from one and only one cluster

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

## Function to generate heatmap
# lineage - lineage for heatmap
# title - name of lineage to save plot
# rows - whether to print rownames

heatMap <- function(lineage, title, rows=F){
  yhatSmooth <- predictSmooth(gam, gene=lineage, nPoints=50, tidy=F)
  
  heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))), cluster_cols=F, 
                         show_rownames=rows,show_colnames=F, 
                         filename = paste0("/home/tr455/project/final/subset/plots/",
                                           data_prefix, "_" ,title, ".png"))
 
}

print(seurat_sling@lineages)

### Association Test

rowData(gam)$assocRes <- associationTest(gam, l2fc = log2(2), lineages=T)

assocRes <- rowData(gam)$assocRes

## Code adapted from tradeSeq downstream analysis vignette

if(conds == F){
  # Number of lineages will vary by dataset
  l1 <- rownames(assocRes)[
    which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.05)
  ] 
  
  l2 <- rownames(assocRes)[
    which(p.adjust(assocRes$pvalue_2, "fdr") <= 0.05)
  ] 
  
} else {
  # If the dataset has conditions
  if(data_prefix == "Kriegstein"){
    l1_cond <- rownames(assocRes)[
      which(p.adjust(assocRes$pvalue_lineage1_conditionControl, "fdr") <= 0.05)
    ]
    
    l2_cond <- rownames(assocRes)[
      which(p.adjust(assocRes$pvalue_lineage2_conditionControl, "fdr") <= 0.05)
    ]
    
    # Example if you need to check other conditions
    # l1_cond <- rownames(assocRes)[
    #   which(p.adjust(assocRes$pvalue_lineage1_conditionASD, "fdr") <= 0.05)
    # ]
    
  } else {
    l1_cond <- rownames(assocRes)[
      which(p.adjust(assocRes$`pvalue_lineage1_conditionNot Applicable`, "fdr") <= 0.05)
    ]
    
    l2_cond <- rownames(assocRes)[
      which(p.adjust(assocRes$`pvalue_lineage2_conditionNot Applicable`, "fdr") <= 0.05)
    ]
  }
}

### function to generate heatmaps and save lineages


getLins <- function() {
  
  ## The folling code block finds problem genes and removes them from the GAM.
  ## Once this is run, we will have to rerun the association test if 
  ## any genes were found.
  
  if(conds == F) {
    to_check <- union(l1, l2)
  } else {
    to_check <- union(l1_cond, l2_cond)
  }
 
  yhatSmooth <- predictSmooth(gam, gene=to_check, nPoints=50, tidy=F)
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
  
  print(unique(infs))
  print(unique(NAs))
  print(unique(zeros))
  
  to_remove <- union(unique(infs), union(unique(NAs), unique(zeros)))
  to_remove <- unique(to_remove)
  sub <- rownames(gam) %in% to_remove
  gam <- gam[which(sub==F),]
  
  ## If we had to remove genes from the GAM we rerun the test here
  
  rowData(gam)$assocRes <- associationTest(gam, l2fc = log2(2), lineages=T)
  
  assocRes <- rowData(gam)$assocRes
  
  if(conds == F){
    l1 <- rownames(assocRes)[
      which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.05)
    ] 
    
    l2 <- rownames(assocRes)[
      which(p.adjust(assocRes$pvalue_2, "fdr") <= 0.05)
    ] 
  } else {
    # If the dataset has conditions
    if(data_prefix == "Kriegstein"){
      l1_cond <- rownames(assocRes)[
        which(p.adjust(assocRes$pvalue_lineage1_conditionControl, "fdr") <= 0.05)
      ]
      
      l2_cond <- rownames(assocRes)[
        which(p.adjust(assocRes$pvalue_lineage2_conditionControl, "fdr") <= 0.05)
      ]
    } else {
      l1_cond <- rownames(assocRes)[
        which(p.adjust(assocRes$`pvalue_lineage1_conditionNot Applicable`, "fdr") <= 0.05)
      ]
      
      l2_cond <- rownames(assocRes)[
        which(p.adjust(assocRes$`pvalue_lineage2_conditionNot Applicable`, "fdr") <= 0.05)
      ]
    }
  }
  
  if(conds == F){
    heatMap(l1, "l1")
    if(length(l2) > 0) heatMap(l2, "l2")
    
    saveRDS(l1, paste0("/home/tr455/project/final/subset/", data_prefix, "_l1"))
    if(length(l2) > 0) saveRDS(l2, paste0("/home/tr455/project/final/subset/", data_prefix, "_l2"))
    
  } else{
    heatMap(l1_cond, "l1_cond")
    if(length(l2_cond) > 0) heatMap(l2_cond, "l2_cond")
    
    saveRDS(l1_cond, paste0("/home/tr455/project/final/subset/", data_prefix, "_l1_cond"))
    if(length(l2_cond) > 0) saveRDS(l2_cond, paste0("/home/tr455/project/final/subset/", data_prefix, "_l2_cond"))
  }
}

# Might have to run this line several times if problem genes remain after 
# initial removals

getLins()