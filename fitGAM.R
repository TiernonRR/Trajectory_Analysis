library(tradeSeq)
library(SingleCellExperiment)
library(reticulate)
library(BiocParallel)
set.seed(3)

# Read in commands fro command line 

args = commandArgs(trailingOnly = TRUE)


data_prefix <- args[1]
use_conds = args[2]
i = strtoi(args[3])
j = strtoi(args[4])

# Prints left in as they're helpful for debugging

print(j)
print(data_prefix)

seurat_sling <-
  readRDS(
    paste0(
      "/home/tr455/project/final/subset/",
      data_prefix,
      "_umap_slingResult_temp.rds"
    )
  )

U_final <-
  readRDS(paste0(
    "/home/tr455/project/final/subset/",
    data_prefix,
    "_U_it.rds"
  ))

counts <-
  readRDS(paste0(
    "/home/tr455/project/final/subset/",
    data_prefix,
    "_counts_it.rds"
  ))

print(dim(U_final))
print(dim(counts))

if (use_conds == "TRUE") {
  conds <-
    readRDS(paste0(
      "/home/tr455/project/final/subset/",
      data_prefix,
      "_conds_it.rds"
    ))
  
}

if (use_conds == "FALSE") {
  gam <- fitGAM(as.matrix(counts),
                U = U_final,
                sds = seurat_sling,
                genes = i:j)
} else {
  gam <- fitGAM(
    as.matrix(counts),
    conditions = factor(conds),
    U = U_final,
    sds = seurat_sling,
    genes = i:j
  )
}

# Write GAM to folder of respective dataset

saveRDS(
  gam,
  file = paste0(
    "/home/tr455/project/final/subset_tmp/",
    data_prefix,
    "_gam/",
    data_prefix,
    args[3],
    "_gam_it2.rds"
  )
)
