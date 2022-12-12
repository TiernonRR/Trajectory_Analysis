library(SeuratDisk)
library(Seurat)
library(reticulate)
library(harmony)
library(stringr)
source("/home/tr455/project/scripts/getMetadata.R") # Link to getMetadata file
set.seed(42)


data_prefix = "DevBrain"

# Initial preparations for seurat object

seurat_obj_query <-
  LoadH5Seurat(
    file = paste0(
      "/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",
      data_prefix,
      "_raw_combined_UMIs_AnnData_res6_filtered.h5seurat"
    ),
    meta.data = FALSE,
    misc = FALSE
  )

scipy <- import("scipy")
np <- import("numpy")
raw_mat <-
  scipy$sparse$load_npz(
    paste0(
      "/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",
      data_prefix,
      "_raw_counts.npz"
    )
  )

cell_IDs <-
  scan(
    paste0(
      "/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",
      data_prefix,
      "_cell_IDs.txt"
    ),
    character()
  )

gene_names <-
  scan(
    paste0(
      "/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",
      data_prefix,
      "_gene_names.txt"
    ),
    character()
  )

dimnames(raw_mat) <- list(gene_names, cell_IDs)


seurat_obj_query <-
  SetAssayData(object = seurat_obj_query,
               slot = "counts",
               new.data = raw_mat)

#### Read in the previously predicted cell type labels and make them the primary 'identities' of the cells

predictions <-
  read.csv(
    paste0(
      "/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",
      data_prefix,
      "_Azimuth_predictions.csv"
    ),
    row.names = 1,
    stringsAsFactors = F
  )

seurat_obj_query <-
  AddMetaData(seurat_obj_query, metadata = predictions, col.name = 'predictions')

Idents(seurat_obj_query) <- seurat_obj_query$predictions


#### Carry out normalization of the counts

# We do not run this for Ma_Sestan
# Can remove condition once Barcodes for Ma_Sestan are obtained.

if (data_prefix != "Ma_Sestan") {
  to_filter <-
    readLines(
      paste0(
        "/home/tr455/project/filter/",
        data_prefix,
        "_tobefiltered_barcodes.txt"
      )
    )
  sub <- colnames(seurat_obj_query) %in% to_filter
  seurat_obj_query <- seurat_obj_query[, which(sub == F)]
}


#### Isolate the excitatory neurons/layer-aware neurons and then run the PCA and UMAP
seurat_obj_query <-
  subset(x = seurat_obj_query,
         idents = c("L2/3 IT", "L6 IT", "L4 IT", "L5 IT"))

seurat_obj_query <-
  subset(x = seurat_obj_query, subset = nCount_RNA > 1000 &
           nFeature_RNA > 500)

seurat_obj_query <-
  SCTransform(
    seurat_obj_query,
    variable.features.n = NULL,
    variable.features.rv.th = 1.3,
    conserve.memory = T
  )

seurat_obj_query <- subset(x = seurat_obj_query, features=VariableFeatures(object=seurat_obj_query))



# Age Filter
if (data_prefix == "Kriegstein" || data_prefix == "Ma_Sestan") {
  metadata <- get_metadata(data_prefix)
  metadata <- data.frame(metadata)
  if (data_prefix == "Kriegstein")
    lb <- 13
  if (data_prefix == "Ma_Sestan")
    lb <- 19
  to_keep <- which(metadata$Age > lb)
  seurat_obj_query <- seurat_obj_query[, to_keep]
}

#### Run dimensionality reduction methods

# Harmony batch correction

to_split <- "-"
seurat_obj_query@meta.data$samples <-
  str_split(rownames(seurat_obj_query@meta.data), to_split, simplify = T)[, 1]

seurat_obj_query <- RunPCA(seurat_obj_query)
seurat_obj_query <-
  RunHarmony(seurat_obj_query,
             "samples",
             plot.convergences = T,
             project.dim = F)

# UMAP

seurat_obj_query <-
  RunUMAP(seurat_obj_query, reduction = "harmony", dims = 1:30)

png(file = paste0(
  "/home/tr455/project/subset/",
  data_prefix,
  "_Test_UMAP_tmp.png"
))

DimPlot(
  seurat_obj_query,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.5
) + NoLegend()
dev.off()



SaveH5Seurat(
  seurat_obj_query,
  filename = paste0("/home/tr455/project/subset/", data_prefix, "_seurat_temp"),
  overwrite = T
)