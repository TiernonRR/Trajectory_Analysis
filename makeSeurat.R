library(SeuratDisk)
library(Seurat)
library(reticulate)
library(harmony)
library(stringr)
set.seed(42)


data_prefix = "CMC" # "IsoHuB" # replace this with IsoHub
coding <- read.delim("/home/tr455/project/Protein_coding_genes_biomart.txt")
#print("Before load...")
seurat_obj_query <- LoadH5Seurat(file = paste0("/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/", data_prefix, "_raw_combined_UMIs_AnnData_res6_filtered.h5seurat"), meta.data=FALSE, misc=FALSE)
#print("After load...")

scipy <- import("scipy")
np <- import("numpy")
raw_mat <- scipy$sparse$load_npz(paste0("/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",data_prefix,"_raw_counts.npz"))
cell_IDs <- scan(paste0("/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",data_prefix,"_cell_IDs.txt"),character())
gene_names <- scan(paste0("/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",data_prefix,"_gene_names.txt"),character())
dimnames(raw_mat) <- list(gene_names, cell_IDs)


print("before assay set.....")
seurat_obj_query <- SetAssayData(object = seurat_obj_query, slot = "counts", new.data = raw_mat)
print("before pred...")
#### Read in the previously predicted cell type labels and make them the primary 'identities' of the cells
predictions <- read.csv(paste0("/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",data_prefix,"_Azimuth_predictions.csv"),row.names=1,stringsAsFactors = F)

seurat_obj_query <- AddMetaData(seurat_obj_query, metadata = predictions, col.name='predictions')
print("after meta...")
Idents(seurat_obj_query) <- seurat_obj_query$predictions
print(unique(Idents(seurat_obj_query)))
#### Carry out normalization of the counts
seurat_obj_query <- SCTransform(seurat_obj_query, variable.features.n = NULL, variable.features.rv.th = 1.3,conserve.memory=T)
print(unique(Idents(seurat_obj_query)))
print(str(seurat_obj_query))
print("after sct...")

to_filter <- readLines(paste0("/home/tr455/project/filter/",data_prefix, "_tobefiltered_barcodes.txt"))
sub <- colnames(seurat_obj_query) %in% to_filter
seurat_obj_query <- seurat_obj_query[,which(sub==F)]

#### Isolate the excitatory neurons/layer-aware neurons and then run the PCA and UMAP
seurat_obj_query <- subset(x = seurat_obj_query, idents = c("L2/3 IT", "L6 IT", "L4 IT", "L5 IT"))  #, "L6 IT Car3", 
                                                           # "L5 ET", "L6 CT", "L5/6 NP", "L6b"))
seurat_obj_query <- subset(x = seurat_obj_query, subset = nCount_RNA > 1000 &
                                                          nFeature_RNA > 500)

 keep <- which(rownames(seurat_obj_query %in% coding$hgnc_symbol))
seurat_obj_query <- seurat_obj_query[keep, ]
print("after subset")
#### Run dimensionality reduction methods
# Harmony batch correction

## Below is for CMC only
sample_conv <- read.table("/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/CMC_sampleID_to_individualID.tsv")
temp <- sapply(seurat_obj_query$donors, function(t){
  row <- which(sample_conv == t)
  sample_conv[row, 2]
})
seurat_obj_query@meta.data$samples <- temp
## End for CMC only section

#to_split <- "-"
#seurat_obj_query@meta.data$samples <- str_split(rownames(seurat_obj_query@meta.data), to_split, simplify=T)[,1]
#seurat_obj_query <- RunPCA(seurat_obj_query)
seurat_obj_query <- RunHarmony(seurat_obj_query, "samples", plot.convergences=T, project.dim=F)

#Find Clusters Find Neighbors

#seurat_obj_query <- FindNeighbors(seurat_obj_query, dims = 1:30)
#seurat_obj_query <- FindClusters(seurat_obj_query, resolution = 0.5)

seurat_obj_query <- RunUMAP(seurat_obj_query, reduction="harmony", dims = 1:30)

png(file=paste0("/home/tr455/project/subset/", data_prefix,"_Test_UMAP_sub.png"))
DimPlot(seurat_obj_query, reduction = "umap", label = TRUE,pt.size = 0.5) + NoLegend()
dev.off()

# png(file=paste0(data_prefix,"_vlnplot_analysis.png"))
# VlnPlot(seurat_obj_query, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
# dev.off()

#print("after pca")
#seurat_obj_query <- RunUMAP(seurat_obj_query, reduction="harmony",dims = 1:30)


#png(file=paste0(data_prefix,"_Test_UMAP_.png"))
#DimPlot(seurat_obj_query, reduction = "umap", group.by = "samples",label = TRUE, pt.size = 0.5) + NoLegend()
#dev.off()

print("after umap")
seurat_obj_query <- subset(x = seurat_obj_query, features=VariableFeatures(object=seurat_obj_query))
print("after subset final")

SaveH5Seurat(seurat_obj_query, filename=paste0("/home/tr455/project/subset/", data_prefix, "_seurat_sub"), overwrite = T)