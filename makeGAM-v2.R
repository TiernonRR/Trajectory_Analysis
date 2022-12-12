library(tradeSeq)
library(SingleCellExperiment)
library(reticulate)
library(BiocParallel)
#library(batchtools)

#tmpl <- "/home/tr455/project/slurm.tmpl" # system.file(package="batchtools", "templates", "slurm-simple.tmpl")
#noquote(readLines((tmpl)))

#param <- BatchtoolsParam(workers=4, cluster="slurm", template=tmpl)
#register(param)

args = commandArgs(trailingOnly=TRUE)#print("registered param...")

set.seed(3)
data_prefix <- args[1] #"IsoHuB"
use_conds = args[2]
i = strtoi(args[3])
j = strtoi(args[4]) 
print(j)
print(data_prefix)
seurat_sling <- readRDS(paste0("/home/tr455/project/subset/", data_prefix,"_umap_slingResult_sub.rds"))
U_final <- readRDS(paste0("/home/tr455/project/subset/", data_prefix, "_U_it.rds"))
counts <- readRDS(paste0("/home/tr455/project/subset/", data_prefix, "_counts_it.rds"))
print(dim(U_final))
print(dim(counts))
if (use_conds == "TRUE"){

conds <- readRDS(paste0("/home/tr455/project/subset/", data_prefix, "_conds_it.rds"))

}

if (use_conds == "FALSE"){
  gam <- fitGAM(as.matrix(counts), 
                U=U_final, 
 #               parallel=T,
  #              BPPARAM=param,
                sds=seurat_sling,
		genes=i:j) 
} else {
  gam <- fitGAM(as.matrix(counts), 
              conditions = factor(conds),
              U=U_final, 
   #           parallel=T,
    #          BPPARAM=param,
              sds=seurat_sling,
	      genes=i:j)
}
saveRDS(gam, file=paste0("/home/tr455/project/subset/", data_prefix,"_gam/",data_prefix, args[3],"_gam_it.rds"))

