library(SeuratDisk)
library(Seurat)
library(reticulate)

#### Use this file to construct your joblist.txt file for the dead simple queue


data_prefix <- "DevBrain"

# How many genes are in its
seurat_obj_query <-
  LoadH5Seurat(paste0(
    "/home/tr455/project/final/subset/",
    data_prefix,
    "_seurat_temp.h5seurat"
  ))

ngenes <- nrow(seurat_obj_query)

# How many genes do you want to consider at a each job
num_split <- 500

to_split <- 1:ngenes
chunks <- split(to_split, ceiling(seq_along(to_split)/num_split))

vect <- rep("", length(chunks))
conditions <- " TRUE "
index <- 1

#### Iterate over chunks to make each line of jobslist file
for(chunk in chunks) {
  c_min <- min(chunk)
  c_max <- max(chunk)
  vect[index] <- 
    paste0("module load miniconda; conda init bash; source activate myenv; Rscript fitGAM.R ", data_prefix, conditions,  c_min, " ", c_max)
  index <- index + 1
}

#### Write jobslist vector to your desired file

fileConn <- file(paste0("/home/tr455/project/final/scripts/", data_prefix, "_joblist_tmp2.txt"))
writeLines(vect, fileConn)
close(fileConn)
