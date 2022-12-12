library(SeuratDisk)
library(Seurat)
library(reticulate)
library(tradeSeq)
library(SingleCellExperiment, quietly = T)
library(stringr)
library(BiocParallel)
library(ggplot2)

# Set conditions = T if using a dataset with control and conditon samples
get_metadata <- function(data_prefix, conditions = F) {
  seurat_obj_query <-
    LoadH5Seurat(paste0(
      "/home/tr455/project/subset/",
      data_prefix,
      "_seurat_temp.h5seurat"
    ))
  
  if (is.null(seurat_obj_query@meta.data$donors)) {
    seurat_obj_query@meta.data$donors <-
      seurat_obj_query@meta.data$samples
  }
  
  if (data_prefix == "Ma_Sestan" || data_prefix == "DevBrain") {
    seurat_obj_query@meta.data$donors <- str_split(seurat_obj_query@meta.data$donors, "_", simplify = TRUE)[, 1]
  }
  
  n <- length(seurat_obj_query@meta.data$donors)
  m <- 4
  lst <- rep(0, m)
  
  if (data_prefix == "IsoHuB") {
    metadata <-
      read.csv(
        paste0(
          "/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",
          data_prefix,
          "_metadata_reduced.csv"
        )
      )
  } else{
    metadata <-
      read.delim(
        paste0(
          "/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",
          data_prefix,
          "_reduced_metadata.tsv"
        )
      )
    if (data_prefix == "Kriegstein")
      metadata$individualID <- metadata$Sample.Library.No.
  }
  
  # Get a table of donor counts so we can calculate proportions
  donor_counts <- table(seurat_obj_query@meta.data$donors)
  
  donors <- unique(seurat_obj_query@meta.data$donors)
  
  # Get proportions of donors on a cell identity basis
  M <-
    t(sapply(donors, function(x)
      table(Idents(seurat_obj_query)[which(seurat_obj_query@meta.data$donors == x)])))
  M_proportions <- prop.table(M, 2)
  M_df <- data.frame(M_proportions)
  M_df$donor <- donors
  
  ids <-
    sapply(as.character(Idents(seurat_obj_query)), function(x)
      str_replace_all(x, " ", "."))
  ids <- sapply(ids, function(x)
    str_replace_all(x, "/", "."))
  tmp <- cbind(seurat_obj_query@meta.data$donors, ids)
  
  # Build U matrix that will be passed to the fitGAM function to correct
  # for covariates
  
  U <-
    t(sapply(seurat_obj_query@meta.data$donors, function(x)
      metadata[which(metadata$individualID == x)[1], c("Sex", "Age", "PMI")]))
  U[, "Sex"] <-
    lapply(U[, "Sex"], function(x)
      as.integer(tolower(x) == "male"))
  U[, "Age"] <- lapply(U[, "Age"], function(x)
    as.integer(x))
  U[, "PMI"] <- lapply(U[, "PMI"], function(x)
    as.double(x))
  U <-
    cbind(U,
          sapply(seurat_obj_query@meta.data$donors, function(x)
            donor_counts[x]))
  lst <- rep(0, length(ids))
  
  for (i in 1:length(ids)) {
    lst[i] <- M_df[tmp[i, ][[1]], tmp[i, ][[2]]]
  }
  
  U <- cbind(U, lst)
  colnames(U) <- c("Sex", "Age", "PMI", "dcounts", "props")
  
  # Convert to a dataframe and finall a model matrix
  U_df <- as.data.frame(U)
  
  sex <- as.factor(unlist(U_df$Sex))
  pmi <- unlist(U_df$PMI)
  pmi[is.na(pmi)] <- mean(pmi[!is.na(pmi)])
  age <- unlist(U_df$Age)
  prop <- unlist(U_df$props)
  dcount <- unlist(U_df$dcounts)
  options(na.action = "na.pass")
  U_final <- model.matrix( ~ age + pmi + sex + prop + dcount)
  counts <- seurat_obj_query@assays$SCT@counts
 
  # Save R data that is required for fitGAM function
  saveRDS(U_final,
          file = paste0("/home/tr455/project/subset/", data_prefix, "_U_it.rds"))
  
  if (conditions == T) {
    conds <-
      sapply(seurat_obj_query@meta.data$donors, function(x)
        metadata[which(metadata$individualID == x)[1], "Diagnosis"])
    saveRDS(conds,
            file = paste0(
              "/home/tr455/project/subset/",
              data_prefix,
              "_conds_it.rds"
            ))
  }
  saveRDS(counts,
          file = paste0(
            "/home/tr455/project/subset/",
            data_prefix,
            "_counts_it.rds"
          ))
  
  return(U)
}

temp <- get_metadata("Kriegstein")
temp <- get_metadata("IsoHuB")
temp <- get_metadata("Ma_Sestan")
temp <- get_metadata("DevBrain")
