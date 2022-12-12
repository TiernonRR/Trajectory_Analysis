library(SingleCellExperiment)
library(tradeSeq)

## Use this file to combine result of dead simple queue job

data_prefix <- "DevBrain"

## Path to directory containing result of dead simple queue jobs
path <- paste0("/home/tr455/project/final/subset_tmp/", data_prefix, "_gam/")

files <- dir(path = path)
combined_GAM <- readRDS(paste0(path, files[1]))

## Iterate over files in directory and combine
for(i in 2:length(files)) {
  gam <- readRDS(paste0(path, files[i]))
  combined_GAM <- rbind(combined_GAM, gam)
}

## Save GAM to desired directory
saveRDS(combined_GAM, paste0("/home/tr455/project/final/subset/", data_prefix, "_gam_final.rds"))

