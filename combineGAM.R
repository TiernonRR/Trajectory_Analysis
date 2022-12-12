library(SingleCellExperiment)
data_prefix <- "IsoHuB"
path <- paste0("/home/tr455/project/temp/", data_prefix, "_gam/", data_prefix)

# Missing 1501 and 501 - CMC
# missing 1 = krie

lst <- c("501", "1001", "1501", "2001", "2501", "3001", "3501", "4001", 
         "4501", "5001", "5501", "6001", "6501", "7001", "7039")

combined_GAM <- readRDS(paste0(path, "1_gam_it.rds"))
for (i in 1:length(lst)){
    l <- lst[i]
    print(l)
    
    gam <- readRDS(paste0(path, l, "_gam_it.rds"))
    print(dim(gam))
    print(dim(combined_GAM))
    combined_GAM <- rbind(combined_GAM, gam)
}

saveRDS(combined_GAM, paste0("/home/tr455/project/", data_prefix, "_gam_it.rds"))
temp <- readRDS(paste0("/home/tr455/project/", data_prefix, "_gam_it.rds"))

gam1 <- readRDS(paste0(path, "501", "_gam_it.rds"))

final <- rbind(temp, gam1)
