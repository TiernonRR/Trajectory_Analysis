library(SeuratDisk)
library(Seurat)
library(reticulate)
library(tradeSeq)
library(SingleCellExperiment, quietly = T)
library(stringr)
library(BiocParallel)
library(ggplot2)


set.seed(3)
data_prefix <- "CMC"

seurat_obj_query <- LoadH5Seurat(paste0("/home/tr455/project/subset/", data_prefix, "_seurat_sub.h5seurat"))
print(seurat_obj_query)
#seurat_obj_query_it <- subset(x = seurat_obj_query, idents = c("L2/3 IT", "L6 IT", "L4 IT", "L5 IT"))
#SaveH5Seurat(seurat_obj_query_it, filename=paste0("/home/tr455/project/", data_prefix, "_seurat_it"))
#print(colnames(seurat_obj_query))
# For DevBrain
#seurat_obj_query$donors <- str_split(colnames(seurat_obj_query), "_", simplify=T)[,1]
# Otherwise
splits <- str_split(colnames(seurat_obj_query), "-", simplify=T)
seurat_obj_query$donors <- paste(splits[,1], splits[,2], sep = "-")

n <- length(seurat_obj_query$donors)
m <- 4 
lst <- rep(0, m)
  
seurat_sling <- readRDS(paste0("/home/tr455/project/subset/", data_prefix,"_umap_slingResult_sub.rds"))
if(data_prefix == "IsoHuB"){
metadata <- read.csv(paste0("/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",data_prefix, "_metadata_reduced.csv"))#read.table("/home/tr455/project/Kriegstein_reduced_metadata.tsv")
#colnames(metadata) <- c("X", "V2", "V4", "V5", "V3")
# rename columns to match rest of code
} else{
metadata <- read.delim(paste0("/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/",data_prefix, "_reduced_individual_metadata.tsv"))
#metadata$individualID <- metadata$Sample.Library.No.
}

sample_conv <- read.table("/gpfs/gibbs/pi/gerstein/pse5/Kriegstein_data/CMC_sampleID_to_individualID.tsv")
temp <- sapply(seurat_obj_query$donors, function(t){
  row <- which(sample_conv == t)
  sample_conv[row, 2]
})
seurat_obj_query$donors <- temp

donor_counts <- table(seurat_obj_query$donors)

png(file=paste0("/home/tr455/project/plots/", data_prefix,"_hist.png"))
plot(donor_counts, type="h", las=2, cex.axis=0.60, ylab="Cell Counts")
dev.off()

donors <- unique(seurat_obj_query$donors)
M <- t(sapply(donors, function(x) table(Idents(seurat_obj_query)[which(seurat_obj_query$donors == x)])))
M_proportions <- prop.table(M, 2)
M_df <- data.frame(M_proportions)
M_df$donor <- donors

ids <- sapply(as.character(Idents(seurat_obj_query)), function(x) str_replace_all(x, " ", "."))
ids <- sapply(ids, function(x) str_replace_all(x, "/", "."))
tmp <- cbind(seurat_obj_query$donors, ids)
  
U <- t(sapply(seurat_obj_query$donors, function(x) metadata[which(metadata$individualID == x)[1],c("Gender", "Age", "PMI")]))
U[,"Gender"] <- sapply(U[,"Gender"], function(x) as.integer(tolower(x) == "male") )
U[,"Age"] <- lapply(U[,"Age"], function(x){
  if (x == "89+") return(89)
  else return(as.integer(x))
} )
U[,"PMI"] <- lapply(U[,"PMI"], function(x) as.double(x) )
U <- cbind(U, sapply(seurat_obj_query$donors, function(x) donor_counts[x]))
lst <- rep(0, length(ids))

for(i in 1:length(ids)){
  lst[i] <- M_df[tmp[i,][[1]], tmp[i,][[2]]]
}



# reminder for kriegstein - Sample.Library.No.
# else - individualID
U <- cbind(U, lst)
colnames(U) <- c("Sex", "Age", "PMI", "dcounts", "props")


U_df <- as.data.frame(U)

sex <- as.factor(unlist(U_df$Sex))
pmi <- unlist(U_df$PMI)
pmi[is.na(pmi)] <- mean(pmi[!is.na(pmi)])
age <- unlist(U_df$Age)
prop <- unlist(U_df$props)
dcount <- unlist(U_df$dcounts)
options(na.action="na.pass")
U_final <- model.matrix(~ age + pmi + sex + prop + dcount)
#U_final[is.na(U_final), "pmi"] <- mean(U_final[!is.na(U_final), "pmi"])
counts <- seurat_obj_query@assays$SCT@counts
conditions <- T
saveRDS(U_final,file=paste0("/home/tr455/project/subset/", data_prefix, "_U_it.rds"))
if (conditions == T){
    conds <- sapply(seurat_obj_query$donors, function(x) metadata[which(metadata$individualID == x)[1],"Diagnosis"])
    saveRDS(conds,file=paste0("/home/tr455/project/subset/", data_prefix, "_conds_it.rds"))
}
saveRDS(counts,file=paste0("/home/tr455/project/subset/", data_prefix, "_counts_it.rds"))

#print("before gam...")
#if (data_prefix=="IsoHuB"){
#  gam <- fitGAM(as.matrix(counts), 
#                U=U_final, 
#                parallel=T,
  #              BPPARAM=param,
#                sds=seurat_sling, 
#                nknots=9, 
#                genes=1:2000)
#} else {
#  conds <- sapply(seurat_obj_query$donors, function(x) metadata[which(metadata$V2 == x)[1],"V6"])
#  gam <- fitGAM(as.matrix(counts), 
#              conditions = factor(conds),
#              U=U_final, 
#              parallel=T,
 #             BPPARAM=param,
#              sds=seurat_sling, 
#              nknots=9, 
#              genes=1:2000)
#}
#saveRDS(gam, file=paste0("/home/tr455/project/",data_prefix, "_gam.rds"))
