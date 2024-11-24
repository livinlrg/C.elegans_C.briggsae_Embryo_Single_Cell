library(monocle3)
library(Seurat)
library(VisCello)
library(MatrixExtra)

# library(ggplot2)
# library(tibble) # for add_column
# library(viridis)
# library(reshape2)
# library(ggrepel)
# library(ggExtra) # For ggmarginal
# library(dplyr)
# library(SeuratWrappers)

dir = "/kimdata/livinlrg/scAnalysis/WS290/"
cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240829.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

# joint object with global projection
seurat_object <- readRDS(paste0(dir, "Objects/seurat_objects/seurat_rpca_integration_joined_20240823.rds"))

# tissue subsets
seurat_tissue_list <- readRDS(paste0(dir, "Objects/seurat_objects/so_cell_class_list_20240827.rds"))

# time subsets
seurat_time_list <- readRDS(paste0(dir, "Objects/seurat_objects/so_time_bin_list_20241010.rds"))

## Generate 3D umaps
for(cell_class in names(seurat_tissue_list)) {
  print(cell_class)
  seurat_tissue_list[[cell_class]] <- RunUMAP(seurat_tissue_list[[cell_class]],
                                              reduction = "integrated.rpca", dims = 1:50, seed = 100,
                                              n.components = 3L, reduction.name = "umap.rpca.3d")
}

for(time_bin in names(seurat_time_list)) {
  print(time_bin)
  seurat_time_list[[time_bin]] <- RunUMAP(seurat_time_list[[time_bin]],
                                          reduction = "integrated.rpca", dims = 1:50, seed = 100,
                                          n.components = 3L, reduction.name = "umap.rpca.3d")
}

seurat_object <- RunUMAP(seurat_object, reduction = "integrated.rpca", dims = 1:300, seed = 100,
                         n.components = 3L, reduction.name = "umap.rpca.3d")

joint_list <- append(seurat_tissue_list, seurat_time_list)

# Make cello object
pData <- data.frame(colData(cds))
fData <- data.frame(rowData(cds))
umap_proj <- Embeddings(seurat_object, reduction = "umap.rpca")
matrix <- exprs(cds)

saveRDS(pData, file = paste0(dir, "Objects/cello_creation/pdata.rds"))
saveRDS(fData, file = paste0(dir, "Objects/cello_creation/fdata.rds"))
saveRDS(umap_proj, file = paste0(dir, "Objects/cello_creation/global_embedding.rds"))
saveRDS(matrix, file = paste0(dir, "Objects/cello_creation/matrix.rds"))
saveRDS(matrix, file = paste0(dir, "Objects/cello_creation/matrix.rds"))
matrix <- readRDS(file = paste0(dir, "Objects/cello_creation/matrix.rds"))

for(object in names(joint_list)) {
  print(object)
  umap_proj <- Embeddings(joint_list[[object]], reduction = "umap.rpca")
  saveRDS(umap_proj, file = paste0(dir, "Objects/cello_creation/", object, "_embeddings.rds"))
}

# Cello  object
fData$symbol <- rownames(matrix)
fData <- fData %>% select("symbol", everything())
pData <- pData %>% select("dataset3", everything())

eset <- new("ExpressionSet",
            assayData = assayDataNew("environment", exprs = matrix, 
                                     norm_exprs = log1p(matrix)),
            phenoData =  new("AnnotatedDataFrame", data = pData),
            featureData = new("AnnotatedDataFrame", data = fData))

clist <- list()
for(object in names(joint_list)) {
  print(object)
  umap_proj_2d <- Embeddings(joint_list[[object]], reduction = "umap.rpca")
  colnames(umap_proj_2d) <- c("UMAP_1", "UMAP_2")
  
  umap_proj_3d <- Embeddings(joint_list[[object]], reduction = "umap.rpca.3d")
  colnames(umap_proj_3d) <- c("UMAP_1", "UMAP_2", "UMAP_3")
  
  cur_idx = match(rownames(umap_proj_2d), colnames(eset))
  cello <- new("Cello", name = object, idx = cur_idx)

  cello@proj <- list(assign(paste0(object, " UMAP [2D]"), umap_proj_2d),
                        assign(paste0(object, " UMAP [3D]"), umap_proj_3d))
  names(cello@proj) <- c(paste0(object, " UMAP [2D]"), paste0(object, " UMAP [3D]"))

  clist[[object]] <- cello
}

umap_proj_2d <- Embeddings(seurat_object, reduction = "umap.rpca")
colnames(umap_proj_2d) <- c("UMAP_1", "UMAP_2")

umap_proj_3d <- Embeddings(seurat_object, reduction = "umap.rpca.3d")
colnames(umap_proj_3d) <- c("UMAP_1", "UMAP_2", "UMAP_3")

cur_idx = match(rownames(umap_proj_2d), colnames(eset))
cello <- new("Cello", name = "Global", idx = cur_idx)

cello@proj <- list(assign(paste0("Global", " UMAP [2D]"), umap_proj_2d),
                   assign(paste0("Global", " UMAP [3D]"), umap_proj_3d))
names(cello@proj) <- c(paste0("Global", " UMAP [2D]"), paste0("Global", " UMAP [3D]"))

clist[["Global"]] <- cello
clist <- clist[c("Global", names(seurat_time_list), names(seurat_tissue_list))]

names(clist) <- c("Global",
                  "0 to 150 embryo time", "0 to 200 embryo time", "0 to 250 embryo time", "0 to 300 embryo time",
                  names(seurat_tissue_list))

dir.create(paste0(dir, "Objects/cello_creation/cello_internal/"))
saveRDS(eset, paste0(dir, "Objects/cello_creation/cello_internal/eset.rds"))
saveRDS(clist, paste0(dir, "Objects/cello_creation/cello_internal/clist.rds"))

for(umap in names(clist)) {
  print(umap)
  
  clist[[paste0("C. elegans ", umap)]] <- new("Cello",
                                              name = eval(paste0("C. elegans ", umap)),
                                              idx = clist[[umap]]@idx[which(eset[,clist[[umap]]@idx]$species == "C.elegans")])
  
  temp_list <- list(data.frame(clist[[umap]]@proj[[1]])[which(eset[,clist[[umap]]@idx]$species == "C.elegans"),],
                    data.frame(clist[[umap]]@proj[[2]])[which(eset[,clist[[umap]]@idx]$species == "C.elegans"),])
  names(temp_list) <- c(paste0("C. elegans ", umap, " UMAP [2D]"),
                        paste0("C. elegans ", umap, " UMAP [3D]"))
  
  clist[[paste0("C. elegans ", umap)]]@proj <- temp_list
  clist[[paste0("C. briggsae ", umap)]] <- new("Cello",
                                               name = eval(paste0("C. briggsae", umap)),
                                               idx = clist[[umap]]@idx[which(eset[,clist[[umap]]@idx]$species == "C.briggsae")])
  
  temp_list <- list(data.frame(clist[[umap]]@proj[[1]])[which(eset[,clist[[umap]]@idx]$species == "C.briggsae"),],
                    data.frame(clist[[umap]]@proj[[2]])[which(eset[,clist[[umap]]@idx]$species == "C.briggsae"),])
  names(temp_list) <- c(paste0("C. briggsae ", umap, " UMAP [2D]"),
                        paste0("C. briggsae ", umap, " UMAP [3D]"))
  clist[[paste0("C. briggsae ", umap)]]@proj <- temp_list
}

saveRDS(clist, paste0(dir, "Objects/cello_creation/cello_internal/clist.rds"))

############################
# External VisCello object #
############################

# just clean the fdata and the pdata
# then update the names of the genes to:
# gene_short_name WBGene Cbr_short_name CbrWBGene

matrix <- readRDS(file = paste0(dir, "Objects/cello_creation/matrix.rds"))
pData_internal <- readRDS(file = paste0(dir, "Objects/cello_creation/pdata.rds"))
fData_internal <- readRDS(file = paste0(dir, "Objects/cello_creation/fdata.rds"))
clist_internal <- readRDS(paste0(dir, "Objects/cello_creation/cello_internal/clist.rds"))

fData_external <- fData_internal[,c("id", "gene_short_name", 
                                    "gene_long_name", "elegans_id", 
                                    "elegans_gene_short_name", "elegans_gene_long_name", 
                                    "briggsae_id", "briggsae_gene_short_name",
                                    "briggsae_gene_long_name", "gene.type", 
                                    "orthology_conf", "WormCategory.1",
                                    "WormCategory.2", "WormCategory.3")]

colnames(fData_external) <- c("id", "gene_short_name", 
                              "gene_long_name", "elegans_id", 
                              "elegans_gene_short_name", "elegans_gene_long_name", 
                              "briggsae_id", "briggsae_gene_short_name",
                              "briggsae_gene_long_name", "gene.type", 
                              "orthology_confidence", "WormCategory.1",
                              "WormCategory.2", "WormCategory.3")

pData_external <- pData_internal[,c("cell_type", "lineage_broad", "species", "n.umi", "barcode",
                                    "cell.type", "cell.subtype", "plot.cell.type", "lineage_packer", 
                                    "smoothed.embryo.time", "embryo.time.bin", "germline_pseudotime",
                                    "time.point", "dataset3", "genotype", "filter", "mt.frac", 
                                    "cel_lineage_specific", "cbr_lineage_specific", "lineage_specific",
                                    "cel_lineage_broad", "cbr_lineage_broad",
                                    "cel_lineage_complete", "cbr_lineage_complete", "lineage_complete")]

colnames(pData_external) <- c("cell_type", "lineage_broad", "species", "n.umi", "barcode",
                              "packer_cell.type", "packer_cell.subtype", "packer_plot.cell.type", "packer_lineage", 
                              "smoothed.embryo.time", "embryo.time.bin", "germline_pseudotime",
                              "time.point", "batch", "genotype", "filter", "mt.frac", 
                              "cel_lineage_specific", "cbr_lineage_specific", "lineage_specific",
                              "cel_lineage_broad", "cbr_lineage_broad",
                              "cel_lineage_complete", "cbr_lineage_complete", "lineage_complete")

colSums(apply(pData_external, 2, is.na))
cols_to_update <- c("packer_cell.type",
                    "packer_cell.subtype",
                    "packer_plot.cell.type",
                    "packer_lineage")
for(cur_col in cols_to_update) {
  pData_external[,cur_col] <- ifelse(is.na(pData_external[,cur_col]), "unannotated", pData_external[,cur_col])
}
pData_external[,"smoothed.embryo.time"] <- ifelse(is.na(pData_external[,"smoothed.embryo.time"]),
                                                  -1, pData_external[,"smoothed.embryo.time"])
pData_external[,"embryo.time.bin"] <- ifelse(is.na(pData_external[,"embryo.time.bin"]),
                                                  "unannotated", pData_external[,"embryo.time.bin"])
pData_external[,"germline_pseudotime"] <- ifelse(is.na(pData_external[,"germline_pseudotime"]),
                                                  -1, pData_external[,"germline_pseudotime"])
pData_external[,"time.point"] <- ifelse(is.na(pData_external[,"time.point"]),
                                        "unannotated", pData_external[,"time.point"])
pData_external[,"filter"] <- ifelse(is.na(pData_external[,"filter"]),
                                        "none", pData_external[,"filter"])

matrix_external <- matrix

# gene_short_name WBGene Cbr_short_name CbrWBGene
common_rownames <- paste0(fData_external[fData_external$gene.type == "common",]$gene_short_name, " ",
                          fData_external[fData_external$gene.type == "common",]$id, " ", 
                          fData_external[fData_external$gene.type == "common",]$briggsae_gene_short_name, " ",
                          fData_external[fData_external$gene.type == "common",]$briggsae_id)
cel_rownames <- paste0(fData_external[fData_external$gene.type == "ce.unique",]$gene_short_name, " ",
                       fData_external[fData_external$gene.type == "ce.unique",]$id)
cbr_rownames <- paste0(fData_external[fData_external$gene.type == "cb.unique",]$briggsae_gene_short_name, " ",
                       fData_external[fData_external$gene.type == "cb.unique",]$briggsae_id)

rownames(fData_external) <- c(common_rownames, cel_rownames, cbr_rownames)
rownames(matrix_external) <- c(common_rownames, cel_rownames, cbr_rownames)

fData_external$gene_short_name <- rownames(fData_external)

norm_matrix <- Matrix(log1p(matrix_external/pData_internal$Size_Factor), sparse = T)

eset <- new("ExpressionSet",
            assayData = assayDataNew("environment",
                                     exprs = matrix_external, 
                                     norm_exprs = norm_matrix),
            phenoData =  new("AnnotatedDataFrame", data = pData_external),
            featureData = new("AnnotatedDataFrame", data = fData_external))

saveRDS(eset, paste0(dir, "Objects/cello_creation/cell_external/eset.rds"))
saveRDS(clist_internal, paste0(dir, "Objects/cello_creation/cell_external/clist.rds"))







