## Integrate the elegans and briggsae global dataset using Seurat RPCA
# This script will also generate a joint cds using the orthologous gene names
# Aug 22nd, 2024

# Load libraries
library(Seurat)
library(monocle3)
library(viridis)
library(future)
library(SeuratWrappers)

plan("multicore", workers = 8)
options(future.globals.maxSize = 45000 * 1024^2)

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

# load gene annoitation files
elegans_gtf <- data.frame(rtracklayer::import("/kimdata/livinlrg/scAnalysis/Reference/elegans/C_elegans.ws290.extendedUTRv2.format.gtf"))
briggsae_gtf <- data.frame(rtracklayer::import("/kimdata/livinlrg/scAnalysis/Reference/briggsae/C_briggsae.ws290.extendedUTRv3.format.gtf"))

# Load the data
# With background correction
cds_cel <- readRDS(paste0(dir, "Upload_20240819/cds.Celegans_WS260_WS290.bg_corr2.repeat.rds"))
cds_cbr <- readRDS(paste0(dir, "Upload_20240819/cds.briggsae_ws290_ws260_supp.cells.3.bg_corr.lt75.rds"))

# Without background correction
cds_cel <- readRDS(paste0(dir, "Upload_20240819/cds.Celegans_WS260_WS290.repeat.rds"))
cds_cbr <- readRDS(paste0(dir, "Upload_20240819/cds.briggsae_ws290_ws260_supp.cells.3.rds"))

# add embryo time estimation for the briggsae data
time_est_cds <- readRDS(paste0(dir, "Upload_20240822/cds.briggsae_ws290_ws260_supp.cells.3.rds"))

colData(cds_cbr)$embryo.time <- colData(time_est_cds)[colnames(cds_cbr),]$embryo.time.new
colData(cds_cbr)$smoothed.embryo.time <-colData(time_est_cds)[colnames(cds_cbr),]$smoothed.embryo.time.unified

rm(time_est_cds)

# generate a joint cds using the orthologous gene names
cel_rowdata <- rowData(cds_cel)
cbr_rowdata <- rowData(cds_cbr)

# record the original gene names in a new column
cbr_rowdata$briggsae_id <- cbr_rowdata$id
cbr_rowdata$briggsae_gene_short_name <- cbr_rowdata$gene_short_name
cbr_rowdata$briggsae_gene_long_name <- cbr_rowdata$gene_long_name

# take the elegans name if it exists, then the briggsae name if not
cbr_rowdata$id <- ifelse(is.na(cbr_rowdata$elegans_id), cbr_rowdata$briggsae_id, cbr_rowdata$elegans_id)
cbr_rowdata$gene_short_name <- ifelse(is.na(cbr_rowdata$elegans_gene_short_name), cbr_rowdata$briggsae_gene_short_name, cbr_rowdata$elegans_gene_short_name)
cbr_rowdata$gene_long_name <- ifelse(is.na(cbr_rowdata$elegans_gene_long_name), cbr_rowdata$briggsae_gene_long_name, cbr_rowdata$elegans_gene_long_name)

rownames(cbr_rowdata) <- cbr_rowdata$gene_short_name

# store the elegans data in a new column to match with the briggsae object
cel_rowdata$elegans_id <- cel_rowdata$id
cel_rowdata$elegans_gene_short_name <- cel_rowdata$gene_short_name
cel_rowdata$elegans_gene_long_name <- cel_rowdata$gene_long_name

rownames(cel_rowdata) <- cel_rowdata$gene_short_name

# need to fix the mtDNA genes not being labeled as common
cel_rowdata[cel_rowdata$id %in% unique(elegans_gtf[elegans_gtf$seqnames == "MtDNA",]$gene_id),]$gene.type <- "common"

# use the matching common gene names to update the briggsae names
cel_rowdata$briggsae_id <- NA
cel_rowdata$briggsae_gene_short_name <- NA
cel_rowdata$briggsae_gene_long_name <- NA

cel_rowdata[cel_rowdata$gene.type == "common",]$briggsae_id <- cbr_rowdata[cel_rowdata[cel_rowdata$gene.type == "common",]$gene_short_name,]$briggsae_id
cel_rowdata[cel_rowdata$gene.type == "common",]$briggsae_gene_short_name <- cbr_rowdata[cel_rowdata[cel_rowdata$gene.type == "common",]$gene_short_name,]$briggsae_gene_short_name
cel_rowdata[cel_rowdata$gene.type == "common",]$briggsae_gene_long_name <- cbr_rowdata[cel_rowdata[cel_rowdata$gene.type == "common",]$gene_short_name,]$briggsae_gene_long_name

# sort the rowdata to have common genes on top and in the same order between the two rowData's
cel_rowdata <- cel_rowdata[order(cel_rowdata$gene.type, decreasing = TRUE),]

# need to homogenize colnames
'''
colnames(cel_rowdata)
 [1] "id"                       "gene_short_name"         
 [3] "gene_long_name"           "WormCategory.1"          
 [5] "WormCategory.2"           "WormCategory.3"          
 [7] "ele_operon"               "TFclass"                 
 [9] "imaged_dataset"           "ws260.orthology_conf"    
[11] "use_for_ordering"         "orthology_conf"          
[13] "gene.type"                "elegans_id"              
[15] "elegans_gene_short_name"  "elegans_gene_long_name"  
[17] "briggsae_id"              "briggsae_gene_short_name"
[19] "briggsae_gene_long_name"
'''
cel_rowdata <- cel_rowdata[,c("id", "gene_short_name", "gene_long_name", "elegans_id", "elegans_gene_short_name", "elegans_gene_long_name", "briggsae_id", "briggsae_gene_short_name", "briggsae_gene_long_name", "gene.type", "orthology_conf", "ws260.orthology_conf", "WormCategory.1", "WormCategory.2", "WormCategory.3", "ele_operon", "TFclass", "imaged_dataset", "use_for_ordering")]
colnames(cel_rowdata)[! colnames(cel_rowdata) %in% colnames(cbr_rowdata)]
'''
[1] "ws260.orthology_conf" "WormCategory.1"       "WormCategory.2"      
[4] "WormCategory.3"       "ele_operon"           "TFclass"             
[7] "imaged_dataset"       "use_for_ordering" 
'''
# need to homogenize colnames
cbr_rowdata[,c("ws260.orthology_conf", "WormCategory.1", "WormCategory.2", "WormCategory.3", "ele_operon", "TFclass", "imaged_dataset", "use_for_ordering")] <- NA
cbr_rowdata <- cbr_rowdata[,c("id", "gene_short_name", "gene_long_name", "elegans_id", "elegans_gene_short_name", "elegans_gene_long_name", "briggsae_id", "briggsae_gene_short_name", "briggsae_gene_long_name", "gene.type", "orthology_conf", "ws260.orthology_conf", "WormCategory.1", "WormCategory.2", "WormCategory.3", "ele_operon", "TFclass", "imaged_dataset", "use_for_ordering")]

# sort the rowdata to have common genes on top and in the same order between the two rowData's
# then add a column to the briggsae one that adds the c. elegans name
cel_rowdata <- cel_rowdata[order(cel_rowdata$gene.type, decreasing = TRUE),]
cbr_rowdata <- cbr_rowdata[order(cbr_rowdata$gene.type, decreasing = TRUE),]

joint_rowdata <- rbind(cel_rowdata, cbr_rowdata[cbr_rowdata$gene.type == "cb.unique",])

# now need to homogenize the colData columns
cel_coldata <- colData(cds_cel)
cbr_coldata <- colData(cds_cbr)

colnames(cel_coldata)[! colnames(cel_coldata) %in% colnames(cbr_coldata)]
'''
[1] "dataset"  "dataset2"
'''
cbr_coldata$dataset <- NA
cbr_coldata$dataset2 <- NA

colnames(cbr_coldata)[! colnames(cbr_coldata) %in% colnames(cel_coldata)]
'''
 [1] "batch4"          "batch3"          "batch_300"       "batch_400"      
 [5] "batch_500"       "batch_120"       "batch_180"       "batch_240"      
 [9] "batch_360"       "cell_short"      "source"          "ws260_cell"     
[13] "global_cluster"  "supp.cell.type"  "supp.cell.type2" "supp.cell.type3"
'''
cel_coldata[,c("batch4", "batch3", "batch_300", "batch_400", "batch_500", "batch_120", "batch_180", "batch_240", "batch_360", "cell_short", "source", "ws260_cell", "global_cluster", "supp.cell.type", "supp.cell.type2", "supp.cell.type3")] <- NA

cel_coldata <- cel_coldata[,match(colnames(cbr_coldata), colnames(cel_coldata))]

joint_coldata <- rbind(cel_coldata, cbr_coldata)

# now need to rearrange the count matrices
cel_counts <- counts(cds_cel)
cbr_counts <- counts(cds_cbr)

# first construct cel portion
cel_counts <- cel_counts[joint_rowdata[joint_rowdata$gene.type %in% c("common", "ce.unique"),]$id,]

# create empty cbr_counts sparse matrix
cel_empty_matrix <- as(matrix(0, nrow = nrow(joint_rowdata[joint_rowdata$gene.type %in% c("cb.unique"),]), ncol = ncol(cel_counts)), "dgCMatrix")

rownames(cel_empty_matrix) <- joint_rowdata[joint_rowdata$gene.type %in% c("cb.unique"),]$id
colnames(cel_empty_matrix) <- colnames(cel_counts)

cel_joint_counts <- rbind(cel_counts, cel_empty_matrix)

# now construct cbr portion
cbr_empty_matrix <- as(matrix(0, nrow = nrow(joint_rowdata[joint_rowdata$gene.type %in% c("ce.unique"),]), ncol = ncol(cbr_counts)), "dgCMatrix")

rownames(cbr_empty_matrix) <- joint_rowdata[joint_rowdata$gene.type %in% c("ce.unique"),]$id
colnames(cbr_empty_matrix) <- colnames(cbr_counts)

# match briggsae_id in the joint_rowdata with cbr_counts rownames, then select the corresponding id and replace the rownames
joint_rowdata_cbr_subset <- joint_rowdata[joint_rowdata$gene.type %in% c("common", "cb.unique"),]
rownames(joint_rowdata_cbr_subset) <- joint_rowdata_cbr_subset$briggsae_id

rownames(cbr_counts) <- joint_rowdata_cbr_subset[rownames(cbr_counts),]$id

cbr_joint_counts <- rbind(cbr_empty_matrix, cbr_counts)

joint_counts <- cbind(cel_joint_counts[joint_rowdata$id,], cbr_joint_counts[joint_rowdata$id,])

rownames(joint_counts) <- joint_rowdata$gene_short_name

# create a monocle3 object
cds <- new_cell_data_set(joint_counts, cell_metadata = joint_coldata, gene_metadata = joint_rowdata)

colData(cds)$dataset3 <- ifelse(is.na(colData(cds)$dataset), colData(cds)$batch4, colData(cds)$dataset)

# pop over to CellQualityControl.R to add filter labels and adjust some cell annotations

# saveRDS(cds, paste0(dir, "Objects/cds_20240822.rds"))
cds <- readRDS(paste0(dir, "Objects/cds_20240822.rds"))

# originally forgot on 20240808 to subset to just the common genes
common_so <- CreateSeuratObject(counts = counts(cds)[rowData(cds)[rowData(cds)$gene.type %in% "common",]$gene_short_name,
                                                     ! colData(cds)$filter %in% c("doublet", "mutant", "damaged", "stressed")],
                                project = "Aug_2024",
                                min.features = 0,
                                min.cells = 0)

common_so <- AddMetaData(common_so, data.frame(colData(cds[rowData(cds)[rowData(cds)$gene.type %in% "common",]$gene_short_name,
                                                           ! colData(cds)$filter %in% c("doublet", "mutant", "damaged", "stressed")])))

common_so.list <- split(common_so, f = common_so$dataset3)

common_so.list <- NormalizeData(common_so.list)
common_so.list <- FindVariableFeatures(common_so.list)
common_so.list <- ScaleData(common_so.list)
common_so.list <- RunPCA(common_so.list, npcs = 300)

saveRDS(common_so.list, paste0(dir, "Objects/seurat_pre_integration_20240823.rds"))
common_so.list <- readRDS(paste0(dir, "Objects/seurat_pre_integration_20240823.rds"))

# mainstream integration
common_so.list <- IntegrateLayers(object = common_so.list,
                                  method = RPCAIntegration,
                                  orig.reduction = "pca",
                                  new.reduction = "integrated.rpca",
                                  dims = 1:300,
                                  verbose = TRUE)

print("Integration Complete")

common_so.list <- FindNeighbors(common_so.list, reduction = "integrated.rpca", dims = 1:300)
common_so.list <- FindClusters(common_so.list, resolution = 2, cluster.name = "rpca_clusters")

common_so.list <- RunUMAP(common_so.list, reduction = "integrated.rpca", dims = 1:300, reduction.name = "umap.rpca", seed.use = 20240826)

# Visualization
p1 <- DimPlot(common_so.list, reduction = "umap.rpca", group.by = "dataset3", raster = FALSE)
p2 <- DimPlot(common_so.list, reduction = "umap.rpca", group.by = "cell_type", label = TRUE, raster = FALSE, repel = TRUE, pt.size = 0.1) + NoLegend()
p3 <- DimPlot(common_so.list, reduction = "umap.rpca", group.by = "species", raster = FALSE)
p4 <- FeaturePlot(common_so.list,
                  reduction = "umap.rpca",
                  features = "smoothed.embryo.time",
                  raster = FALSE,
                  label = FALSE) +
  scale_color_viridis()
p5 <- DimPlot(common_so.list, reduction = "umap.rpca", group.by = "plot.cell.type", label = TRUE, repel = TRUE, pt.size = 0.1, raster = FALSE) + NoLegend()
p6 <- DimPlot(common_so.list, reduction = "umap.rpca", group.by = "rpca_clusters", label = TRUE, raster = FALSE, repel = TRUE, pt.size = 0.1) + NoLegend()

pdf(paste0(dir, "Global_RPCA_UMAP_20240823.pdf"), width = 60, height = 40)
print(p1 + p2 + p3 + p4 + p5 + p6)
dev.off()

pdf(paste0(dir, "Global_RPCA_UMAP_cell_labels_20240823.pdf"), width = 20, height = 20)
print(p2)
dev.off()

pdf(paste0(dir, "Global_RPCA_UMAP_clusters_20240823.pdf"), width = 20, height = 20)
print(p6)
dev.off()

# assign clusters to cell classes
cell_by_cluster <- table(common_so.list[["seurat_clusters"]][,1], common_so.list[["cell_type"]][,1])
total_cell_counts <- colSums(cell_by_cluster)

cell_by_cluster_percentages <- apply(cell_by_cluster, 1, function(x) (x / total_cell_counts) * 100)
cell_by_cluster_percentages_df <- data.frame(cell_by_cluster_percentages)
colnames(cell_by_cluster_percentages_df) <- colnames(cell_by_cluster_percentages)

fwd_append_df <- data.frame(cell_type = rownames(cell_by_cluster_percentages_df))
rownames(fwd_append_df) <- fwd_append_df$cell_type

fwd_append_df$cell_class <- NA
for(cell_class in names(ListOfCellTypes)) {
  fwd_append_df[fwd_append_df$cell_type %in% ListOfCellTypes[[cell_class]],]$cell_class <- cell_class
}
cell_by_cluster_percentages_df <- cbind(fwd_append_df, cell_by_cluster_percentages_df)


write.csv(cell_by_cluster_percentages_df,
          paste0(dir, "CellQC/global_cluster_cell_type_percentages.csv"), row.names = FALSE)
write.csv(cell_by_cluster_percentages_df %>% group_by(cell_class) %>% select(-cell_type) %>% summarise_all(mean),
          paste0(dir, "CellQC/global_cluster_cell_class_percentages.csv"), row.names = FALSE)

Idents(common_so.list) <- common_so.list[["seurat_clusters"]][,1]

new.cluster.ids <- c("Muscle", "Muscle", "Hypodermis and seam", "Muscle", "Muscle", "Muscle", "Hypodermis and seam", "Muscle", "Ciliated neurons", "Hypodermis and seam", "Early", "Muscle", "Muscle", "Pharynx and rectal", "Intestine", "Pharynx and rectal", "Non-ciliated neurons", "Early", "Hypodermis and seam", "Muscle", "Early", "Hypodermis and seam", "Intestine", "Non-ciliated neurons", "Early", "Hypodermis and seam", "Hypodermis and seam", "Muscle", "Hypodermis and seam", "Pharynx and rectal", "Glia and excretory", "Non-ciliated neurons", "Non-ciliated neurons", "Early", "Hypodermis and seam", "Non-ciliated neurons", "Germline", "Hypodermis and seam", "Mesoderm", "Early", "Muscle", "Muscle", "Muscle", "Mesoderm", "Early", "Hypodermis and seam", "Ciliated neurons", "Non-ciliated neurons", "Hypodermis and seam", "Non-ciliated neurons", "Early", "Non-ciliated neurons", "Pharynx and rectal", "Non-ciliated neurons", "Early", "Muscle", "Pharynx and rectal", "Mesoderm", "Ciliated neurons", "Mesoderm", "Ciliated neurons", "Ciliated neurons", "Non-ciliated neurons", "Non-ciliated neurons", "Glia and excretory", "Early", "Pharynx and rectal", "Ciliated neurons", "Muscle", "Mesoderm", "Early", "Pharynx and rectal", "Early", "Intestine", "Non-ciliated neurons", "Early", "Non-ciliated neurons", "Early", "Non-ciliated neurons", "Ciliated neurons", "Glia and excretory", "Early", "Glia and excretory", "Early", "Non-ciliated neurons", "Ciliated neurons", "Pharynx and rectal", "Non-ciliated neurons", "Pharynx and rectal", "Mesoderm", "Muscle", "Hypodermis and seam", "Non-ciliated neurons", "Non-ciliated neurons", "Non-ciliated neurons", "Early", "Early", "Glia and excretory", "Early", "Non-ciliated neurons", "Non-ciliated neurons", "Non-ciliated neurons", "Pharynx and rectal", "Non-ciliated neurons", "Non-ciliated neurons", "Ciliated neurons", "Ciliated neurons", "Early", "Non-ciliated neurons", "Non-ciliated neurons", "Early", "Non-ciliated neurons", "Pharynx and rectal", "Glia and excretory", "Ciliated neurons", "Glia and excretory", "Early", "Muscle", "Pharynx and rectal", "Ciliated neurons", "Pharynx and rectal", "Ciliated neurons", "Non-ciliated neurons", "Early", "Ciliated neurons", "Ciliated neurons", "Non-ciliated neurons", "Non-ciliated neurons", "Hypodermis and seam", "Ciliated neurons", "Pharynx and rectal", "Ciliated neurons", "Pharynx and rectal", "Non-ciliated neurons", "Non-ciliated neurons", "Non-ciliated neurons", "Non-ciliated neurons", "Muscle", "Early", "Pharynx and rectal", "Early", "Early", "Early", "Pharynx and rectal", "Early", "Glia and excretory", "Non-ciliated neurons", "Non-ciliated neurons", "Non-ciliated neurons", "Early", "Non-ciliated neurons", "Glia and excretory", "Early", "Mesoderm", "Early", "Early")
names(new.cluster.ids) <- levels(common_so.list)
common_so.list <- RenameIdents(common_so.list, new.cluster.ids)

pdf(paste0(dir, "Global_RPCA_UMAP_tissues_20240823.pdf"), width = 20, height = 20)
print(DimPlot(common_so.list, reduction = "umap.rpca", label = TRUE, raster = FALSE, repel = TRUE, pt.size = 0.1) + NoLegend())
dev.off()

saveRDS(common_so.list, paste0(dir, "Objects/seurat_rpca_integration_20240823.rds"))

common_so.rpca <- JoinLayers(common_so.list)

saveRDS(common_so.rpca, paste0(dir, "Objects/seurat_rpca_integration_joined_20240823.rds"))

combine_datasets_list <- list()
combine_datasets_list[["Germline"]] <- data.frame(to = c("Murray_b02_Murray_r17", "Murray_b02_Murray_r17",
                                                         "Waterston_400_minutes_Waterston_500_1_minutes_Waterston_500_2_minutes", "Waterston_400_minutes_Waterston_500_1_minutes_Waterston_500_2_minutes", "Waterston_400_minutes_Waterston_500_1_minutes_Waterston_500_2_minutes",
                                                         "Ce_wt_600_minutes_Ce_mec3_600_minutes", "Ce_wt_600_minutes_Ce_mec3_600_minutes",
                                                         "batch_400_batch_500", "batch_400_batch_500"),
                                                  from = c("Murray_b02", "Murray_r17",
                                                           "Waterston_400_minutes", "Waterston_500_1_minutes", "Waterston_500_2_minutes",
                                                           "Ce_wt_600_minutes", "Ce_mec3_600_minutes",
                                                           "batch_500", "batch_400"),
                                                  dimensions = 30,
                                                  k.weight = 100)
combine_datasets_list[["Glia and excretory"]] <- data.frame(to = c(),
                                                  from = c(),
                                                  dimensions = 50,
                                                  k.weight = 80)

combine_datasets_list[["Intestine"]] <- data.frame(to = c("Waterston_500_1_minutes_Waterston_500_2_minutes", "Waterston_500_1_minutes_Waterston_500_2_minutes",
                                                          "Murray_b01_Murray_b02", "Murray_b01_Murray_b02"),
                                                            from = c("Waterston_500_1_minutes", "Waterston_500_2_minutes",
                                                                     "Murray_b01", "Murray_b02"),
                                                            dimensions = 50,
                                                            k.weight = 70)

combine_datasets_list[["Mesoderm"]] <- data.frame(to = c(NA),
                                                  from = c(NA),
                                                  dimensions = 50,
                                                  k.weight = 85)

so_cell_class_list <- list()
for(cell_class in c("Muscle", "Hypodermis and seam", "Ciliated neurons", "Non-ciliated neurons", "Pharynx and rectal", "Glia and excretory", "Intestine", "Early", "Mesoderm", "Germline")) {
  print(cell_class)
  
  temp_so <- subset(common_so.rpca, cells = WhichCells(common_so.rpca, idents = cell_class))
  temp_so[["integrate_by"]] <- temp_so[["dataset3"]]
  
  # filter for dataset3 which have number of cells > 50
  # or combine like datasets with low batch effect
  if(cell_class %in% names(combine_datasets_list)) {
    for(i in 1:nrow(combine_datasets_list[[cell_class]])) {
      temp_so[["integrate_by"]][temp_so[["integrate_by"]][, 1] %in% combine_datasets_list[[cell_class]][i,"from"], "integrate_by"] <- combine_datasets_list[[cell_class]][i,"to"]
    }
  }
  
  # then integrate between species in these subsets
  so_cell_class_list[[cell_class]]  <- split(temp_so, f = temp_so$integrate_by)
  
  so_cell_class_list[[cell_class]] <- NormalizeData(so_cell_class_list[[cell_class]])
  so_cell_class_list[[cell_class]] <- FindVariableFeatures(so_cell_class_list[[cell_class]])
  so_cell_class_list[[cell_class]] <- ScaleData(so_cell_class_list[[cell_class]])
  so_cell_class_list[[cell_class]] <- RunPCA(so_cell_class_list[[cell_class]], npcs = 50)
}

saveRDS(so_cell_class_list, paste0(dir, "Objects/so_cell_class_list_20240827.rds"))

BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))

colData(cds)$cell_type_bin <- ""
for(cur_species in c("C.elegans", "C.briggsae")) {
  for(cell_type in names(BarcodeBinsList[[cur_species]])) {
    cell_names <- names(BarcodeBinsList[[cur_species]][[cell_type]])
    for(i in 1:length(BarcodeBinsList[[cur_species]][[cell_type]])) {
      colData(cds)[BarcodeBinsList[[cur_species]][[cell_type]][[i]],]$cell_type_bin <- paste0(colData(cds)[BarcodeBinsList[[cur_species]][[cell_type]][[i]],]$cell_type_bin, "_", cell_names[i])
    }
  }
}

colData(cds)$cell_type_bin <- gsub("^_", "", colData(cds)$cell_type_bin)
colData(cds)$cell_type_bin[colData(cds)$cell_type_bin == ""] <- "unassigned"

for(cell_class in c("Muscle", "Hypodermis and seam", "Ciliated neurons", "Non-ciliated neurons", "Pharynx and rectal", "Glia and excretory", "Intestine", "Early", "Mesoderm", "Germline")) {
  print(cell_class)
  if(cell_class %in% names(combine_datasets_list)) {
    n_dims = unique(combine_datasets_list[[cell_class]][1,"dimensions"])
    k.weight = unique(combine_datasets_list[[cell_class]][1,"k.weight"])
  } else {
    n_dims = 50
    k.weight = 100
  }
  so_cell_class_list[[cell_class]] <- IntegrateLayers(object = so_cell_class_list[[cell_class]],
                                                      method = RPCAIntegration,
                                                      orig.reduction = "pca",
                                                      new.reduction = "integrated.rpca",
                                                      dims = 1:n_dims,
                                                      k.weight = k.weight,
                                                      verbose = TRUE)

  so_cell_class_list[[cell_class]] <- FindNeighbors(so_cell_class_list[[cell_class]], reduction = "integrated.rpca", dims = 1:n_dims)
  so_cell_class_list[[cell_class]] <- FindClusters(so_cell_class_list[[cell_class]], resolution = 2, cluster.name = "rpca_clusters")

  so_cell_class_list[[cell_class]] <- RunUMAP(so_cell_class_list[[cell_class]],
                                              reduction = "integrated.rpca", dims = 1:n_dims,
                                              reduction.name = "umap.rpca", seed = 100)
  
  so_cell_class_list[[cell_class]][["plot_cell_type"]] <- so_cell_class_list[[cell_class]][["cell_type"]]
  
  subset_cell_types <- unique(so_cell_class_list[[cell_class]][["plot_cell_type"]])[,1]
  subset_cell_types <- unlist(sapply(subset_cell_types, function(cur_cell_type) {
    temp_length <- sum(so_cell_class_list[[cell_class]][["plot_cell_type"]][,1] %in% cur_cell_type)
    if(temp_length > 10) {
      return(cur_cell_type)
    }
  }))
  
  so_cell_class_list[[cell_class]][["plot_cell_type"]][! so_cell_class_list[[cell_class]][["plot_cell_type"]][,1] %in% subset_cell_types, "plot_cell_type"] <- "other"
  
  so_cell_class_list[[cell_class]] <- AddMetaData(so_cell_class_list[[cell_class]], data.frame(colData(cds)[colnames(so_cell_class_list[[cell_class]]),"cell_type_bin", drop = FALSE]))
  
  p1 <- DimPlot(so_cell_class_list[[cell_class]], reduction = "umap.rpca", group.by = "dataset3", raster = FALSE)
  p2 <- DimPlot(so_cell_class_list[[cell_class]], reduction = "umap.rpca", group.by = "plot_cell_type", label = TRUE, raster = FALSE, repel = TRUE, pt.size = 0.1) + NoLegend()
  p3 <- DimPlot(so_cell_class_list[[cell_class]], reduction = "umap.rpca", group.by = "species", raster = FALSE)
  p4 <- FeaturePlot(so_cell_class_list[[cell_class]],
                    reduction = "umap.rpca",
                    features = "smoothed.embryo.time",
                    raster = FALSE,
                    label = FALSE) +
    scale_color_viridis()
  p5 <- DimPlot(so_cell_class_list[[cell_class]], reduction = "umap.rpca", group.by = "cell_type_bin", label = TRUE, raster = FALSE, repel = TRUE, pt.size = 0.1) + NoLegend()
  p6 <- DimPlot(so_cell_class_list[[cell_class]], reduction = "umap.rpca", group.by = "rpca_clusters", label = TRUE, raster = FALSE, repel = TRUE, pt.size = 0.1) + NoLegend()
  
  pdf(paste0(dir, cell_class, "_RPCA_UMAP_20240823.pdf"), width = 45, height = 30)
  print(p1 + p2 + p3 + p4 + p5 + p6)
  dev.off()
}

saveRDS(so_cell_class_list, paste0(dir, "Objects/so_cell_class_list_20240827.rds"))

cds_cell_class_list <- lapply(so_cell_class_list, function(so) {
  as.cell_data_set(JoinLayers(so))
})

cds_cell_class_list <- lapply(cds_cell_class_list, function(cur_cds) {
  rowData(cur_cds) <- rowData(cds)[rownames(rowData(cur_cds)),]
  return(cur_cds)
})

cds_cell_class_list <- lapply(cds_cell_class_list, function(cur_cds) {
  cur_cds@clusters@listData[["UMAP.RPCA"]]$cluster_result <- NULL
  return(cur_cds)
})

saveRDS(cds_cell_class_list, paste0(dir, "Objects/cds_cell_class_list.clean.rds"))

library(SeuratWrappers)

cds.rpca <- as.cell_data_set(common_so.rpca)
rowData(cds.rpca) <- rowData(cds)[rownames(rowData(cds.rpca)),]

# Bob might have an issue with igraph
# I can try to get rid of that part of the cds
cds.rpca@clusters@listData[["UMAP.RPCA"]]$cluster_result <- NULL

saveRDS(cds.rpca, paste0(dir, "Objects/cds.rpca.clean.rds"))
cds.rpca <- readRDS(paste0(dir, "Objects/cds.rpca.clean.rds"))

saveRDS(cds.rpca@int_colData$reducedDims$UMAP.RPCA, paste0(dir, "Objects/UMAP_rpca_coordinates_global.rds"))
coord_rpca <- readRDS(paste0(dir, "Objects/UMAP_rpca_coordinates_global.rds"))

cds_filt <- readRDS(paste0(dir, "Objects/cds_bg_20240829.rds"))

# add the rpca clusters to the cds_filt
# global cluster then _ then sub_umap_cluster
colData(cds_filt)$global_cluster <- NA
colData(cds_filt)[rownames(common_so.rpca[["seurat_clusters"]]),]$global_cluster <- common_so.rpca[["seurat_clusters"]][,1]

colData(cds_filt)$global_local_cluster <- NA
for(cell_class in c("Muscle", "Hypodermis and seam", "Ciliated neurons", "Non-ciliated neurons", "Pharynx and rectal", "Glia and excretory", "Intestine", "Early", "Mesoderm", "Germline")) {
  colData(cds_filt)[rownames(so_cell_class_list[[cell_class]][["rpca_clusters"]]),]$global_local_cluster <- paste0(colData(cds_filt)[rownames(so_cell_class_list[[cell_class]][["rpca_clusters"]]),]$global_cluster, "_", so_cell_class_list[[cell_class]][["rpca_clusters"]][,1])
}

saveRDS(cds_filt, paste0(dir, "Objects/cds_bg_20240903.rds"))

saveRDS(tissue_barcodes, paste0(dir, "Objects/tissue_barcodes_20240903.rds"))
