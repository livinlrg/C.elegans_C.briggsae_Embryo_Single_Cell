library(presto) # for faster seurat wilcox test
library(Seurat)
library(monocle3)
library(parallel)
library(dplyr)

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240829.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

rm(cds)

common_cds <- cds[rowData(cds)$gene.type == "common",]

so_common <- readRDS(paste0(dir, "Objects/seurat_rpca_integration_joined_20240823.rds"))

BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))
ProBarcodeBinsList <- readRDS(paste0(dir, "Objects/ProBarcodeBinsList.rds"))
JointBarcodeBinsList <- readRDS(paste0(dir, "Objects/JointBarcodeBinsList.rds"))

# DEG genes
# Create Ident.1 being eleangs and Ident.2 being briggsae
# barcodes [[species]][[cell_type_time_bin]][[time_bin]]
DEG_list <- mclapply(names(JointBarcodeBinsList[["C.elegans"]]), function(cur_cell_type) {
  print(cur_cell_type)
  DEG_temp <- list()
  for(i in 1:length(JointBarcodeBinsList[["C.elegans"]][[cur_cell_type]])) {
    cel_barcodes <- JointBarcodeBinsList[["C.elegans"]][[cur_cell_type]][[i]]
    cbr_barcodes <- JointBarcodeBinsList[["C.briggsae"]][[cur_cell_type]][[i]]
    
    so_common_temp <- subset(so_common, cells = c(cel_barcodes, cbr_barcodes))
    Idents(so_common_temp) <- "species"
    
    cur_cell_type_bin <- names(JointBarcodeBinsList[["C.elegans"]][[cur_cell_type]])[i]
    
    DEG_temp[[cur_cell_type_bin]] <- FindMarkers(object = so_common_temp, ident.1 = "C.elegans", ident.2 = "C.briggsae",
                            min.pct = 0, logfc.threshold = 0, test.use = "wilcox", only.pos = FALSE)
    DEG_temp[[cur_cell_type_bin]]$cell_type <- cur_cell_type_bin
    DEG_temp[[cur_cell_type_bin]]$gene <- rownames(DEG_temp[[cur_cell_type_bin]])
  }
  return(DEG_temp)
}, mc.cores = 32)

DEG_df <- do.call(rbind, lapply(DEG_list, function(DEG) {
  return(do.call(rbind, DEG))
}))

saveRDS(DEG_df, paste0(dir, "Objects/DEG_df.rds"))
# add in a lot of metadata columns later in the script


##############################
# Markers
##############################

# Marker calling for terminal cell types
# Going to want to bulk cell barcodes that are from the same time together
TermMarker_list <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  print(cur_species)
  so_temp <- CreateSeuratObject(counts = counts(cds_filt[rowData(cds_filt)$gene.type %in% c("common", ifelse(cur_species == "C.elegans",
                                                                                                             "ce.unique",
                                                                                                             "cb.unique")),
                                                         colData(cds_filt)$species == cur_species]),
                                       project = "Markers",
                                       min.features = 0,
                                       min.cells = 0)
  so_temp <- NormalizeData(object = so_temp)
  
  TermMarker_list[[cur_species]] <- mclapply(names(BarcodeBinsList[[cur_species]]), function(cur_cell_type) {
    print(cur_cell_type)
    so_temp@meta.data$cell_type <- NA
    
    # pull barcodes of this cell type
    barcodes <- unlist(BarcodeBinsList[[cur_species]][[cur_cell_type]])
    so_temp@meta.data[barcodes,]$cell_type <- cur_cell_type
    
    Idents(so_temp) <- "cell_type"
    temp_markers <- FindMarkers(object = so_temp, ident.1 = cur_cell_type,
                                min.pct = 0, logfc.threshold = 0, test.use = "wilcox", only.pos = TRUE)
    temp_markers$cell_type <- cur_cell_type
    return(temp_markers)
  }, mc.cores = 4)
}

# need to add gene name to the markers
for(cur_species in c("C.elegans", "C.briggsae")) {
  TermMarker_list[[cur_species]] <- lapply(TermMarker_list[[cur_species]], function(cur_cell_type_markers) {
    cur_cell_type_markers$gene <- rownames(cur_cell_type_markers)
    return(cur_cell_type_markers)
  })
}

TemMarker_list_stored <- list()
TemMarker_list_stored[["C.elegans"]] <- TermMarker_list[["C.elegans"]]
TemMarker_list_stored[["C.briggsae"]] <- TermMarker_list[["C.briggsae"]]

# TermMarker_list[["C.elegans"]] <- TemMarker_list_stored[["C.elegans"]]
# TermMarker_list[["C.briggsae"]] <- TemMarker_list_stored[["C.briggsae"]]

TermMarker_list[["C.elegans"]] <- do.call(rbind, TermMarker_list[["C.elegans"]])
TermMarker_list[["C.briggsae"]] <- do.call(rbind, TermMarker_list[["C.briggsae"]])

saveRDS(TermMarker_list, paste0(dir, "Objects/TermMarker_list.rds"))

# Progenitor marker calling
ProMarker_list <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  pro_cds_temp <- cds_filt[rowData(cds_filt)$gene.type %in% c("common", ifelse(cur_species == "C.elegans",
                                                                          "ce.unique",
                                                                          "cb.unique")),
                           colData(cds_filt)$species == cur_species & ! is.na(colData(cds_filt)$lineage_broad)]
  so_temp <- CreateSeuratObject(counts = counts(pro_cds_temp),
                                project = "Markers",
                                min.features = 0,
                                min.cells = 0)
  so_temp <- NormalizeData(object = so_temp)
  so_temp <- AddMetaData(so_temp, data.frame(colData(pro_cds_temp)))
  Idents(so_temp) <- "lineage_broad"
  
  ProMarker_list[[cur_species]] <- FindAllMarkers(object = so_temp,
                                                  min.pct = 0, logfc.threshold = 0, test.use = "wilcox", only.pos = TRUE)
}

saveRDS(ProMarker_list, paste0(dir, "Objects/ProMarker_list.rds"))

'
shared_marker <- Marker in both and orthologous

cel_marker_not_shared <- Marker in cel, not expressed in that cell in the other
cbr_marker_not_shared <- Marker in cbr, not expressed in that cell in the other

cel_marker_still_detected <- Marker in cel, expressed more highly in other cells in the other dataset or poorly detected in that cell, but still detected in other species
cbr_marker_still_detected <- Marker in cbr, expressed more highly in other cells in the other dataset or poorly detected in that cell, but still detected in other species

not_marker <- Not detected in either

non_ortho_marker <- Marker in one species, but not a 1:1 ortholog
'

################################
# Add TerminalMarkers metadata #
################################

ProMarker_list <- readRDS(paste0(dir, "Objects/ProMarker_list.rds"))
TermMarker_list <- readRDS(paste0(dir, "Objects/TermMarker_list.rds"))

TPMListBootstrap_term <- readRDS(paste0(dir, "Objects/TPMTimeBinsListBootstrap_CellCorrection.rds"))
TPMListBootstrapMean_term <- readRDS(paste0(dir, "Objects/TPMListBootstrapMean_CellCorrection.rds"))
TPMListBootstrap_pro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))

tpm_matrix_list <- readRDS(paste0(dir, "Objects/tpm_matrix_list_cell_bg.rds"))
tpm_matrix_list_filt <- readRDS(paste0(dir, "Objects/tpm_matrix_list_filt_cell_bg.rds"))

CellTable <- readRDS(paste0(dir, "Objects/CellTable_Names_20240901.rds"))
gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))

#####################
# Value Calculation #
#####################

# Calculate Tau using the bootstrap TPM values for all genes (non-orthologous as well) for markers (mean cell types)
# Tau

TauList_pro <- list()
TauList_term <- list()
TauList_joint <- list()

#Tau = Sum(1-expNormalizedByMax)/1-numberOfTissues
for(cur_species in c("C.elegans", "C.briggsae")) {
  TauList_joint[[cur_species]] <- apply(log(tpm_matrix_list_filt[[cur_species]] + 1),
                                        1, function(row) {
                                          return(sum(1 - (row/max(row)))/(length(row) - 1))
                                        })
  
  TauList_pro[[cur_species]] <- apply(log(TPMListBootstrap_pro[[cur_species]] + 1),
                                      1, function(row) {
                                        return(sum(1 - (row/max(row)))/(length(row) - 1))
                                      })
  
  TauList_term[[cur_species]] <- apply(log(TPMListBootstrapMean_term[[cur_species]] + 1),
                                       1, function(row) {
                                         return(sum(1 - (row/max(row)))/(length(row) - 1))
                                       })
}

saveRDS(TauList_pro, paste0(dir, "Objects/TauList_pro_bootstrapTPM_cell_bg.rds"))
saveRDS(TauList_term, paste0(dir, "Objects/TauList_term_bootstrapTPM_cell_bg.rds"))
saveRDS(TauList_joint, paste0(dir, "Objects/TauList_joint_bootstrapTPM_cell_bg.rds"))

TauList_pro <- readRDS(paste0(dir, "Objects/TauList_pro_bootstrapTPM_cell_bg.rds"))
TauList_term <- readRDS(paste0(dir, "Objects/TauList_term_bootstrapTPM_cell_bg.rds"))
TauList_joint <- readRDS(paste0(dir, "Objects/TauList_joint_bootstrapTPM_cell_bg.rds"))

# make a matrix of log2fc values for each cell type in each species
log2fc_matrix_list <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  log2fc_matrix_list[[cur_species]] <- matrix(NA,
                                              nrow = nrow(tpm_matrix_list[[cur_species]]),
                                              ncol = length(unique(TermMarker_list[[cur_species]]$cell_type)))
  colnames(log2fc_matrix_list[[cur_species]]) <- unique(TermMarker_list[[cur_species]]$cell_type)
  rownames(log2fc_matrix_list[[cur_species]]) <- rownames(tpm_matrix_list[[cur_species]])

  for(cur_cell_type in unique(TermMarker_list[[cur_species]]$cell_type)) {
    cur_cell_type_lineages <- CellTable[CellTable$MergedDatasetName %in% cur_cell_type,]$Lineage
    cur_cell_type_sans_lineages <- CellTable[! CellTable$MergedDatasetName %in% cur_cell_type,]$Lineage
    log2fc_matrix_list[[cur_species]][,cur_cell_type] <- log2(rowMeans(tpm_matrix_list[[cur_species]][, cur_cell_type_lineages, drop = FALSE]) + 1) -
      log2(rowMeans(tpm_matrix_list[[cur_species]][, cur_cell_type_sans_lineages, drop = FALSE], na.rm = TRUE) + 1)
  }
}

max_tpm_list_term <- lapply(TPMListBootstrapMean_term, function(tpm_table) apply(tpm_table, 1, max))
max_tpm_list_pro <- lapply(TPMListBootstrap_pro, function(tpm_table) apply(tpm_table, 1, max))

mean_tpm_list_term <- lapply(TPMListBootstrapMean_term, function(tpm_table) apply(tpm_table, 1, mean))
mean_tpm_list_pro <- lapply(TPMListBootstrap_pro, function(tpm_table) apply(tpm_table, 1, mean))

# add the cell_gene column
TermMarker_list[["C.elegans"]]$cell_gene <- paste0(TermMarker_list[["C.elegans"]]$cell_type, "_", TermMarker_list[["C.elegans"]]$gene)
rownames(TermMarker_list[["C.elegans"]]) <- TermMarker_list[["C.elegans"]]$cell_gene
TermMarker_list[["C.briggsae"]]$cell_gene <- paste0(TermMarker_list[["C.briggsae"]]$cell_type, "_", TermMarker_list[["C.briggsae"]]$gene)
rownames(TermMarker_list[["C.briggsae"]]) <- TermMarker_list[["C.briggsae"]]$cell_gene

TermMarker_list_match <- list()
# reciprocal species hits
TermMarker_list_match[["C.elegans"]] <- left_join(TermMarker_list[["C.elegans"]], TermMarker_list[["C.briggsae"]],
                                            by = c("cell_gene", "cell_type", "gene"),
                                            suffix = c(".cel", ".cbr"))
TermMarker_list_match[["C.briggsae"]] <- left_join(TermMarker_list[["C.briggsae"]], TermMarker_list[["C.elegans"]],
                                             by = c("cell_gene", "cell_type", "gene"),
                                             suffix = c(".cbr", ".cel"))

TermMarker_list <- TermMarker_list_match

for(cur_species in c("C.elegans", "C.briggsae")) {
  if("E-BE45912.2" %in% TermMarker_list[[cur_species]]$gene) {
    TermMarker_list[[cur_species]][TermMarker_list[[cur_species]]$gene %in% "E-BE45912.2",]$gene <- "E_BE45912.2"
  }
  
  # add log2fc values to the markers
  TermMarker_list[[cur_species]]$cel_tpm_log2fc <- apply(TermMarker_list[[cur_species]], 1, function(row) {
    if(row["gene"] %in% rownames(log2fc_matrix_list[["C.elegans"]])) {
      log2fc_matrix_list[["C.elegans"]][row["gene"], row["cell_type"]]
    } else {
      NA
    }
  })
  
  TermMarker_list[[cur_species]]$cbr_tpm_log2fc <- apply(TermMarker_list[[cur_species]], 1, function(row) {
    if(row["gene"] %in% rownames(log2fc_matrix_list[["C.briggsae"]])) {
      log2fc_matrix_list[["C.briggsae"]][row["gene"], row["cell_type"]]
    } else {
      NA
    }
  })

  # focal cell type tpm
  TermMarker_list[[cur_species]]$cel_tpm <- apply(TermMarker_list[[cur_species]], 1, function(row) {
    if(row["gene"] %in% rownames(TPMListBootstrapMean_term[["C.elegans"]])) {
      TPMListBootstrapMean_term[["C.elegans"]][row["gene"], row["cell_type"]]
    } else {
      NA
    }
  })

  TermMarker_list[[cur_species]]$cbr_tpm <- apply(TermMarker_list[[cur_species]], 1, function(row) {
    if(row["gene"] %in% rownames(TPMListBootstrapMean_term[["C.briggsae"]])) {
      TPMListBootstrapMean_term[["C.briggsae"]][row["gene"], row["cell_type"]]
    } else {
      NA
    }
  })

  # max tpm
  TermMarker_list[[cur_species]]$cel_max_tpm_term <- max_tpm_list_term[["C.elegans"]][TermMarker_list[[cur_species]]$gene]
  TermMarker_list[[cur_species]]$cbr_max_tpm_term <- max_tpm_list_term[["C.briggsae"]][TermMarker_list[[cur_species]]$gene]
  
  TermMarker_list[[cur_species]]$cel_max_tpm_pro <- max_tpm_list_pro[["C.elegans"]][TermMarker_list[[cur_species]]$gene]
  TermMarker_list[[cur_species]]$cbr_max_tpm_pro <- max_tpm_list_pro[["C.briggsae"]][TermMarker_list[[cur_species]]$gene]
  
  # mean tpm
  TermMarker_list[[cur_species]]$cel_mean_tpm_term <- mean_tpm_list_term[["C.elegans"]][TermMarker_list[[cur_species]]$gene]
  TermMarker_list[[cur_species]]$cbr_mean_tpm_term <- mean_tpm_list_term[["C.briggsae"]][TermMarker_list[[cur_species]]$gene]
  
  TermMarker_list[[cur_species]]$cel_mean_tpm_pro <- mean_tpm_list_pro[["C.elegans"]][TermMarker_list[[cur_species]]$gene]
  TermMarker_list[[cur_species]]$cbr_mean_tpm_pro <- mean_tpm_list_pro[["C.briggsae"]][TermMarker_list[[cur_species]]$gene]

  # tau
  TermMarker_list[[cur_species]]$cel_tau_pro <- TauList_pro[["C.elegans"]][TermMarker_list[[cur_species]]$gene]
  TermMarker_list[[cur_species]]$cbr_tau_pro <- TauList_pro[["C.briggsae"]][TermMarker_list[[cur_species]]$gene]
  
  TermMarker_list[[cur_species]]$cel_tau_term <- TauList_term[["C.elegans"]][TermMarker_list[[cur_species]]$gene]
  TermMarker_list[[cur_species]]$cbr_tau_term <- TauList_term[["C.briggsae"]][TermMarker_list[[cur_species]]$gene]
  
  TermMarker_list[[cur_species]]$cel_tau_joint <- TauList_joint[["C.elegans"]][TermMarker_list[[cur_species]]$gene]
  TermMarker_list[[cur_species]]$cbr_tau_joint <- TauList_joint[["C.briggsae"]][TermMarker_list[[cur_species]]$gene]
}

# add orthogroup number, status, and worm.cat
# size: 994956
TermMarker_list[["C.elegans"]] <- left_join(TermMarker_list[["C.elegans"]],
                                            gff_list_mrna[["elegans"]][,c("cds_gene_name", "gene.type", "orthology_conf",
                                                                          "OG", "cel_OG_count", "cbr_OG_count",
                                                                          "WormCat.1", "WormCat.2", "WormCat.3")],
                                            by = c("gene" = "cds_gene_name"), na_matches = "never")

# size: 1178644
TermMarker_list[["C.briggsae"]] <- left_join(TermMarker_list[["C.briggsae"]],
                                             gff_list_mrna[["briggsae"]][,c("cds_gene_name", "gene.type", "orthology_conf",
                                                                            "OG", "cel_OG_count", "cbr_OG_count",
                                                                            "WormCat.1", "WormCat.2", "WormCat.3")],
                                             by = c("gene" = "cds_gene_name"), na_matches = "never")

rownames(TermMarker_list[["C.elegans"]]) <- TermMarker_list[["C.elegans"]]$cell_gene
rownames(TermMarker_list[["C.briggsae"]]) <- TermMarker_list[["C.briggsae"]]$cell_gene

col_order_cel <- c("cell_type", "gene", "cell_gene",
                   "p_val.cel", "avg_log2FC.cel", "pct.1.cel", "pct.2.cel", "p_val_adj.cel",
                   "p_val.cbr", "avg_log2FC.cbr", "pct.1.cbr", "pct.2.cbr", "p_val_adj.cbr",
                   "cel_tpm", "cbr_tpm", "cel_tpm_log2fc", "cbr_tpm_log2fc",
                   "cel_max_tpm_term", "cbr_max_tpm_term", "cel_max_tpm_pro", "cbr_max_tpm_pro",
                   "cel_mean_tpm_term", "cbr_mean_tpm_term", "cel_mean_tpm_pro", "cbr_mean_tpm_pro",
                   "cel_tau_term", "cbr_tau_term", "cel_tau_pro", "cbr_tau_pro", "cel_tau_joint", "cbr_tau_joint",
                   "gene.type", "orthology_conf", "OG", "cel_OG_count", "cbr_OG_count", "WormCat.1", "WormCat.2", "WormCat.3")

col_order_cbr <- c("cell_type", "gene", "cell_gene",
                   "p_val.cbr", "avg_log2FC.cbr", "pct.1.cbr", "pct.2.cbr", "p_val_adj.cbr",
                   "p_val.cel", "avg_log2FC.cel", "pct.1.cel", "pct.2.cel", "p_val_adj.cel",
                   "cbr_tpm", "cel_tpm", "cbr_tpm_log2fc", "cel_tpm_log2fc",
                   "cbr_max_tpm_term", "cel_max_tpm_term", "cbr_max_tpm_pro", "cel_max_tpm_pro",
                   "cbr_mean_tpm_term", "cel_mean_tpm_term", "cbr_mean_tpm_pro", "cel_mean_tpm_pro",
                   "cbr_tau_term", "cel_tau_term", "cbr_tau_pro", "cel_tau_pro", "cbr_tau_joint", "cel_tau_joint",
                   "gene.type", "orthology_conf", "OG", "cel_OG_count", "cbr_OG_count", "WormCat.1", "WormCat.2", "WormCat.3")

TermMarker_list[["C.elegans"]] <- TermMarker_list[["C.elegans"]][,col_order_cel]
TermMarker_list[["C.briggsae"]] <- TermMarker_list[["C.briggsae"]][,col_order_cbr]

saveRDS(TermMarker_list, paste0(dir, "Objects/TermMarker_list_metadata_cell_bg.rds"))

####################################
# Add ProgenitorsMarkers metadata
####################################

colnames(ProMarker_list[["C.elegans"]])[6] <- "cell_type"
colnames(ProMarker_list[["C.briggsae"]])[6] <- "cell_type"

ProMarker_list[["C.elegans"]]$cell_type <- as.character(ProMarker_list[["C.elegans"]]$cell_type)
ProMarker_list[["C.briggsae"]]$cell_type <- as.character(ProMarker_list[["C.briggsae"]]$cell_type)

ProMarker_list[["C.elegans"]] <- ProMarker_list[["C.elegans"]][ProMarker_list[["C.elegans"]]$cell_type != "unassigned",]
ProMarker_list[["C.briggsae"]] <- ProMarker_list[["C.briggsae"]][ProMarker_list[["C.briggsae"]]$cell_type != "unassigned",]

ProMarker_list_stored <- list()
ProMarker_list_stored[["C.elegans"]] <- ProMarker_list[["C.elegans"]]
ProMarker_list_stored[["C.briggsae"]] <- ProMarker_list[["C.briggsae"]]

# ProMarker_list[["C.elegans"]] <- ProMarker_list_stored[["C.elegans"]]
# ProMarker_list[["C.briggsae"]] <- ProMarker_list_stored[["C.briggsae"]]

# make a matrix of log2fc values for each cell type in each species
log2fc_matrix_list_pro <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  log2fc_matrix_list_pro[[cur_species]] <- matrix(NA,
                                              nrow = nrow(tpm_matrix_list[[cur_species]]),
                                              ncol = length(unique(ProMarker_list[[cur_species]]$cell_type)))
  colnames(log2fc_matrix_list_pro[[cur_species]]) <- sort(unique(ProMarker_list[[cur_species]]$cell_type))
  rownames(log2fc_matrix_list_pro[[cur_species]]) <- rownames(tpm_matrix_list[[cur_species]])
  
  for(cur_cell_type in unique(ProMarker_list[[cur_species]]$cell_type)) {
    cur_cell_type_lineages <- CellTable[CellTable$MergedDatasetName %in% cur_cell_type,]$Lineage
    cur_cell_type_sans_lineages <- CellTable[! CellTable$MergedDatasetName %in% cur_cell_type,]$Lineage
    log2fc_matrix_list_pro[[cur_species]][,cur_cell_type] <- log2(rowMeans(tpm_matrix_list[[cur_species]][, cur_cell_type_lineages, drop = FALSE], na.rm = TRUE) + 1) -
      log2(rowMeans(tpm_matrix_list[[cur_species]][, cur_cell_type_sans_lineages, drop = FALSE], na.rm = TRUE) + 1)
  }
}

# make a matrix of log2fc values for each cell type in each species without the terminal cell types
log2fc_matrix_list_just_pro <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  log2fc_matrix_list_just_pro[[cur_species]] <- matrix(NA,
                                                  nrow = nrow(tpm_matrix_list[[cur_species]]),
                                                  ncol = length(unique(ProMarker_list[[cur_species]]$cell_type)))
  colnames(log2fc_matrix_list_just_pro[[cur_species]]) <- sort(unique(ProMarker_list[[cur_species]]$cell_type))
  rownames(log2fc_matrix_list_just_pro[[cur_species]]) <- rownames(tpm_matrix_list[[cur_species]])
  
  for(cur_cell_type in unique(ProMarker_list[[cur_species]]$cell_type)) {
    cur_cell_type_lineages <- CellTable[CellTable$MergedDatasetName %in% cur_cell_type,]$Lineage
    cur_cell_type_sans_lineages <- CellTable[(! CellTable$MergedDatasetName %in% cur_cell_type) & CellTable$Cell == "progenitor",]$Lineage
    log2fc_matrix_list_just_pro[[cur_species]][,cur_cell_type] <- log2(rowMeans(tpm_matrix_list[[cur_species]][, cur_cell_type_lineages, drop = FALSE], na.rm = TRUE) + 1) -
      log2(rowMeans(tpm_matrix_list[[cur_species]][, cur_cell_type_sans_lineages, drop = FALSE], na.rm = TRUE) + 1)
  }
}

# add the cell_gene column
ProMarker_list[["C.elegans"]]$cell_gene <- paste0(ProMarker_list[["C.elegans"]]$cell_type, "_", ProMarker_list[["C.elegans"]]$gene)
rownames(ProMarker_list[["C.elegans"]]) <- ProMarker_list[["C.elegans"]]$cell_gene
ProMarker_list[["C.briggsae"]]$cell_gene <- paste0(ProMarker_list[["C.briggsae"]]$cell_type, "_", ProMarker_list[["C.briggsae"]]$gene)
rownames(ProMarker_list[["C.briggsae"]]) <- ProMarker_list[["C.briggsae"]]$cell_gene

# reciprocal species hits
ProMarker_list_match <- list()
ProMarker_list_match[["C.elegans"]] <- left_join(ProMarker_list[["C.elegans"]], ProMarker_list[["C.briggsae"]],
                                            by = c("cell_gene", "cell_type", "gene"),
                                            suffix = c(".cel", ".cbr"))
ProMarker_list_match[["C.briggsae"]] <- left_join(ProMarker_list[["C.briggsae"]], ProMarker_list[["C.elegans"]],
                                             by = c("cell_gene", "cell_type", "gene"),
                                             suffix = c(".cbr", ".cel"))

ProMarker_list <- ProMarker_list_match

for(cur_species in c("C.elegans", "C.briggsae")) {
  if("E-BE45912.2" %in% ProMarker_list[[cur_species]]$gene) {
    ProMarker_list[[cur_species]][ProMarker_list[[cur_species]]$gene %in% "E-BE45912.2",]$gene <- "E_BE45912.2"
  }
  
  # add log2fc values to the markers
  ProMarker_list[[cur_species]]$cel_tpm_log2fc <- apply(ProMarker_list[[cur_species]], 1, function(row) {
    if(row["gene"] %in% rownames(log2fc_matrix_list_pro[["C.elegans"]])) {
      log2fc_matrix_list_pro[["C.elegans"]][row["gene"], row["cell_type"]]
    } else {
      NA
    }
  })
  
  ProMarker_list[[cur_species]]$cbr_tpm_log2fc <- apply(ProMarker_list[[cur_species]], 1, function(row) {
    if(row["gene"] %in% rownames(log2fc_matrix_list_pro[["C.briggsae"]])) {
      log2fc_matrix_list_pro[["C.briggsae"]][row["gene"], row["cell_type"]]
    } else {
      NA
    }
  })
  
  ProMarker_list[[cur_species]]$cel_tpm_log2fc_just_pro <- apply(ProMarker_list[[cur_species]], 1, function(row) {
    if(row["gene"] %in% rownames(log2fc_matrix_list_just_pro[["C.elegans"]])) {
      log2fc_matrix_list_just_pro[["C.elegans"]][row["gene"], row["cell_type"]]
    } else {
      NA
    }
  })
  
  ProMarker_list[[cur_species]]$cbr_tpm_log2fc_just_pro <- apply(ProMarker_list[[cur_species]], 1, function(row) {
    if(row["gene"] %in% rownames(log2fc_matrix_list_just_pro[["C.briggsae"]])) {
      log2fc_matrix_list_just_pro[["C.briggsae"]][row["gene"], row["cell_type"]]
    } else {
      NA
    }
  })
  
  # focal cell type tpm
  ProMarker_list[[cur_species]]$cel_tpm <- apply(ProMarker_list[[cur_species]], 1, function(row) {
    if(row["gene"] %in% rownames(TPMListBootstrap_pro[["C.elegans"]])) {
      TPMListBootstrap_pro[["C.elegans"]][row["gene"], row["cell_type"]]
    } else {
      NA
    }
  })
  
  ProMarker_list[[cur_species]]$cbr_tpm <- apply(ProMarker_list[[cur_species]], 1, function(row) {
    if(row["gene"] %in% rownames(TPMListBootstrap_pro[["C.briggsae"]])) {
      TPMListBootstrap_pro[["C.briggsae"]][row["gene"], row["cell_type"]]
    } else {
      NA
    }
  })

  # max tpm
  ProMarker_list[[cur_species]]$cel_max_tpm_term <- max_tpm_list_term[["C.elegans"]][ProMarker_list[[cur_species]]$gene]
  ProMarker_list[[cur_species]]$cbr_max_tpm_term <- max_tpm_list_term[["C.briggsae"]][ProMarker_list[[cur_species]]$gene]
  
  ProMarker_list[[cur_species]]$cel_max_tpm_pro <- max_tpm_list_pro[["C.elegans"]][ProMarker_list[[cur_species]]$gene]
  ProMarker_list[[cur_species]]$cbr_max_tpm_pro <- max_tpm_list_pro[["C.briggsae"]][ProMarker_list[[cur_species]]$gene]
  
  # mean tpm
  ProMarker_list[[cur_species]]$cel_mean_tpm_term <- mean_tpm_list_term[["C.elegans"]][ProMarker_list[[cur_species]]$gene]
  ProMarker_list[[cur_species]]$cbr_mean_tpm_term <- mean_tpm_list_term[["C.briggsae"]][ProMarker_list[[cur_species]]$gene]
  
  ProMarker_list[[cur_species]]$cel_mean_tpm_pro <- mean_tpm_list_pro[["C.elegans"]][ProMarker_list[[cur_species]]$gene]
  ProMarker_list[[cur_species]]$cbr_mean_tpm_pro <- mean_tpm_list_pro[["C.briggsae"]][ProMarker_list[[cur_species]]$gene]
  
  # tau
  ProMarker_list[[cur_species]]$cel_tau_pro <- TauList_pro[["C.elegans"]][ProMarker_list[[cur_species]]$gene]
  ProMarker_list[[cur_species]]$cbr_tau_pro <- TauList_pro[["C.briggsae"]][ProMarker_list[[cur_species]]$gene]
  
  ProMarker_list[[cur_species]]$cel_tau_term <- TauList_term[["C.elegans"]][ProMarker_list[[cur_species]]$gene]
  ProMarker_list[[cur_species]]$cbr_tau_term <- TauList_term[["C.briggsae"]][ProMarker_list[[cur_species]]$gene]
  
  ProMarker_list[[cur_species]]$cel_tau_joint <- TauList_joint[["C.elegans"]][ProMarker_list[[cur_species]]$gene]
  ProMarker_list[[cur_species]]$cbr_tau_joint <- TauList_joint[["C.briggsae"]][ProMarker_list[[cur_species]]$gene]
}

# ProMarker_list[["C.elegans"]] <- ProMarker_list[["C.elegans"]][,! colnames(ProMarker_list[["C.elegans"]]) %in% c("gene.type", "orthology_conf",
#                                                                                 "OG", "cel_OG_count", "cbr_OG_count",
#                                                                                 "WormCat.1", "WormCat.2", "WormCat.3")]

# add orthogroup number, status, and worm.cat
# size: 678650
ProMarker_list[["C.elegans"]] <- left_join(ProMarker_list[["C.elegans"]],
                                            gff_list_mrna[["elegans"]][,c("cds_gene_name", "gene.type", "orthology_conf",
                                                                          "OG", "cel_OG_count", "cbr_OG_count",
                                                                          "WormCat.1", "WormCat.2", "WormCat.3")],
                                            by = c("gene" = "cds_gene_name"), na_matches = "never")

# size: 572811
ProMarker_list[["C.briggsae"]] <- left_join(ProMarker_list[["C.briggsae"]],
                                             gff_list_mrna[["briggsae"]][,c("cds_gene_name", "gene.type", "orthology_conf",
                                                                            "OG", "cel_OG_count", "cbr_OG_count",
                                                                            "WormCat.1", "WormCat.2", "WormCat.3")],
                                             by = c("gene" = "cds_gene_name"), na_matches = "never")

rownames(ProMarker_list[["C.elegans"]]) <- ProMarker_list[["C.elegans"]]$cell_gene
rownames(ProMarker_list[["C.briggsae"]]) <- ProMarker_list[["C.briggsae"]]$cell_gene

saveRDS(ProMarker_list, paste0(dir, "Objects/ProMarker_list_metadata_cell_bg.rds"))

TermMarker_list <- readRDS(paste0(dir, "Objects/TermMarker_list_metadata_cell_bg.rds"))
ProMarker_list <- readRDS(paste0(dir, "Objects/ProMarker_list_metadata_cell_bg.rds"))

#### Filter markers list

# Number of markers
term_markers_filt <- list()
term_markers_filt[["C.elegans"]] <- TermMarker_list[["C.elegans"]] %>%
  filter(cel_tpm >= 80 & cel_tpm_log2fc >= 1 & p_val_adj.cel < 0.05 & avg_log2FC.cel >= 1)
term_markers_filt[["C.briggsae"]] <- TermMarker_list[["C.briggsae"]] %>%
  filter(cbr_tpm >= 80 & cbr_tpm_log2fc >= 1 & p_val_adj.cbr < 0.05 & avg_log2FC.cbr >= 1)

pro_markers_filt <- list()
pro_markers_filt[["C.elegans"]] <- ProMarker_list[["C.elegans"]] %>%
  filter(cel_tpm >= 80 & cel_tpm_log2fc_just_pro >= 1 & p_val_adj.cel < 0.05 & avg_log2FC.cel >= 1)
pro_markers_filt[["C.briggsae"]] <- ProMarker_list[["C.briggsae"]] %>%
  filter(cbr_tpm >= 80 & cbr_tpm_log2fc_just_pro >= 1 & p_val_adj.cbr < 0.05 & avg_log2FC.cbr >= 1)

saveRDS(term_markers_filt, paste0(dir, "Objects/TermMarker_list_filt_metadata_cell_bg.rds"))
saveRDS(pro_markers_filt, paste0(dir, "Objects/ProMarker_list_filt_metadata_cell_bg.rds"))


# Add in DEG_df metadata columns

DEG_df <- readRDS(paste0(dir, "Objects/DEG_df.rds"))
cell_data_term <- readRDS(paste0(dir, "Objects/cell_data_term_cell_bg.rds"))

# size: 8439943

DEG_df$cell_type_bin <- DEG_df$cell_type
DEG_df$cell_type <- ifelse(DEG_df$cell_type_bin %in% cell_data_term$cell_type_bin, 
                           cell_data_term[DEG_df$cell_type_bin,]$cell_type,
                           DEG_df$cell_type_bin)

TPMListComb <- list()
TPMListComb[["C.elegans"]] <- cbind(TPMListBootstrap_term[["C.elegans"]], TPMListBootstrap_pro[["C.elegans"]])
TPMListComb[["C.briggsae"]] <- cbind(TPMListBootstrap_term[["C.briggsae"]], TPMListBootstrap_pro[["C.briggsae"]])

# focal cell type tpm
DEG_df$cel_tpm <- apply(DEG_df, 1, function(row) {
  return(TPMListComb[["C.elegans"]][row["gene"], row["cell_type_bin"]])
})

DEG_df$cbr_tpm <- apply(DEG_df, 1, function(row) {
  return(TPMListComb[["C.briggsae"]][row["gene"], row["cell_type_bin"]])
})

# max tpm
max_tpm_list_term <- lapply(TPMListBootstrap_term, function(tpm_table) apply(tpm_table, 1, max))
max_tpm_list_pro <- lapply(TPMListBootstrap_pro, function(tpm_table) apply(tpm_table, 1, max))

DEG_df$cel_max_tpm_term <- max_tpm_list_term[["C.elegans"]][DEG_df$gene]
DEG_df$cbr_max_tpm_term <- max_tpm_list_term[["C.briggsae"]][DEG_df$gene]

DEG_df$cel_max_tpm_pro <- max_tpm_list_pro[["C.elegans"]][DEG_df$gene]
DEG_df$cbr_max_tpm_pro <- max_tpm_list_pro[["C.briggsae"]][DEG_df$gene]

# mean tpm
mean_tpm_list_term <- lapply(TPMListBootstrap_term, function(tpm_table) apply(tpm_table, 1, mean))
mean_tpm_list_pro <- lapply(TPMListBootstrap_pro, function(tpm_table) apply(tpm_table, 1, mean))

DEG_df$cel_mean_tpm_term <- mean_tpm_list_term[["C.elegans"]][DEG_df$gene]
DEG_df$cbr_mean_tpm_term <- mean_tpm_list_term[["C.briggsae"]][DEG_df$gene]

DEG_df$cel_mean_tpm_pro <- mean_tpm_list_pro[["C.elegans"]][DEG_df$gene]
DEG_df$cbr_mean_tpm_pro <- mean_tpm_list_pro[["C.briggsae"]][DEG_df$gene]

TauList_pro <- readRDS(paste0(dir, "Objects/TauList_pro_bootstrapTPM.rds"))
TauList_term <- readRDS(paste0(dir, "Objects/TauList_term_bootstrapTPM.rds"))
TauList_joint <- readRDS(paste0(dir, "Objects/TauList_joint_bootstrapTPM.rds"))

# tau
DEG_df$cel_tau_pro <- TauList_pro[["C.elegans"]][DEG_df$gene]
DEG_df$cbr_tau_pro <- TauList_pro[["C.briggsae"]][DEG_df$gene]
DEG_df$cel_tau_term <- TauList_term[["C.elegans"]][DEG_df$gene]
DEG_df$cbr_tau_term <- TauList_term[["C.briggsae"]][DEG_df$gene]
DEG_df$cel_tau_joint <- TauList_joint[["C.elegans"]][DEG_df$gene]
DEG_df$cbr_tau_joint <- TauList_joint[["C.briggsae"]][DEG_df$gene]

DEG_df <- left_join(DEG_df,
                    gff_list_mrna[["elegans"]][,c("cds_gene_name", "gene.type", "orthology_conf",
                                                  "OG", "cel_OG_count", "cbr_OG_count",
                                                  "WormCat.1", "WormCat.2", "WormCat.3")],
                    by = c("gene" = "cds_gene_name"), na_matches = "never")

saveRDS(DEG_df, paste0(dir, "Objects/DEG_df_metadata_cell_bg.rds"))

deg_filt <- DEG_df %>% filter((cel_tpm > 80 | cbr_tpm > 80) & (avg_log2FC > 1 | avg_log2FC < -1) & p_val_adj < 0.05)
saveRDS(deg_filt, paste0(dir, "Objects/DEG_df_filt_metadata_cell_bg.rds"))


