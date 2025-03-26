
library(Matrix)
library(monocle3)
library(Seurat)

## Create expression matrix

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

## Need 3 things for the loom file creation using python (RunScenic.py)
# WS290_elegans_expression.mtx
# WS290_elegans_barcodes.csv
# WS290_elegans_genes.csv

# source 
# source(paste0(dir, "../Scripts/CalcFunctions.R"))
cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240903.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

rm(cds)

# C. elegans matrix
comb_seurat_cel <- CreateSeuratObject(counts(cds_filt[rowData(cds_filt)$gene.type %in% c("common", "ce.unique"),
                                                  colData(cds_filt)$species == "C.elegans"]),
                                  project = "WS290_elegans",
                                  assay = "RNA")

# setwd(paste0(dir, "SCENIC/"))
setwd("/scratch/livinlrg/WS290/")

# writeMM(comb_seurat_cel@assays$RNA@layers$counts, paste0(dir, "Objects/WS290_elegans_expression.mtx"))
writeMM(comb_seurat_cel@assays$RNA@layers$counts, paste0("/scratch/livinlrg/WS290/WS290_elegans_expression.mtx"))

barcodes <- data.frame(colnames(comb_seurat_cel))
colnames(barcodes) <- 'Barcode'
write.csv(barcodes, paste0("/scratch/livinlrg/WS290/WS290_elegans_barcodes.csv"),
          quote = FALSE, row.names = FALSE)

genes <- data.frame(rownames(comb_seurat_cel))
colnames(genes) <- 'Gene'
write.csv(genes, paste0("/scratch/livinlrg/WS290/WS290_elegans_genes.csv"),
          quote = FALSE, row.names = FALSE)

# C. briggsae matrix
comb_seurat_cbr <- CreateSeuratObject(counts(cds_filt[rowData(cds_filt)$gene.type %in% c("common", "cb.unique"),
                                                      colData(cds_filt)$species == "C.briggsae"]),
                                      project = "WS290_briggsae",
                                      assay = "RNA")

writeMM(comb_seurat_cbr@assays$RNA@layers$counts, paste0("/scratch/livinlrg/WS290/WS290_briggsae_expression.mtx"))

barcodes <- data.frame(colnames(comb_seurat_cbr))
colnames(barcodes) <- 'Barcode'
write.csv(barcodes, paste0("/scratch/livinlrg/WS290/WS290_briggsae_barcodes.csv"),
          quote = FALSE, row.names = FALSE)

genes <- data.frame(rownames(comb_seurat_cbr))
colnames(genes) <- 'Gene'
write.csv(genes,  paste0("/scratch/livinlrg/WS290/WS290_briggsae_genes.csv"),
          quote = FALSE, row.names = FALSE)

## Out for TPM input to SCENIC
cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint.rds"))
cell_data$min_cell_count <- apply(cell_data[,c("cel_cell_count", "cbr_cell_count")], 1, function(x) min(x, na.rm = TRUE))

TPMListBootstrap_term <- readRDS(paste0(dir, "Objects/TPMListBootstrap_CellCorrection.rds"))
TPMListBootstrap_pro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))

cell_type_time_bins <- readRDS(paste0(dir, "Objects/cell_type_time_bins.rds"))

out_ele_tpm <- cbind(TPMListBootstrap_term[["C.elegans"]],
                     TPMListBootstrap_pro[["C.elegans"]])[,cell_data_joint[cell_data_joint$min_cell_count >= 15,]$cell_type_bin]
out_cbr_tpm <- cbind(TPMListBootstrap_term[["C.briggsae"]],
                     TPMListBootstrap_pro[["C.briggsae"]])[,cell_data_joint[cell_data_joint$min_cell_count >= 15,]$cell_type_bin]

write.csv(out_ele_tpm, paste0(dir, "Objects/scenic/cel_expression.csv"))
write.csv(out_cbr_tpm, paste0(dir, "Objects/scenic/cbr_expression.csv"))

write.csv(log(out_ele_tpm + 1, 2), paste0(dir, "Objects/scenic/cel_expression_log.csv"))
write.csv(log(out_cbr_tpm + 1, 2), paste0(dir, "Objects/scenic/cbr_expression_log.csv"))
