
library(monocle3)
library(Seurat)
library(ggplot2)
library(viridis)
library(reshape2)
library(ggrepel)
library(ggExtra) # For ggmarginal
library(dplyr)
library(cowplot) # for plot_grid
library(ggsankey)
library(tidyr)
library(purrr)

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"
source(paste0(dir, "../Scripts/WS290/CalcFunctions.R"))

# monocle objects
cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240903.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

# gene and cell objects
gene_data <- readRDS(paste0(dir, "Objects/gene_data_joint_cell_bg.rds"))
cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))
cell_data_bins <- readRDS(paste0(dir, "Objects/cell_data_joint_cell_bg.rds"))

gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))

TPMListBootstrap_term <- readRDS(paste0(dir, "Objects/TPMTimeBinsListBootstrap_CellCorrection.rds"))
TPMListBootstrapMean_term <- readRDS(paste0(dir, "Objects/TPMListBootstrapMean_CellCorrection.rds"))
BinarizeBootstrapList_term <- readRDS(paste0(dir, "Objects/BinarizeBootstrapTimeBinsList_CellCorrection.rds"))
BinarizeBootstrapListMean_term <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListMean_CellCorrection.rds"))

TPMListBootstrap_pro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))
BinarizeBootstrapList_pro <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListPro_CellCorrection.rds"))

# John handled
# term_markers <- readRDS(paste0(dir, "Objects/TermMarker_list_filt_metadata_cell_bg.rds"))
# pro_markers <- readRDS(paste0(dir, "Objects/ProMarker_list_filt_metadata_cell_bg.rds"))
# joint_markers <- readRDS(paste0(dir, "Objects/joint_markers_cell_bg.rds"))

deg <- readRDS(paste0(dir, "Objects/DEG_df_filt_metadata_cell_bg.rds"))

BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))
ProBarcodeBinsList <- readRDS(paste0(dir, "Objects/ProBarcodeBinsList.rds"))
JointBarcodeBinsList <- readRDS(paste0(dir, "Objects/JointBarcodeBinsList.rds"))

NewBobPath <- "/kimdata/livinlrg/scAnalysis/BobDataComb/"
og_data <- readRDS(paste0(NewBobPath, "Objects/WS290/og_data.rds"))

# cell_data
out_cell_data <- cell_data
out_cell_data$cell_type_bin <- unlist(sapply(out_cell_data$cell_type_bin, function(bins) {
  paste(bins, collapse = ", ")
}))

write.table(out_cell_data,
            paste0(dir, "Tables/cell_data_joint_mean_cell_bg.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# gene_data

# colnames to keep
gene_data_colnames_to_keep <- c("gene",
  "jsd_median_term",
  "jsd_lower_term",
  "jsd_upper_term",
  "pcor_median_term",
  "pcor_lower_term",
  "pcor_upper_term",
  "scor_median_term",
  "scor_lower_term",
  "scor_upper_term",
  "cos_median_term",
  "cos_lower_term",
  "cos_upper_term",
  "cel_tau_median_term",
  "cel_tau_lower_term",
  "cel_tau_upper_term",
  "cbr_tau_median_term",
  "cbr_tau_lower_term",
  "cbr_tau_upper_term",
  "jsd_median_pro",
  "jsd_lower_pro",
  "jsd_upper_pro",
  "pcor_median_pro",
  "pcor_lower_pro",
  "pcor_upper_pro",
  "scor_median_pro",
  "scor_lower_pro",
  "scor_upper_pro",
  "cos_median_pro",
  "cos_lower_pro",
  "cos_upper_pro",
  "cel_tau_median_pro",
  "cel_tau_lower_pro",
  "cel_tau_upper_pro",
  "cbr_tau_median_pro",
  "cbr_tau_lower_pro",
  "cbr_tau_upper_pro",
  "jsd_median_joint",
  "jsd_lower_joint",
  "jsd_upper_joint",
  "pcor_median_joint",
  "pcor_lower_joint",
  "pcor_upper_joint",
  "scor_median_joint",
  "scor_lower_joint",
  "scor_upper_joint",
  "cos_median_joint",
  "cos_lower_joint",
  "cos_upper_joint",
  "cel_tau_median_joint",
  "cel_tau_lower_joint",
  "cel_tau_upper_joint",
  "cbr_tau_median_joint",
  "cbr_tau_lower_joint",
  "cbr_tau_upper_joint",
  "cel_max_tpm_term",
  "cbr_max_tpm_term",
  "max_tpm_term",
  "cel_max_tpm_pro",
  "cbr_max_tpm_pro",
  "max_tpm_pro",
  "elegans_id",
  "elegans_gene_long_name",
  "briggsae_id",
  "briggsae_gene_short_name",
  "briggsae_gene_long_name",
  "orthology_conf",
  "percent_identity",
  "percent_similarity",
  "syntenic",
  "OG",
  "cel_OG_count",
  "cbr_OG_count",
  "WormCat.1",
  "WormCat.2",
  "WormCat.3",
  "omega",
  # "Ma_omega",
  "Cutter_Ka",
  "Cutter_Ks",
  "Cutter_Ka.Ks",
  "PS.Value",
  "PS.Name",
  # "operon",
  "maternal",
  "RecombRegion",
  "RNAi_or_Allele"
  # "embryonic_lethal_direct",
  # "embryonic_lethal_direct_inferred",
  # "larval_lethal_direct",
  # "larval_lethal_direct_inferred",
  # "lethal_direct",
  # "lethal_direct_inferred"
  )

write.table(gene_data[, gene_data_colnames_to_keep],
            paste0(dir, "Tables/gene_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Orthogroup table
out_og_data <- og_data

colnames(out_og_data) <- c("orthogroup", "cel_genes", "cbr_genes",
                           "total_gene_count", "WBGene_names", "transcript_names",
                           "gene_names", "omega", "saturated_branch_count", "tree")

write.table(out_og_data,
            paste0(dir, "Tables/og_data.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

##############
# TPM Tables #
##############

JointTPM <- list()
JointBinarize <- list()
JointTPM[["C.elegans"]] <- cbind(TPMListBootstrap_pro[["C.elegans"]], TPMListBootstrap_term[["C.elegans"]])
JointBinarize[["C.elegans"]] <- cbind(BinarizeBootstrapList_pro[["C.elegans"]], BinarizeBootstrapList_term[["C.elegans"]])

JointTPM[["C.briggsae"]] <- cbind(TPMListBootstrap_pro[["C.briggsae"]], TPMListBootstrap_term[["C.briggsae"]])
JointBinarize[["C.briggsae"]] <- cbind(BinarizeBootstrapList_pro[["C.briggsae"]], BinarizeBootstrapList_term[["C.briggsae"]])

# Now calculate % in for each bin
# Need to use the barcodes and the cds_filt
percent_in_list <- list()
for(cur_species in unique(cds_filt$species)) {
  matrix <- counts(cds_filt[rowData(cds_filt)$gene.type %in% c("common", 
                                                               ifelse(cur_species == "C.elegans",
                                                                      "ce.unique",
                                                                      "cb.unique")),
                            colData(cds_filt)$species == cur_species])
  percent_in_list[[cur_species]] <- list()
  
  for(cell_type in names(JointBarcodeBinsList[[cur_species]])) {
    for(cell_type_bin in names(JointBarcodeBinsList[[cur_species]][[cell_type]])) {
      barcodes <- JointBarcodeBinsList[[cur_species]][[cell_type]][[cell_type_bin]]
      percent_in_list[[cur_species]][[cell_type_bin]] <-sapply(rowSums(matrix[, barcodes] > 0), function(x) x / length(barcodes))
    }
  }
  percent_in_list[[cur_species]] <- do.call(cbind, percent_in_list[[cur_species]])
}

percent_in_list[["C.elegans"]] <- percent_in_list[["C.elegans"]] * 100
percent_in_list[["C.briggsae"]] <- percent_in_list[["C.briggsae"]] * 100

# gene cell_type cell_type_bin median_tpm lower_tpm percent_in
expression_table <- list()
for(cur_species in names(percent_in_list)) {
  for(cell_type in names(JointBarcodeBinsList[[cur_species]])) {
    if(cell_type %in% c("hyp3", "ABalaapppp/ABalapaapp")) {
      next
    }
    for(cell_type_bin in names(JointBarcodeBinsList[[cur_species]][[cell_type]])) {
      expression_table[[cur_species]][[cell_type_bin]] <- data.frame(cell_type = cell_type,
                                                                     cell_type_bin = cell_type_bin,
                                                                     gene = rownames(JointTPM[[cur_species]]),
                                                                     median_tpm = JointTPM[[cur_species]][,cell_type_bin],
                                                                     lower_tpm = JointBinarize[[cur_species]][,cell_type_bin],
                                                                     percent_in = percent_in_list[[cur_species]][,cell_type_bin])
    }
  }
  expression_table[[cur_species]] <- do.call(rbind, expression_table[[cur_species]])
}

write.table(expression_table[["C.elegans"]],
            paste0(dir, "Tables/cel_tpm_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(expression_table[["C.briggsae"]],
            paste0(dir, "Tables/cbr_tpm_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)





