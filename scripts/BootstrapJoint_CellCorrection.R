

# Joining of the terminal and progenitor data into the same object
library(monocle3)
library(Seurat)
library(boot)
library(parallel)
library(dplyr)

set.seed(2022)

# source 
dir <- "/kimdata/livinlrg/scAnalysis/WS290/"
NewBobPath <- "/kimdata/livinlrg/scAnalysis/BobDataComb/"
source(paste0(dir, "../Scripts/CalcFunctions.R"))

# Terminal data
cell_data_term <- readRDS(paste0(dir, "Objects/cell_data_term_cell_bg.rds"))
cell_data_term_mean <- readRDS(paste0(dir, "Objects/cell_data_term_mean_cell_bg.rds"))
gene_data_term <- readRDS(paste0(dir, "Objects/gene_data_term_cell_bg.rds"))
TPMListBootstrap_term <- readRDS(paste0(dir, "Objects/TPMTimeBinsListBootstrap_CellCorrection.rds"))
TPMListBootstrapMean_term <- readRDS(paste0(dir, "Objects/TPMListBootstrapMean_CellCorrection.rds"))

cel_cell_gene_list_term <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cel_cell_gene_list.rds"))
cbr_cell_gene_list_term <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cbr_cell_gene_list.rds"))

# Progenitor Data
cell_data_pro <- readRDS(paste0(dir, "Objects/cell_data_pro_cell_bg.rds"))
gene_data_pro <- readRDS(paste0(dir, "Objects/gene_data_pro_cell_bg.rds"))
TPMListBootstrap_pro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))

cel_cell_gene_list_pro <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cel_cell_gene_list_pro.rds"))
cbr_cell_gene_list_pro <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cbr_cell_gene_list_pro.rds"))

old_cell_data_pro <- readRDS(paste0(NewBobPath, "Objects/cell_data_pro_20231122.rds"))

old_cell_data_pro_subset <- old_cell_data_pro %>%
  select("cell_type", "lineage.group", "div.stage") %>%
  rename("lineage_group" = "lineage.group") %>%
  rename("div_stage" = "div.stage")

old_cell_data_pro_subset$div_stage <- gsub(" celled", "", old_cell_data_pro_subset$div_stage)
old_cell_data_pro_subset$lineage_group <- gsub("\\(", " \\(", old_cell_data_pro_subset$lineage_group)
old_cell_data_pro_subset$lineage_group <- gsub("28_cell_or_earlier", "28 cell or earlier", old_cell_data_pro_subset$lineage_group)

CellTable <- readRDS(paste0(dir, "Objects/CellTable_Names_20240901.rds"))
cell_type_time_bins <- readRDS(paste0(dir, "Objects/cell_type_time_bins.rds"))
CellTableTimeBins <- readRDS(paste0(dir, "Objects/CellTableTimeBins.rds"))

# load barcode bins from EmbryoTimeBins
BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))
ProBarcodeBinsList <- readRDS(paste0(dir, "Objects/ProBarcodeBinsList.rds"))
JointBarcodeBinsList <- readRDS(paste0(dir, "Objects/JointBarcodeBinsList.rds"))

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_no_bg_20240828.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]
rm(cds)

gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))
synteny <- readRDS(paste0(dir, "Objects/synteny_filt.rds"))

#############
# TPM joint #
#############

# Cell_type x gene x species
tpm_matrix_list <- list()
tpm_matrix_t_list <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  temp_gene_tpm <- rbind(t(TPMListBootstrapMean_term[[cur_species]]), t(TPMListBootstrap_pro[[cur_species]]))
  temp_gene_tpm <- data.frame(temp_gene_tpm)
  temp_gene_tpm$cell <- rownames(temp_gene_tpm)
  
  tpm_matrix_list[[cur_species]] <- left_join(CellTable[,c("Lineage", "MergedDatasetName")], temp_gene_tpm, by = c("MergedDatasetName" = "cell"))
  rownames(tpm_matrix_list[[cur_species]]) <- tpm_matrix_list[[cur_species]]$Lineage
  
  tpm_matrix_t_list[[cur_species]] <- t(tpm_matrix_list[[cur_species]][, -c(1, 2)])
  rownames(tpm_matrix_t_list[[cur_species]]) <- rownames(TPMListBootstrapMean_term[[cur_species]])
}
tpm_matrix_list <- tpm_matrix_t_list
saveRDS(tpm_matrix_list, paste0(dir, "Objects/tpm_matrix_list_cell_bg.rds"))

tpm_matrix_list_filt <- list()
tpm_matrix_list_filt[["C.elegans"]] <- tpm_matrix_list[["C.elegans"]][, which(colSums(tpm_matrix_list[["C.elegans"]] > 0) > 0)]
tpm_matrix_list_filt[["C.briggsae"]] <- tpm_matrix_list[["C.briggsae"]][, which(colSums(tpm_matrix_list[["C.briggsae"]] > 0) > 0)]

saveRDS(tpm_matrix_list_filt, paste0(dir, "Objects/tpm_matrix_list_filt_cell_bg.rds"))

# Cell_type x gene x species
tpm_matrix_time_list <- list()
tpm_matrix_time_t_list <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  temp_gene_tpm <- rbind(t(TPMListBootstrap_term[[cur_species]]), t(TPMListBootstrap_pro[[cur_species]]))
  temp_gene_tpm <- data.frame(temp_gene_tpm)
  temp_gene_tpm$cell <- rownames(temp_gene_tpm)
  
  tpm_matrix_time_list[[cur_species]] <- left_join(CellTableTimeBins[,c("Lineage", "TPMName")], temp_gene_tpm, by = c("TPMName" = "cell"))
  rownames(tpm_matrix_time_list[[cur_species]]) <- tpm_matrix_time_list[[cur_species]]$Lineage
  
  tpm_matrix_time_t_list[[cur_species]] <- t(tpm_matrix_time_list[[cur_species]][, -c(1, 2)])
  rownames(tpm_matrix_time_t_list[[cur_species]]) <- rownames(TPMListBootstrap_term[[cur_species]])
}
tpm_matrix_time_list <- tpm_matrix_time_t_list
saveRDS(tpm_matrix_time_list, paste0(dir, "Objects/tpm_matrix_time_list_cell_bg.rds"))

tpm_matrix_time_list_filt <- list()
tpm_matrix_time_list_filt[["C.elegans"]] <- tpm_matrix_time_list[["C.elegans"]][, which(colSums(tpm_matrix_time_list[["C.elegans"]] > 0) > 0)]
tpm_matrix_time_list_filt[["C.briggsae"]] <- tpm_matrix_time_list[["C.briggsae"]][, which(colSums(tpm_matrix_time_list[["C.briggsae"]] > 0) > 0)]

saveRDS(tpm_matrix_time_list_filt, paste0(dir, "Objects/tpm_matrix_time_list_filt_cell_bg.rds"))

##########
# Gene data
##########

gene_data <- left_join(gene_data_term, gene_data_pro, by = "gene", suffix = c("_term", "_pro"))

# per gene matrix
distance_list <- mclapply(names(cel_cell_gene_list_pro), function(this.gene) {
  print(paste0(this.gene, ": ", grep(this.gene, names(cel_cell_gene_list_pro)), " of ", length(names(cel_cell_gene_list_pro))))
  cel_tpm_matrix <- cbind(cel_cell_gene_list_pro[[this.gene]], cel_cell_gene_list_term[[this.gene]])
  cbr_tpm_matrix <- cbind(cbr_cell_gene_list_pro[[this.gene]], cbr_cell_gene_list_term[[this.gene]])
  
  cel_temp_gene_tpm <- data.frame(t(cel_tpm_matrix))
  cel_temp_gene_tpm$cell <- rownames(cel_temp_gene_tpm)
  
  cbr_temp_gene_tpm <- data.frame(t(cbr_tpm_matrix))
  cbr_temp_gene_tpm$cell <- rownames(cbr_temp_gene_tpm)
  
  cel_tpm_matrix <- left_join(CellTableTimeBins[,c("Lineage", "TPMName")], cel_temp_gene_tpm, by = c("TPMName" = "cell"))
  rownames(cel_tpm_matrix) <- cel_tpm_matrix$Lineage
  
  cbr_tpm_matrix <- left_join(CellTableTimeBins[,c("Lineage", "TPMName")], cbr_temp_gene_tpm, by = c("TPMName" = "cell"))
  rownames(cbr_tpm_matrix) <- cbr_tpm_matrix$Lineage
  
  tpm_matrix_list <- list()
  tpm_matrix_list[["C.elegans"]] <- t(cel_tpm_matrix[, -c(1, 2)])
  tpm_matrix_list[["C.briggsae"]] <- t(cbr_tpm_matrix[, -c(1, 2)])
  
  pseudocount <- 1/sum(! is.na(tpm_matrix_list[["C.elegans"]][1,]))
  
  jsd_df <- data.frame(i = seq(1, 1000), jsd = NA, gene = this.gene)
  pcor_df <- data.frame(i = seq(1, 1000), pcor = NA, gene = this.gene)
  scor_df <- data.frame(i = seq(1, 1000), scor = NA, gene = this.gene)
  cos_df <- data.frame(i = seq(1, 1000), cos = NA, gene = this.gene)
  cel_tau_df <- data.frame(i = seq(1, 1000), cel_tau = NA, gene = this.gene)
  cbr_tau_df <- data.frame(i = seq(1, 1000), cbr_tau = NA, gene = this.gene)
  
  for(i in seq(1, 1000)) {
    p = tpm_matrix_list[["C.elegans"]][i,] + pseudocount
    q = tpm_matrix_list[["C.briggsae"]][i,] + pseudocount
    
    # Filter NA's out
    p_filt = p[! (is.na(p) | is.na(q))]
    q_filt = q[! (is.na(p) | is.na(q))]
    
    jsd_df[i,"jsd"] <- sqrt(js_divg(p = p_filt/sum(p_filt), q = q_filt/sum(q_filt)))
    cos_df[i,"cos"] <- cos_dist(p_filt, q_filt)
    
    p = log2(tpm_matrix_list[["C.elegans"]][i,] + 1)
    q = log2(tpm_matrix_list[["C.briggsae"]][i,] + 1)
    
    # Filter NA's outq
    p_filt = p[! (is.na(p) | is.na(q))]
    q_filt = q[! (is.na(p) | is.na(q))]
    
    pcor_df[i, "pcor"] <- cor(p_filt, q_filt)
    scor_df[i, "scor"] <- cor(p_filt, q_filt, method = "spearman")
    
    # calculate tau in matrix format
    #Tau = Sum(1-expNormalizedByMax)/1-numberOfTissues
    cel_tau_df[i, "cel_tau"] <- sum(1 - (p_filt/max(p_filt)))/(length(p_filt) - 1)
    cbr_tau_df[i, "cbr_tau"] <- sum(1 - (q_filt/max(q_filt)))/(length(q_filt) - 1)
  }
  return(left_join(left_join(left_join(left_join(left_join(jsd_df, pcor_df, by = c("i" = "i", "gene" = "gene")),
                                                 scor_df, by = c("i" = "i", "gene" = "gene")),
                                       cos_df, by = c("i" = "i", "gene" = "gene")),
                             cel_tau_df, by = c("i" = "i", "gene" = "gene")),
                   cbr_tau_df, by = c("i" = "i", "gene" = "gene")))
}, mc.cores = 20)

distance_df <- do.call(rbind, distance_list)

gene_data <- left_join(gene_data, distance_df %>%
            group_by(gene) %>%
            summarise(jsd_median_joint = median(jsd, na.rm = TRUE),
                      jsd_lower_joint = quantile(jsd, na.rm = TRUE, probs = c(0.025)), 
                      jsd_upper_joint = quantile(jsd, na.rm = TRUE, probs = c(0.975)),
                      pcor_median_joint = median(pcor, na.rm = TRUE),
                      pcor_lower_joint = quantile(pcor, na.rm = TRUE, probs = c(0.025)), 
                      pcor_upper_joint = quantile(pcor, na.rm = TRUE, probs = c(0.975)),
                      scor_median_joint = median(scor, na.rm = TRUE),
                      scor_lower_joint = quantile(scor, na.rm = TRUE, probs = c(0.025)), 
                      scor_upper_joint = quantile(scor, na.rm = TRUE, probs = c(0.975)),
                      cos_median_joint = median(cos, na.rm = TRUE),
                      cos_lower_joint = quantile(cos, na.rm = TRUE, probs = c(0.025)), 
                      cos_upper_joint = quantile(cos, na.rm = TRUE, probs = c(0.975)),
                      cel_tau_median_joint = median(cel_tau, na.rm = TRUE),
                      cel_tau_lower_joint = quantile(cel_tau, na.rm = TRUE, probs = c(0.025)),
                      cel_tau_upper_joint = quantile(cel_tau, na.rm = TRUE, probs = c(0.975)),
                      cbr_tau_median_joint = median(cbr_tau, na.rm = TRUE),
                      cbr_tau_lower_joint = quantile(cbr_tau, na.rm = TRUE, probs = c(0.025)),
                      cbr_tau_upper_joint = quantile(cbr_tau, na.rm = TRUE, probs = c(0.975))),
            by = "gene",
            suffix = c("", ""))

rownames(gene_data) <- gene_data$gene
saveRDS(gene_data, paste0(dir, "Objects/gene_data_joint_cell_bg.rds"))

# Max TPM elegans terminal
gene_data$cel_max_tpm_term <- apply(TPMListBootstrap_term[["C.elegans"]][gene_data$gene,], 1, max)
# Max TPM briggsae terminal
gene_data$cbr_max_tpm_term <- apply(TPMListBootstrap_term[["C.briggsae"]][gene_data$gene,], 1, max)
# Max TPM terminal
gene_data$max_tpm_term <- apply(gene_data[, c("cel_max_tpm_term", "cbr_max_tpm_term")], 1, max)

# Max TPM elegans progenitor
gene_data$cel_max_tpm_pro <- apply(TPMListBootstrap_pro[["C.elegans"]][gene_data$gene,], 1, max)
# Max TPM briggsae progenitor
gene_data$cbr_max_tpm_pro <- apply(TPMListBootstrap_pro[["C.briggsae"]][gene_data$gene,], 1, max)
# Max TPM progenitor
gene_data$max_tpm_pro <- apply(gene_data[, c("cel_max_tpm_pro", "cbr_max_tpm_pro")], 1, max)

# add some mutant metadata to the gff_list
gff_list_mrna[["elegans"]]$RNAi_or_Allele <- NA
gff_list_mrna[["elegans"]][which(gff_list_mrna[["elegans"]]$RNAi_Phenotype == "N.A." & gff_list_mrna[["elegans"]]$Allele_Phenotype == "N.A."),]$RNAi_or_Allele <- FALSE # None annotated
gff_list_mrna[["elegans"]][which(gff_list_mrna[["elegans"]]$RNAi_Phenotype != "N.A." | gff_list_mrna[["elegans"]]$Allele_Phenotype != "N.A."),]$RNAi_or_Allele <- TRUE # Has/Doesn't have

gene_data <- left_join(gene_data, gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene.type == "common",
                                                             c("cds_gene_name", "elegans_id", "elegans_gene_long_name",
                                                               "briggsae_id", "briggsae_gene_short_name", "briggsae_gene_long_name",
                                                               "orthology_conf", "percent_identity", "percent_similarity", "syntenic",
                                                               "OG", "cel_OG_count", "cbr_OG_count", "WormCat.1", "WormCat.2", "WormCat.3",
                                                               "omega", "Ma_omega", "Cutter_Ka", "Cutter_Ks", "Cutter_Ka.Ks", 
                                                               "PS.Value", "PS.Name",
                                                               "operon",  "maternal", "RecombRegion",
                                                               "RNAi_or_Allele", "embryonic_lethal_direct", "embryonic_lethal_direct_inferred",
                                                               "larval_lethal_direct", "larval_lethal_direct_inferred", "lethal_direct", "lethal_direct_inferred")],
                       by = c("gene" = "cds_gene_name"))

rownames(gene_data) <- gene_data$gene
saveRDS(gene_data, paste0(dir, "Objects/gene_data_joint_cell_bg.rds"))

##########
# Cell data
##########

# Time binned
cell_data_joint <- rbind(cell_data_term, cell_data_pro)

# Time mean'd 
cell_data_joint_mean <- rbind(cell_data_term_mean, cell_data_pro)

saveRDS(cell_data_joint, paste0(dir, "Objects/cell_data_joint_cell_bg.rds"))

term_markers <- readRDS(paste0(dir, "Objects/TermMarker_list_metadata_cell_bg.rds"))
term_markers_filt <- readRDS(paste0(dir, "Objects/TermMarker_list_filt_metadata_cell_bg.rds"))

pro_markers <- readRDS(paste0(dir, "Objects/ProMarker_list_metadata_cell_bg.rds"))
pro_markers_filt <- readRDS(paste0(dir, "Objects/ProMarker_list_filt_metadata_cell_bg.rds"))

deg <- readRDS(paste0(dir, "Objects/DEG_df_metadata_cell_bg.rds"))
deg_filt <- readRDS(paste0(dir, "Objects/DEG_df_filt_metadata_cell_bg.rds"))

BinarizeBootstrapList <- readRDS(paste0(dir, "Objects/BinarizeBootstrapTimeBinsList_CellCorrection.rds"))
BinarizeBootstrapListPro <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListPro_CellCorrection.rds"))

# cel_markers	Number of elegans markers
cel_markers_total <- data.frame(term_markers_filt[["C.elegans"]] %>%
                                  group_by(cell_type) %>%
                                  summarise(count = n()))

cel_markers_total <- left_join(cel_markers_total, cell_data_joint[,c("cell_type", "cell_type_bin")])
pro_markers_total <- data.frame(pro_markers_filt[["C.elegans"]] %>%
                                  group_by(cell_type) %>%
                                  summarise(count = n()))
pro_markers_total$cell_type_bin <- pro_markers_total$cell_type

cel_markers_total <- rbind(cel_markers_total,
                           pro_markers_total)

cell_data_joint$cel_markers <- NA
cell_data_joint[cel_markers_total$cell_type_bin,]$cel_markers <- cel_markers_total$count

# cbr_markers	Number of briggsae markers
cbr_markers_total <- data.frame(term_markers_filt[["C.briggsae"]] %>%
                                  group_by(cell_type) %>%
                                  summarise(count = n()))

cbr_markers_total <- left_join(cbr_markers_total, cell_data_joint[,c("cell_type", "cell_type_bin")])
pro_markers_total <- data.frame(pro_markers_filt[["C.briggsae"]] %>%
                                  group_by(cell_type) %>%
                                  summarise(count = n()))
pro_markers_total$cell_type_bin <- pro_markers_total$cell_type

cbr_markers_total <- rbind(cbr_markers_total,
                           pro_markers_total)

cell_data_joint$cbr_markers <- NA
cell_data_joint[cbr_markers_total$cell_type_bin,]$cbr_markers <- cbr_markers_total$count

# shared_markers	Number of shared markers
# either_markers	Union of number of markers

# cel_markers_common	Number of markers in 1:1 orthology set
cel_markers_total <- data.frame(term_markers_filt[["C.elegans"]] %>%
                                  filter(gene.type %in% "common") %>%
                                  group_by(cell_type) %>%
                                  summarise(count = n()))

cel_markers_total <- left_join(cel_markers_total, cell_data_joint[,c("cell_type", "cell_type_bin")])
pro_markers_total <- data.frame(pro_markers_filt[["C.elegans"]] %>%
                                  filter(gene.type %in% "common") %>%
                                  group_by(cell_type) %>%
                                  summarise(count = n()))
pro_markers_total$cell_type_bin <- pro_markers_total$cell_type

cel_markers_total <- rbind(cel_markers_total,
                           pro_markers_total)

cell_data_joint$cel_markers_common <- NA
cell_data_joint[cel_markers_total$cell_type_bin,]$cel_markers_common <- cel_markers_total$count

# cbr_markers_common	Number of markers in 1:1 orthology set
cbr_markers_total <- data.frame(term_markers_filt[["C.briggsae"]] %>%
                                  group_by(cell_type) %>%
                                  filter(gene.type %in% "common") %>%
                                  summarise(count = n()))

cbr_markers_total <- left_join(cbr_markers_total, cell_data_joint[,c("cell_type", "cell_type_bin")])
pro_markers_total <- data.frame(pro_markers_filt[["C.briggsae"]] %>%
                                  filter(gene.type %in% "common") %>%
                                  group_by(cell_type) %>%
                                  summarise(count = n()))
pro_markers_total$cell_type_bin <- pro_markers_total$cell_type

cbr_markers_total <- rbind(cbr_markers_total,
                           pro_markers_total)

cell_data_joint$cbr_markers_common <- NA
cell_data_joint[cbr_markers_total$cell_type_bin,]$cbr_markers_common <- cbr_markers_total$count

# cel_markers_private	Number of markers that are private to elegans
cel_markers_total <- data.frame(term_markers_filt[["C.elegans"]] %>%
                                  filter(! gene.type %in% "common") %>%
                                  group_by(cell_type) %>%
                                  summarise(count = n()))

cel_markers_total <- left_join(cel_markers_total, cell_data_joint[,c("cell_type", "cell_type_bin")])
pro_markers_total <- data.frame(pro_markers_filt[["C.elegans"]] %>%
                                  filter(! gene.type %in% "common") %>%
                                  group_by(cell_type) %>%
                                  summarise(count = n()))
pro_markers_total$cell_type_bin <- pro_markers_total$cell_type

cel_markers_total <- rbind(cel_markers_total,
                           pro_markers_total)

cell_data_joint$cel_markers_non_one_to_one <- NA
cell_data_joint[cel_markers_total$cell_type_bin,]$cel_markers_non_one_to_one <- cel_markers_total$count

# cbr_markers_private	Number of markers that are private to briggsae
cbr_markers_total <- data.frame(term_markers_filt[["C.briggsae"]] %>%
                                  filter(! gene.type %in% "common") %>%
                                  group_by(cell_type) %>%
                                  summarise(count = n()))

cbr_markers_total <- left_join(cbr_markers_total, cell_data_joint[,c("cell_type", "cell_type_bin")])
pro_markers_total <- data.frame(pro_markers_filt[["C.briggsae"]] %>%
                                  filter(! gene.type %in% "common") %>%
                                  group_by(cell_type) %>%
                                  summarise(count = n()))
pro_markers_total$cell_type_bin <- pro_markers_total$cell_type

cbr_markers_total <- rbind(cbr_markers_total,
                           pro_markers_total)

cell_data_joint$cbr_markers_non_one_to_one <- NA
cell_data_joint[cbr_markers_total$cell_type_bin,]$cbr_markers_non_one_to_one <- cbr_markers_total$count

# Genes detected
genes_detected_term_cel <- colSums(BinarizeBootstrapList[["C.elegans"]] > 0)
genes_detected_pro_cel <- colSums(BinarizeBootstrapListPro[["C.elegans"]] > 0)

genes_detected_term_cbr <- colSums(BinarizeBootstrapList[["C.briggsae"]] > 0)
genes_detected_pro_cbr <- colSums(BinarizeBootstrapListPro[["C.briggsae"]] > 0)

genes_detected_cel <- c(genes_detected_term_cel, genes_detected_pro_cel)
genes_detected_cbr <- c(genes_detected_term_cbr, genes_detected_pro_cbr)
cell_data_joint <- left_join(cell_data_joint, data.frame(genes_detected_bootstrap_cel = genes_detected_cel,
                                                         genes_detected_bootstrap_cbr = genes_detected_cbr,
                                                         cell_type_bin = names(genes_detected_cel)))

# Cell count
rownames(cell_data_joint) <- cell_data_joint$cell_type_bin
cell_data_joint$cel_cell_count <- NA
cell_data_joint$cbr_cell_count <- NA
for(cell_type_bin in cell_data_joint$cell_type_bin) {
  cur_cell_type <- unique(cell_data_joint[cell_type_bin,]$cell_type)
  cell_data_joint[cell_type_bin,]$cel_cell_count <- length(JointBarcodeBinsList[["C.elegans"]][[cur_cell_type]][[cell_type_bin]])
  cell_data_joint[cell_type_bin,]$cbr_cell_count <- length(JointBarcodeBinsList[["C.briggsae"]][[cur_cell_type]][[cell_type_bin]])
}
cell_data_joint$min_cell_count <- apply(cell_data_joint[,c("cel_cell_count", "cbr_cell_count")], 1, function(x) min(x))

# Mean embryo time
rownames(cell_data_joint) <- cell_data_joint$cell_type_bin
cell_data_joint$embryo_time <- NA
for(cell_type_bin in cell_data_joint$cell_type_bin) {
  cur_cell_type <- unique(cell_data_joint[cell_type_bin,]$cell_type)
  cell_data_joint[cell_type_bin,]$embryo_time <- median(colData(cds_filt)[JointBarcodeBinsList[["C.elegans"]][[cur_cell_type]][[cell_type_bin]],]$smoothed.embryo.time)
}

# Mean number of UMI
cell_data_joint$cel_median_umi <- NA
cell_data_joint$cbr_median_umi <- NA
for(cell_type_bin in cell_data_joint$cell_type_bin) {
  cur_cell_type <- unique(cell_data_joint[cell_type_bin,]$cell_type)
  cell_data_joint[cell_type_bin,]$cel_median_umi <- median(colData(cds_filt)[JointBarcodeBinsList[["C.elegans"]][[cur_cell_type]][[cell_type_bin]],]$n.umi)
  cell_data_joint[cell_type_bin,]$cbr_median_umi <- median(colData(cds_filt)[JointBarcodeBinsList[["C.briggsae"]][[cur_cell_type]][[cell_type_bin]],]$n.umi)
}

# Between species DEG
deg_count <- deg_filt %>% group_by(cell_type_bin) %>% distinct(gene) %>% summarise(deg = n())

cell_data_joint <- left_join(cell_data_joint, deg_count, by = c("cell_type_bin" = "cell_type_bin"))

cell_data_joint <- left_join(cell_data_joint, old_cell_data_pro_subset)
cell_data_joint[cell_data_joint$cell_class == "progenitor" & is.na(cell_data_joint$lineage_group),]$lineage_group <- c("Other AB", "Other AB", "28 cell or earlier")
cell_data_joint[cell_data_joint$cell_class == "progenitor" & is.na(cell_data_joint$div_stage),]$div_stage <- c("350", "350", "15")

cell_data_joint[is.na(cell_data_joint$div_stage),]$div_stage <- "600"
cell_data_joint$div_stage <- factor(cell_data_joint$div_stage, levels = c("15", "28", "50", "100", "200", "350", "600"))

rownames(cell_data_joint) <- cell_data_joint$cell_type_bin
saveRDS(cell_data_joint, paste0(dir, "Objects/cell_data_joint_cell_bg.rds"))

###############
# Time meaned #
###############

cell_data_joint_mean <- rbind(cell_data_term_mean, cell_data_pro)

# [1] "cel_markers"                "cbr_markers"               
# [3] "cel_markers_common"         "cbr_markers_common"        
# [5] "cel_markers_non_one_to_one" "cbr_markers_non_one_to_one"
cell_data_joint_mean <- left_join(cell_data_joint_mean, cell_data_joint[! duplicated(cell_data_joint$cell_type),
                                                                        c("cell_type", "cel_markers", "cbr_markers",
                                                                          "cel_markers_common", "cbr_markers_common",
                                                                          "cel_markers_non_one_to_one", "cbr_markers_non_one_to_one")])

# [7] "genes_detected_bootstrap"
cell_data_joint_mean$genes_detected_bootstrap_cel <- NA
cell_data_joint_mean$genes_detected_bootstrap_cbr <- NA
rownames(cell_data_joint_mean) <- cell_data_joint_mean$cell_type
for(cur_cell_type in cell_data_joint_mean$cell_type) {
  cell_data_joint_mean[cur_cell_type,]$genes_detected_bootstrap_cel <- sum(rowSums(cbind(BinarizeBootstrapList[["C.elegans"]], BinarizeBootstrapListPro[["C.elegans"]])[,unlist(cell_data_joint_mean[cur_cell_type,]$cell_type_bin), drop = FALSE]) > 0)
  cell_data_joint_mean[cur_cell_type,]$genes_detected_bootstrap_cbr <- sum(rowSums(cbind(BinarizeBootstrapList[["C.briggsae"]], BinarizeBootstrapListPro[["C.briggsae"]])[,unlist(cell_data_joint_mean[cur_cell_type,]$cell_type_bin), drop = FALSE]) > 0)
}

# "cel_cell_count" "cbr_cell_count"
rownames(cell_data_joint_mean) <- cell_data_joint_mean$cell_type
cell_data_joint_mean$cel_cell_count <- NA
cell_data_joint_mean$cbr_cell_count <- NA
for(cur_cell_type in cell_data_joint_mean$cell_type) {
  cell_data_joint_mean[cur_cell_type,]$cel_cell_count <- length(unlist(JointBarcodeBinsList[["C.elegans"]][[cur_cell_type]]))
  cell_data_joint_mean[cur_cell_type,]$cbr_cell_count <- length(unlist(JointBarcodeBinsList[["C.briggsae"]][[cur_cell_type]]))
}

cell_data_joint_mean$min_cell_count <- pmin(cell_data_joint_mean$cel_cell_count, cell_data_joint_mean$cbr_cell_count)

# "embryo_time"
cell_data_joint_mean$embryo_time <- NA
for(cur_cell_type in cell_data_joint_mean$cell_type) {
  cell_data_joint_mean[cur_cell_type,]$embryo_time <- median(colData(cds_filt)[unlist(JointBarcodeBinsList[["C.elegans"]][[cur_cell_type]]),]$smoothed.embryo.time)
}

# "cel_median_umi" "cbr_median_umi"
cell_data_joint_mean$cel_median_umi <- NA
cell_data_joint_mean$cbr_median_umi <- NA
for(cur_cell_type in cell_data_joint_mean$cell_type) {
  cell_data_joint_mean[cur_cell_type,]$cel_median_umi <- median(colData(cds_filt)[unlist(JointBarcodeBinsList[["C.elegans"]][[cur_cell_type]]),]$n.umi)
  cell_data_joint_mean[cur_cell_type,]$cbr_median_umi <- median(colData(cds_filt)[unlist(JointBarcodeBinsList[["C.briggsae"]][[cur_cell_type]]),]$n.umi)
}

# "deg"
cell_data_joint_mean <- left_join(cell_data_joint_mean,
                                  deg_filt %>% group_by(cell_type) %>% distinct(gene) %>% summarise(deg = n()),
                                  by = c("cell_type" = "cell_type"))

cell_data_joint_mean$deg <- NULL
for(cur_cell_type in cell_data_joint_mean$cell_type) {
  num_unique_deg <- length(unique(deg_filt[deg_filt %in% cur_cell_type,]$gene))
  cell_data_joint_mean[cell_data_joint_mean$cell_type == cur_cell_type,]$deg <- num_unique_deg
}

# cell_data_joint_mean <- left_join(cell_data_joint_mean, deg_filt, by = c("cell_type" = "cell_type"))

cell_data_joint_mean <- left_join(cell_data_joint_mean, old_cell_data_pro_subset)
cell_data_joint_mean[cell_data_joint_mean$cell_class == "progenitor" & is.na(cell_data_joint_mean$lineage_group),]$lineage_group <- c("Other AB", "Other AB", "28 cell or earlier")
cell_data_joint_mean[cell_data_joint_mean$cell_class == "progenitor" & is.na(cell_data_joint_mean$div_stage),]$div_stage <- c("350", "350", "15")

cell_data_joint_mean[is.na(cell_data_joint_mean$div_stage),]$div_stage <- "600"
cell_data_joint_mean$div_stage <- factor(cell_data_joint_mean$div_stage, levels = c("15", "28", "50", "100", "200", "350", "600"))

rownames(cell_data_joint_mean) <- cell_data_joint_mean$cell_type
saveRDS(cell_data_joint_mean, paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))
cell_data_joint_mean <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))

#########################################
# Remove hyp3 and ABalaapppp/ABalapaapp #
#########################################

cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_cell_bg.rds"))
cell_data <- cell_data[cell_data$cell_type != "hyp3" & cell_data$cell_type != "ABalaapppp/ABalapaapp",]
cell_data["T_270_330:330_390","cell_class"] <- "Hypodermis and seam"
saveRDS(cell_data, paste0(dir, "Objects/cell_data_joint_cell_bg.rds"))

cell_data_mean <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))
cell_data_mean <- cell_data_mean[cell_data_mean$cell_type != "hyp3" & cell_data_mean$cell_type != "ABalaapppp/ABalapaapp",]
cell_data_mean["T","cell_class"] <- "Hypodermis and seam"
saveRDS(cell_data_mean, paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))

term_markers <- readRDS(paste0(dir, "Objects/TermMarker_list_filt_metadata_cell_bg.rds"))
term_markers <- lapply(term_markers, function(x) {
  return(x[! x$cell_type %in% c("hyp3", "ABalaapppp/ABalapaapp"),])
})
saveRDS(term_markers, paste0(dir, "Objects/TermMarker_list_filt_metadata_cell_bg.rds"))

pro_markers <- readRDS(paste0(dir, "Objects/ProMarker_list_filt_metadata_cell_bg.rds"))
pro_markers <- lapply(pro_markers, function(x) {
  return(x[! x$cell_type %in% c("hyp3", "ABalaapppp/ABalapaapp"),])
})
saveRDS(pro_markers, paste0(dir, "Objects/ProMarker_list_filt_metadata_cell_bg.rds"))

TPMListBootstrap_term <- readRDS(paste0(dir, "Objects/TPMTimeBinsListBootstrap_CellCorrection.rds"))
TPMListBootstrap_term <- lapply(TPMListBootstrap_term, function(x) {
  return(x[,! colnames(x) %in% c("hyp3_390_450:580_650")])
})
saveRDS(TPMListBootstrap_term, paste0(dir, "Objects/TPMTimeBinsListBootstrap_CellCorrection.rds"))

TPMListBootstrapMean_term <- readRDS(paste0(dir, "Objects/TPMListBootstrapMean_CellCorrection.rds"))
TPMListBootstrapMean_term <- lapply(TPMListBootstrapMean_term, function(x) {
  return(x[,! colnames(x) %in% c("hyp3")])
})
saveRDS(TPMListBootstrapMean_term, paste0(dir, "Objects/TPMListBootstrapMean_CellCorrection.rds"))

BinarizeBootstrapList_term <- readRDS(paste0(dir, "Objects/BinarizeBootstrapList_CellCorrection.rds"))
BinarizeBootstrapList_term <- lapply(BinarizeBootstrapList_term, function(x) {
  return(x[,! colnames(x) %in% c("hyp3_390_450:580_650")])
})
saveRDS(BinarizeBootstrapList_term, paste0(dir, "Objects/BinarizeBootstrapList_CellCorrection.rds"))

BinarizeBootstrapListMean_term <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListMean_CellCorrection.rds"))
BinarizeBootstrapListMean_term <- lapply(BinarizeBootstrapListMean_term, function(x) {
  return(x[,! colnames(x) %in% c("hyp3")])
})
saveRDS(BinarizeBootstrapListMean_term, paste0(dir, "Objects/BinarizeBootstrapListMean_CellCorrection.rds"))

TPMListBootstrap_pro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))
TPMListBootstrap_pro <- lapply(TPMListBootstrap_pro, function(x) {
  return(x[,! colnames(x) %in% c("ABalaapppp/ABalapaapp")])
})
saveRDS(TPMListBootstrap_pro, paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))

BinarizeBootstrapListPro <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListPro_CellCorrection.rds"))
BinarizeBootstrapListPro <- lapply(BinarizeBootstrapListPro, function(x) {
  return(x[,! colnames(x) %in% c("ABalaapppp/ABalapaapp")])
})
saveRDS(BinarizeBootstrapListPro, paste0(dir, "Objects/BinarizeBootstrapListPro_CellCorrection.rds"))

BootstrapListUpper <- readRDS(paste0(dir, "Objects/BinarizeBootstrapTimeBinsListUpper_CellCorrection.rds"))
BootstrapListProUpper <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListUpperPro_CellCorrection.rds"))

BootstrapListUpperJoint <- list()
BootstrapListUpperJoint[["C.elegans"]] <- cbind(BootstrapListUpper[["C.elegans"]], BootstrapListProUpper[["C.elegans"]])
BootstrapListUpperJoint[["C.briggsae"]] <- cbind(BootstrapListUpper[["C.briggsae"]], BootstrapListProUpper[["C.briggsae"]])

BootstrapListUpperJoint <- lapply(BootstrapListUpperJoint, function(x) {
  return(x[,! colnames(x) %in% c("hyp3_390_450:580_650", "ABalaapppp/ABalapaapp")])
})

saveRDS(BootstrapListUpperJoint, paste0(dir, "Objects/BootstrapListUpperJoint_CellCorrection.rds"))

######################
# Bootstrap matrices #
######################

BootstrapListUpperJoint <- readRDS(paste0(dir, "Objects/BootstrapListUpperJoint_CellCorrection.rds"))

BinarizeBootstrapListPro <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListPro_CellCorrection.rds"))
BinarizeBootstrapList_term <- readRDS(paste0(dir, "Objects/BinarizeBootstrapList_CellCorrection.rds"))

BootstrapListLowerJoint <- list()
BootstrapListLowerJoint[["C.elegans"]] <- cbind(BinarizeBootstrapList_term[["C.elegans"]], BinarizeBootstrapListPro[["C.elegans"]])
BootstrapListLowerJoint[["C.briggsae"]] <- cbind(BinarizeBootstrapList_term[["C.briggsae"]], BinarizeBootstrapListPro[["C.briggsae"]])

# Cell_type x gene x species
tpm_matrix_time_list <- list()
tpm_matrix_time_t_list <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  temp_gene_tpm <- t(BootstrapListLowerJoint[[cur_species]])
  cell_names <- rownames(temp_gene_tpm)
  temp_gene_tpm <- data.frame(temp_gene_tpm)
  temp_gene_tpm$cell <- cell_names
  
  tpm_matrix_time_list[[cur_species]] <- left_join(CellTableTimeBins[,c("Lineage", "TPMName")], temp_gene_tpm, by = c("TPMName" = "cell"))
  rownames(tpm_matrix_time_list[[cur_species]]) <- tpm_matrix_time_list[[cur_species]]$Lineage
  
  tpm_matrix_time_t_list[[cur_species]] <- t(tpm_matrix_time_list[[cur_species]][, -c(1, 2)])
  rownames(tpm_matrix_time_t_list[[cur_species]]) <- rownames(BootstrapListLowerJoint[[cur_species]])
}
tpm_matrix_time_list <- tpm_matrix_time_t_list

saveRDS(tpm_matrix_time_list, paste0(dir, "Objects/tpm_matrix_lower_ci_time_list_cell_bg.rds"))

tpm_matrix_time_list_filt <- list()
tpm_matrix_time_list_filt[["C.elegans"]] <- tpm_matrix_time_list[["C.elegans"]][, which(colSums(tpm_matrix_time_list[["C.elegans"]] > 0) > 0)]
tpm_matrix_time_list_filt[["C.briggsae"]] <- tpm_matrix_time_list[["C.briggsae"]][, which(colSums(tpm_matrix_time_list[["C.briggsae"]] > 0) > 0)]

saveRDS(tpm_matrix_time_list_filt, paste0(dir, "Objects/tpm_matrix_lower_ci_time_list_filt_cell_bg.rds"))

# Cell_type x gene x species
tpm_matrix_time_list <- list()
tpm_matrix_time_t_list <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  temp_gene_tpm <- t(BootstrapListUpperJoint[[cur_species]])
  cell_names <- rownames(temp_gene_tpm)
  temp_gene_tpm <- data.frame(temp_gene_tpm)
  temp_gene_tpm$cell <- cell_names
  
  tpm_matrix_time_list[[cur_species]] <- left_join(CellTableTimeBins[,c("Lineage", "TPMName")], temp_gene_tpm, by = c("TPMName" = "cell"))
  rownames(tpm_matrix_time_list[[cur_species]]) <- tpm_matrix_time_list[[cur_species]]$Lineage
  
  tpm_matrix_time_t_list[[cur_species]] <- t(tpm_matrix_time_list[[cur_species]][, -c(1, 2)])
  rownames(tpm_matrix_time_t_list[[cur_species]]) <- rownames(BootstrapListUpperJoint[[cur_species]])
}
tpm_matrix_time_list <- tpm_matrix_time_t_list

saveRDS(tpm_matrix_time_list, paste0(dir, "Objects/tpm_matrix_upper_ci_time_list_cell_bg.rds"))

tpm_matrix_time_list_filt <- list()
tpm_matrix_time_list_filt[["C.elegans"]] <- tpm_matrix_time_list[["C.elegans"]][, which(colSums(tpm_matrix_time_list[["C.elegans"]] > 0) > 0)]
tpm_matrix_time_list_filt[["C.briggsae"]] <- tpm_matrix_time_list[["C.briggsae"]][, which(colSums(tpm_matrix_time_list[["C.briggsae"]] > 0) > 0)]

saveRDS(tpm_matrix_time_list_filt, paste0(dir, "Objects/tpm_matrix_upper_ci_time_list_filt_cell_bg.rds"))

tpm_matrix_lower <- readRDS(paste0(dir, "Objects/tpm_matrix_lower_ci_time_list_filt_cell_bg.rds"))
tpm_matrix_upper <- readRDS(paste0(dir, "Objects/tpm_matrix_upper_ci_time_list_filt_cell_bg.rds"))