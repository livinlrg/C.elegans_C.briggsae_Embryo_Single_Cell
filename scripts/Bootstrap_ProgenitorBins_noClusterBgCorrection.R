# for each set of barcodes coming from the time bin barcode list, I will take a sample of it with replacement.
# inside of this set of barcodes, the percentage of the background will be estimated and multiplied by those
# experiments background tpm.

library(monocle3)
library(Seurat)
library(boot)
library(parallel)
library(dplyr)

set.seed(2022)

# source 
dir <- "/kimdata/livinlrg/scAnalysis/WS290/"
source(paste0(dir, "../Scripts/CalcFunctions.R"))
cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240903.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

rm(cds)

# load barcode bins from EmbryoTimeBins
BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))
ProBarcodeBinsList <- readRDS(paste0(dir, "Objects/ProBarcodeBinsList.rds"))
JointBarcodeBinsList <- readRDS(paste0(dir, "Objects/JointBarcodeBinsList.rds"))

dir.create(paste0(dir, "Objects/BootObjects_CellCorrection/"))

##########################################################################
# Load a bunch of functions
##########################################################################

# Normalize expression matrix by Size factor
get.norm.expr.matrix <- function(cds) {
  mat = counts(cds)
  mat@x = mat@x / rep.int(pData(cds)$Size_Factor, diff(mat@p))
  return(mat)
}

## Simple bootstrap TPM calculations w/o background subtraction
BinTPM <- function(matrix, i) {
  matrix_sub <- matrix[i,]
  matrix_colMean <- Matrix::colMeans(matrix_sub)
  colSum <- sum(matrix_colMean)
  return(matrix_colMean/colSum * 1000000)
}

##########################################################################
# Bootstrap TPM for Progenitor Cell Types
##########################################################################
bootstrapTPMProList <- list()
column = "cell_type"
for(cur_species in c("C.elegans", "C.briggsae")) {
  print(cur_species)
  cur_cds <- cds_filt[rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                            "ce.unique",
                                                            "cb.unique") |
                        rowData(cds_filt)$gene.type == "common", #genes
                      colData(cds_filt)$species == cur_species &
                        rownames(colData(cds_filt)) %in% unlist(ProBarcodeBinsList[[cur_species]])]
  
  cur_cds$batch <-  cur_cds$dataset3
  
  cds_fData <- data.frame(fData(cur_cds))
  norm.matrix = get.norm.expr.matrix(cur_cds) # normalize by size factor
  
  cell.bins <- names(ProBarcodeBinsList[[cur_species]])
  
  bootstrapTPMProList[[cur_species]] <- list()
  for(this.bin in cell.bins) {
    print(paste0("Starting: ", this.bin))
    # gene
    # cell.type
    # bootstrapped median TPM
    # 97.5% CI
    # 95% CI
    # 5% CI
    # 2.5% CI
    
    cell_type <- this.bin
    # cds_pData <- data.frame(pData(cur_cds)[ProBarcodeBinsList[[cur_species]][[cell_type]][[this.bin]],])
    
    bootObject <- boot(t(norm.matrix[,ProBarcodeBinsList[[cur_species]][[cell_type]][[this.bin]]]),
                       BinTPM,
                       parallel = "multicore",
                       ncpus = 64,
                       R = 1000)
    
    saveRDS(bootObject, paste0(dir, "Objects/BootObjects_CellCorrection/bootObject_", gsub("/", ".", this.bin), "_", cur_species, ".rds"))
    
    temp_ci <- matrix(unlist(lapply(1:dim(bootObject$t)[2], function (i) {
      # [,1] 0.95 [4:5,1]
      # [,2] 0.975 [4:5,2]
      temp.ci <- boot.ci(bootObject,
                         type = "perc",
                         conf = c(0.8, 0.95),
                         index = i)$percent
      return(c(temp.ci[2,4], temp.ci[1,4], temp.ci[1,5], temp.ci[2,5]))
    }
    )), ncol = 4, byrow = TRUE)
    
    bootstrapTPMProList[[cur_species]][[this.bin]] <- data.frame(gene = names(bootObject$t0), # gene (is actually barcode here?)
                                                                 cell.type = rep(this.bin, dim(bootObject$t)[2]), # cell.type
                                                                 median.tpm = apply(bootObject$t, 2, median), # median boot.strapped TPM
                                                                 ci_5 = temp_ci[,1],
                                                                 ci_20 = temp_ci[,2],
                                                                 ci_80 = temp_ci[,3],
                                                                 ci_95 = temp_ci[,4])
  }
}

saveRDS(bootstrapTPMProList, paste0(dir, "Objects/bootstrapTPMProList_CellCorrection.rds"))
bootstrapTPMProList <- readRDS(paste0(dir, "Objects/bootstrapTPMProList_CellCorrection.rds"))

TPMListBootstrapPro <- list()
for(cur_species in names(bootstrapTPMProList)) {
  
  TPMListBootstrapPro[[cur_species]] <- sapply(bootstrapTPMProList[[cur_species]], function(x) {
    return(x$median.tpm)
  })
  
  TPMListBootstrapPro[[cur_species]] <- data.frame(TPMListBootstrapPro[[cur_species]])
  colnames(TPMListBootstrapPro[[cur_species]]) <- sapply(bootstrapTPMProList[[cur_species]], function(x) {head(x$cell.type, n = 1)})
  rownames(TPMListBootstrapPro[[cur_species]]) <- bootstrapTPMProList[[cur_species]][[1]]$gene
}

names(TPMListBootstrapPro) <- c("C.elegans", "C.briggsae")

TPMListBootstrapProSubset <- list()
TPMListBootstrapProSubset[["C.elegans"]] <- TPMListBootstrapPro[["C.elegans"]][rownames(TPMListBootstrapPro[["C.elegans"]]) %in% rownames(TPMListBootstrapPro[["C.briggsae"]]),]
TPMListBootstrapProSubset[["C.briggsae"]] <- TPMListBootstrapPro[["C.briggsae"]][rownames(TPMListBootstrapPro[["C.briggsae"]]) %in% rownames(TPMListBootstrapPro[["C.elegans"]]),]

# Reordered
TPMListBootstrapProSubset[["C.elegans"]] <- TPMListBootstrapProSubset[["C.elegans"]][match(rownames(TPMListBootstrapProSubset[["C.briggsae"]]),
                                                                                           rownames(TPMListBootstrapProSubset[["C.elegans"]])),]

saveRDS(TPMListBootstrapPro, paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))
saveRDS(TPMListBootstrapProSubset, paste0(dir, "Objects/TPMListBootstrapProSubset_CellCorrection.rds"))


## Lower quantile bootstrap cutoff
TPMDataframeBootstrapPro <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(TPMDataframeBootstrapPro) <- c(colnames(bootstrapTPMProList[[1]][[1]]), "species")
for(cur_species in names(bootstrapTPMProList)) {
  # print(paste0("Starting: ", cur_species))
  for(dataframe in bootstrapTPMProList[[cur_species]]) {
    # print(dim(dataframe))
    if(is.null(dim(dataframe))) {
      print(head(dataframe))
    } else {
      dataframe$species = cur_species
      TPMDataframeBootstrapPro <- rbind(TPMDataframeBootstrapPro, dataframe)
    }
  }
}

BinarizeBootstrapListPro <- list()
for(cur_species in levels(factor(TPMDataframeBootstrapPro$species))) {
  BinarizeBootstrapListPro[[cur_species]] <- dcast(TPMDataframeBootstrapPro[TPMDataframeBootstrapPro$species == cur_species, c(1,2,4)], gene ~ cell.type, value.var = "ci_5")
  rownames(BinarizeBootstrapListPro[[cur_species]]) <- BinarizeBootstrapListPro[[cur_species]]$gene
  BinarizeBootstrapListPro[[cur_species]]$gene <- NULL
  BinarizeBootstrapListPro[[cur_species]] <- BinarizeBootstrapListPro[[cur_species]][match(rownames(TPMListBootstrapPro[[cur_species]]), rownames(BinarizeBootstrapListPro[[cur_species]])),]
}

BinarizeBootstrapListPro[["C.elegans"]] <- BinarizeBootstrapListPro[["C.elegans"]][,match(colnames(BinarizeBootstrapListPro[["C.elegans"]]), colnames(BinarizeBootstrapListPro[["C.elegans"]]))]
BinarizeBootstrapListPro[["C.briggsae"]] <- BinarizeBootstrapListPro[["C.briggsae"]][,match(colnames(BinarizeBootstrapListPro[["C.briggsae"]]), colnames(BinarizeBootstrapListPro[["C.briggsae"]]))]

saveRDS(BinarizeBootstrapListPro, paste0(dir, "Objects/BinarizeBootstrapListPro_CellCorrection.rds"))

BinarizeBootstrapListProSubset <- list()
BinarizeBootstrapListProSubset[["C.elegans"]] <- BinarizeBootstrapListPro[["C.elegans"]][rownames(BinarizeBootstrapListPro[["C.elegans"]]) %in% rownames(BinarizeBootstrapListPro[["C.briggsae"]]),]
BinarizeBootstrapListProSubset[["C.briggsae"]] <- BinarizeBootstrapListPro[["C.briggsae"]][rownames(BinarizeBootstrapListPro[["C.briggsae"]]) %in% rownames(BinarizeBootstrapListPro[["C.elegans"]]),]

# Reordered
BinarizeBootstrapListProSubset[["C.elegans"]] <- BinarizeBootstrapListProSubset[["C.elegans"]][match(rownames(BinarizeBootstrapListProSubset[["C.briggsae"]]),
                                                                                                     rownames(BinarizeBootstrapListProSubset[["C.elegans"]])),]

saveRDS(BinarizeBootstrapListProSubset, paste0(dir, "Objects/BinarizeBootstrapListProSubset_CellCorrection.rds"))


################################################
# Calculate cell distances from boot objects
################################################

# iteratively loads the boot objects, calculates the cell distances, and saves the iterations in an object to calculate the gene distances on
library(reldist) #for distributions and gini
library(reshape2)
source(paste0(dir, "../Scripts/CalcFunctions.R"))

TPMListBootstrapPro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))

cell.bins <- colnames(TPMListBootstrapPro[["C.elegans"]])

cel_gene_in <- rownames(TPMListBootstrapPro[["C.elegans"]]) %in% rownames(TPMListBootstrapPro[["C.briggsae"]])
cbr_gene_in <- rownames(TPMListBootstrapPro[["C.briggsae"]]) %in% rownames(TPMListBootstrapPro[["C.elegans"]])
gene_in_names <- rownames(TPMListBootstrapPro[["C.elegans"]])[cel_gene_in]

cel_cell_gene_list <- list()
cbr_cell_gene_list <- list()

# Want to pre-allocate the memory for each of the matrices
# Empty matrix at each iteration that is i by cell.type
temp_matrix <- matrix(nrow = 1000, ncol = length(cell.bins), dimnames = list(NULL, cell.bins))

for(this.gene in gene_in_names) {
  cel_cell_gene_list[[this.gene]] <- temp_matrix
  cbr_cell_gene_list[[this.gene]] <- temp_matrix
}

jsd_list <- list()
cor_list <- list()
cos_list <- list()
gini_list <- list()

# The gene distance calculation is constructed here
# and the cell distance is calculated afterwards
for(this.bin in cell.bins) {
  bin.number <- grep(paste0("^",this.bin,"$"), cell.bins)
  print(paste0(this.bin, ": ", bin.number, " of 304"), )
  
  # Read in each cell one at a time
  bootObject_cel <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection/bootObject_", gsub("/", ".", this.bin), "_", "C.elegans", ".rds"))
  bootObject_cbr <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection/bootObject_", gsub("/", ".", this.bin), "_", "C.briggsae", ".rds"))
  
  print(paste0("loaded ", this.bin, ": ", bin.number, " of 304"), )
  
  temp_cel <- bootObject_cel$t[,cel_gene_in]
  temp_cbr <- bootObject_cbr$t[,cbr_gene_in]
  
  colnames(temp_cel) <- gene_in_names
  colnames(temp_cbr) <- gene_in_names
  
  temp_cel_melt <- melt(temp_cel) # Var1 is iteration, Var2 is gene?
  temp_cbr_melt <- melt(temp_cbr) # Var1 is iteration, Var2 is gene?
  
  print(paste0("melted ", this.bin, ": ", bin.number, " of 304"), )
  
  # store the value of each gene from each cell as a column
  for(i in seq(0, (length(gene_in_names) - 1) * 1000, by = 1000)) { 
    this.gene <- gene_in_names[(i/1000) + 1]
    cel_cell_gene_list[[this.gene]][, bin.number] <- temp_cel_melt[(i + 1):(i + 1000),]$value
    cbr_cell_gene_list[[this.gene]][, bin.number] <- temp_cbr_melt[(i + 1):(i + 1000),]$value
  }
  
  print(paste0("Calc. distances on ", this.bin, ": ", bin.number, " of 304"), )
  
  # calculate cell distances
  temp_jsd <- list()
  temp_cor <- list()
  temp_cos <- list()
  temp_cel_gini <- list()
  temp_cbr_gini <- list()
  
  pseudocount = 1 / length(gene_in_names)
  for(i in seq(1, 1000)) {
    cel_tpm = temp_cel[i,] + pseudocount
    cbr_tpm = temp_cbr[i,] + pseudocount
    temp_jsd[i] <- sqrt(js_divg(p = cel_tpm/sum(cel_tpm), q = cbr_tpm/sum(cbr_tpm)))
    
    temp_cor[i] <- cor(log2(temp_cel[i,] + 1),
                       log2(temp_cbr[i,] + 1))
    
    temp_cos[i] <- cos_dist(cel_tpm,
                            cbr_tpm)
    
    temp_cel_gini[i] <- gini(log2(temp_cel[i,] + 1))
    temp_cbr_gini[i] <- gini(log2(temp_cbr[i,] + 1))
  }
  
  jsd_list[[this.bin]] <- data.frame(cell = this.bin,
                                     iteration = seq(1, 1000),
                                     jsd = unlist(temp_jsd))
  
  cor_list[[this.bin]] <- data.frame(cell = this.bin,
                                     iteration = seq(1, 1000),
                                     cor = unlist(temp_cor))
  
  cos_list[[this.bin]] <- data.frame(cell = this.bin,
                                     iteration = seq(1, 1000),
                                     cos = unlist(temp_cos))
  
  gini_list[[this.bin]] <- data.frame(cell = this.bin,
                                      iteration = seq(1, 1000),
                                      cel_gini = unlist(temp_cel_gini),
                                      cbr_gini = unlist(temp_cbr_gini))
}

saveRDS(cel_cell_gene_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/cel_cell_gene_list_pro.rds"))
saveRDS(cbr_cell_gene_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/cbr_cell_gene_list_pro.rds"))

jsd_df_bins <- do.call(rbind, jsd_list)
cor_df_bins <- do.call(rbind, cor_list)
cos_df_bins <- do.call(rbind, cos_list)
gini_df_bins <- do.call(rbind, gini_list)

saveRDS(jsd_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/jsd_df_bins_pro.rds"))
saveRDS(cor_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/cor_df_bins_pro.rds"))
saveRDS(cos_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/cos_df_bins_pro.rds"))
saveRDS(gini_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/gini_df_bins_pro.rds"))

################################################
# Calculate gene distances
################################################

library(parallel)

cel_cell_gene_list <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cel_cell_gene_list_pro.rds"))
cbr_cell_gene_list <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cbr_cell_gene_list_pro.rds"))

pseudocount = 1/length(cel_cell_gene_list[[1]][1,])

## Jensen-Shannon distance
pro_gene_jsd_list <- mclapply(names(cel_cell_gene_list), function(this.gene) {
  temp_jsd <- list()
  for(i in seq(1, 1000)) {
    p = cel_cell_gene_list[[this.gene]][i,] + pseudocount
    q = cbr_cell_gene_list[[this.gene]][i,] + pseudocount
    
    # Filter NA's out
    p = p[! (is.na(p) | is.na(q))]
    q = q[! (is.na(p) | is.na(q))]
    
    temp_jsd[i] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
  return(unlist(temp_jsd))
}, mc.cores = 12)

saveRDS(pro_gene_jsd_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_jsd_list.rds"))

temp_gene_names <- lapply(names(cel_cell_gene_list), function(this.gene) {
  return(rep(this.gene, 1000))
})

pro_gene_jsd_df <- data.frame(gene = unlist(temp_gene_names),
                               i = rep(seq(1, 1000), length(pro_gene_jsd_list)),
                               jsd = unlist(pro_gene_jsd_list))

saveRDS(pro_gene_jsd_df, paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_jsd_df.rds"))

## Pearson correlation
pro_gene_cor_list <- mclapply(names(cel_cell_gene_list), function(this.gene) {
  temp_cor <- list()
  for(i in seq(1, 1000)) {
    p = log2(cel_cell_gene_list[[this.gene]][i,] + 1)
    q = log2(cbr_cell_gene_list[[this.gene]][i,] + 1)
    
    # Filter NA's out
    p = p[! (is.na(p) | is.na(q))]
    q = q[! (is.na(p) | is.na(q))]
    
    temp_cor[i] <- cor(p, q)
  }
  return(unlist(temp_cor))
}, mc.cores = 12)

saveRDS(pro_gene_cor_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_cor_list.rds"))

pro_gene_cor_df <- data.frame(gene = unlist(temp_gene_names),
                               i = rep(seq(1, 1000), length(pro_gene_cor_list)),
                               cor = unlist(pro_gene_cor_list))

saveRDS(pro_gene_cor_df, paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_cor_df.rds"))

## Spearman correlation
pro_gene_cor_spearman_list <- mclapply(names(cel_cell_gene_list), function(this.gene) {
  temp_cor <- list()
  for(i in seq(1, 1000)) {
    p = log2(cel_cell_gene_list[[this.gene]][i,] + 1)
    q = log2(cbr_cell_gene_list[[this.gene]][i,] + 1)
    
    # Filter NA's out
    p = p[! (is.na(p) | is.na(q))]
    q = q[! (is.na(p) | is.na(q))]
    
    temp_cor[i] <- cor(p, q, method = "spearman")
  }
  return(unlist(temp_cor))
}, mc.cores = 12)

saveRDS(pro_gene_cor_spearman_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_cor_spearman_list.rds"))

pro_gene_cor_spearman_df <- data.frame(gene = unlist(temp_gene_names),
                                        i = rep(seq(1, 1000), length(pro_gene_cor_spearman_list)),
                                        cor = unlist(pro_gene_cor_spearman_list))

saveRDS(pro_gene_cor_spearman_df, paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_cor_spearman_df.rds"))

## Cosine angle
pro_gene_cos_list <- mclapply(names(cel_cell_gene_list), function(this.gene) {
  temp_cos <- list()
  for(i in seq(1, 1000)) {
    p = cel_cell_gene_list[[this.gene]][i,] + pseudocount
    q = cbr_cell_gene_list[[this.gene]][i,] + pseudocount
    
    # Filter NA's out
    p = p[! (is.na(p) | is.na(q))]
    q = q[! (is.na(p) | is.na(q))]
    
    temp_cos[i] <- cos_dist(p, q)
  }
  return(unlist(temp_cos))
}, mc.cores = 12)

saveRDS(pro_gene_cos_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_cos_list.rds"))

pro_gene_cos_df <- data.frame(gene = unlist(temp_gene_names),
                               i = rep(seq(1, 1000), length(pro_gene_cos_list)),
                               cos = unlist(pro_gene_cos_list))

saveRDS(pro_gene_cos_df, paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_cos_df.rds"))

## C. elegans / C. briggsae tau
pro_tau <- mclapply(names(cel_cell_gene_list), function(this.gene) {
  temp_cel_tau <- list()
  temp_cbr_tau <- list()
  
  for(i in seq(1, 1000)) {
    p = log2(cel_cell_gene_list[[this.gene]][i,] + 1)
    q = log2(cbr_cell_gene_list[[this.gene]][i,] + 1)
    
    # Filter NA's out
    p_filt = p[! (is.na(p) | is.na(q))]
    q_filt = q[! (is.na(p) | is.na(q))]
    
    temp_cel_tau[i] <- sum(1 - (p_filt/max(p_filt)))/(length(p_filt) - 1)
    temp_cbr_tau[i] <- sum(1 - (q_filt/max(q_filt)))/(length(q_filt) - 1)
  }
  return(data.frame(gene = this.gene,
                    i = seq(1, 1000),
                    cel_tau = unlist(temp_cel_tau),
                    cbr_tau = unlist(temp_cbr_tau)))
}, mc.cores = 20)

pro_tau_df <- do.call(rbind, pro_tau)

saveRDS(pro_tau_df, paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_tau_df.rds"))

#############
# Cell Data #
#############


bg_perc <- readRDS(paste0(dir, "Objects/bg_perc.rds"))
bg_perc_pro <- bg_perc[bg_perc$cell_type %in% names(ProBarcodeBinsList[["C.elegans"]]),]

jsd_df_bins <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/jsd_df_bins_pro.rds"))
cor_df_bins <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cor_df_bins_pro.rds"))
cos_df_bins <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cos_df_bins_pro.rds"))
gini_df_bins <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/gini_df_bins_pro.rds"))

# Have the same bins as the terminal cell types
jsd_df_bins$cell_type <- jsd_df_bins$cell
jsd_df_bins$cell_class <- "progenitor"
cor_df_bins$cell_type <- cor_df_bins$cell
cor_df_bins$cell_class <- "progenitor"
cos_df_bins$cell_type <- cos_df_bins$cell
cos_df_bins$cell_class <- "progenitor"
gini_df_bins$cell_type <- gini_df_bins$cell
gini_df_bins$cell_class <- "progenitor"

saveRDS(jsd_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/bootstrap_cell_jsd_df_pro.rds"))
saveRDS(cor_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/bootstrap_cell_cor_df_pro.rds"))
saveRDS(cos_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/bootstrap_cell_cos_df_pro.rds"))
saveRDS(gini_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/bootstrap_cell_gini_df_pro.rds"))

# 95% confidence interval and median
cell_data <- data.frame(cell_type_bin = unique(jsd_df_bins$cell))

cell_data$cell_class <- "progenitor"
cell_data$cell_type <- cell_data$cell_type_bin

# Jensen-Shannon Distance
cell_data <- left_join(cell_data, jsd_df_bins %>%
                         group_by(cell) %>%
                         summarise(jsd_median = median(jsd)), by = c('cell_type_bin' = 'cell')) # median
cell_data <- left_join(cell_data, (jsd_df_bins %>%
                                     group_by(cell) %>%
                                     arrange(desc(jsd), .by_group = TRUE) %>%
                                     filter(row_number() == 25) %>%
                                     rename("jsd_upper" = "jsd"))[, c("cell", "jsd_upper")], by = c('cell_type_bin' = 'cell')) # Upper
cell_data <- left_join(cell_data, (jsd_df_bins %>%
                                     group_by(cell) %>%
                                     arrange(desc(jsd), .by_group = TRUE) %>%
                                     filter(row_number() == 975) %>%
                                     rename("jsd_lower" = "jsd"))[, c("cell", "jsd_lower")], by = c('cell_type_bin' = 'cell')) # Lower

# Pearson correlation
cell_data <- left_join(cell_data, cor_df_bins %>%
                         group_by(cell) %>%
                         summarise(cor_median = median(cor)), by = c('cell_type_bin' = 'cell')) # median
cell_data <- left_join(cell_data, (cor_df_bins %>%
                                     group_by(cell) %>%
                                     arrange(desc(cor), .by_group = TRUE) %>%
                                     filter(row_number() == 25) %>%
                                     rename("cor_upper" = "cor"))[, c("cell", "cor_upper")], by = c('cell_type_bin' = 'cell')) # Upper
cell_data <- left_join(cell_data, (cor_df_bins %>%
                                     group_by(cell) %>%
                                     arrange(desc(cor), .by_group = TRUE) %>%
                                     filter(row_number() == 975) %>%
                                     rename("cor_lower" = "cor"))[, c("cell", "cor_lower")], by = c('cell_type_bin' = 'cell')) # Lower

# Cosine distance
cell_data <- left_join(cell_data, cos_df_bins %>%
                         group_by(cell) %>%
                         summarise(cos_median = median(cos)), by = c('cell_type_bin' = 'cell')) # median
cell_data <- left_join(cell_data, (cos_df_bins %>%
                                     group_by(cell) %>%
                                     arrange(desc(cos), .by_group = TRUE) %>%
                                     filter(row_number() == 25) %>%
                                     rename("cos_upper" = "cos"))[, c("cell", "cos_upper")], by = c('cell_type_bin' = 'cell')) # Upper
cell_data <- left_join(cell_data, (cos_df_bins %>% 
                                     group_by(cell) %>%
                                     arrange(desc(cos), .by_group = TRUE) %>%
                                     filter(row_number() == 975) %>%
                                     rename("cos_lower" = "cos"))[, c("cell", "cos_lower")], by = c('cell_type_bin' = 'cell')) # Lower

# Gini coefficient
cell_data <- left_join(cell_data, gini_df_bins %>%
                         group_by(cell) %>%
                         summarise(cel_gini_median = median(cel_gini)), by = c('cell_type_bin' = 'cell')) # median
cell_data <- left_join(cell_data, gini_df_bins %>%
                         group_by(cell) %>%
                         summarise(cbr_gini_median = median(cbr_gini)), by = c('cell_type_bin' = 'cell')) # median
cell_data <- left_join(cell_data, (gini_df_bins %>%
                                     group_by(cell) %>%
                                     arrange(desc(cel_gini), .by_group = TRUE) %>%
                                     filter(row_number() == 25) %>%
                                     rename("cel_gini_upper" = "cel_gini"))[, c("cell", "cel_gini_upper")], by = c('cell_type_bin' = 'cell')) # Upper
cell_data <- left_join(cell_data, (gini_df_bins %>% 
                                     group_by(cell) %>%
                                     arrange(desc(cel_gini), .by_group = TRUE) %>%
                                     filter(row_number() == 975) %>%
                                     rename("cel_gini_lower" = "cel_gini"))[, c("cell", "cel_gini_lower")], by = c('cell_type_bin' = 'cell')) # Lower
cell_data <- left_join(cell_data, (gini_df_bins %>% 
                                     group_by(cell) %>%
                                     arrange(desc(cbr_gini), .by_group = TRUE) %>%
                                     filter(row_number() == 25) %>%
                                     rename("cbr_gini_upper" = "cbr_gini"))[, c("cell", "cbr_gini_upper")], by = c('cell_type_bin' = 'cell')) # Upper
cell_data <- left_join(cell_data, (gini_df_bins %>% 
                                     group_by(cell) %>%
                                     arrange(desc(cbr_gini), .by_group = TRUE) %>%
                                     filter(row_number() == 975) %>%
                                     rename("cbr_gini_lower" = "cbr_gini"))[, c("cell", "cbr_gini_lower")], by = c('cell_type_bin' = 'cell')) # Lower

rownames(cell_data) <- cell_data$cell_type_bin

cell_data <- cell_data[,c("cell_type_bin", "cell_type", "cell_class", "jsd_median", "jsd_upper", 
                          "jsd_lower", "cor_median", "cor_upper", "cor_lower", 
                          "cos_median", "cos_upper", "cos_lower", "cel_gini_median",
                          "cbr_gini_median", "cel_gini_upper", "cel_gini_lower", "cbr_gini_upper", 
                          "cbr_gini_lower")]

saveRDS(cell_data, paste0(dir, "Objects/cell_data_pro_cell_bg.rds"))

##############
# Gene Data  #
##############

pro_gene_jsd_df <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_jsd_df.rds"))
pro_gene_cor_df <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_cor_df.rds"))
pro_gene_cor_spearman_df <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_cor_spearman_df.rds"))
pro_gene_cos_df <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_cos_df.rds"))
pro_tau_df <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/pro_gene_tau_df.rds"))

# 95% confidence interval and median
gene_data <- data.frame(gene = unique(pro_gene_jsd_df$gene))

distance_df <- left_join(left_join(left_join(left_join(pro_gene_jsd_df, pro_gene_cor_df, by = c("i" = "i", "gene" = "gene")),
                                             pro_gene_cor_spearman_df, by = c("i" = "i", "gene" = "gene"), suffix = c("_pearson", "_spearman")),
                                   pro_gene_cos_df, by = c("i" = "i", "gene" = "gene")),
                         pro_tau_df, by = c("i" = "i", "gene" = "gene"))

gene_data <- left_join(gene_data, distance_df %>%
                         group_by(gene) %>%
                         summarise(jsd_median = median(jsd, na.rm = TRUE),
                                   jsd_lower = quantile(jsd, na.rm = TRUE, probs = c(0.025)), 
                                   jsd_upper = quantile(jsd, na.rm = TRUE, probs = c(0.975)),
                                   pcor_median = median(cor_pearson, na.rm = TRUE),
                                   pcor_lower = quantile(cor_pearson, na.rm = TRUE, probs = c(0.025)), 
                                   pcor_upper = quantile(cor_pearson, na.rm = TRUE, probs = c(0.975)),
                                   scor_median = median(cor_spearman, na.rm = TRUE),
                                   scor_lower = quantile(cor_spearman, na.rm = TRUE, probs = c(0.025)), 
                                   scor_upper = quantile(cor_spearman, na.rm = TRUE, probs = c(0.975)),
                                   cos_median = median(cos, na.rm = TRUE),
                                   cos_lower = quantile(cos, na.rm = TRUE, probs = c(0.025)), 
                                   cos_upper = quantile(cos, na.rm = TRUE, probs = c(0.975)),
                                   cel_tau_median = median(cel_tau, na.rm = TRUE),
                                   cel_tau_lower = quantile(cel_tau, na.rm = TRUE, probs = c(0.025)),
                                   cel_tau_upper = quantile(cel_tau, na.rm = TRUE, probs = c(0.975)),
                                   cbr_tau_median = median(cbr_tau, na.rm = TRUE),
                                   cbr_tau_lower = quantile(cbr_tau, na.rm = TRUE, probs = c(0.025)),
                                   cbr_tau_upper = quantile(cbr_tau, na.rm = TRUE, probs = c(0.975))),
                       by = "gene",
                       suffix = c("", ""))

rownames(gene_data) <- gene_data$gene
saveRDS(gene_data, paste0(dir, "Objects/gene_data_pro_cell_bg.rds"))


