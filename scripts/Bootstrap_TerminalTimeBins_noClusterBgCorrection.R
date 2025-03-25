# for each set of barcodes coming from the time bin barcode list, I will take a sample of it with replacement.
# inside of this set of barcodes, the percentage of the background will be estimated and multiplied by those
# experiments background tpm.

library(monocle3)
library(Seurat)
library(boot)
library(parallel)
library(dplyr)

set.seed(2022)

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

# source 
# source(paste0(dir, "../Scripts/CalcFunctions.R"))
cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240903.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

rm(cds)

# load barcode bins from EmbryoTimeBins
BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))
ProBarcodeBinsList <- readRDS(paste0(dir, "Objects/ProBarcodeBinsList.rds"))
JointBarcodeBinsList <- readRDS(paste0(dir, "Objects/JointBarcodeBinsList.rds"))

dir.create(paste0(dir, "Objects/BootObjects_CellCorrection/"))

lookup_terminal <- data.frame(cell_type_bins = unlist(lapply(BarcodeBinsList[["C.elegans"]], names)),
                              cell_type = NA)

for(cell_type in names(BarcodeBinsList[["C.elegans"]])) {
  for(cell_type_bin in names(BarcodeBinsList[["C.elegans"]][[cell_type]])) {
    lookup_terminal[lookup_terminal$cell_type_bins == cell_type_bin,]$cell_type <- cell_type
  }
}

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
# Bootstrap TPM for Terminal Cell Types
##########################################################################

bootstrapTPMList <- list()
column = "cell_type"
for(cur_species in c("C.elegans", "C.briggsae")) {
  print(cur_species)
  cur_cds <- cds_filt[rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                            "ce.unique",
                                                            "cb.unique") |
                        rowData(cds_filt)$gene.type == "common", #genes
                      colData(cds_filt)$species == cur_species &
                        rownames(colData(cds_filt)) %in% unlist(BarcodeBinsList[[cur_species]])]
  
  cur_cds$batch <-  cur_cds$dataset3
  
  cds_fData <- data.frame(fData(cur_cds))
  norm.matrix = get.norm.expr.matrix(cur_cds) # normalize by size factor
  
  cell.bins <- lookup_terminal$cell_type_bins
  
  bootstrapTPMList[[cur_species]] <- list()
  for(this.bin in cell.bins) {
    print(paste0("Starting: ", this.bin))
    # gene
    # cell.type
    # bootstrapped median TPM
    # 97.5% CI
    # 95% CI
    # 5% CI
    # 2.5% CI
    
    cell_type <- lookup_terminal[lookup_terminal$cell_type_bin == this.bin,]$cell_type
    
    # norm.expr.temp <- norm.expr[BarcodeBinsList[[cur_species]][[cell_type]][[this.bin]],]
    # cds_pData <- data.frame(pData(cur_cds)[BarcodeBinsList[[cur_species]][[cell_type]][[this.bin]],])
    
    bootObject <- boot(t(norm.matrix[,BarcodeBinsList[[cur_species]][[cell_type]][[this.bin]]]),
                       BinTPM,
                       parallel = "multicore",
                       ncpus = 64,
                       R = 1000)
    
    saveRDS(bootObject, paste0(dir, "Objects/BootObjects_CellCorrection/bootObject_", this.bin, "_", cur_species, ".rds"))
    
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
    
    bootstrapTPMList[[cur_species]][[this.bin]] <- data.frame(gene = names(bootObject$t0), # gene 
                                                              cell.type = rep(this.bin, dim(bootObject$t)[2]), # cell.type
                                                              median.tpm = apply(bootObject$t, 2, median), # median boot.strapped TPM
                                                              ci_5 = temp_ci[,1],
                                                              ci_20 = temp_ci[,2],
                                                              ci_80 = temp_ci[,3],
                                                              ci_95 = temp_ci[,4])
  }
}

saveRDS(bootstrapTPMList, paste0(dir, "Objects/bootstrapTPMTimeBinsList_CellCorrection.rds"))
bootstrapTPMList <- readRDS(paste0(dir, "Objects/bootstrapTPMTimeBinsList_CellCorrection.rds"))

TPMListBootstrap <- list()
for(cur_species in names(bootstrapTPMList)) {
  TPMListBootstrap[[cur_species]] <- sapply(bootstrapTPMList[[cur_species]], function(x) {
    return(x$median.tpm)
  })
  
  TPMListBootstrap[[cur_species]] <- data.frame(TPMListBootstrap[[cur_species]])
  colnames(TPMListBootstrap[[cur_species]]) <- sapply(bootstrapTPMList[[cur_species]], function(x) {head(x$cell.type, n = 1)})
  rownames(TPMListBootstrap[[cur_species]]) <- bootstrapTPMList[[cur_species]][[1]]$gene
}

names(TPMListBootstrap) <- c("C.elegans", "C.briggsae")

TPMListBootstrapSubset <- list()
TPMListBootstrapSubset[["C.elegans"]] <- TPMListBootstrap[["C.elegans"]][rownames(TPMListBootstrap[["C.elegans"]]) %in% rownames(TPMListBootstrap[["C.briggsae"]]),]
TPMListBootstrapSubset[["C.briggsae"]] <- TPMListBootstrap[["C.briggsae"]][rownames(TPMListBootstrap[["C.briggsae"]]) %in% rownames(TPMListBootstrap[["C.elegans"]]),]

# Reordered
TPMListBootstrapSubset[["C.elegans"]] <- TPMListBootstrapSubset[["C.elegans"]][match(rownames(TPMListBootstrapSubset[["C.briggsae"]]),
                                                                                     rownames(TPMListBootstrapSubset[["C.elegans"]])),]

saveRDS(TPMListBootstrap, paste0(dir, "Objects/TPMTimeBinsListBootstrap_CellCorrection.rds"))
saveRDS(TPMListBootstrapSubset, paste0(dir, "Objects/TPMTimeBinsListBootstrapSubset_CellCorrection.rds"))

TPMListBootstrap <- readRDS(paste0(dir, "Objects/TPMTimeBinsListBootstrap_CellCorrection.rds"))

## Lower quantile bootstrap cutoff
TPMDataframeBootstrap <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(TPMDataframeBootstrap) <- c(colnames(bootstrapTPMList[[1]][[1]]), "species")
for(cur_species in names(bootstrapTPMList)) {
  # print(paste0("Starting: ", cur_species))
  for(dataframe in bootstrapTPMList[[cur_species]]) {
    print(dim(dataframe))
    if(is.null(dim(dataframe))) {
      print(head(dataframe))
    } else {
      dataframe$species = cur_species
      TPMDataframeBootstrap <- rbind(TPMDataframeBootstrap, dataframe)
    }
  }
}

saveRDS(TPMDataframeBootstrap, paste0(dir, "Objects/TPMDataframeBootstrap_CellCorrection.rds"))
TPMTimeBinsList <- readRDS(paste0(dir, "Objects/TPMTimeBinsList.rds"))

library(reshape2)

BinarizeBootstrapList <- list()
for(cur_species in levels(factor(TPMDataframeBootstrap$species))) {
  BinarizeBootstrapList[[cur_species]] <- dcast(TPMDataframeBootstrap[TPMDataframeBootstrap$species == cur_species, c(1,2,4)], gene ~ cell.type, value.var = "ci_5")
  rownames(BinarizeBootstrapList[[cur_species]]) <- BinarizeBootstrapList[[cur_species]]$gene
  BinarizeBootstrapList[[cur_species]]$gene <- NULL
  BinarizeBootstrapList[[cur_species]] <- BinarizeBootstrapList[[cur_species]][match(rownames(TPMTimeBinsList[[cur_species]]), rownames(BinarizeBootstrapList[[cur_species]])),]
}

BinarizeBootstrapList[["C.elegans"]] <- BinarizeBootstrapList[["C.elegans"]][,match(colnames(BinarizeBootstrapList[["C.briggsae"]]), colnames(BinarizeBootstrapList[["C.elegans"]]))]

saveRDS(BinarizeBootstrapList, paste0(dir, "Objects/BinarizeBootstrapTimeBinsList_CellCorrection.rds"))

BinarizeBootstrapListSubset <- list()
BinarizeBootstrapListSubset[["C.elegans"]] <- BinarizeBootstrapList[["C.elegans"]][rownames(BinarizeBootstrapList[["C.elegans"]]) %in% rownames(BinarizeBootstrapList[["C.briggsae"]]),]
BinarizeBootstrapListSubset[["C.briggsae"]] <- BinarizeBootstrapList[["C.briggsae"]][rownames(BinarizeBootstrapList[["C.briggsae"]]) %in% rownames(BinarizeBootstrapList[["C.elegans"]]),]

# Reordered
BinarizeBootstrapListSubset[["C.elegans"]] <- BinarizeBootstrapListSubset[["C.elegans"]][match(rownames(BinarizeBootstrapListSubset[["C.briggsae"]]),
                                                                                               rownames(BinarizeBootstrapListSubset[["C.elegans"]])),]

saveRDS(BinarizeBootstrapListSubset, paste0(dir, "Objects/BinarizeBootstrapTimeBinsListSubset_CellCorrection.rds"))

## Upper quantile bootstrap cutoff
BinarizeBootstrapList_Upper <- list()

for(cur_species in levels(factor(TPMDataframeBootstrap$species))) {
  BinarizeBootstrapList_Upper[[cur_species]] <- dcast(TPMDataframeBootstrap[TPMDataframeBootstrap$species == cur_species, c(1,2,7)], gene ~ cell.type, value.var = "ci_95")
  rownames(BinarizeBootstrapList_Upper[[cur_species]]) <- BinarizeBootstrapList_Upper[[cur_species]]$gene
  BinarizeBootstrapList_Upper[[cur_species]]$gene <- NULL
  BinarizeBootstrapList_Upper[[cur_species]] <- BinarizeBootstrapList_Upper[[cur_species]][match(rownames(TPMTimeBinsList[[cur_species]]), rownames(BinarizeBootstrapList_Upper[[cur_species]])),]
}

BinarizeBootstrapList_Upper[["C.elegans"]] <- BinarizeBootstrapList_Upper[["C.elegans"]][,match(colnames(BinarizeBootstrapList_Upper[["C.briggsae"]]), colnames(BinarizeBootstrapList_Upper[["C.elegans"]]))]

saveRDS(BinarizeBootstrapList_Upper, paste0(dir, "Objects/BinarizeBootstrapTimeBinsListUpper_CellCorrection.rds"))

BinarizeBootstrapListSubset_Upper <- list()
BinarizeBootstrapListSubset_Upper[["C.elegans"]] <- BinarizeBootstrapList_Upper[["C.elegans"]][rownames(BinarizeBootstrapList_Upper[["C.elegans"]]) %in% rownames(BinarizeBootstrapList_Upper[["C.briggsae"]]),]
BinarizeBootstrapListSubset_Upper[["C.briggsae"]] <- BinarizeBootstrapList_Upper[["C.briggsae"]][rownames(BinarizeBootstrapList_Upper[["C.briggsae"]]) %in% rownames(BinarizeBootstrapList_Upper[["C.elegans"]]),]

# Reordered
BinarizeBootstrapListSubset_Upper[["C.elegans"]] <- BinarizeBootstrapListSubset_Upper[["C.elegans"]][match(rownames(BinarizeBootstrapListSubset_Upper[["C.briggsae"]]),
                                                                                               rownames(BinarizeBootstrapListSubset_Upper[["C.elegans"]])),]

saveRDS(BinarizeBootstrapListSubset_Upper, paste0(dir, "Objects/BinarizeBootstrapTimeBinsListUpperSubset_CellCorrection.rds"))


################################################
# Calculate cell distances from boot objects
################################################

# iteratively loads the boot objects, calculates the cell distances, and saves the iterations in an object to calculate the gene distances on
library(reldist) #for distributions and gini
source(paste0(dir, "../Scripts/CalcFunctions.R"))

TPMListBootstrap <- readRDS(paste0(dir, "Objects/TPMTimeBinsListBootstrap_CellCorrection.rds"))

cell.bins <- colnames(TPMListBootstrap[["C.elegans"]])

cel_gene_in <- rownames(TPMListBootstrap[["C.elegans"]]) %in% rownames(TPMListBootstrap[["C.briggsae"]])
cbr_gene_in <- rownames(TPMListBootstrap[["C.briggsae"]]) %in% rownames(TPMListBootstrap[["C.elegans"]])
gene_in_names <- rownames(TPMListBootstrap[["C.elegans"]])[cel_gene_in]

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
  print(paste0(this.bin, ": ", bin.number, " of 313"), )
  
  # Read in each cell one at a time
  bootObject_cel <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection/bootObject_", this.bin, "_", "C.elegans", ".rds"))
  bootObject_cbr <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection/bootObject_", this.bin, "_", "C.briggsae", ".rds"))
  
  print(paste0("loaded ", this.bin, ": ", bin.number, " of 313"), )
  
  temp_cel <- bootObject_cel$t[,cel_gene_in]
  temp_cbr <- bootObject_cbr$t[,cbr_gene_in]
  
  colnames(temp_cel) <- gene_in_names
  colnames(temp_cbr) <- gene_in_names
  
  temp_cel_melt <- melt(temp_cel) # Var1 is iteration, Var2 is gene?
  temp_cbr_melt <- melt(temp_cbr) # Var1 is iteration, Var2 is gene?
  
  print(paste0("melted ", this.bin, ": ", bin.number, " of 313"), )
  
  # store the value of each gene from each cell as a column
  for(i in seq(0, (length(gene_in_names) - 1) * 1000, by = 1000)) { 
    this.gene <- gene_in_names[(i/1000) + 1]
    cel_cell_gene_list[[this.gene]][, bin.number] <- temp_cel_melt[(i + 1):(i + 1000),]$value
    cbr_cell_gene_list[[this.gene]][, bin.number] <- temp_cbr_melt[(i + 1):(i + 1000),]$value
  }
  
  print(paste0("Calc. distances on ", this.bin, ": ", bin.number, " of 313"), )
  
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

saveRDS(cel_cell_gene_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/cel_cell_gene_list.rds"))
saveRDS(cbr_cell_gene_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/cbr_cell_gene_list.rds"))

jsd_df_bins <- do.call(rbind, jsd_list)
cor_df_bins <- do.call(rbind, cor_list)
cos_df_bins <- do.call(rbind, cos_list)
gini_df_bins <- do.call(rbind, gini_list)

saveRDS(jsd_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/jsd_df_bins.rds"))
saveRDS(cor_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/cor_df_bins.rds"))
saveRDS(cos_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/cos_df_bins.rds"))
saveRDS(gini_df_bins, paste0(dir, "Objects/BootObjects_CellCorrection_out/gini_df_bins.rds"))

################################################
# Calculate gene distances
################################################

library(parallel)

cel_cell_gene_list <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cel_cell_gene_list.rds"))
cbr_cell_gene_list <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cbr_cell_gene_list.rds"))

pseudocount = 1/length(cel_cell_gene_list[[1]][1,])

## Jensen-Shannon distance
term_gene_jsd_list <- mclapply(names(cel_cell_gene_list), function(this.gene) {
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
}, mc.cores = 24)

saveRDS(term_gene_jsd_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/term_gene_jsd_list.rds"))

temp_gene_names <- lapply(names(cel_cell_gene_list), function(this.gene) {
  return(rep(this.gene, 1000))
})

term_gene_jsd_df <- data.frame(gene = unlist(temp_gene_names),
                               i = rep(seq(1, 1000), length(term_gene_jsd_list)),
                               jsd = unlist(term_gene_jsd_list))

saveRDS(term_gene_jsd_df, paste0(dir, "Objects/BootObjects_CellCorrection_out/term_gene_jsd_df.rds"))

## Pearson correlation
term_gene_cor_list <- mclapply(names(cel_cell_gene_list), function(this.gene) {
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
}, mc.cores = 24)

saveRDS(term_gene_cor_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/term_gene_cor_list.rds"))

term_gene_cor_df <- data.frame(gene = unlist(temp_gene_names),
                               i = rep(seq(1, 1000), length(term_gene_cor_list)),
                               cor = unlist(term_gene_cor_list))

saveRDS(term_gene_cor_df, paste0(dir, "Objects/BootObjects_CellCorrection_out/term_gene_cor_df.rds"))

## Spearman correlation
term_gene_cor_spearman_list <- mclapply(names(cel_cell_gene_list), function(this.gene) {
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

saveRDS(term_gene_cor_spearman_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/term_gene_cor_spearman_list.rds"))

term_gene_cor_spearman_df <- data.frame(gene = unlist(temp_gene_names),
                                        i = rep(seq(1, 1000), length(term_gene_cor_spearman_list)),
                                        cor = unlist(term_gene_cor_spearman_list))

saveRDS(term_gene_cor_spearman_df, paste0(dir, "Objects/BootObjects_CellCorrection_out/term_gene_cor_spearman_df.rds"))

## Cosine angle
term_gene_cos_list <- mclapply(names(cel_cell_gene_list), function(this.gene) {
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
}, mc.cores = 8)

saveRDS(term_gene_cos_list, paste0(dir, "Objects/BootObjects_CellCorrection_out/term_gene_cos_list.rds"))

term_gene_cos_df <- data.frame(gene = unlist(temp_gene_names),
                               i = rep(seq(1, 1000), length(term_gene_cos_list)),
                               cos = unlist(term_gene_cos_list))

saveRDS(term_gene_cos_df, paste0(dir, "Objects/BootObjects_CellCorrection_out/term_gene_cos_df.rds"))

## C. elegans / C. briggsae tau
term_tau <- mclapply(names(cel_cell_gene_list), function(this.gene) {
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

term_tau_df <- do.call(rbind, term_tau)

saveRDS(term_tau_df, paste0(dir, "Objects/BootObjects_CellCorrection_out/term_gene_tau_df.rds"))

#############
# Cell Data #
#############

bg_perc <- readRDS(paste0(dir, "Objects/bg_perc.rds"))
ListOfCellTypes <- readRDS(paste0(dir, "Objects/ListOfCellTypes.rds"))

ListOfCellTypes[["Muscle"]] <- c(ListOfCellTypes[["Muscle"]], "BWM_anterior", "BWM_posterior")
ListOfCellTypes[["Pharynx and rectal"]] <- c(ListOfCellTypes[["Pharynx and rectal"]], "mc2", "pm3_pm4_pm5")
ListOfCellTypes[["Mesoderm"]] <- c(ListOfCellTypes[["Mesoderm"]], "GLR")

bg_perc_term <- bg_perc[bg_perc$cell_type %in% names(BarcodeBinsList[["C.elegans"]]),]

jsd_df_bins <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/jsd_df_bins.rds"))
cor_df_bins <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cor_df_bins.rds"))
cos_df_bins <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cos_df_bins.rds"))
gini_df_bins <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/gini_df_bins.rds"))

cell_df_bins <- left_join(left_join(left_join(jsd_df_bins, cor_df_bins, by = c("iteration" = "iteration", "cell" = "cell")),
                                    cos_df_bins, by = c("iteration" = "iteration", "cell" = "cell")),
                          gini_df_bins, by = c("iteration" = "iteration", "cell" = "cell"))

cell_data <- data.frame(cell_type_bin = unique(cell_df_bins$cell))

cell_data$cell_type <- NA
cell_data$cell_class <- NA

for(cur_cell_type in bg_perc_term$cell_type_bin) {
  print(cur_cell_type)
  cell_type <- unique(bg_perc[bg_perc$cell_type_bin %in% cur_cell_type,]$cell_type)
  cell_class <- names(ListOfCellTypes)[grepl(paste0("\\b", cell_type,"\\b"), ListOfCellTypes)]
  
  cell_data[cell_data$cell_type_bin %in% cur_cell_type,]$cell_type <- cell_type
  cell_data[cell_data$cell_type_bin %in% cur_cell_type,]$cell_class <- cell_class
}

cell_data <- left_join(cell_data, cell_df_bins %>%
                         group_by(cell) %>%
                         summarise(jsd_median = median(jsd, na.rm = TRUE),
                                   jsd_lower = quantile(jsd, na.rm = TRUE, probs = c(0.025)), 
                                   jsd_upper = quantile(jsd, na.rm = TRUE, probs = c(0.975)),
                                   cor_median = median(cor, na.rm = TRUE),
                                   cor_lower = quantile(cor, na.rm = TRUE, probs = c(0.025)), 
                                   cor_upper = quantile(cor, na.rm = TRUE, probs = c(0.975)),
                                   cos_median = median(cos, na.rm = TRUE),
                                   cos_lower = quantile(cos, na.rm = TRUE, probs = c(0.025)), 
                                   cos_upper = quantile(cos, na.rm = TRUE, probs = c(0.975)),
                                   cel_gini_median = median(cel_gini, na.rm = TRUE),
                                   cel_gini_lower = quantile(cel_gini, na.rm = TRUE, probs = c(0.025)), 
                                   cel_gini_upper = quantile(cel_gini, na.rm = TRUE, probs = c(0.975)),
                                   cbr_gini_median = median(cbr_gini, na.rm = TRUE),
                                   cbr_gini_lower = quantile(cbr_gini, na.rm = TRUE, probs = c(0.025)), 
                                   cbr_gini_upper = quantile(cbr_gini, na.rm = TRUE, probs = c(0.975))),
                       by = c("cell_type_bin" = "cell"),
                       suffix = c("", ""))

rownames(cell_data) <- cell_data$cell_type_bin
saveRDS(cell_data, paste0(dir, "Objects/cell_data_term_cell_bg.rds"))

##############
# Gene Data  #
##############

# 95% confidence interval and median
gene_data <- data.frame(gene = unique(term_gene_jsd_df$gene))

distance_df <- left_join(left_join(left_join(left_join(term_gene_jsd_df, term_gene_cor_df, by = c("i" = "i", "gene" = "gene")),
                                             term_gene_cor_spearman_df, by = c("i" = "i", "gene" = "gene"), suffix = c("_pearson", "_spearman")),
                                   term_gene_cos_df, by = c("i" = "i", "gene" = "gene")),
                         term_tau_df, by = c("i" = "i", "gene" = "gene"))

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
saveRDS(gene_data, paste0(dir, "Objects/gene_data_term_cell_bg.rds"))

##############
# Mean Data  #
##############

cell_type_time_bins <- readRDS(paste0(dir, "Objects/cell_type_time_bins.rds"))

TPMListBootstrapMean <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  TPMListBootstrapMean[[cur_species]] <- matrix(nrow = nrow(TPMListBootstrap[[cur_species]]), ncol = length(unique(cell_type_time_bins$cell_type)))
  colnames(TPMListBootstrapMean[[cur_species]]) <- unique(cell_type_time_bins$cell_type)
  rownames(TPMListBootstrapMean[[cur_species]]) <- rownames(TPMListBootstrap[[cur_species]])
  for(cell_type in unique(cell_type_time_bins$cell_type)) {
    if(nrow(cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]) > 1) {
      TPMListBootstrapMean[[cur_species]][, cell_type] <- rowMeans(TPMListBootstrap[[cur_species]][, paste0(cell_type, "_", cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$start_bin, ":",
                                                                                                  cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$end_bin)])
    } else {
      TPMListBootstrapMean[[cur_species]][, cell_type] <- TPMListBootstrap[[cur_species]][, paste0(cell_type, "_", cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$start_bin, ":",
                                                                                         cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$end_bin)]
    }
  }
  TPMListBootstrapMean[[cur_species]] <- TPMListBootstrapMean[[cur_species]][, unique(cell_type_time_bins$cell_type)]
}

saveRDS(TPMListBootstrapMean, paste0(dir, "Objects/TPMListBootstrapMean_CellCorrection.rds"))

BinarizeBootstrapListMean <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  BinarizeBootstrapListMean[[cur_species]] <- matrix(nrow = nrow(BinarizeBootstrapList[[cur_species]]), ncol = length(unique(cell_type_time_bins$cell_type)))
  colnames(BinarizeBootstrapListMean[[cur_species]]) <- unique(cell_type_time_bins$cell_type)
  rownames(BinarizeBootstrapListMean[[cur_species]]) <- rownames(BinarizeBootstrapList[[cur_species]])
  for(cell_type in unique(cell_type_time_bins$cell_type)) {
    if(nrow(cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]) > 1) {
      BinarizeBootstrapListMean[[cur_species]][, cell_type] <- rowMeans(BinarizeBootstrapList[[cur_species]][, paste0(cell_type, "_", cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$start_bin, ":",
                                                                                                            cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$end_bin)])
    } else {
      BinarizeBootstrapListMean[[cur_species]][, cell_type] <- BinarizeBootstrapList[[cur_species]][, paste0(cell_type, "_", cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$start_bin, ":",
                                                                                                   cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$end_bin)]
    }
  }
  BinarizeBootstrapListMean[[cur_species]] <- BinarizeBootstrapListMean[[cur_species]][, unique(cell_type_time_bins$cell_type)]
}

saveRDS(BinarizeBootstrapListMean, paste0(dir, "Objects/BinarizeBootstrapListMean_CellCorrection.rds"))

# CellData that uses these time bins

# -> Mean of the distances based off of the mean values
cell_data_mean <- data.frame(cell_data %>% group_by(cell_type) %>% summarise_all(mean))
rownames(cell_data_mean) <- cell_data_mean$cell_type
cell_data_mean <- cell_data_mean[cell_data[! duplicated(cell_data$cell_type),]$cell_type,]
cell_data_mean$cell_class <- cell_data[! duplicated(cell_data$cell_type),]$cell_class

for(cur_cell_type in cell_data_mean$cell_type) {
  cell_data_mean[cur_cell_type,]$cell_type_bin <- list(cell_data[cell_data$cell_type == cur_cell_type,]$cell_type_bin)
}

saveRDS(cell_data_mean, paste0(dir, "Objects/cell_data_term_mean_cell_bg.rds"))





