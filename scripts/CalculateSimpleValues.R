

library(presto) # for faster seurat wilcox test
library(Seurat)
library(monocle3)
library(parallel)

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

# load cds_bg_corrected
cds <- readRDS(paste0(dir, "Objects/cds_bg_20240829.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

rm(cds)

# load barcodes for terminal time bins
BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))
ProBarcodeBinsList <- readRDS(paste0(dir, "Objects/ProBarcodeBinsList.rds"))
JointBarcodeBinsList <- readRDS(paste0(dir, "Objects/JointBarcodeBinsList.rds"))

cell_type_time_bins <- readRDS(paste0(dir, "Objects/cell_type_time_bins.rds"))
CellTable <- readRDS(paste0(dir, "Objects/CellTable_20240826.rds"))

time_bin_vector <- c("lt_100", "100_130", "130_170", "170_210", "210_270", "270_330", "330_390", "390_450", "450_510", "510_580", "580_650", "650_710", "gt_710")

time_bin_df <- data.frame(bins = time_bin_vector,
                          start = c(0, 100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710),
                          end = c(100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710, 1000))

# Normalize expression matrix by Size factor
get.norm.expr.matrix <- function(cds) {
  mat = counts(cds)
  mat@x = mat@x / rep.int(pData(cds)$Size_Factor, diff(mat@p))
  return(mat)
}

## TPM calculations
BinTPM <- function(matrix) {
  matrix_rowMean <- Matrix::rowMeans(matrix)
  colSum <- sum(matrix_rowMean)
  return(matrix_rowMean/colSum * 1000000)
}

## TPM calculations
BinRobustTPM <- function(matrix) {
  matrix_rowMean <- apply(matrix, 1, function(x) {
    x = sort(x);
    mean(x[2:(length(x)-1)]) # try to make estimates more robust against outliers
  })
  colSum <- sum(matrix_rowMean)
  return(matrix_rowMean/colSum * 1000000)
}

# Basic TPM Calculation
CalcTPM <- function(cur_cds, column) {
  norm.expr = get.norm.expr.matrix(cur_cds) # normalize by size factor
  cell.bins <- unique(pData(cur_cds)[[column]])[! is.na(unique(pData(cur_cds)[[column]]))]
  
  bin.norm.means = sapply(cell.bins, function(this.bin) {
    if(sum(pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]) < 2) {
      print("Yes")
      norm.expr[,pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]]
    } else {
      Matrix::rowMeans(norm.expr[,pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]])
    }
  }) #for each cell.bin do the rowMean
  
  bin.tpm = base::sweep(
    bin.norm.means, 2,
    Matrix::colSums(bin.norm.means), "/") * 1000000
  
  return(bin.tpm)
}

TPMList <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  TPMList[[cur_species]] <- CalcTPM(cds_filt[rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                                                   "ce.unique",
                                                                                   "cb.unique") |
                                               rowData(cds_filt)$gene.type == "common",
                                             colData(cds_filt)$species == cur_species], "cell_type")
}

# tpm terminal
TPMBinsList <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  TPMBinsList[[cur_species]] <- mclapply(names(BarcodeBinsList[[cur_species]]), function(cell_type) { 
    print(cell_type)
    temp_tpm_out <- data.frame(matrix(nrow = sum(rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                                                       "ce.unique",
                                                                                       "cb.unique") |
                                                   rowData(cds_filt)$gene.type == "common"), ncol = 0))
    rownames(temp_tpm_out) <- rownames(cds_filt)[rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                                                       "ce.unique",
                                                                                       "cb.unique") |
                                                   rowData(cds_filt)$gene.type == "common"]
    for(i in seq(1, length(BarcodeBinsList[["C.elegans"]][[cell_type]]))) {
      temp_tpm <- data.frame(BinTPM(get.norm.expr.matrix(cds_filt[rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                                                                        "ce.unique",
                                                                                                        "cb.unique") |
                                                                    rowData(cds_filt)$gene.type == "common",
                                                                  BarcodeBinsList[[cur_species]][[cell_type]][[i]]])))
        temp_tpm_out <- cbind(temp_tpm_out, temp_tpm)
    }
    colnames(temp_tpm_out) <- names(BarcodeBinsList[["C.elegans"]][[cell_type]])
    return(temp_tpm_out)
  }, mc.cores = 32)
  TPMBinsList[[cur_species]] <- do.call("cbind", TPMBinsList[[cur_species]])
}

## Create a matched version
TPMBinsListSubset <- list()
TPMBinsListSubset[["C.elegans"]] <- TPMBinsList[["C.elegans"]][,colnames(TPMBinsList[["C.elegans"]]) %in% colnames(TPMBinsList[["C.briggsae"]])]
TPMBinsListSubset[["C.briggsae"]] <- TPMBinsList[["C.briggsae"]][,colnames(TPMBinsList[["C.briggsae"]]) %in% colnames(TPMBinsList[["C.elegans"]])]
TPMBinsListSubset[["C.elegans"]] <- TPMBinsListSubset[["C.elegans"]][rownames(TPMBinsListSubset[["C.elegans"]]) %in% rownames(TPMBinsListSubset[["C.briggsae"]]),]
TPMBinsListSubset[["C.briggsae"]] <- TPMBinsListSubset[["C.briggsae"]][rownames(TPMBinsListSubset[["C.briggsae"]]) %in% rownames(TPMBinsListSubset[["C.elegans"]]),]

# tpm progenitor
TPMBinsListPro <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  TPMBinsListPro[[cur_species]] <- CalcTPM(cds_filt[rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                                                          "ce.unique",
                                                                                          "cb.unique") |
                                                      rowData(cds_filt)$gene.type == "common",
                                                    colData(cds_filt)$species == cur_species], "lineage_broad")
}

TPMBinsListPro[["C.elegans"]] <- TPMBinsListPro[["C.elegans"]][,sort(colnames(TPMBinsListPro[["C.briggsae"]]))]
TPMBinsListPro[["C.briggsae"]] <- TPMBinsListPro[["C.briggsae"]][,sort(colnames(TPMBinsListPro[["C.elegans"]]))]
TPMBinsListProSubset <- list()
TPMBinsListProSubset[["C.elegans"]] <- TPMBinsListPro[["C.elegans"]][rownames(TPMBinsListPro[["C.elegans"]]) %in% rownames(TPMBinsListPro[["C.briggsae"]]),]
TPMBinsListProSubset[["C.briggsae"]] <- TPMBinsListPro[["C.briggsae"]][rownames(TPMBinsListPro[["C.briggsae"]]) %in% rownames(TPMBinsListPro[["C.elegans"]]),]

TPMBinsListComp <- list()
TPMBinsListComp[["C.elegans"]] <- cbind(TPMBinsList[["C.elegans"]], TPMBinsListPro[["C.elegans"]])
TPMBinsListComp[["C.briggsae"]] <- cbind(TPMBinsList[["C.briggsae"]], TPMBinsListPro[["C.briggsae"]])

TPMBinsListComp[["C.elegans"]] <- TPMBinsListComp[["C.elegans"]][-which(names(TPMBinsListComp[["C.elegans"]]) == "unassigned"),]
TPMBinsListComp[["C.briggsae"]] <- TPMBinsListComp[["C.briggsae"]][-which(names(TPMBinsListComp[["C.briggsae"]]) == "unassigned"),]

TPMBinsListCompSubset <- TPMBinsListComp
TPMBinsListCompSubset[["C.elegans"]] <- TPMBinsListCompSubset[["C.elegans"]][rownames(TPMBinsListCompSubset[["C.elegans"]]) %in% rownames(TPMBinsListCompSubset[["C.briggsae"]]),]
TPMBinsListCompSubset[["C.briggsae"]] <- TPMBinsListCompSubset[["C.briggsae"]][rownames(TPMBinsListCompSubset[["C.briggsae"]]) %in% rownames(TPMBinsListCompSubset[["C.elegans"]]),]

saveRDS(TPMBinsList, paste0(dir, "Objects/TPMTimeBinsList.rds"))
saveRDS(TPMBinsListSubset, paste0(dir, "Objects/TPMTimeBinsListSubset.rds"))
saveRDS(TPMBinsListPro, paste0(dir, "Objects/TPMBinsListPro.rds"))
saveRDS(TPMBinsListProSubset, paste0(dir, "Objects/TPMBinsListProSubset.rds"))
saveRDS(TPMBinsListComp, paste0(dir, "Objects/TPMTimeBinsListComp.rds"))
saveRDS(TPMBinsListCompSubset, paste0(dir, "Objects/TPMTimeBinsListCompSubset.rds"))

# create a version of the tpm time bins that is just one value per bin based on taking the mean of the values
TPMBinsListMean <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  TPMBinsListMean[[cur_species]] <- matrix(nrow = nrow(TPMBinsList[[cur_species]]), ncol = length(unique(cell_type_time_bins$cell_type)))
  colnames(TPMBinsListMean[[cur_species]]) <- unique(cell_type_time_bins$cell_type)
  rownames(TPMBinsListMean[[cur_species]]) <- rownames(TPMBinsList[[cur_species]])
  for(cell_type in unique(cell_type_time_bins$cell_type)) {
    if(nrow(cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]) > 1) {
      TPMBinsListMean[[cur_species]][, cell_type] <- rowMeans(TPMBinsList[[cur_species]][, paste0(cell_type, "_", cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$start_bin, ":",
                                                                                                  cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$end_bin)])
    } else {
      TPMBinsListMean[[cur_species]][, cell_type] <- TPMBinsList[[cur_species]][, paste0(cell_type, "_", cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$start_bin, ":",
                                                                                         cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$end_bin)]
    }
  }
}

TPMBinsListMeanSubset <- list()
TPMBinsListMeanSubset[["C.elegans"]] <- TPMBinsListMean[["C.elegans"]][,colnames(TPMBinsListMean[["C.elegans"]]) %in% colnames(TPMBinsListMean[["C.briggsae"]])]
TPMBinsListMeanSubset[["C.briggsae"]] <- TPMBinsListMean[["C.briggsae"]][,colnames(TPMBinsListMean[["C.briggsae"]]) %in% colnames(TPMBinsListMean[["C.elegans"]])]
TPMBinsListMeanSubset[["C.elegans"]] <- TPMBinsListMeanSubset[["C.elegans"]][rownames(TPMBinsListMeanSubset[["C.elegans"]]) %in% rownames(TPMBinsListMeanSubset[["C.briggsae"]]),]
TPMBinsListMeanSubset[["C.briggsae"]] <- TPMBinsListMeanSubset[["C.briggsae"]][rownames(TPMBinsListMeanSubset[["C.briggsae"]]) %in% rownames(TPMBinsListMeanSubset[["C.elegans"]]),]

TPMBinsListMeanComp <- list()
TPMBinsListMeanComp[["C.elegans"]] <- cbind(TPMBinsListMean[["C.elegans"]], TPMBinsListPro[["C.elegans"]])
TPMBinsListMeanComp[["C.briggsae"]] <- cbind(TPMBinsListMean[["C.briggsae"]], TPMBinsListPro[["C.briggsae"]])

TPMBinsListMeanComp[["C.elegans"]] <- TPMBinsListMeanComp[["C.elegans"]][,which(colnames(TPMBinsListMeanComp[["C.elegans"]]) != "unassigned")]
TPMBinsListMeanComp[["C.briggsae"]] <- TPMBinsListMeanComp[["C.briggsae"]][,which(colnames(TPMBinsListMeanComp[["C.briggsae"]]) != "unassigned")]

TPMBinsListMeanCompSubset <- TPMBinsListMeanComp
TPMBinsListMeanCompSubset[["C.elegans"]] <- TPMBinsListMeanCompSubset[["C.elegans"]][rownames(TPMBinsListMeanCompSubset[["C.elegans"]]) %in% rownames(TPMBinsListMeanCompSubset[["C.briggsae"]]),]
TPMBinsListMeanCompSubset[["C.briggsae"]] <- TPMBinsListMeanCompSubset[["C.briggsae"]][rownames(TPMBinsListMeanCompSubset[["C.briggsae"]]) %in% rownames(TPMBinsListMeanCompSubset[["C.elegans"]]),]

saveRDS(TPMBinsListMean, paste0(dir, "Objects/TPMTimeMeanBinsList.rds"))
saveRDS(TPMBinsListMeanSubset, paste0(dir, "Objects/TPMTimeMeanBinsListSubset.rds"))

# make into joint with the circle plot CellTable thing
# First we need to reform the CellTable to include the extra terminal time/pseudotime bins
# cell table
CellTable$MergedDatasetName <- ifelse(CellTable$TerminalDatasetName != "", CellTable$TerminalDatasetName, CellTable$cel_lineage_broad)
rownames(CellTable) <- CellTable$Lineage

# Body wall muscle
CellTable[CellTable$Lineage %in% c("MSapappp", "MSappapp", "MSpppapp", "MSppappp"),]$MergedDatasetName <- "BWM_headrow1_in"
CellTable[CellTable$Lineage %in% c("MSapaaap", "MSapapap", "MSppapap", "MSppaaap"),]$MergedDatasetName <- "BWM_headrow1_out"
CellTable[CellTable$Lineage %in% c("MSappppa", "Daapa", "Dpapa", "MSpppppa"),]$MergedDatasetName <- "BWM_headrow2_in"
CellTable[CellTable$Lineage %in% c("MSapappa", "MSapppaa", "MSppppaa", "MSppappa"),]$MergedDatasetName <- "BWM_headrow2_out"
CellTable[CellTable$Lineage %in% c("Daaaa","Daaap","Daapp","Dpaaa","Dpaap","Dpapp"),]$MergedDatasetName <- "BWM_anterior"
CellTable[CellTable$Lineage %in% c("Dapaa","Dapap","Dappaa","Dappap","Dapppa","Dapppp","Dppaa","Dppap","Dpppaa","Dpppap","Dppppa","Dppppp","MSaapppaa","MSaapppap","MSaappppa","MSaappppp","MSpapppaa","MSpapppap","MSpappppa","MSpappppp", "ABprpppppaa", "MSapppap", "MSappppp", "MSpappaa", "MSpappap", "MSppppap", "MSpppppp"),]$MergedDatasetName <- "BWM_middle"
CellTable[CellTable$Lineage %in% c("Capaaaa","Capaaap","Capaapa","Capaapp","Capapaa","Capapap","Capappa","Capappp","Cappaaa","Cappaap","Cappapa","Capppaa","Capppap","Cppaaaa","Cppaaap","Cppaapa","Cppaapp","Cppapaa","Cppapap","Cppappa","Cppappp","Cpppaaa","Cpppaap","Cpppapa","Cppppaa","Cppppap"),]$MergedDatasetName <- "BWM_posterior" 
CellTable[CellTable$Lineage %in% c("Cappapp","Cappppd","Cappppv","Cpppapp","Cpppppd","Cpppppv"),]$MergedDatasetName <- "BWM_far_posterior"

## Intestine
CellTable[CellTable$Lineage %in% c("Ealaad", "Ealaav", "Earaad", "Earaav"),]$MergedDatasetName <- "Intestine_anterior"
CellTable[CellTable$Lineage %in% c("Eplppa", "Eplppp", "Eprppa", "Eprppp"),]$MergedDatasetName <- "Intestine_far_posterior"
CellTable[CellTable$Lineage %in% c("Ealpp", "Earpp", "Eplap", "Eprap", "Eplpa", "Eprpa"),]$MergedDatasetName <- "Intestine_middle_and_posterior"
CellTable[CellTable$Lineage %in% c("Ealpa", "Earpa", "Ealap", "Earap", "Eplaa", "Epraa"),]$MergedDatasetName <- "Intestine_middle"

# hyp1V
CellTable[CellTable$Lineage %in% c("ABarappaapa"),]$MergedDatasetName <- "hyp1V"
CellTable[CellTable$MergedDatasetName %in% c("MC?"),]$MergedDatasetName <- "MC"
CellTable[CellTable$Lineage %in% c("MSpppaaa"),]$MergedDatasetName <- "hmc_homolog"

CellTable[CellTable$MergedDatasetName %in% c("GLR_1/GLR_2"),]$MergedDatasetName <- "GLR"
CellTable[CellTable$MergedDatasetName %in% c("pm3_pm4_pm5a/b/c"),]$MergedDatasetName <- "pm3_pm4_pm5"
CellTable[CellTable$MergedDatasetName %in% c("mc2a/mc2b"),]$MergedDatasetName <- "mc2"
CellTable[CellTable$MergedDatasetName %in% c("URB_and_possibly_URA"),]$MergedDatasetName <- "URB_and_URA"

saveRDS(CellTable, paste0(dir, "Objects/CellTable_Names_20240901.rds"))

CellTableTimeBins <- CellTable[! (CellTable$MergedDatasetName %in% cell_type_time_bins$cell_type),]

# check everything is good
cell_type_time_bins[! cell_type_time_bins$cell_type %in% CellTable$MergedDatasetName,]

CellTableTimeBins$TPMName <- CellTableTimeBins$MergedDatasetName
cell_type_counter = 0
cur_cell_type = cell_type_time_bins[1,]$cell_type
## Add a TPM name column that includes the time information
# cut this line [cell_type_time_bins$elegans_count >= 20 & cell_type_time_bins$briggsae_count >= 20,]
for(i in 1:nrow(cell_type_time_bins)) {
  print(cell_type_time_bins[i,"cell_type"])
  if(min(cell_type_time_bins[i,c("elegans_count", "briggsae_count")]) < 20) {
    print(paste("Skip: ", cell_type_time_bins[i,"cell_type"]))
  } else {
    if(cell_type_time_bins[i,]$cell_type != cur_cell_type) {
      cell_type_counter = 1 # reset counter for iterating over cell types
      cur_cell_type = cell_type_time_bins[i,]$cell_type
    } else {
      cell_type_counter = cell_type_counter + 1
    }
    
    temp_CellTableTimeBins <- CellTable[which(CellTable$MergedDatasetName == cell_type_time_bins[i,]$cell_type),]
    nhits = nrow(temp_CellTableTimeBins)
    temp_CellTableTimeBins$Lineage <- paste0(temp_CellTableTimeBins$Lineage, "_", cell_type_time_bins[i,]$start_bin, ":", cell_type_time_bins[i,]$end_bin)
    temp_CellTableTimeBins$TPMName <- paste0(temp_CellTableTimeBins$MergedDatasetName, "_", cell_type_time_bins[i,]$start_bin, ":", cell_type_time_bins[i,]$end_bin)
    if(cell_type_counter == 1) {
      CellTableTimeBins <- rbind(CellTableTimeBins, temp_CellTableTimeBins)
    } else {
      temp_CellTableTimeBins$Parent <- CellTableTimeBins[(nrow(CellTableTimeBins) - nhits + 1):nrow(CellTableTimeBins),]$Lineage
      temp_CellTableTimeBins$ParentDatasetName <- CellTableTimeBins[(nrow(CellTableTimeBins) - nhits + 1):nrow(CellTableTimeBins),]$TPMName
      CellTableTimeBins <- rbind(CellTableTimeBins, temp_CellTableTimeBins)
    }
  }
}

# Remove because postembryonic
CellTableTimeBins <- CellTableTimeBins[! rownames(CellTableTimeBins) %in% c("ABpraappaaa", "ABpraappaap", "ABplappaapa", "ABplappaapp", "ABprappappa", "ABprappappp"),]

saveRDS(CellTableTimeBins, paste0(dir, "Objects/CellTableTimeBins.rds"))
CellTableTimeBins <- readRDS(paste0(dir, "Objects/CellTableTimeBins.rds"))

# Cell_type x gene x species
tpm_matrix_list <- list()
tpm_matrix_t_list <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  temp_gene_tpm <- rbind(t(TPMBinsListMean[[cur_species]]), t(TPMBinsListPro[[cur_species]]))
  temp_gene_tpm <- data.frame(temp_gene_tpm)
  temp_gene_tpm$cell <- rownames(temp_gene_tpm)
  
  tpm_matrix_list[[cur_species]] <- left_join(CellTable[,c("Lineage", "MergedDatasetName")], temp_gene_tpm, by = c("MergedDatasetName" = "cell"))
  rownames(tpm_matrix_list[[cur_species]]) <- tpm_matrix_list[[cur_species]]$Lineage
  
  tpm_matrix_t_list[[cur_species]] <- t(tpm_matrix_list[[cur_species]][, -c(1, 2)])
  rownames(tpm_matrix_t_list[[cur_species]]) <- rownames(TPMBinsListMean[[cur_species]])
}
tpm_matrix_list <- tpm_matrix_t_list
saveRDS(tpm_matrix_list, paste0(dir, "Objects/tpm_matrix_list.rds"))

tpm_matrix_list_filt <- list()
tpm_matrix_list_filt[["C.elegans"]] <- tpm_matrix_list[["C.elegans"]][, which(colSums(tpm_matrix_list[["C.elegans"]] > 0) > 0)]
tpm_matrix_list_filt[["C.briggsae"]] <- tpm_matrix_list[["C.briggsae"]][, which(colSums(tpm_matrix_list[["C.briggsae"]] > 0) > 0)]

saveRDS(tpm_matrix_list_filt, paste0(dir, "Objects/tpm_matrix_list_filt.rds"))

# Cell_type x gene x species
tpm_matrix_time_list <- list()
tpm_matrix_time_t_list <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  temp_gene_tpm <- rbind(t(TPMBinsList[[cur_species]]), t(TPMBinsListPro[[cur_species]]))
  temp_gene_tpm <- data.frame(temp_gene_tpm)
  temp_gene_tpm$cell <- rownames(temp_gene_tpm)
  
  tpm_matrix_time_list[[cur_species]] <- left_join(CellTableTimeBins[,c("Lineage", "TPMName")], temp_gene_tpm, by = c("TPMName" = "cell"))
  rownames(tpm_matrix_time_list[[cur_species]]) <- tpm_matrix_time_list[[cur_species]]$Lineage
  
  tpm_matrix_time_t_list[[cur_species]] <- t(tpm_matrix_time_list[[cur_species]][, -c(1, 2)])
  rownames(tpm_matrix_time_t_list[[cur_species]]) <- rownames(TPMBinsList[[cur_species]])
}
tpm_matrix_time_list <- tpm_matrix_time_t_list
saveRDS(tpm_matrix_time_list, paste0(dir, "Objects/tpm_matrix_time_list.rds"))

tpm_matrix_time_list_filt <- list()
tpm_matrix_time_list_filt[["C.elegans"]] <- tpm_matrix_time_list[["C.elegans"]][, which(colSums(tpm_matrix_time_list[["C.elegans"]] > 0) > 0)]
tpm_matrix_time_list_filt[["C.briggsae"]] <- tpm_matrix_time_list[["C.briggsae"]][, which(colSums(tpm_matrix_time_list[["C.briggsae"]] > 0) > 0)]

# calculate tau in matrix format
TauList <- list()
#Tau = Sum(1-expNormalizedByMax)/1-numberOfTissues
for(cur_species in c("C.elegans", "C.briggsae")) {
  TauList[[cur_species]] <- apply(log(tpm_matrix_list_filt[[cur_species]] + 1),
                                  1, function(row) {
                                    return(sum(1 - (row/max(row)))/(length(row) - 1))
                                  })
}

saveRDS(TauList, paste0(dir, "Objects/TauList.rds"))

###################################
# trinarize
###################################

trinarize <- function(a = 1.5, b = 2, f  = 0.2, n, k) {
  # Calculate confidence in measurements
  # We calculate probability p that at least f of the cells detect (in each group),
  # and compare with pep, setting the binary pattern to 1 if p > pep,
  # -1 if p < (1 - pep) and 0 otherwise.
  # Args:
  #   k (int):	Number of observed positive cells
  #   n (int):	Total number of cells
  #   a (int):  Hyperparameter for beta(a,b)
  #   b (int):  Hyperparameter for beta(a,b)
  
  incb = pbeta(f, a + k, b - k + n, lower.tail = TRUE)
  if(incb == 0) {
    p = 1
  } else {
    p = 1 - exp(log(incb) + lbeta(a + k, b - k + n) + lgamma(a + b + n) - lgamma(a + k) - lgamma(b - k + n))
  }
  return(p)
}

trinarize_cds <- function(cur_cds, barcode_file, PEP = 0.05, f = 0.10) {
  # Wraps trinarize on a cds per cell.type based on a cell.type column
  # Depending on posterior probability, counts as detected, not detected, 
  # or intermediate as 1, -1, or 0 respectively
  #
  # Args:
  #   cur_cds (cds):	Cds to be analyzed
  #   column (char):	Column of cell.type info in cds pData
  #
  cell.bins <- names(barcode_file)
  
  gene_list <- rownames(cur_cds)
  bin.binarize <- data.frame(matrix(ncol = length(cell.bins), nrow = length(gene_list)))
  
  colnames(bin.binarize) <- cell.bins
  rownames(bin.binarize) <- gene_list
  
  for(cell_type in cell.bins) {
    print(paste0("Starting: ", cell_type))
    
    temp_counts = counts(cur_cds[,which(colData(cur_cds)[[column]] == cell_type)])
    total_cells = length(temp_counts[1,])
    
    bin.binarize[[cell_type]] <- sapply(rowSums(temp_counts > 0), function(non_zero_counts) {
      p = trinarize(n = total_cells, k = non_zero_counts, f = f)
      
      if(p > 1 - PEP) { # Detected
        return(1)
      } else if(p < PEP) { # Not detected
        return(-1)
      } else { # Else is indeterminate
        return(0)
      }
    })
  }
  return(bin.binarize)
}

trinarize_cds(cur_cds = cds_filt, column = "cell_type")

# from pyscenic_motif_read_in
lineage_list <- readRDS(paste0(dir, "Objects/lineage_list.rds"))
cell_class_barcodes <- readRDS(paste0(dir, "Objects/cell_class_barcodes.rds"))
time_bin_vector <- c("lt_100", "100_130", "130_170", "170_210", "210_270", "270_330", "330_390", "390_450", "450_510", "510_580", "580_650", "650_710", "gt_710")

TPMLineageList <- list()
time_class_metadata_list <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  print(cur_species)
  species_cds <- cds_filt[rowData(cds_filt)$gene.type %in% c("common", 
                                                             ifelse(cur_species == "C.elegans",
                                                                    "ce.unique",
                                                                    "cb.unique")),
                          colData(cds_filt)$species == cur_species]
  TPMLineageList[[cur_species]] <- data.frame(matrix(ncol = 0, nrow = nrow(rowData(species_cds))))
  time_class_metadata_list[[cur_species]] <- data.frame(matrix(ncol = 4, nrow = 0))
  colnames(time_class_metadata_list[[cur_species]]) <- c("cell_class", "time_bin", "species", "n_cells")
  for(cell_class in names(cell_class_barcodes)) {
    print(cell_class)
    cell_class_cds <- species_cds[,colData(species_cds)$barcode %in% cell_class_barcodes[[cell_class]]]
    for(time_bin in time_bin_vector) {
      if(! time_bin %in% unique(colData(cell_class_cds)$embryo.time.bin)) {
        next
      } else {
        class_time_cds <- cell_class_cds[,colData(cell_class_cds)$embryo.time.bin == time_bin]
        temp_df <- data.frame(BinTPM(get.norm.expr.matrix(class_time_cds)))
        colnames(temp_df) <- paste0(cell_class, "_", time_bin)
        
        TPMLineageList[[cur_species]] <- cbind(TPMLineageList[[cur_species]], 
                                               temp_df)
        
        time_class_metadata_list[[cur_species]] <- rbind(time_class_metadata_list[[cur_species]],
                                                         data.frame(cell_class = cell_class,
                                                              time_bin = time_bin,
                                                              species = cur_species,
                                                              n_cells = ncol(class_time_cds)))
      }
    }
  }
}
time_class_metadata_df <- do.call(rbind, time_class_metadata_list)

time_class_metadata_df <- time_class_metadata_df %>%
  filter(n_cells > 20) %>%
  complete(time_bin, cell_class) %>%
  filter(! is.na(n_cells)) %>% arrange(cell_class, time_bin)

saveRDS(TPMLineageList, paste0(dir, "Objects/TPMLineageList.rds"))
saveRDS(time_class_metadata_df, paste0(dir, "Objects/time_class_metadata_df.rds"))

write.table(TPMLineageList[["C.elegans"]], paste0(dir, "Tables/TPMLineageList_elegans.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(TPMLineageList[["C.briggsae"]], paste0(dir, "Tables/TPMLineageList_briggsae.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

# early_cell_tpm calculation
early_cell_tpm_list <- list()
early_cell_tpm_list[["C.elegans"]] <- CalcTPM(cds_filt[rowData(cds_filt)$gene.type %in% c("common", "ce.unique"),
                 which(colData(cds_filt)$species == "C.elegans" & colData(cds_filt)$lineage_packer == "28_cell_or_earlier")], "lineage_packer")

early_cell_tpm_list[["C.briggsae"]] <- CalcTPM(cds_filt[rowData(cds_filt)$gene.type %in% c("common", "cb.unique"),
                 which(colData(cds_filt)$species == "C.briggsae" & colData(cds_filt)$manual_lineage2 == "28_cell_or_earlier")], "manual_lineage2")

saveRDS(early_cell_tpm_list, paste0(dir, "Objects/early_cell_tpm_list.rds"))
