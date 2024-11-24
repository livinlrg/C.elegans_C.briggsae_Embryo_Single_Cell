library(data.table)
library(dplyr)
library(RcisTarget)
library(parallel)
library(ggplot2)
library(tidyr)
library(reshape2)
library(monocle3)
library(ggrepel)
library(ggExtra)
library(cowplot)
library(viridis)
library(pheatmap)

OldPath <- "/kimdata/livinlrg/scAnalysis/BobDataComb/"
dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

gene_data <- readRDS(paste0(dir, "Objects/gene_data_joint_cell_bg.rds"))
cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))
cell_data_bins <- readRDS(paste0(dir, "Objects/cell_data_joint_cell_bg.rds"))

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240903.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]
rm(cds)

cds_split <- readRDS(paste0(dir, "Objects/cds_objects/cds_split.rds"))

OGList <- readRDS(paste0(OldPath, "Objects/WS290/OGList.rds"))
OGNames <- readRDS(paste0(OldPath, "Objects/WS290/OGNames.rds"))
OGList_transcript_name <- readRDS(paste0(OldPath, "Objects/WS290/OGList_transcript_name.rds"))
OGList_short_name <- readRDS(paste0(OldPath, "Objects/WS290/OGList_short_name.rds"))
OGFullSpeciesList <- readRDS(paste0(OldPath, "Objects/WS290/OGFullSpeciesList.rds"))
og_data <- readRDS(paste0(OldPath, "Objects/WS290/og_data.rds"))

gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))

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


## Basically need to cycle through the orthogroups, and grab those rows from the matrix and assign them to these new orthogroups by adding them up

## Need to match the og information to the cds gene names
# short_name doesn't work
# OG_list is the wbgene names
# at the moment, cds name is the cds_gene_name
# could convert to the wbgene name where c.elegans takes dominance

# Either create a cds object that has the two WBG gene names separated
# Or generate two matrices initially and then combine them

matrix <- counts(cds_split)

OGVectorList <- lapply(OGList[lapply(OGList, function(x) sum(x > 0)) > 0], function(OG) {
  if(length(OG) == 1) {
    return(matrix[OG,])
  } else {
    return(colSums(matrix[OG,], sparseResult = TRUE))
  }
})

# Convert to sparseMatrix
OGMatrixList <- lapply(OGVectorList, as, "sparseMatrix")

# Convert to spareMatrix in unlisted format
OGMatrix <- do.call(cbind, OGMatrixList)

OGMatrix <- t(OGMatrix)

rownames(OGMatrix) <- names(OGList[lapply(OGList, function(x) sum(x > 0)) > 0])

saveRDS(OGMatrix, paste0(dir, "Objects/OGMatrix.rds"))

## Need to sort through the OG's and pull out the number of elegans and briggsae genes in each
## Number of species from each OG
cel_length <- lapply(OGFullSpeciesList, function(OG) {
  return(sum(OG == "elegans", na.rm = TRUE))
})
names(cel_length) <- OGNames
cel_length <- unlist(cel_length)

cbr_length <- lapply(OGFullSpeciesList, function(OG) {
  return(sum(OG == "briggsae", na.rm = TRUE))
})
names(cbr_length) <- OGNames
cbr_length <- unlist(cbr_length)

# Metadata for object
getMode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

### WormCat
OGList_WormCat1_fData <- mclapply(OGList_Raw, function(OG) {
  temp <- data.frame(strsplit(OG, " "))
  colnames(temp) <- c("OG")
  
  temp_cel <- temp[grep("CELEG.", temp$OG),]
  temp_cel <- gsub("CELEG.", "", temp_cel)
  temp_cel <- gsub("[a-z]$", "", temp_cel)
  
  temp_cel_wormcat <- gff_list_mrna[["elegans"]][which(gsub("[a-z]$", "", gff_list_mrna[["elegans"]]$short_name) %in% temp_cel),]$WormCat.1
  
  return(getMode(temp_cel_wormcat))
}, mc.cores = 16)

OGList_WormCat2_fData <- mclapply(OGList_Raw, function(OG) {
  temp <- data.frame(strsplit(OG, " "))
  colnames(temp) <- c("OG")
  
  temp_cel <- temp[grep("CELEG.", temp$OG),]
  temp_cel <- gsub("CELEG.", "", temp_cel)
  temp_cel <- gsub("[a-z]$", "", temp_cel)

  temp_cel_wormcat <- gff_list_mrna[["elegans"]][which(gsub("[a-z]$", "", gff_list_mrna[["elegans"]]$short_name) %in% temp_cel),]$WormCat.2
  
  return(getMode(temp_cel_wormcat))
}, mc.cores = 16)

OGList_WormCat3_fData <- mclapply(OGList_Raw, function(OG) {
  temp <- data.frame(strsplit(OG, " "))
  colnames(temp) <- c("OG")
  
  temp_cel <- temp[grep("CELEG.", temp$OG),]
  temp_cel <- gsub("CELEG.", "", temp_cel)
  temp_cel <- gsub("[a-z]$", "", temp_cel)

  temp_cel_wormcat <- gff_list_mrna[["elegans"]][which(gsub("[a-z]$", "", gff_list_mrna[["elegans"]]$short_name) %in% temp_cel),]$WormCat.3
  
  return(getMode(temp_cel_wormcat))
}, mc.cores = 16)

## Make a monocle object from the orthogroup thing
f_data <- data.frame(names(OGList[lapply(OGList, function(x) sum(x > 0)) > 0]))
colnames(f_data) <- "Orthogroup"
rownames(f_data) <- f_data$Orthogroup

f_data$cel_gene_count <- unlist(cel_length[lapply(OGList, function(x) sum(x > 0)) > 0])
f_data$cbr_gene_count <- unlist(cbr_length[lapply(OGList, function(x) sum(x > 0)) > 0])

f_data$WormCat1 <- unlist(OGList_WormCat1_fData[lapply(OGList, function(x) sum(x > 0)) > 0])
f_data$WormCat2 <- unlist(OGList_WormCat2_fData[lapply(OGList, function(x) sum(x > 0)) > 0])
f_data$WormCat3 <- unlist(OGList_WormCat3_fData[lapply(OGList, function(x) sum(x > 0)) > 0])

f_data$OGgenes <- OGList[lapply(OGList, function(x) sum(x > 0)) > 0]

og_cds <- new_cell_data_set(as(OGMatrix, "dgTMatrix"),
                            cell_metadata = colData(cds_split),
                            gene_metadata = f_data)

og_cds$gene_short_name <- og_cds$Orthogroup

rowData(og_cds)$species_gene <- NA
rowData(og_cds)[rowData(og_cds)$cel_gene_count > 0 & rowData(og_cds)$cbr_gene_count > 0,]$species_gene <- "common"
rowData(og_cds)[rowData(og_cds)$cel_gene_count > 0 & rowData(og_cds)$cbr_gene_count == 0,]$species_gene <- "ce.unique"
rowData(og_cds)[rowData(og_cds)$cel_gene_count == 0 & rowData(og_cds)$cbr_gene_count > 0,]$species_gene <- "cb.unique"

saveRDS(og_cds, paste0(dir, "Objects/cds_objects/og_cds.rds"))

############################
# Calculate orthogroup TPM #
############################

og_cds <- readRDS(paste0(dir, "Objects/cds_objects/og_cds.rds"))

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


# tpm terminal
TPMBinsList <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  TPMBinsList[[cur_species]] <- mclapply(names(BarcodeBinsList[[cur_species]]), function(cell_type) { 
    
    print(cell_type)
    temp_tpm_out <- data.frame(matrix(nrow = sum(rowData(og_cds)$species_gene == ifelse(cur_species == "C.elegans",
                                                                                        "ce.unique",
                                                                                        "cb.unique") |
                                                   rowData(og_cds)$species_gene == "common"), ncol = 0))
    rownames(temp_tpm_out) <- rownames(og_cds[rowData(og_cds)$species_gene == ifelse(cur_species == "C.elegans",
                                                                                     "ce.unique",
                                                                                     "cb.unique") |
                                                rowData(og_cds)$species_gene == "common",])
    
    for(i in seq(1, length(BarcodeBinsList[["C.elegans"]][[cell_type]]))) {
      temp_tpm <- data.frame(BinTPM(get.norm.expr.matrix(og_cds[rowData(og_cds)$species_gene == ifelse(cur_species == "C.elegans",
                                                                                                       "ce.unique",
                                                                                                       "cb.unique") |
                                                                  rowData(og_cds)$species_gene == "common",
                                                                BarcodeBinsList[[cur_species]][[cell_type]][[i]]])))
      temp_tpm_out <- cbind(temp_tpm_out, temp_tpm)
    }
    colnames(temp_tpm_out) <- names(BarcodeBinsList[["C.elegans"]][[cell_type]])
    return(temp_tpm_out)
  }, mc.cores = 42)
  TPMBinsList[[cur_species]] <- do.call("cbind", TPMBinsList[[cur_species]])
}

# tpm progenitor
TPMBinsListPro <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  TPMBinsListPro[[cur_species]] <- CalcTPM(og_cds[rowData(og_cds)$species_gene == ifelse(cur_species == "C.elegans",
                                                                                         "ce.unique",
                                                                                         "cb.unique") |
                                                    rowData(og_cds)$species_gene == "common",
                                                  colData(og_cds)$species == cur_species], "lineage_broad")
}

TPMBinsListPro[["C.elegans"]] <- TPMBinsListPro[["C.elegans"]][,sort(colnames(TPMBinsListPro[["C.briggsae"]]))]
TPMBinsListPro[["C.briggsae"]] <- TPMBinsListPro[["C.briggsae"]][,sort(colnames(TPMBinsListPro[["C.elegans"]]))]

TPMBinsListComp <- list()
TPMBinsListComp[["C.elegans"]] <- cbind(TPMBinsList[["C.elegans"]], TPMBinsListPro[["C.elegans"]])
TPMBinsListComp[["C.briggsae"]] <- cbind(TPMBinsList[["C.briggsae"]], TPMBinsListPro[["C.briggsae"]])

TPMBinsListComp[["C.elegans"]] <- TPMBinsListComp[["C.elegans"]][-which(names(TPMBinsListComp[["C.elegans"]]) == "unassigned"),]
TPMBinsListComp[["C.briggsae"]] <- TPMBinsListComp[["C.briggsae"]][-which(names(TPMBinsListComp[["C.briggsae"]]) == "unassigned"),]

saveRDS(TPMBinsList, paste0(dir, "Objects/TPMTimeBinsList_og.rds"))
saveRDS(TPMBinsListPro, paste0(dir, "Objects/TPMBinsListPro_og.rds"))
saveRDS(TPMBinsListComp, paste0(dir, "Objects/TPMTimeBinsListComp_og.rds"))

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

TPMBinsListMeanComp <- list()
TPMBinsListMeanComp[["C.elegans"]] <- cbind(TPMBinsListMean[["C.elegans"]], TPMBinsListPro[["C.elegans"]])
TPMBinsListMeanComp[["C.briggsae"]] <- cbind(TPMBinsListMean[["C.briggsae"]], TPMBinsListPro[["C.briggsae"]])

TPMBinsListMeanComp[["C.elegans"]] <- TPMBinsListMeanComp[["C.elegans"]][,which(colnames(TPMBinsListMeanComp[["C.elegans"]]) != "unassigned")]
TPMBinsListMeanComp[["C.briggsae"]] <- TPMBinsListMeanComp[["C.briggsae"]][,which(colnames(TPMBinsListMeanComp[["C.briggsae"]]) != "unassigned")]

saveRDS(TPMBinsListMean, paste0(dir, "Objects/TPMTimeMeanBinsList_og.rds"))
saveRDS(TPMBinsListMeanComp, paste0(dir, "Objects/TPMBinsListMeanComp_og.rds"))

TPMJoint <- TPMBinsListMeanComp
TPMJoint[["C.elegans"]] <- TPMJoint[["C.elegans"]][rownames(TPMJoint[["C.elegans"]]) %in% rownames(TPMJoint[["C.briggsae"]]),]
TPMJoint[["C.briggsae"]] <- TPMJoint[["C.briggsae"]][rownames(TPMJoint[["C.briggsae"]]) %in% rownames(TPMJoint[["C.elegans"]]),]

pseudocount = 1/dim(TPMJoint[["C.elegans"]])[2]
jsd_matrix_joint <- matrix(nrow = dim(TPMJoint[["C.elegans"]])[2],
                           ncol = dim(TPMJoint[["C.briggsae"]])[2])
rownames(jsd_matrix_joint) <- colnames(TPMJoint[["C.elegans"]])
colnames(jsd_matrix_joint) <- colnames(TPMJoint[["C.briggsae"]])
for(ele_cell in rownames(jsd_matrix_joint)) {
  for(bri_cell in colnames(jsd_matrix_joint)) {
    p = TPMJoint[["C.elegans"]][,ele_cell] + pseudocount
    q = TPMJoint[["C.briggsae"]][,bri_cell] + pseudocount
    jsd_matrix_joint[ele_cell, bri_cell] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
}

cell_data <- left_join(cell_data, data.frame(jsd_og = diag(jsd_matrix_joint),
                                cell_type = colnames(jsd_matrix_joint)))

# Load gene based TPM
TPMListBootstrap_term <- readRDS(paste0(dir, "Objects/TPMTimeBinsListBootstrap_CellCorrection.rds"))

TPMListBootstrap_term_subset <- list()
TPMListBootstrap_term_subset[["C.elegans"]] <- TPMListBootstrap_term[["C.elegans"]][gene_data$gene,]
TPMListBootstrap_term_subset[["C.briggsae"]] <- TPMListBootstrap_term[["C.briggsae"]][gene_data$gene,]

TPMBinsListMean <- readRDS(paste0(dir, "Objects/TPMTimeMeanBinsList.rds"))
TPMBinsListMeanSubset <- readRDS(paste0(dir, "Objects/TPMTimeMeanBinsListSubset.rds"))

TPMListBootstrap_pro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))
BinarizeBootstrapListPro <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListPro_CellCorrection.rds"))

TPMListBootstrap_pro_subset <- list()
TPMListBootstrap_pro_subset[["C.elegans"]] <- TPMListBootstrap_pro[["C.elegans"]][gene_data$gene,]
TPMListBootstrap_pro_subset[["C.briggsae"]] <- TPMListBootstrap_pro[["C.briggsae"]][gene_data$gene,]

TPMJoint <- list()
TPMJoint[["C.elegans"]] <- cbind(TPMBinsListMeanSubset[["C.elegans"]], TPMListBootstrap_pro_subset[["C.elegans"]])
TPMJoint[["C.briggsae"]] <- cbind(TPMBinsListMeanSubset[["C.briggsae"]], TPMListBootstrap_pro_subset[["C.briggsae"]])

pseudocount = 1/dim(TPMJoint[["C.elegans"]])[2]
jsd_matrix <- matrix(nrow = dim(TPMJoint[["C.elegans"]])[2],
                           ncol = dim(TPMJoint[["C.briggsae"]])[2])
rownames(jsd_matrix) <- colnames(TPMJoint[["C.elegans"]])
colnames(jsd_matrix) <- colnames(TPMJoint[["C.briggsae"]])
for(ele_cell in rownames(jsd_matrix)) {
  for(bri_cell in colnames(jsd_matrix)) {
    p = TPMJoint[["C.elegans"]][,ele_cell] + pseudocount
    q = TPMJoint[["C.briggsae"]][,bri_cell] + pseudocount
    jsd_matrix[ele_cell, bri_cell] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
}

cell_data <- left_join(cell_data, data.frame(jsd_raw = diag(jsd_matrix),
                                             cell_type = colnames(jsd_matrix)))

reg <- lm(formula = jsd_og ~ jsd_raw,
          data = cell_data)

#get intercept and slope value
coeff <- coefficients(reg)          
intercept <-coeff[1]
slope <- coeff[2]

pdf(paste0(dir, "Plots/ortho_plots/jsd.pdf"), width = 5, height = 5)
cell_data %>%
  mutate(plot = ifelse(cell_class == "progenitor", as.character(div_stage), cell_class)) %>%
  ggplot(aes(x = jsd_raw, y = jsd_og, color = plot, fill = plot, label = cell_type)) +
  geom_point(pch = 21, alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  annotate("text", x = c(-Inf), y = c(Inf), hjust = -0.1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype="dashed") +
  scale_x_continuous(name = "Gene based cell distance") +
  scale_y_continuous(name = "Orthogroup based cell distance") +
  scale_fill_manual(values = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                                    '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF",
                                                    'Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                                    'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                                    'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1")) + 
  scale_color_manual(limits = c('15', '28', '50', '100', '200', '350', "600",
                                'Ciliated neurons', 'Germline', 'Glia and excretory',
                                'Hypodermis and seam', 'Intestine', 'Muscle', 'Mesoderm',
                                'Non-ciliated neurons', 'Pharynx and rectal'), 
                     values = colorspace::darken(c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                                   '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF",
                                                   'Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                                   'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                                   'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1"), 0.1)) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.title.x = element_text(color = "#009E73", size = 18),
        axis.title.y = element_text(color = "#56B4E9", size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

# How many more genes incorporated in the orthogroup based cell distance
# 13679 total in the comparison set in gene_data

# 33433 total genes in the orthogroups

# 13559 in gene_data
# 19874 not in gene_data$elegans_id
unlist(lapply(rowData(og_cds)[rowData(og_cds)$species_gene == "common",]$OGgenes, function(x) {
  ! unlist(strsplit(x, ",")) %in% gene_data$elegans_id
})) %>% sum()

# 13539 in gene_data
# 19894 not in gene_data$briggsae_id
unlist(lapply(rowData(og_cds)[rowData(og_cds)$species_gene == "common",]$OGgenes, function(x) {
  ! unlist(strsplit(x, ",")) %in% gene_data$briggsae_id
})) %>% sum()

# 6335 new genes incorporated in the orthogroup based cell distance
unlist(lapply(rowData(og_cds)[rowData(og_cds)$species_gene == "common",]$OGgenes, function(x) {
  ! (unlist(strsplit(x, ",")) %in% gene_data$elegans_id | unlist(strsplit(x, ",")) %in% gene_data$briggsae_id)
})) %>% sum()

# 3382 new C. elegans genes incorporated in the orthogroup based cell distance
sum(unlist(lapply(rowData(og_cds)[rowData(og_cds)$species_gene == "common",]$OGgenes, function(x) {
  unlist(strsplit(x, ","))[! (unlist(strsplit(x, ",")) %in% gene_data$elegans_id | unlist(strsplit(x, ",")) %in% gene_data$briggsae_id)]
})) %in% gff_list_mrna[["elegans"]]$elegans_id)

sum(unlist(lapply(rowData(og_cds)[rowData(og_cds)$species_gene == "common",]$OGgenes, function(x) {
  unlist(strsplit(x, ","))[! (unlist(strsplit(x, ",")) %in% gene_data$elegans_id | unlist(strsplit(x, ",")) %in% gene_data$briggsae_id)]
})) %in% gff_list_mrna[["briggsae"]]$gene_name)
