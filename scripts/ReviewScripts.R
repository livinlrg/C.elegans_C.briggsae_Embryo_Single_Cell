

library(Seurat)
library(monocle3)
library(parallel)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(cowplot)
library(gridExtra)
library(ggtext)

source("/kimdata/livinlrg/scAnalysis/Scripts/WS290/CalcFunctions.R")

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"
cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_cell_bg.rds"))
cell_data_mean <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))

# load cds_bg_corrected
cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240829.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

rm(cds)

tpm_matrix_list <- readRDS(paste0(dir, "Objects/tpm_matrix_time_list_filt_cell_bg.rds"))

# load barcodes for terminal time bins
BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))
ProBarcodeBinsList <- readRDS(paste0(dir, "Objects/ProBarcodeBinsList.rds"))
JointBarcodeBinsList <- readRDS(paste0(dir, "Objects/JointBarcodeBinsList.rds"))

TPMListBootstrap_term <- readRDS(paste0(dir, "Objects/TPMTimeBinsListBootstrap_CellCorrection.rds"))
TPMListBootstrapMean_term <- readRDS(paste0(dir, "Objects/TPMListBootstrapMean_CellCorrection.rds"))
TPMListBootstrap_pro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))

TPMListJoint <- list()
TPMListJoint[["C.elegans"]] <- cbind(TPMListBootstrapMean_term[["C.elegans"]], TPMListBootstrap_pro[["C.elegans"]])
TPMListJoint[["C.briggsae"]] <- cbind(TPMListBootstrapMean_term[["C.briggsae"]], TPMListBootstrap_pro[["C.briggsae"]])

cell_type_time_bins <- readRDS(paste0(dir, "Objects/cell_type_time_bins.rds"))
CellTable <- readRDS(paste0(dir, "Objects/CellTable_20240826.rds"))

time_bin_vector <- c("lt_100", "100_130", "130_170", "170_210", "210_270", "270_330", "330_390", "390_450", "450_510", "510_580", "580_650", "650_710", "gt_710")

time_bin_df <- data.frame(bins = time_bin_vector,
                          start = c(0, 100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710),
                          end = c(100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710, 1000))

# tpm terminal
CompBinsList <- list()
for(cur_species in c("C.elegans")) {
  CompBinsList[[cur_species]] <- mclapply(names(BarcodeBinsList[[cur_species]]), function(cell_type) { 
    print(cell_type)
    
    temp_comp_out <- data.frame(matrix(nrow = 0, ncol = 3))
    colnames(temp_comp_out) <- c("Var1", "Freq", "cell_type_bin")
    
    for(i in seq(1, length(BarcodeBinsList[["C.elegans"]][[cell_type]]))) {
      temp_comp <- data.frame(table(colData(cds_filt[rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                                                            "ce.unique",
                                                                                            "cb.unique") |
                                                        rowData(cds_filt)$gene.type == "common",
                                                      BarcodeBinsList[[cur_species]][[cell_type]][[i]]])$genotype) / length(BarcodeBinsList[[cur_species]][[cell_type]][[i]]) * 100)
      temp_comp$cell_type_bin <- names(BarcodeBinsList[["C.elegans"]][[cell_type]])[i]
      
      temp_comp_out <- rbind(temp_comp_out, temp_comp)
    }
    return(temp_comp_out)
  }, mc.cores = 32)
  CompBinsList[[cur_species]] <- do.call("rbind", CompBinsList[[cur_species]])
}

CompBinsList <- lapply(CompBinsList, function(x) {
  return(x[! x$cell_type_bin %in% c("hyp3_390_450:580_650"),])
})

pdf(paste0(dir, "Plots/test.pdf"), width = 30, height = 5)
CompBinsList[["C.elegans"]] %>%
  mutate(cell_class = factor(cell_data[CompBinsList[["C.elegans"]]$cell_type_bin, "cell_class"], levels = unique(cell_data$cell_class))) %>%
  arrange(cell_class, cell_type_bin) %>%
  mutate(cell_type_bin = factor(cell_type_bin, levels = unique(cell_type_bin))) %>%
ggplot() +
  geom_bar(aes(x = cell_type_bin, y = Freq, fill = Var1, group = Var1), stat = "identity") +
  scale_y_continuous(name = "Percent of cell type count", expand = c(0, 0)) +
  scale_fill_manual(breaks = c("wt", "mec-3", "M03D44", "ceh-9"), labels = c("Wildtype", "mec-3", "M03D4.4", "ceh-9"),
                    values = c("#56B4E9", "#D55E00", "#CC79A7", "#009E73")) +
  theme(panel.background = element_rect(fill = 'white'),
        plot.background = element_rect(fill = 'white', color = NA),
        axis.line = element_line(color = "grey80", size = 1),
        axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = 6.5, angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.y = element_line(color = "grey80", size = 0.25),
        panel.grid.minor.y = element_line(color = "grey80", size = 0.05),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.title = element_blank(),
        legend.box.background = element_rect(colour = "white", fill = "white"),
        legend.key = element_rect(colour = "white", fill = "white"),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal")
dev.off()

################# Per Genotype TPM #################

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
TPMGenoList <- list()
for(cur_geno in c("ceh-9", "M03D44", "mec-3", "wt", "all")) {
  if(cur_geno == "all") {
    TPMGenoList[[cur_geno]] <- mclapply(names(BarcodeBinsList[["C.elegans"]]), function(cell_type) { 
      print(cell_type)
      temp_tpm_out <- data.frame(matrix(nrow = sum(rowData(cds_filt)$gene.type == "ce.unique" | rowData(cds_filt)$gene.type == "common"), ncol = 0))
      rownames(temp_tpm_out) <- rownames(cds_filt)[rowData(cds_filt)$gene.type == "ce.unique" | rowData(cds_filt)$gene.type == "common"]
      
      for(i in seq(1, length(BarcodeBinsList[["C.elegans"]][[cell_type]]))) {
        temp_tpm <- data.frame(BinTPM(get.norm.expr.matrix(cds_filt[rowData(cds_filt)$gene.type == "ce.unique" | rowData(cds_filt)$gene.type == "common",
                                                                    colnames(cds_filt) %in% BarcodeBinsList[["C.elegans"]][[cell_type]][[i]]])))
        temp_tpm_out <- cbind(temp_tpm_out, temp_tpm)
      }
      colnames(temp_tpm_out) <- names(BarcodeBinsList[["C.elegans"]][[cell_type]])
      return(temp_tpm_out)
    }, mc.cores = 32)
    TPMGenoList[[cur_geno]] <- do.call("cbind", TPMGenoList[[cur_geno]])
  } else {
    TPMGenoList[[cur_geno]] <- mclapply(names(BarcodeBinsList[["C.elegans"]]), function(cell_type) { 
      print(cell_type)
      temp_tpm_out <- data.frame(matrix(nrow = sum(rowData(cds_filt)$gene.type == "ce.unique" | rowData(cds_filt)$gene.type == "common"), ncol = 0))
      rownames(temp_tpm_out) <- rownames(cds_filt)[rowData(cds_filt)$gene.type == "ce.unique" | rowData(cds_filt)$gene.type == "common"]
      
      for(i in seq(1, length(BarcodeBinsList[["C.elegans"]][[cell_type]]))) {
        temp_tpm <- data.frame(BinTPM(get.norm.expr.matrix(cds_filt[rowData(cds_filt)$gene.type == "ce.unique" | rowData(cds_filt)$gene.type == "common",
                                                                    colnames(cds_filt) %in% BarcodeBinsList[["C.elegans"]][[cell_type]][[i]] &
                                                                      colData(cds_filt)$genotype %in% cur_geno])))
        temp_tpm_out <- cbind(temp_tpm_out, temp_tpm)
      }
      colnames(temp_tpm_out) <- names(BarcodeBinsList[["C.elegans"]][[cell_type]])
      return(temp_tpm_out)
    }, mc.cores = 32)
    TPMGenoList[[cur_geno]] <- do.call("cbind", TPMGenoList[[cur_geno]])
  }
}

# create a version of the tpm time bins that is just one value per bin based on taking the mean of the values
TPMGenoListMean <- list()
for(cur_geno in c("ceh-9", "M03D44", "mec-3", "wt", "all")) {
  TPMGenoListMean[[cur_geno]] <- matrix(nrow = nrow(TPMGenoList[[cur_geno]]), ncol = length(unique(cell_type_time_bins$cell_type)))
  colnames(TPMGenoListMean[[cur_geno]]) <- unique(cell_type_time_bins$cell_type)
  rownames(TPMGenoListMean[[cur_geno]]) <- rownames(TPMGenoList[[cur_geno]])
  for(cell_type in unique(cell_type_time_bins[cell_type_time_bins$cell_type != "hyp3",]$cell_type)) {
    if(nrow(cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]) > 1) {
      TPMGenoListMean[[cur_geno]][, cell_type] <- rowMeans(TPMGenoList[[cur_geno]][, paste0(cell_type, "_", cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$start_bin, ":",
                                                                                                  cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$end_bin)], na.rm = TRUE)
    } else {
      TPMGenoListMean[[cur_geno]][, cell_type] <- TPMGenoList[[cur_geno]][, paste0(cell_type, "_", cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$start_bin, ":",
                                                                                         cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$end_bin)]
    }
  }
}

TPMGenoList <- lapply(TPMGenoList, function(x) {
  return(x[,! colnames(x) %in% c("hyp3_390_450:580_650")])
})

TPMGenoListMean <- lapply(TPMGenoListMean, function(x) {
  return(x[,! colnames(x) %in% c("hyp3")])
})

shared_genes <- intersect(rownames(TPMGenoList[["wt"]]), rownames(TPMListBootstrap_term[["C.briggsae"]]))

geno_jsd_list <- list()
for(cur_geno in c("ceh-9", "M03D44", "mec-3", "wt", "all")) {
  print(cur_geno)
  
  pseudocount = 1/dim(TPMGenoList[[cur_geno]][, ! is.na(colSums(TPMGenoList[[cur_geno]]))])[1]
  geno_jsd_list[[cur_geno]] <- data.frame(matrix(nrow = dim(TPMGenoList[[cur_geno]][, ! is.na(colSums(TPMGenoList[[cur_geno]]))])[2],
                                                 ncol = 2))
  colnames(geno_jsd_list[[cur_geno]]) <- c("cell_type", "jsd")
  geno_jsd_list[[cur_geno]]$cell_type <- colnames(TPMGenoList[[cur_geno]][, ! is.na(colSums(TPMGenoList[[cur_geno]]))])
  rownames(geno_jsd_list[[cur_geno]]) <- geno_jsd_list[[cur_geno]]$cell_type
  
  for(cur_cell in rownames(geno_jsd_list[[cur_geno]])) {
    p = TPMGenoList[[cur_geno]][shared_genes, cur_cell] + pseudocount
    q = TPMListBootstrap_term[["C.briggsae"]][shared_genes, cur_cell] + pseudocount
    geno_jsd_list[[cur_geno]][cur_cell,"jsd"] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
}

cell_geno_data <- data.frame(cell_type = colnames(TPMListBootstrap_term[["C.briggsae"]]))
rownames(cell_geno_data) <- cell_geno_data$cell_type
for(cur_geno in c("ceh-9", "M03D44", "mec-3", "wt", "all")) {
  print(cur_geno)
  cell_geno_data[,cur_geno] <- NA
  cell_geno_data[geno_jsd_list[[cur_geno]]$cell_type, cur_geno] <- geno_jsd_list[[cur_geno]]$jsd
}

cell_geno_data$cell_class <- cell_data[cell_geno_data$cell_type, "cell_class"]
cell_geno_data$cell_type_mean <- cell_data[cell_geno_data$cell_type, "cell_type"]

# cell count terminal
for(cur_geno in c("ceh-9", "M03D44", "mec-3", "wt")) {
  cell_geno_data[,paste0(cur_geno, "_cell_count")] <- NA
  for(cell_type in names(BarcodeBinsList[["C.elegans"]])) { 
    for(i in seq(1, length(BarcodeBinsList[["C.elegans"]][[cell_type]]))) {
      cell_geno_data[names(BarcodeBinsList[["C.elegans"]][[cell_type]])[i], paste0(cur_geno, "_cell_count")] <- sum(colnames(cds_filt) %in% BarcodeBinsList[["C.elegans"]][[cell_type]][[i]] &
                                                                                                                      colData(cds_filt)$genotype == cur_geno)
    }
  }
}

cell_geno_data$all_cell_count <- rowSums(cell_geno_data[,c("ceh-9_cell_count", "M03D44_cell_count", "mec-3_cell_count", "wt_cell_count")])


cols = list(cell_class = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                           'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                           'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1"))

data_columns = c("ceh-9", "M03D44", "mec-3", "wt", "all")

cell_geno_data <- cell_geno_data[! is.na(cell_geno_data$cell_type),]

plots <- list()
plot_index <- 1

for (i in 1:length(data_columns)) {
  for (j in 1:length(data_columns)) {
    # Dynamically reference column names
    x_col <- sym(data_columns[j])
    y_col <- sym(data_columns[i])
    
    #ceh-9_cell_count
    cur_cell_geno_data <- cell_geno_data[cell_geno_data[,paste0(data_columns[j], "_cell_count")] > 20 & cell_geno_data[,paste0(data_columns[i], "_cell_count"),] > 20,]
    
    rank_sum_test <- p.adjust(wilcox.test(cur_cell_geno_data[,data_columns[j]],
                                          cur_cell_geno_data[,data_columns[i]],
                                          alternative = "two.sided")$p.value, method = "BH", n = 12)
    lm_results <- summary(lm(cur_cell_geno_data[,data_columns[j]] ~ cur_cell_geno_data[,data_columns[i]]))$r.squared
    print(paste0("x: ", data_columns[j], " y: ", data_columns[i], " p-value: ", rank_sum_test, " lm: ", lm_results))

    if(i == j) {
      p <- ggplot(cur_cell_geno_data, aes(x = forcats::fct_reorder(cell_class, !!y_col), y = !!y_col, color = cell_class)) +
        geom_boxplot() +
        scale_color_manual(values = cols$cell_class) +
        theme_minimal() +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
    } else {
      p <- ggplot(cur_cell_geno_data, aes(x = !!x_col, y = !!y_col, color = cell_class, label = cell_type)) +
        geom_point(alpha = 0.7) +
        geom_text_repel() +
        scale_color_manual(values = cols$cell_class) +
        theme_minimal() +
        theme(legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1))
    }
    
    
    plots[[plot_index]] <- p
    plot_index <- plot_index + 1
  }
}

# Arrange plots in a grid using cowplot
pdf(paste0(dir, "test.pdf"), width = 12, height = 12)
plot_grid(plotlist = plots, ncol = 5)
dev.off()


## Mean version ##
cell_count_limit <- 25

cell_geno_data_mean <- data.frame(cell_type = unique(cell_geno_data$cell_type_mean),
                                  cell_class = cell_data_mean[unique(cell_geno_data$cell_type_mean), "cell_class"])
rownames(cell_geno_data_mean) <- cell_geno_data_mean$cell_type
for(cur_cell_type in unique(cell_geno_data$cell_type_mean)) {
  for(cur_geno in c("ceh-9", "M03D44", "mec-3", "wt", "all")) {
    temp_geno_data <- cell_geno_data[cell_geno_data$cell_type_mean == cur_cell_type &
                                       cell_geno_data[,paste0(cur_geno, "_cell_count")] > cell_count_limit,
                                     cur_geno]
    
    cell_geno_data_mean[cur_cell_type, cur_geno] <- mean(temp_geno_data, na.rm = TRUE)
    
    cell_geno_data_mean[cur_cell_type, paste0(cur_geno, "_cell_count")] <- min(cell_geno_data[cell_geno_data$cell_type_mean == cur_cell_type &
                                                                                                cell_geno_data[,paste0(cur_geno, "_cell_count")] > cell_count_limit,
                                                                                              paste0(cur_geno, "_cell_count")],
                                                                               na.rm = TRUE)
    cell_geno_data_mean[cur_cell_type, paste0(cur_geno, "_cell_count")] <- ifelse(is.infinite(cell_geno_data_mean[cur_cell_type, paste0(cur_geno, "_cell_count")]), NA, cell_geno_data_mean[cur_cell_type, paste0(cur_geno, "_cell_count")])
  }
}

plots <- list()
plot_index <- 1

for (i in 1:length(data_columns)) {
  for (j in 1:length(data_columns)) {
    # Dynamically reference column names
    x_col <- sym(data_columns[j])
    y_col <- sym(data_columns[i])
    
    cur_cell_geno_data <- cell_geno_data_mean
    cur_cell_geno_data$min_cell_count <- apply(cur_cell_geno_data[,c(paste0(x_col, "_cell_count"), paste0(y_col, "_cell_count"))], 1, min)
    
    cor_test <- cor(cur_cell_geno_data[,data_columns[j]], cur_cell_geno_data[,data_columns[i]], use = "pairwise.complete.obs")
    lm_results <- summary(lm(cur_cell_geno_data[,data_columns[j]] ~ cur_cell_geno_data[,data_columns[i]]))$r.squared
    print(paste0("x: ", data_columns[j], " y: ", data_columns[i], " p-value: ", rank_sum_test, " lm: ", lm_results, " cor: ", cor_test))
    
    if(i == j) {
      cur_cell_geno_data <- cur_cell_geno_data[! is.na(cur_cell_geno_data[[as.character(y_col)]]), ]
      p <- ggplot(cur_cell_geno_data, aes(x = forcats::fct_reorder(cell_class, !!y_col, .na_rm = TRUE), y = !!y_col, color = cell_class)) +
        geom_boxplot() +
        scale_color_manual(values = cols$cell_class) +
        theme_minimal() +
        theme(legend.position = "none",
              axis.title.x = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1))
    } else if (i > j) {
      p <- ggplot(cur_cell_geno_data, aes(x = !!x_col, y = !!y_col, color = cell_class, label = cell_type)) +
        geom_point(alpha = 0.7) +
        geom_text_repel() +
        scale_color_manual(values = cols$cell_class) +
        theme_minimal() +
        theme(legend.position = "none")
    } else if (i < j) {
      ## annotate text from tests
      text_label <- paste0(x_col, " vs. ", y_col,
                           # "\n mean cell count = ", round(mean(cur_cell_geno_data$min_cell_count, na.rm = TRUE), 0),
                           # "\n p = ", formatC(rank_sum_test, format = "e", digits = 2),
                           "\nR² = ", round(lm_results, 2),
                           "\nCor = ", round(cor_test, 2))

      p <- ggplot() +
        annotate("label", x = 0.5, y = 0.5, label = text_label, size = 5) +
        theme_void()
    }
    plots[[plot_index]] <- p
    plot_index <- plot_index + 1
  }
}

# Arrange plots in a grid using cowplot
pdf(paste0(dir, "test2.pdf"), width = 12, height = 12)
plot_grid(plotlist = plots, ncol = 5)
dev.off()

embryo_time_temp <- colData(cds_filt)[colData(cds_filt)$species != "C.briggsae",c("smoothed.embryo.time", "genotype")]
p <- ggplot(embryo_time_temp, aes(x = smoothed.embryo.time, color = genotype)) +
  geom_density() +
  # scale_color_manual() +
  theme_minimal() +
  theme(legend.position = "top")

pdf(paste0(dir, "test3.pdf"), width = 5, height = 2.5)
print(p)
dev.off()

#### Gene based ####
geno_gene_jsd_list <- list()
# use cell_geno_data to help filter some cell types.
cell_count_limit <- 25

for(cur_geno in c("ceh-9", "M03D44", "mec-3", "wt", "all")) {
  print(cur_geno)
  
  cells_to_use <- cell_geno_data[cell_geno_data[,paste0(cur_geno, "_cell_count")] > cell_count_limit, "cell_type"]
  
  pseudocount = 1/length(cells_to_use)
  
  geno_gene_jsd_list[[cur_geno]] <- data.frame(matrix(nrow = length(shared_genes), ncol = 2))
  colnames(geno_gene_jsd_list[[cur_geno]]) <- c("gene", "jsd")
  
  geno_gene_jsd_list[[cur_geno]]$gene <- shared_genes
  rownames(geno_gene_jsd_list[[cur_geno]]) <- geno_gene_jsd_list[[cur_geno]]$gene
  
  for(cur_gene in rownames(geno_gene_jsd_list[[cur_geno]])) {
    p = unlist(TPMGenoList[[cur_geno]][cur_gene, cells_to_use] + pseudocount)
    q = unlist(TPMListBootstrap_term[["C.briggsae"]][cur_gene, cells_to_use] + pseudocount)
    geno_gene_jsd_list[[cur_geno]][cur_gene,"jsd"] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
}

geno_gene_data <- left_join(geno_gene_jsd_list[["ceh-9"]], geno_gene_jsd_list[["M03D44"]], by = "gene") %>%
  left_join(geno_gene_jsd_list[["mec-3"]], by = "gene") %>%
  left_join(geno_gene_jsd_list[["wt"]], by = "gene") %>%
  left_join(geno_gene_jsd_list[["all"]], by = "gene")

colnames(geno_gene_data) <- c("gene", "ceh-9", "M03D44", "mec-3", "wt", "all")
rownames(geno_gene_data) <- geno_gene_data$gene

for(cur_geno in c("ceh-9", "M03D44", "mec-3", "wt", "all")) {
  print(cur_geno)
  
  cells_to_use <- cell_geno_data[cell_geno_data[,paste0(cur_geno, "_cell_count")] > cell_count_limit, "cell_type"]
  for(cur_gene in rownames(geno_gene_data)) {
    geno_gene_data[cur_gene, paste0(cur_geno, "_cel_max_exp")] <- max(TPMGenoList[[cur_geno]][cur_gene, cells_to_use], na.rm = TRUE)
    geno_gene_data[cur_gene, paste0(cur_geno, "_cbr_max_exp")] <- max(TPMListBootstrap_term[["C.briggsae"]][cur_gene, cells_to_use], na.rm = TRUE)
  }
}


plots <- list()
plot_index <- 1

for (i in 1:length(data_columns)) {
  for (j in 1:length(data_columns)) {
    # Dynamically reference column names
    x_col <- sym(data_columns[j])
    y_col <- sym(data_columns[i])
    
    cur_gene_geno_data <- geno_gene_data[geno_gene_data[,paste0(x_col, "_cel_max_exp")] > 80 |
                                           geno_gene_data[,paste0(x_col, "_cbr_max_exp")] > 80 &
                                           geno_gene_data[,paste0(y_col, "_cel_max_exp")] > 80 |
                                           geno_gene_data[,paste0(y_col, "_cbr_max_exp")] > 80,]
    # cur_gene_geno_data <- geno_gene_data
    
    cor_test <- cor(cur_gene_geno_data[,data_columns[j]], cur_gene_geno_data[,data_columns[i]], use = "pairwise.complete.obs")
    lm_results <- summary(lm(cur_gene_geno_data[,data_columns[j]] ~ cur_gene_geno_data[,data_columns[i]]))$r.squared
    print(paste0("x: ", data_columns[j], " y: ", data_columns[i], " p-value: ", rank_sum_test, " lm: ", lm_results, " cor: ", cor_test))
    
    if(i == j) {
      p <- ggplot(cur_gene_geno_data, aes(x = !!y_col)) +
        geom_density() +
        theme_minimal() +
        theme(legend.position = "none")
    } else if (i > j) {
      p <- ggplot(cur_gene_geno_data, aes(x = !!x_col, y = !!y_col)) +
        geom_point(alpha = 0.3, size = 0.25) +
        theme_minimal() +
        theme(legend.position = "none")
    } else if (i < j) {
      ## annotate text from tests
      text_label <- paste0(x_col, " vs. ", y_col,
                           "\nR² = ", round(lm_results, 2),
                           "\nCor = ", round(cor_test, 2))
      
      p <- ggplot() +
        annotate("label", x = 0.5, y = 0.5, label = text_label, size = 5) +
        theme_void()
    }
    plots[[plot_index]] <- p
    plot_index <- plot_index + 1
  }
}

# Arrange plots in a grid using cowplot
pdf(paste0(dir, "test4.pdf"), width = 12, height = 12)
plot_grid(plotlist = plots, ncol = 5)
dev.off()



gene_data <- readRDS(paste0(dir, "Objects/gene_data_joint_cell_bg.rds"))
geno_gene_data <- left_join(geno_gene_data, gene_data[, c("gene", "cel_tau_median_term", "cbr_tau_median_term")], by = "gene")

geno_gene_data$tau <- apply(geno_gene_data[,c("cel_tau_median_term", "cbr_tau_median_term")], 1, max, na.rm = TRUE)

geno_gene_data %>% filter(all_cel_max_exp > 100 | all_cbr_max_exp > 100) %>%
  mutate(selector = log2(all/M03D44)) %>%
  filter(selector > 1 | selector < -1)


###### 2x gene based ######
geno_gene_jsd_filt_list <- list()
# use cell_geno_data to help filter some cell types.
cell_count_limit <- 25

for(cur_geno_x in c("ceh-9", "M03D44", "mec-3", "wt", "all")) {
  for(cur_geno_y in c("ceh-9", "M03D44", "mec-3", "wt", "all")) {
    print(paste0(cur_geno_x, " vs. ", cur_geno_y))
    
    cells_to_use_x <- cell_geno_data[cell_geno_data[,paste0(cur_geno_x, "_cell_count")] > cell_count_limit, "cell_type"]
    cells_to_use_y <- cell_geno_data[cell_geno_data[,paste0(cur_geno_y, "_cell_count")] > cell_count_limit, "cell_type"]
    
    cells_to_use <- intersect(cells_to_use_x, cells_to_use_y)
    pseudocount = 1/length(cells_to_use)
    
    geno_gene_jsd_filt_list[[paste0(cur_geno_x, "_", cur_geno_y)]] <- data.frame(matrix(nrow = length(shared_genes), ncol = 1))
    colnames(geno_gene_jsd_filt_list[[paste0(cur_geno_x, "_", cur_geno_y)]]) <- c(paste0(cur_geno_x, "_", cur_geno_y))
    rownames(geno_gene_jsd_filt_list[[paste0(cur_geno_x, "_", cur_geno_y)]]) <- shared_genes
    
    geno_gene_jsd_filt_list[[paste0(cur_geno_x, "_", cur_geno_y)]][, paste0(cur_geno_x, "_", cur_geno_y)] <- unlist(lapply(shared_genes, function(cur_gene) {
      p = unlist(TPMGenoList[[cur_geno_x]][cur_gene, cells_to_use] + pseudocount)
      q = unlist(TPMListBootstrap_term[["C.briggsae"]][cur_gene, cells_to_use] + pseudocount)
      return(sqrt(js_divg(p = p/sum(p), q = q/sum(q))))
    }))
  }
}

geno_gene_filt_data <- do.call(cbind, geno_gene_jsd_filt_list)
geno_gene_filt_data$gene <- rownames(geno_gene_filt_data)

for(cur_geno_x in c("ceh-9", "M03D44", "mec-3", "wt", "all")) {
  for(cur_geno_y in c("ceh-9", "M03D44", "mec-3", "wt", "all")) {
    print(paste0(cur_geno_x, " vs. ", cur_geno_y))
    
    cells_to_use_x <- cell_geno_data[cell_geno_data[,paste0(cur_geno_x, "_cell_count")] > cell_count_limit, "cell_type"]
    cells_to_use_y <- cell_geno_data[cell_geno_data[,paste0(cur_geno_y, "_cell_count")] > cell_count_limit, "cell_type"]
    
    cells_to_use <- intersect(cells_to_use_x, cells_to_use_y)
    
    for(cur_gene in rownames(geno_gene_filt_data)) {
      geno_gene_filt_data[cur_gene, paste0(cur_geno_x, "_", cur_geno_y , "_cel_max_exp")] <- min(max(TPMGenoList[[cur_geno_x]][cur_gene, cells_to_use], na.rm = TRUE),
                                                                                                 max(TPMGenoList[[cur_geno_y]][cur_gene, cells_to_use], na.rm = TRUE), na.rm = TRUE)
      geno_gene_filt_data[cur_gene, paste0(cur_geno_x, "_", cur_geno_y , "_cbr_max_exp")] <- max(TPMListBootstrap_term[["C.briggsae"]][cur_gene, cells_to_use], na.rm = TRUE)
    }
  }
}

plots <- list()
plot_index <- 1

for (i in 1:length(data_columns)) {
  for (j in 1:length(data_columns)) {
    # Dynamically reference column names
    x_col <- sym(paste0(data_columns[j], "_", data_columns[i]))
    y_col <- sym(paste0(data_columns[i], "_", data_columns[j]))
    
    cur_gene_geno_data <- geno_gene_filt_data[geno_gene_filt_data[,paste0(x_col, "_cel_max_exp")] > 80 |
                                                geno_gene_filt_data[,paste0(x_col, "_cbr_max_exp")] > 80 &
                                                geno_gene_filt_data[,paste0(y_col, "_cel_max_exp")] > 80 |
                                                geno_gene_filt_data[,paste0(y_col, "_cbr_max_exp")] > 80,]

    cor_test <- cor(cur_gene_geno_data[,as.character(y_col)], cur_gene_geno_data[,as.character(x_col)], use = "pairwise.complete.obs")
    lm_results <- summary(lm(cur_gene_geno_data[,as.character(y_col)] ~ cur_gene_geno_data[,as.character(x_col)]))$r.squared
    print(paste0("x: ", data_columns[j], " y: ", data_columns[i], " p-value: ", rank_sum_test, " lm: ", lm_results, " cor: ", cor_test))
    
    if(i == j) {
      p <- ggplot(cur_gene_geno_data, aes(x = !!y_col)) +
        geom_density() +
        scale_x_continuous(name = data_columns[j]) +
        theme_minimal() +
        theme(legend.position = "none")
    } else if (i > j) {
      p <- ggplot(cur_gene_geno_data, aes(x = !!x_col, y = !!y_col)) +
        geom_point(alpha = 0.3, size = 0.25) +
        scale_x_continuous(name = data_columns[j]) +
        scale_y_continuous(name = data_columns[i]) +
        theme_minimal() +
        theme(legend.position = "none")
    } else if (i < j) {
      ## annotate text from tests
      text_label <- paste0(data_columns[j], " vs. ", data_columns[i],
                           "\nR² = ", round(lm_results, 2),
                           "\nCor = ", round(cor_test, 2))
      
      p <- ggplot() +
        annotate("label", x = 0.5, y = 0.5, label = text_label, size = 5) +
        theme_void()
    }
    plots[[plot_index]] <- p
    plot_index <- plot_index + 1
  }
}

# Arrange plots in a grid using cowplot
pdf(paste0(dir, "test5.pdf"), width = 12, height = 12)
plot_grid(plotlist = plots, ncol = 5)
dev.off()

geno_gene_filt_data[,"wt_all_cel_max_exp"] > 80 | geno_gene_filt_data[,"wt_all_cbr_max_exp"] > 80
geno_gene_filt_data[,"all_wt_cel_max_exp"] > 80 | geno_gene_filt_data[,"all_wt_cbr_max_exp"] > 80


rownames(geno_gene_data) <- geno_gene_data$gene

geno_gene_data %>%
  filter(all_cel_max_exp > 80 | all_cbr_max_exp > 80) %>%
  filter(wt_cel_max_exp > 80 | wt_cbr_max_exp > 80) %>%
  mutate(diff = log2(all / wt)) %>%
  filter(diff > 0.75 | diff < -0.75) %>%
  arrange(diff) %>%
  select(gene, diff, all, wt, tau, all_cel_max_exp, wt_cel_max_exp, all_cbr_max_exp, wt_cbr_max_exp)

# Wild type versus all
p1 <- geno_gene_data %>%
  filter(all_cel_max_exp > 80 | all_cbr_max_exp > 80) %>%
  filter(wt_cel_max_exp > 80 | wt_cbr_max_exp > 80) %>%
  ggplot(aes(x = all, y = wt)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_x_continuous(name = "Gene Distance (All datasets)") +
  scale_y_continuous(name = "Gene Distance (Wild type)") +
  theme(legend.position = "None",
    axis.text = element_text(size = 12),
    rect = element_rect(fill = "transparent"),
    axis.line = element_line(color="grey80", size=1),
    panel.grid.major = element_line(color="grey80", size=0.25),
    panel.grid.minor = element_line(color="grey80", size=0.05),
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA))

genes_to_plot <- c("pha-4", "hlh-4", "eef-2", "rab-7", "T04A6.1", "vab-7")

p2 <- geno_gene_data %>%
  filter(all_cel_max_exp > 80 | all_cbr_max_exp > 80) %>%
  ggplot(aes(x = all)) +
  geom_density() +
  scale_x_continuous(name = "Gene Distance (All datasets)") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[1], "all"], linetype = "solid", color = "red") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[2], "all"], linetype = "solid", color = "red") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[3], "all"], linetype = "solid", color = "red") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[4], "all"], linetype = "solid", color = "red") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[5], "all"], linetype = "solid", color = "red") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[6], "all"], linetype = "solid", color = "red") +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "none")

p3 <- geno_gene_data %>%
  filter(wt_cel_max_exp > 80 | wt_cbr_max_exp > 80) %>%
  ggplot(aes(x = wt)) +
  geom_density() +
  scale_x_continuous(name = "Gene Distance (Wildtype)") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[1], "wt"], linetype = "solid", color = "red") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[2], "wt"], linetype = "solid", color = "red") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[3], "wt"], linetype = "solid", color = "red") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[4], "wt"], linetype = "solid", color = "red") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[5], "wt"], linetype = "solid", color = "red") +
  geom_vline(xintercept = geno_gene_data[genes_to_plot[6], "wt"], linetype = "solid", color = "red") +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "none")
  
pdf(paste0(dir, "test6.pdf"), width = 12, height = 4)
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

p4_lm <- round(summary(lm(geno_gene_data$all ~ geno_gene_data$wt))$r.squared, 2)
p4_cor <- round(cor(geno_gene_data$all, geno_gene_data$wt), 2)

p4 <- cell_geno_data_mean %>%
  filter(all_cell_count > 20) %>%
  ggplot(aes(x = all, y = wt, color = cell_class, label = cell_type)) +
  geom_point(alpha = 0.7) +
  scale_x_continuous(name = "Cell Distance (All datasets)", limits = c(0.325, 0.58)) +
  scale_y_continuous(name = "Cell Distance (Wild type)", limits = c(0.325, 0.58)) +
  scale_color_manual(values = cols$cell_class) +
  theme(legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))

p5 <- cell_geno_data_mean %>%
  filter(all_cell_count > 20) %>%
  ggplot(aes(x = forcats::fct_reorder(cell_class, all, .na_rm = TRUE), y = all, color = cell_class)) +
  geom_boxplot() +
  scale_y_continuous(name = "Cell Distance (All datasets)",
                     limits = c(0.325, 0.58)) +
  scale_color_manual(values = cols$cell_class) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major.y = element_line(color="grey80", size=0.25),
        panel.grid.minor.y = element_line(color="grey80", size=0.05),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

p6 <- cell_geno_data_mean %>%
  filter(wt_cell_count > 20) %>%
  ggplot(aes(x = forcats::fct_reorder(cell_class, wt, .na_rm = TRUE), y = wt, color = cell_class)) +
  geom_boxplot() +
  scale_y_continuous(name = "Cell Distance (Wild type)",
                     limits = c(0.325, 0.58)) +
  scale_color_manual(values = cols$cell_class) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major.y = element_line(color="grey80", size=0.25),
        panel.grid.minor.y = element_line(color="grey80", size=0.05),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pdf(paste0(dir, "test7.pdf"), width = 12, height = 4)
plot_grid(p4, p5, p6, ncol = 3)
dev.off()

pdf(paste0(dir, "test8.pdf"),  width = 12, height = 8)
plot_grid(p1, p2, p3, p4, p5, p6, ncol = 3)
dev.off()


metdata_for_plot <- data.frame(colData(cds_filt))
# metdata_for_plot$genotype <- ifelse(metdata_for_plot$genotype == "AF16", "Wildtype", metdata_for_plot$genotype)
# metdata_for_plot$genotype <- ifelse(metdata_for_plot$genotype == "N2", "Wildtype", metdata_for_plot$genotype)
metdata_for_plot$genotype <- ifelse(metdata_for_plot$genotype == "wt", "N2", metdata_for_plot$genotype)

BleachTime <- metdata_for_plot %>%
  group_by(genotype, embryo.time.bin) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(genotype) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  tidyr::complete(genotype, embryo.time.bin)

# BleachTime$species <- factor(BleachTime$species, levels = c("C.briggsae", "C.elegans"),
#                              labels = c("C. briggsae", "C. elegans"))

p7 <- ggplot() +
  geom_tile(data = data.frame(BleachTime)[which(BleachTime$genotype != "NA" &
                                                  BleachTime$embryo.time.bin != "NA"),],
            aes(y = genotype,
                x = factor(embryo.time.bin, levels = c("lt_100", "100_130", "130_170",
                                                       "170_210", "210_270", "270_330",
                                                       "330_390", "390_450", "450_510",
                                                       "510_580", "580_650", "650_710", "gt_710")),
                fill = freq), color = "grey80") +
  scale_y_discrete(name = "Genotype") +
  scale_x_discrete(name = "Embryo Time Bin (Minutes)",
                   labels = c("< 100", "100 to 130", "130 to 170", "170 to 210",
                              "210 to 270", "270 to 330", "330 to 390", "390 to 450",
                              "450 to 510", "510 to 580", "580 to 650",
                              "650 to 710", "> 710")) +
  scale_fill_gradientn(name = "Percent",
                       na.value = "#F0F0F0",
                       colors = viridis(4),
                       trans = "log2",
                       breaks = c(10^-5, 10^-4, 10^-3, 10^-2, 10^-1),
                       labels = c("0.001", "0.01", "0.1", "1.0", "10")) +
  theme(rect = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        #axis.line = element_line(color="grey80", size=1),
        axis.line.y = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.ticks = element_line(color = "grey80"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank(),
        strip.placement = "none",
        # plot.margin = margin(0, 5.5, 5.5, 5.5, "pt"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.position = "bottom",
        legend.text = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))

pdf(paste0(dir, "test9.pdf"), height = 2.5, width = 3)
print(p7)
dev.off()

pdf(paste0(dir, "test10.pdf"),  width = 12, height = 8)
plot_grid(p1, plot_grid(p2, p3, ncol = 1), p7, p4, p5, p6, ncol = 3)
dev.off()


######################
# Multiple threshold #
######################

exp_cutoff = 80

gene_data$max_tpm <- apply(gene_data[,c("max_tpm_term", "max_tpm_pro")], 1, max)
gene_data$cel_max_tpm <- apply(gene_data[,c("cel_max_tpm_pro", "cel_max_tpm_term")], 1, max)
gene_data$cbr_max_tpm <- apply(gene_data[,c("cbr_max_tpm_pro", "cbr_max_tpm_term")], 1, max)
gene_data$mean_tau_joint <- (gene_data$cel_tau_median_joint + gene_data$cbr_tau_median_joint) / 2
gene_data$mean_tau_pro <- (gene_data$cel_tau_median_pro + gene_data$cbr_tau_median_pro) / 2
gene_data$mean_tau_term <- (gene_data$cel_tau_median_term + gene_data$cbr_tau_median_joint) / 2

gene_data_sankey <- data.frame(matrix(nrow = dim(gene_data)[1], ncol = 13))
colnames(gene_data_sankey) <- c("gene",
                                "joint_detected", "progenitor_detected", "terminal_detected",
                                "joint_tau", "joint_conservation",
                                "terminal_tau", "terminal_conservation",
                                "progenitor_tau", "progenitor_conservation",
                                "joint_just_conservation", "terminal_just_conservation", "progenitor_just_conservation")
gene_data_sankey$gene <- gene_data$gene

####
# Detected
####

# joint_detected
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm >= exp_cutoff & gene_data$cbr_max_tpm >= exp_cutoff),]$gene,]$joint_detected <- "Both" 
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm < exp_cutoff & gene_data$cbr_max_tpm < exp_cutoff),]$gene,]$joint_detected <- "Neither" 
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm >= exp_cutoff & gene_data$cbr_max_tpm < exp_cutoff),]$gene,]$joint_detected <- "C. elegans" 
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm < exp_cutoff & gene_data$cbr_max_tpm >= exp_cutoff),]$gene,]$joint_detected <- "C. briggsae" 

# progenitor_detected
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm_pro >= exp_cutoff & gene_data$cbr_max_tpm_pro >= exp_cutoff),]$gene,]$progenitor_detected <- "Both" 
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm_pro < exp_cutoff & gene_data$cbr_max_tpm_pro < exp_cutoff),]$gene,]$progenitor_detected <- "Neither" 
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm_pro >= exp_cutoff & gene_data$cbr_max_tpm_pro < exp_cutoff),]$gene,]$progenitor_detected <- "C. elegans" 
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm_pro < exp_cutoff & gene_data$cbr_max_tpm_pro >= exp_cutoff),]$gene,]$progenitor_detected <- "C. briggsae" 

# terminal_detected
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm_term >= exp_cutoff & gene_data$cbr_max_tpm_term >= exp_cutoff),]$gene,]$terminal_detected <- "Both" 
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm_term < exp_cutoff & gene_data$cbr_max_tpm_term < exp_cutoff),]$gene,]$terminal_detected <- "Neither" 
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm_term >= exp_cutoff & gene_data$cbr_max_tpm_term < exp_cutoff),]$gene,]$terminal_detected <- "C. elegans" 
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$cel_max_tpm_term < exp_cutoff & gene_data$cbr_max_tpm_term >= exp_cutoff),]$gene,]$terminal_detected <- "C. briggsae" 



calc_diff_thresholds <- function(gene_data, gene_data_sankey, jsd_perc, tau_perc) {
  exp_cutoff = 80
  jsd_cutoff = c(quantile(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$jsd_median_joint, jsd_perc, na.rm = TRUE),
                 quantile(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$jsd_median_joint, 1 - jsd_perc, na.rm = TRUE))
  names(jsd_cutoff) <- c("lower", "upper")
  tau_cutoff = c(quantile(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$mean_tau_joint, tau_perc, na.rm = TRUE),
                 quantile(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$mean_tau_joint, 1 - tau_perc, na.rm = TRUE))
  names(tau_cutoff) <- c("lower", "upper")
  
  gene_data_sankey$jsd_perc <- jsd_perc
  gene_data_sankey$tau_perc <- tau_perc
  gene_data_sankey$jsd_cutoff_lower <- jsd_cutoff[1]
  gene_data_sankey$jsd_cutoff_upper <- jsd_cutoff[2]
  gene_data_sankey$tau_cutoff_lower <- tau_cutoff[1]
  gene_data_sankey$tau_cutoff_upper <- tau_cutoff[2]
  
  print(paste(jsd_perc, tau_perc))

  ####
  # Tau
  ####
  
  # joint_tau
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                                (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_tau <- "Broad"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_tau <- "Patterned"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_tau <- "Specific"
  
  # progenitor_tau
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_tau <- "Broad"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_tau <- "Patterned"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_tau <- "Specific"
  
  # terminal_tau
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_tau <- "Broad"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_tau <- "Patterned"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_tau <- "Specific"
  
  ####
  # Conservation
  ####
  
  # joint conservation
  gene_data_sankey[gene_data_sankey$joint_detected %in%
                     c("Both", "C. elegans", "C. briggsae"),]$joint_conservation <- paste0(gene_data_sankey[gene_data_sankey$joint_detected %in%
                                                                                                              c("Both", "C. elegans", "C. briggsae"),]$joint_tau, " ", "neutral")
  
  if((length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                                           (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0)) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Broad conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                                            (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Patterned conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Specific conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                                        (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Broad diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Patterned diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Specific diverged"
  }
  
  # just conservation
  gene_data$joint_just_conservation <- "Neutral"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_joint <= jsd_cutoff["lower"]),]$gene,]$joint_just_conservation <- "Conserved"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_joint >= jsd_cutoff["upper"]),]$gene,]$joint_just_conservation <- "Diverged"
  
  
  # terminal conservation
  gene_data_sankey[gene_data_sankey$terminal_detected %in%
                     c("Both", "C. elegans", "C. briggsae"),]$terminal_conservation <- paste0(gene_data_sankey[gene_data_sankey$terminal_detected %in%
                                                                                                                 c("Both", "C. elegans", "C. briggsae"),]$terminal_tau, " ", "neutral")
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Broad conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Patterned conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Specific conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Broad diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Patterned diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Specific diverged"
  }
  
  # just conservation
  gene_data$terminal_just_conservation <- "Neutral"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"]),]$gene,]$terminal_just_conservation <- "Conserved"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"]),]$gene,]$terminal_just_conservation <- "Diverged"
  
  # progenitor conservation
  gene_data_sankey[gene_data_sankey$progenitor_detected %in%
                     c("Both", "C. elegans", "C. briggsae"),]$progenitor_conservation <- paste0(gene_data_sankey[gene_data_sankey$progenitor_detected %in%
                                                                                                                   c("Both", "C. elegans", "C. briggsae"),]$progenitor_tau, " ", "neutral")
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Broad conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Patterned conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Specific conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Broad diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Patterned diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Specific diverged"
  }
  
  # just conservation
  gene_data$progenitor_just_conservation <- "Neutral"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"]),]$gene,]$progenitor_just_conservation <- "Conserved"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"]),]$gene,]$progenitor_just_conservation <- "Diverged"
  
  gene_data_sankey <- left_join(gene_data_sankey, gene_data[,c("gene", "WormCat.1", "WormCat.2", "WormCat.3")])
  rownames(gene_data_sankey) <- gene_data_sankey$gene
  
  gene_data_sankey[which(gene_data_sankey$joint_conservation == "NA neutral"),]$joint_conservation <- NA
  gene_data_sankey[which(gene_data_sankey$terminal_conservation == "NA neutral"),]$terminal_conservation <- NA
  gene_data_sankey[which(gene_data_sankey$progenitor_conservation == "NA neutral"),]$progenitor_conservation <- NA
  
  return(gene_data_sankey)
}

diff_treshold_list <- list()
different_thresholds <- c()
for(tau_perc in seq(0.05, 0.45, by = 0.05)) {
  for(jsd_perc in seq(0.05, 0.45, by = 0.05)) {
    different_thresholds <- c(different_thresholds, paste0(tau_perc, "_", jsd_perc))
  }
}

diff_treshold_list <- lapply(different_thresholds, function(thresh) {
  jsd_perc <- as.numeric(unlist(strsplit(thresh, "_"))[2])
  tau_perc <- as.numeric(unlist(strsplit(thresh, "_"))[1])
  
  calc_diff_thresholds(gene_data, gene_data_sankey, jsd_perc, tau_perc)
})

diff_treshold <- do.call(rbind, diff_treshold_list)

diff_treshold_bin <- diff_treshold %>%
  group_by(jsd_perc, tau_perc) %>%
  summarise(patterned_neutral = sum(joint_conservation == "Patterned neutral", na.rm = TRUE),
            broad_neutral = sum(joint_conservation == "Broad neutral", na.rm = TRUE),
            specific_neutral = sum(joint_conservation == "Specific neutral", na.rm = TRUE),
            patterned_conserved = sum(joint_conservation == "Patterned conserved", na.rm = TRUE),
            broad_conserved = sum(joint_conservation == "Broad conserved", na.rm = TRUE),
            specific_conserved = sum(joint_conservation == "Specific conserved", na.rm = TRUE),
            patterned_diverged = sum(joint_conservation == "Patterned diverged", na.rm = TRUE),
            broad_diverged = sum(joint_conservation == "Broad diverged", na.rm = TRUE),
            specific_diverged = sum(joint_conservation == "Specific diverged", na.rm = TRUE),
            jsd_cutoff_lower = unique(jsd_cutoff_lower),
            jsd_cutoff_upper = unique(jsd_cutoff_upper),
            tau_cutoff_lower = unique(tau_cutoff_lower),
            tau_cutoff_upper = unique(tau_cutoff_upper))

pdf(paste0(dir, "Plots/gene_figure_plots/test11.pdf"),  width = 4, height = 6)
pivot_longer(diff_treshold_bin,
             cols = c("patterned_neutral", "broad_neutral", "specific_neutral",
                      "patterned_conserved", "broad_conserved", "specific_conserved",
                      "patterned_diverged", "broad_diverged", "specific_diverged"),
             names_to = "joint_conservation", values_to = "count") %>%
ggplot() +
  geom_tile(aes(x = tau_perc, y = jsd_perc, fill = count), color = "grey80") +
  scale_x_continuous(name = "Tau Percentile") +
  scale_y_continuous(name = "JSD Percentile") +
  scale_fill_viridis_c(name = "Count", trans = "log10") +
  facet_wrap(joint_conservation ~ ., nrow = 3, ncol = 3) +
  theme(rect = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.line.y = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.ticks = element_line(color = "grey80"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.position = "bottom",
        legend.text = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

exp_cutoff = 80
ecdf(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$jsd_median_joint)(0.45)
ecdf(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$jsd_median_joint)(0.55)

ecdf(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$mean_tau_joint)(0.4)
ecdf(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$mean_tau_joint)(0.7)

########## Try with two thresholds
calc_diff_thresholds_upper_lower <- function(gene_data, gene_data_sankey, jsd_perc, tau_perc) {
  exp_cutoff = 80
  gene_data_subset <- gene_data[(gene_data$cel_max_tpm >= exp_cutoff & gene_data$cbr_max_tpm >= exp_cutoff),]
  
  jsd_cutoff = c(quantile(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$jsd_median_joint, jsd_perc[1], na.rm = TRUE),
                 quantile(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$jsd_median_joint, jsd_perc[2], na.rm = TRUE))
  names(jsd_cutoff) <- c("lower", "upper")
  
  tau_cutoff = c(quantile(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$mean_tau_joint, tau_perc[1], na.rm = TRUE),
                 quantile(gene_data[(gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff),]$mean_tau_joint, tau_perc[2], na.rm = TRUE))
  names(tau_cutoff) <- c("lower", "upper")
  
  gene_data_sankey$jsd_cutoff_lower <- jsd_cutoff[1]
  gene_data_sankey$jsd_cutoff_upper <- jsd_cutoff[2]
  gene_data_sankey$tau_cutoff_lower <- tau_cutoff[1]
  gene_data_sankey$tau_cutoff_upper <- tau_cutoff[2]
  
  print(paste(jsd_perc, tau_perc))
  
  ####
  # Tau
  ####
  
  # joint_tau
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                                (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_tau <- "Broad"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_tau <- "Patterned"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_tau <- "Specific"
  
  # progenitor_tau
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_tau <- "Broad"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_tau <- "Patterned"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_tau <- "Specific"
  
  # terminal_tau
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_tau <- "Broad"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_tau <- "Patterned"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_tau <- "Specific"
  
  ####
  # Conservation
  ####
  
  # joint conservation
  gene_data_sankey[gene_data_sankey$joint_detected %in%
                     c("Both", "C. elegans", "C. briggsae"),]$joint_conservation <- paste0(gene_data_sankey[gene_data_sankey$joint_detected %in%
                                                                                                              c("Both", "C. elegans", "C. briggsae"),]$joint_tau, " ", "neutral")
  
  if((length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                                         (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0)) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Broad conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Patterned conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Specific conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                                        (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Broad diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Patterned diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Specific diverged"
  }
  
  # just conservation
  gene_data$joint_just_conservation <- "Neutral"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_joint <= jsd_cutoff["lower"]),]$gene,]$joint_just_conservation <- "Conserved"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_joint >= jsd_cutoff["upper"]),]$gene,]$joint_just_conservation <- "Diverged"
  
  
  # terminal conservation
  gene_data_sankey[gene_data_sankey$terminal_detected %in%
                     c("Both", "C. elegans", "C. briggsae"),]$terminal_conservation <- paste0(gene_data_sankey[gene_data_sankey$terminal_detected %in%
                                                                                                                 c("Both", "C. elegans", "C. briggsae"),]$terminal_tau, " ", "neutral")
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Broad conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Patterned conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Specific conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Broad diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Patterned diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Specific diverged"
  }
  
  # just conservation
  gene_data$terminal_just_conservation <- "Neutral"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"]),]$gene,]$terminal_just_conservation <- "Conserved"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"]),]$gene,]$terminal_just_conservation <- "Diverged"
  
  # progenitor conservation
  gene_data_sankey[gene_data_sankey$progenitor_detected %in%
                     c("Both", "C. elegans", "C. briggsae"),]$progenitor_conservation <- paste0(gene_data_sankey[gene_data_sankey$progenitor_detected %in%
                                                                                                                   c("Both", "C. elegans", "C. briggsae"),]$progenitor_tau, " ", "neutral")
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Broad conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Patterned conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Specific conserved"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Broad diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Patterned diverged"
  }
  
  if(length(gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                        (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation) > 0) {
    gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                  (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Specific diverged"
  }
  
  # just conservation
  gene_data$progenitor_just_conservation <- "Neutral"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"]),]$gene,]$progenitor_just_conservation <- "Conserved"
  gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"]),]$gene,]$progenitor_just_conservation <- "Diverged"
  
  gene_data_sankey <- left_join(gene_data_sankey, gene_data[,c("gene", "WormCat.1", "WormCat.2", "WormCat.3")])
  rownames(gene_data_sankey) <- gene_data_sankey$gene
  
  gene_data_sankey[which(gene_data_sankey$joint_conservation == "NA neutral"),]$joint_conservation <- NA
  gene_data_sankey[which(gene_data_sankey$terminal_conservation == "NA neutral"),]$terminal_conservation <- NA
  gene_data_sankey[which(gene_data_sankey$progenitor_conservation == "NA neutral"),]$progenitor_conservation <- NA
  
  return(gene_data_sankey)
}


different_thresholds <- c()
for(tau_perc_lower in seq(0.1, 0.4, by = 0.1)) {
  for(tau_perc_upper in seq(0.5, 0.9, by = 0.1)) {
    for(jsd_perc_lower in seq(0.1, 0.4, by = 0.1)) {
      for(jsd_perc_upper in seq(0.6, 0.9, by = 0.1)) {
        different_thresholds <- c(different_thresholds, paste0(tau_perc_lower, "_", tau_perc_upper, "_", jsd_perc_lower, "_", jsd_perc_upper))
      }
    }
  }
}

diff_treshold_list <- lapply(different_thresholds, function(thresh) {
  tau_perc <- c(as.numeric(unlist(strsplit(thresh, "_"))[1]), as.numeric(unlist(strsplit(thresh, "_"))[2]))
  names(tau_perc) <- c("lower", "upper")
  jsd_perc <- c(as.numeric(unlist(strsplit(thresh, "_"))[3]), as.numeric(unlist(strsplit(thresh, "_"))[4]))
  names(jsd_perc) <- c("lower", "upper")
  
  out_df <- calc_diff_thresholds_upper_lower(gene_data, gene_data_sankey, jsd_perc, tau_perc)
  out_df$tau_perc_lower <- tau_perc["lower"]
  out_df$tau_perc_upper <- tau_perc["upper"]
  out_df$jsd_perc_lower <- jsd_perc["lower"]
  out_df$jsd_perc_upper <- jsd_perc["upper"]
  return(out_df)
})

names(diff_treshold_list) <- different_thresholds
diff_treshold <- do.call(rbind, diff_treshold_list)

diff_treshold_bin <- diff_treshold %>%
  group_by(jsd_cutoff_lower, jsd_cutoff_upper, tau_cutoff_lower, tau_cutoff_upper) %>%
  summarise(patterned_neutral = sum(joint_conservation == "Patterned neutral", na.rm = TRUE),
            broad_neutral = sum(joint_conservation == "Broad neutral", na.rm = TRUE),
            specific_neutral = sum(joint_conservation == "Specific neutral", na.rm = TRUE),
            patterned_conserved = sum(joint_conservation == "Patterned conserved", na.rm = TRUE),
            broad_conserved = sum(joint_conservation == "Broad conserved", na.rm = TRUE),
            specific_conserved = sum(joint_conservation == "Specific conserved", na.rm = TRUE),
            patterned_diverged = sum(joint_conservation == "Patterned diverged", na.rm = TRUE),
            broad_diverged = sum(joint_conservation == "Broad diverged", na.rm = TRUE),
            specific_diverged = sum(joint_conservation == "Specific diverged", na.rm = TRUE),
            jsd_cutoff_lower = unique(jsd_cutoff_lower),
            jsd_cutoff_upper = unique(jsd_cutoff_upper),
            tau_cutoff_lower = unique(tau_cutoff_lower),
            tau_cutoff_upper = unique(tau_cutoff_upper),
            tau_perc_lower = unique(tau_perc_lower),
            tau_perc_upper = unique(tau_perc_upper),
            jsd_perc_lower = unique(jsd_perc_lower),
            jsd_perc_uppper = unique(jsd_perc_upper)) %>%
  data.frame()

pdf(paste0(dir, "Plots/gene_figure_plots/test12.pdf"),  width = 16, height = 12)
pivot_longer(diff_treshold_bin,
             cols = c("patterned_neutral", "broad_neutral", "specific_neutral",
                      "patterned_conserved", "broad_conserved", "specific_conserved",
                      "patterned_diverged", "broad_diverged", "specific_diverged"),
             names_to = "joint_conservation", values_to = "count") %>%
  ggplot(aes(x = paste0("Lower: ", round(tau_cutoff_lower, 2), " Upper: ", round(tau_cutoff_upper, 2)),
             y = paste0("Lower: ", round(jsd_cutoff_lower, 2), " Upper: ", round(jsd_cutoff_upper, 2)),
             fill = count, label = count)) +
  geom_tile(color = "grey80") +
  geom_text(size = 2) + 
  scale_x_discrete(name = "Tau Quantile Thresholds (10% to 40% on lower bounds and 50% to 90% on upper)") +
  scale_y_discrete(name = "JSD Quantile Thresholds (10% to 40% on lower bounds and 60% to 90% on upper)") +
  scale_fill_viridis_c(name = "Count", trans = "log10") +
  facet_wrap(joint_conservation ~ ., nrow = 3, ncol = 3) +
  theme(rect = element_rect(fill = "transparent"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.line.y = element_line(color="grey80", size=1),
        axis.line.x = element_line(color="grey80", size=1),
        axis.ticks = element_line(color = "grey80"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.position = "bottom",
        legend.text = element_text(angle = 0, vjust = 1, hjust = 1),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

rownames(diff_treshold_bin) <- paste(diff_treshold_bin$tau_perc_lower, diff_treshold_bin$tau_perc_upper, diff_treshold_bin$jsd_perc_lower, diff_treshold_bin$jsd_perc_uppper, sep = "_")
# 0.1 to 0.4 tau lower
# 0.5 to 0.9 tau upper
# 0.1 to 0.4 jsd lower
# 0.6 to 0.9 jsd upper
diff_treshold_bin_filt <- diff_treshold_bin[c("0.2_0.5_0.4_0.6", # ~main
                                              "0.3_0.6_0.4_0.6", # more broad, less specific
                                              "0.1_0.7_0.4_0.6", # less broad, less specific
                                              "0.2_0.5_0.3_0.7", # less divergent, less conserved
                                              "0.2_0.5_0.4_0.7"),] # less divergent

# Table out


############################
# Wormcat Enrichment Plots #
############################

rownames(gff_list_mrna[["elegans"]]) <- ifelse(! is.na(gff_list_mrna[["elegans"]]$cds_gene_name), gff_list_mrna[["elegans"]]$cds_gene_name, gff_list_mrna[["elegans"]]$gene_name)

# Calculate the enrichment of WormCat terms from Amy Walker
# Code below is derived from their github page:
# https://github.com/dphiggs01/Wormcat
.merger_cats <- function(rgs_annotated_cat, annotated_cat, total_annotations_count, total_rgs_count) {
  merged_cats <- merge(rgs_annotated_cat, annotated_cat, by = "Var1", all.x = TRUE)
  colnames(merged_cats) <- c("Category", "RGS", "AC")
  
  # Step 5: Build contingency table for each category in RGS vs AC
  df <- data.frame(Category = character(),
                   RGS = double(),
                   AC = double(),
                   Fold = double(),
                   PValue = double(),
                   stringsAsFactors = FALSE)
  
  fact_character <- levels(merged_cats$Category)[as.numeric(merged_cats$Category)]
  
  for (i in 1:nrow(merged_cats)) {
    if (is.na(merged_cats$RGS[i]) | is.na(merged_cats$AC[i])) {
      pvalue <- NA
    } else {
      stat <- fisher.test(matrix(c(merged_cats$RGS[i], total_rgs_count,
                                   merged_cats$AC[i],  total_annotations_count),
                                 nrow = 2, ncol = 2),
                          alternative = "greater")
      pvalue <- stat$p.value
    }
    
    df[nrow(df) + 1, ] <- list(Category = fact_character[i],
                               RGS = merged_cats$RGS[i],
                               AC = merged_cats$AC[i],
                               Fold = (merged_cats$RGS[i]/sum(merged_cats$RGS))/(merged_cats$AC[i]/sum(merged_cats$AC)), # added fold enrichment
                               pvalue)
  }
  
  sorted_df <- df[with(df, order(PValue)), ]
  return(sorted_df)
}

.worm_cat_acceptable_pvalues <- function(rgs_fisher_cat) {
  Bonferroni <- NULL
  
  rgs_fisher_cat <- na.omit(rgs_fisher_cat)
  
  rgs_fisher_cat[order(rgs_fisher_cat$PValue), ]
  
  bonferroni <- p.adjust(rgs_fisher_cat$PValue, method = "bonferroni")
  rgs_fisher_cat <- data.frame(rgs_fisher_cat, Bonferroni = bonferroni)
  
  ### Acceptable is to be 0.01 on Bonferroni
  # rgs_fisher_cat <- subset(rgs_fisher_cat, Bonferroni < 0.05)
}

# Create a wrapper function for the enrichment
wormcat_enrichment <- function(rgs, background_gene_set, gff) {
  background_cat_list <- list()
  rgs_cat_list <- list()
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    background_cat_list[[tier]] <- data.frame(table(gff[background_gene_set, tier]))
    rgs_cat_list[[tier]] <- data.frame(table(gff[rgs, tier]))
  }
  
  total_background <- sum(background_cat_list[["WormCat.1"]][,"Freq"])
  total_rgs <- length(rgs)
  
  out_list <- list()
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    out_list[[tier]] <- .merger_cats(rgs_cat_list[[tier]],
                                     background_cat_list[[tier]],
                                     total_background, total_rgs)
    out_list[[tier]] <- .worm_cat_acceptable_pvalues(out_list[[tier]])
  }
  return(out_list)
}



for(filter_set in rownames(diff_treshold_bin_filt)) {
  print(filter_set)
  gene_data_sankey_temp <- diff_treshold_list[[filter_set]]
  background_gene_set <- gene_data_sankey_temp[gene_data_sankey_temp$terminal_detected %in% c("Both", "C. elegans", "C. briggsae"),]$gene
  
  Cat_WormCat_term <- list()
  
  temp_cats <- unique(gene_data_sankey_temp$terminal_conservation)
  temp_cats <- temp_cats[temp_cats %in% c("Specific conserved",
                                          "Patterned conserved",
                                          "Broad conserved",
                                          "Specific diverged",
                                          "Patterned diverged")]
  
  for(cat in temp_cats) {
    print(cat)
    rgs <- gene_data_sankey_temp[which(gene_data_sankey_temp$terminal_conservation == cat),]$gene
    
    Cat_WormCat_term[[cat]] <- wormcat_enrichment(rgs, background_gene_set, gff_list_mrna[["elegans"]])
    
    for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
      Cat_WormCat_term[[cat]][[tier]]$type <- cat
    }
  }
  
  background_gene_set <- gene_data_sankey_temp[gene_data_sankey_temp$progenitor_detected %in% c("Both", "C. elegans", "C. briggsae"),]$gene
  Cat_WormCat_pro <- list()
  
  temp_cats <- unique(gene_data_sankey_temp$progenitor_conservation)
  temp_cats <- temp_cats[temp_cats %in% c("Specific conserved",
                                          "Patterned conserved",
                                          "Broad conserved",
                                          "Specific diverged",
                                          "Patterned diverged")]
  
  for(cat in temp_cats) {
    print(cat)
    rgs <- gene_data_sankey_temp[which(gene_data_sankey_temp$progenitor_conservation == cat),]$gene
    
    Cat_WormCat_pro[[cat]] <- wormcat_enrichment(rgs, background_gene_set, gff_list_mrna[["elegans"]])
    
    for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
      Cat_WormCat_pro[[cat]][[tier]]$type <- cat
    }
  }
  
  Cat_WormCat_term <- transpose(Cat_WormCat_term)
  Cat_WormCat_pro <- transpose(Cat_WormCat_pro)
  
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    Cat_WormCat_term[[tier]] <- do.call("rbind", Cat_WormCat_term[[tier]])
    Cat_WormCat_pro[[tier]] <- do.call("rbind", Cat_WormCat_pro[[tier]])
    
    Cat_WormCat_term[[tier]] <- Cat_WormCat_term[[tier]] %>% complete(Category, type)
    Cat_WormCat_pro[[tier]] <- Cat_WormCat_pro[[tier]] %>% complete(Category, type)
    
    Cat_WormCat_term[[tier]][is.na(Cat_WormCat_term[[tier]]$Fold),]$Bonferroni <- 1
    Cat_WormCat_term[[tier]][is.na(Cat_WormCat_term[[tier]]$Fold),]$Fold <- 0
    
    Cat_WormCat_pro[[tier]][is.na(Cat_WormCat_pro[[tier]]$Fold),]$Bonferroni <- 1
    Cat_WormCat_pro[[tier]][is.na(Cat_WormCat_pro[[tier]]$Fold),]$Fold <- 0
    
    Cat_WormCat_term[[tier]]$signif <- ""
    Cat_WormCat_term[[tier]]$signif <- ifelse(Cat_WormCat_term[[tier]]$Bonferroni < 0.05, "*", Cat_WormCat_term[[tier]]$signif)
    Cat_WormCat_term[[tier]]$signif <- ifelse(Cat_WormCat_term[[tier]]$Bonferroni < 0.005, "**", Cat_WormCat_term[[tier]]$signif)
    Cat_WormCat_term[[tier]]$signif <- ifelse(Cat_WormCat_term[[tier]]$Bonferroni < 0.0005, "***", Cat_WormCat_term[[tier]]$signif)
    
    Cat_WormCat_pro[[tier]]$signif <- ""
    Cat_WormCat_pro[[tier]]$signif <- ifelse(Cat_WormCat_pro[[tier]]$Bonferroni < 0.05, "*", Cat_WormCat_pro[[tier]]$signif)
    Cat_WormCat_pro[[tier]]$signif <- ifelse(Cat_WormCat_pro[[tier]]$Bonferroni < 0.005, "**", Cat_WormCat_pro[[tier]]$signif)
    Cat_WormCat_pro[[tier]]$signif <- ifelse(Cat_WormCat_pro[[tier]]$Bonferroni < 0.0005, "***", Cat_WormCat_pro[[tier]]$signif)
  }
  
  
  ## Simple
  Filtered_categories <- c("Cell cycle", "Development", "DNA: replication", "Cilia", "Neuronal function",
                           "Transcription factor", "Transcription factor: NHR", "Transcription factor: T box",
                           "Transcription factor: bHLH", "Transcription factor: homeodomain", "Transmembrane protein",
                           "Neuronal function: synaptic function: neuropeptide")
  
  # Multilevel
  Filtered_categories <- c("Cell cycle", "Cilia",
                           "Development", "DNA: replication",
                           "Transcription factor: NHR", "Transcription factor: T box",
                           "Signaling: heteromeric G protein: receptor",
                           "Transcription factor: bHLH", "Transcription factor: homeodomain",
                           "Neuronal function: synaptic function: neuropeptide",
                           "Transmembrane protein: seven transmembrane receptor",
                           "Signaling: heteromeric G protein: guanylyl cyclase")
  
  term_wormcat_plot <- do.call("rbind", Cat_WormCat_term) %>% arrange(type) %>% filter(Category %in% Filtered_categories) %>% filter(type != "Broad diverged") %>%
    ggplot(aes(x = factor(Category, levels = c("Cell cycle", 
                                               "Development",
                                               "DNA: replication",
                                               "Cilia",
                                               "Neuronal function: synaptic function: neuropeptide",
                                               "Transmembrane protein: seven transmembrane receptor",
                                               "Signaling: heteromeric G protein: guanylyl cyclase",
                                               "Signaling: heteromeric G protein: receptor",
                                               "Transcription factor: NHR", "Transcription factor: T box",
                                               "Transcription factor: bHLH", "Transcription factor: homeodomain")),
               y = factor(type, levels = c("Patterned diverged", "Specific diverged",
                                           "Broad conserved", "Patterned conserved", "Specific conserved")), fill = Fold)) +
    geom_tile(color="grey80",  linewidth = 0.1) +
    geom_text(aes(label = signif)) +
    geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
    geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
    geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
    geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
    scale_x_discrete(guide = guide_axis(angle = 90),
                     labels = c("Cell cycle", 
                                "Development",
                                "DNA replication",
                                "Cilia",
                                "Neuropeptide",
                                "7TM receptor",
                                "Guanylyl cyclase",
                                "GPCR",
                                "NHR", "T box",
                                "bHLH", "Homeodomain")) +
    scale_fill_gradient(name = "Fold\nEnrich.", low = "white", high="blue",
                        breaks = c(0, 1, 2, 3, 4, 5),
                        labels = c("0", "1", "2", "3", "4", ">5"),
                        limits = c(0, 5),
                        oob = scales::squish) +
    theme(axis.title = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_line(color="grey80", size=0.25),
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.position = "none",
          plot.margin = margin(0, 5.5, 5.5, 5.5, "pt"))
  
  pro_wormcat_plot <- do.call("rbind", Cat_WormCat_pro) %>% arrange(type) %>% filter(Category %in% Filtered_categories) %>% filter(type != "Broad diverged") %>%
    ggplot(aes(x = factor(Category, levels = c("Cell cycle", 
                                               "Development",
                                               "DNA: replication",
                                               "Cilia",
                                               "Neuronal function: synaptic function: neuropeptide",
                                               "Transmembrane protein: seven transmembrane receptor",
                                               "Signaling: heteromeric G protein: guanylyl cyclase",
                                               "Signaling: heteromeric G protein: receptor",
                                               "Transcription factor: NHR", "Transcription factor: T box",
                                               "Transcription factor: bHLH", "Transcription factor: homeodomain")),
               y = factor(type, levels = c("Patterned diverged", "Specific diverged",
                                           "Broad conserved", "Patterned conserved", "Specific conserved")), fill = Fold)) +
    geom_tile(color="grey80",  linewidth = 0.1) +
    geom_text(aes(label = signif)) +
    geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
    geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
    geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
    geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
    scale_x_discrete(guide = guide_axis(angle = 45),
                     labels = c("Cell cycle", 
                                "Development",
                                "DNA replication",
                                "Cilia",
                                "Neuropeptide",
                                "7TM receptor",
                                "Guanylyl cyclase",
                                "NHR",
                                "T box",
                                "bHLH",
                                "Homeodomain")) +
    # scale_fill_distiller(name = "Fold\nEnrich.", palette = "Blues", trans = "reverse") +
    scale_fill_gradient(name = "Fold\nEnrich.", low = "white", high="blue",
                        breaks = c(0, 1, 2, 3, 4, 5),
                        labels = c("0", "1", "2", "3", "4", ">5"),
                        limits = c(0, 5),
                        oob = scales::squish) +  # scale_fill_viridis() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_line(color="grey80", size=0.25),
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.position = "none",
          plot.margin = margin(5.5, 5.5, 0, 5.5, "pt"))
  
  pdf(paste0(dir, "Plots/gene_figure_plots/WormCatEnrichment_", filter_set, ".pdf"), height = 3.75, width = 4)
  print(plot_grid(pro_wormcat_plot, term_wormcat_plot, align = "v", nrow = 2, rel_heights = c(0.6, 1), axis = 'lr'))
  dev.off()
}

############
# TPM Bins #
############

TPMLineageList <- readRDS(paste0(dir, "Objects/TPMLineageList.rds"))

gene_tpm <- rbind(data.frame(tpm = c(unlist(TPMLineageList[["C.elegans"]]["cep-1",]),
                            unlist(TPMLineageList[["C.briggsae"]]["cep-1",])),
                    species = c(rep("C.elegans", length(TPMLineageList[["C.elegans"]]["cep-1",])),
                                rep("C.briggsae", length(TPMLineageList[["C.briggsae"]]["cep-1",]))),
                    gene = "cep-1",
                    cell_bin = c(colnames(TPMLineageList[["C.elegans"]]),
                                 colnames(TPMLineageList[["C.briggsae"]]))),
               data.frame(tpm = c(unlist(TPMLineageList[["C.elegans"]]["gld-1",]),
                                  unlist(TPMLineageList[["C.briggsae"]]["gld-1",])),
                          species = c(rep("C.elegans", length(TPMLineageList[["C.elegans"]]["gld-1",])),
                                      rep("C.briggsae", length(TPMLineageList[["C.briggsae"]]["gld-1",]))),
                          gene = "gld-1",
                          cell_bin = c(colnames(TPMLineageList[["C.elegans"]]),
                                       colnames(TPMLineageList[["C.briggsae"]]))))

gene_tpm <- rbind(gene_tpm,
                  data.frame(tpm = TPMListBootstrapMean_term[["C.elegans"]][c("cep-1","gld-1"),"germline"],
                             species = "C.elegans",
                             gene = c("cep-1", "gld-1"),
                             cell_bin = c("Germline")))
gene_tpm <- rbind(gene_tpm,
                  data.frame(tpm = TPMListBootstrapMean_term[["C.briggsae"]][c("cep-1","gld-1"),"germline"],
                             species = "C.briggsae",
                             gene = c("cep-1", "gld-1"),
                             cell_bin = c("Germline")))

gene_tpm <- gene_tpm[! gene_tpm$cell_bin %in% c("Ciliated neurons_100_130",
                                               "Non-ciliated neurons_lt_100",
                                               "Muscle_lt_100",
                                               "Mesoderm_lt_100",
                                               "Hypodermis and seam_lt_100",
                                               "Muscle_lt_100",
                                               "Pharynx and rectal_lt_100",
                                               "Ciliated neurons_lt_100",
                                               "Ciliated neurons_gt_710",
                                               "Glia and excretory_gt_710",
                                               "Hypodermis and seam_gt_710",
                                               "Pharynx and rectal_580_650",
                                               "Pharynx and rectal_650_710",
                                               "Hypodermis and seam_650_710",
                                               "Pharynx and rectal_gt_710"),]

pdf(paste0(dir, "Plots/cep_1_gld_1_tpm_bins.pdf"), height = 3, width = 16)
ggplot(gene_tpm) +
  geom_tile(aes(x = cell_bin, y = paste(gene, species), fill = tpm)) +
  scale_fill_viridis_c(trans = "log2") +
  theme(#legend.title= element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    rect = element_rect(fill = "transparent"),
    axis.line = element_line(color="grey80", size=1),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  )
dev.off()