library(ggplot2)
library(ggrepel)
library(ggExtra)
library(cowplot)
library(dplyr)
library(viridis)
library(pheatmap)
library(purrr)
library(tidyr)

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"
source(paste0(dir, "../Scripts/WS290/CalcFunctions.R"))

# load barcode bins from EmbryoTimeBins
BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))
ProBarcodeBinsList <- readRDS(paste0(dir, "Objects/ProBarcodeBinsList.rds"))
JointBarcodeBinsList <- readRDS(paste0(dir, "Objects/JointBarcodeBinsList.rds"))

# load cell tables
CellTable <- readRDS(paste0(dir, "Objects/CellTable_Names_20240901.rds"))
cell_type_time_bins <- readRDS(paste0(dir, "Objects/cell_type_time_bins.rds"))
CellTableTimeBins <- readRDS(paste0(dir, "Objects/CellTableTimeBins.rds"))

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240903.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]
rm(cds)

cds_filt_cel <- cds_filt[rowData(cds_filt)$gene.type %in% c("common", "ce.unique"), colData(cds_filt)$species %in% "C.elegans"]
cds_filt_cbr <- cds_filt[rowData(cds_filt)$gene.type %in% c("common", "cb.unique"), colData(cds_filt)$species %in% "C.briggsae"]

gene_data <- readRDS(paste0(dir, "Objects/gene_data_joint_cell_bg.rds"))
cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))

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

library(parallel)
downsample = 20
n_runs = 1000
common_genes <- rowData(cds_filt[rowData(cds_filt)$gene.type %in% c("common"),])$gene_short_name
pseudocount <- 1/length(common_genes)
cell_df_list <- mclapply(names(JointBarcodeBinsList[["C.elegans"]]), function(cell_type_name) {
  temp_cel <- list()
  temp_cbr <- list()
  for(i in 1:n_runs) {
    set.seed(i)
    temp_cel[[i]] <- lapply(JointBarcodeBinsList[["C.elegans"]][[cell_type_name]], function(cell_type_bin) {
      if(length(unlist(cell_type_bin)) > downsample) {
        return(sample(unlist(cell_type_bin), downsample))
      } else {
        return(unlist(cell_type_bin))
      }
    })
    
    temp_cbr[[i]] <- lapply(JointBarcodeBinsList[["C.briggsae"]][[cell_type_name]], function(cell_type_bin) {
      if(length(unlist(cell_type_bin)) > downsample) {
        return(sample(unlist(cell_type_bin), downsample))
      } else {
        return(unlist(cell_type_bin))
      }
    })
  }
  
  TPMList <- list()
  
  temp_cds <- cds_filt_cel[,unlist(JointBarcodeBinsList[["C.elegans"]][[cell_type_name]])]
  temp_matrix <- get.norm.expr.matrix(temp_cds)
  TPMList[["C.elegans"]] <- list()
  for(i in 1:n_runs) {
    TPMList[["C.elegans"]][[i]] <- do.call(cbind, lapply(temp_cel[[i]], function(cell_type_bin) {
      return(BinTPM(temp_matrix[,unlist(cell_type_bin)])[common_genes])
    }))
  }
  
  temp_cds <- cds_filt_cbr[,unlist(JointBarcodeBinsList[["C.briggsae"]][[cell_type_name]])]
  temp_matrix <- get.norm.expr.matrix(temp_cds)
  TPMList[["C.briggsae"]] <- list()
  for(i in 1:n_runs) {
    TPMList[["C.briggsae"]][[i]] <- do.call(cbind, lapply(temp_cbr[[i]], function(cell_type_bin) {
      return(BinTPM(temp_matrix[,unlist(cell_type_bin)])[common_genes])
    }))
  }
  
  cell_df <- data.frame(cell_type = rep(colnames(TPMList[["C.elegans"]][[1]]), each = n_runs),
                        i = rep(1:n_runs, ncol(TPMList[["C.elegans"]][[1]])), jsd = NA, cor = NA)
  rownames(cell_df) <- paste0(cell_df$cell_type, "_", cell_df$i)
  for(cur_cell_bin in colnames(TPMList[["C.elegans"]][[1]])) {
    for(i in 1:n_runs) {
      cell_df[paste0(cur_cell_bin, "_", i),]$cor <- cor(
        log2(TPMList[["C.elegans"]][[i]][, cur_cell_bin] + 1),
        log2(TPMList[["C.briggsae"]][[i]][, cur_cell_bin] + 1),
        method = "pearson",
        use = "na.or.complete")
      
      p = TPMList[["C.elegans"]][[i]][, cur_cell_bin] + pseudocount
      q = TPMList[["C.briggsae"]][[i]][, cur_cell_bin] + pseudocount
      
      cell_df[paste0(cur_cell_bin, "_", i),]$jsd <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
    }
  }
  
  return(cell_df %>% group_by(i) %>% summarise(cell_type = cell_type_name,
                                               jsd = mean(jsd),
                                               cor = mean(cor)))
}, mc.cores = 60)

cell_df <- do.call(rbind, cell_df_list)

# saveRDS(cell_df, paste0(dir, "Objects/cell_df_downsampled_10.rds"))
saveRDS(cell_df, paste0(dir, "Objects/cell_df_downsampled_20.rds"))

cell_df <- readRDS(paste0(dir, "Objects/cell_df_downsampled_20.rds"))


cell_data_downsample <- data.frame(left_join(cell_df %>% group_by(cell_type) %>% summarise(cell_type = unique(cell_type),
                                              median_downsampled_jsd = median(jsd),
                                              mean_downsampled_jsd = mean(jsd),
                                              ci_upper_downsampled_jsd = quantile(jsd, 0.975),
                                              ci_lower_downsampled_jsd = quantile(jsd, 0.025),
                                              median_downsampled_cor = median(cor),
                                              mean_downsampled_cor = mean(cor),
                                              ci_upper_downsampled_cor = quantile(cor, 0.975),
                                              ci_lower_downsampled_cor = quantile(cor, 0.025)),
          cell_data[,c("cell_type", "cell_class", "div_stage", "jsd_median", "cor_median")]))


cor(cell_data_downsample[! cell_data_downsample$cell_type %in% c("ABarppxa", "ABalaapppp/ABalapaapp", "hyp3"),]$jsd_median,
    cell_data_downsample[! cell_data_downsample$cell_type %in% c("ABarppxa", "ABalaapppp/ABalapaapp", "hyp3"),]$mean_downsampled_jsd)

cor(cell_data_downsample[! cell_data_downsample$cell_type %in% c("ABarppxa", "ABalaapppp/ABalapaapp", "hyp3"),]$jsd_median,
    cell_data_downsample[! cell_data_downsample$cell_type %in% c("ABarppxa", "ABalaapppp/ABalapaapp", "hyp3"),]$mean_downsampled_jsd,
    method = "spearman")

reg <- lm(formula = mean_downsampled_jsd ~ jsd_median ,
          data = cell_data_downsample[cell_data_downsample$cell_type != "ABarppxa",])

#get intercept and slope value
coeff <- coefficients(reg)          
intercept <-coeff[1]
slope <- coeff[2]

pdf(file = paste0(dir, "Plots/cell_figure_plots/jsd_downsample_20_to_bootstrap.pdf"), width = 5, height = 5);
cell_data_downsample %>%
  filter(! cell_type %in% c("ABarppxa")) %>%
  ggplot(aes(x = jsd_median, y = mean_downsampled_jsd, color = cell_class, label = cell_type)) +
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype="dashed") +
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  annotate("text", x = c(-Inf), y = c(Inf), hjust = -0.1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  scale_x_continuous(name = "Bootstrapped cell distance (Jensen-Shannon Distance)") +
  scale_y_continuous(name = "Downsampled cell distance (Jensen-Shannon Distance)") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'progenitor' = "grey20")) +
  guides(color = guide_legend(ncol = 2), label = "none") +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()


reg <- lm(formula = mean_downsampled_cor ~ cor_median ,
          data = cell_data_downsample[cell_data_downsample$cell_type != "ABarppxa",])

#get intercept and slope value
coeff <- coefficients(reg)          
intercept <-coeff[1]
slope <- coeff[2]

pdf(file = paste0(dir, "Plots/cell_figure_plots/cor_downsample_20_to_bootstrap.pdf"), width = 5, height = 5);
cell_data_downsample %>%
  filter(! cell_type %in% c("ABarppxa")) %>%
  ggplot(aes(x = cor_median, y = mean_downsampled_cor, color = cell_class, label = cell_type)) +
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype="dashed") +
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  annotate("text", x = c(-Inf), y = c(Inf), hjust = -0.1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  scale_x_continuous(name = "Bootstrapped cell distance (Pearson correlation)") +
  scale_y_continuous(name = "Downsampled cell distance (Pearson correlation)") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'progenitor' = "grey20")) +
  guides(color = guide_legend(ncol = 2), label = "none") +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

#################
# Boxplots plot #
#################
cell_class_plot_order <- c("Germline", "Intestine", "Hypodermis and seam", "Muscle", "Pharynx and rectal", "Mesoderm", "Ciliated neurons", "Glia and excretory", "Non-ciliated neurons", "progenitor")
cell_class_plot_labels <- c("Germline", "Intestine", "Hypodermis and seam", "Muscle", "Pharynx and rectal", "Mesoderm",  "Ciliated neurons", "Glia and excretory", "Non-ciliated neurons", "Progenitors")

jsd_cell_class_pointrange <- cell_data_downsample %>%
  filter(! cell_type %in% "ABarppxa") %>%
  arrange(mean_downsampled_jsd) %>%
  ggplot() + 
  geom_boxplot(aes(x = forcats::fct_reorder(cell_class, mean_downsampled_jsd), y = mean_downsampled_jsd), outlier.shape=NA) +
  geom_pointrange(aes(ymin = ci_lower_downsampled_jsd, ymax = ci_upper_downsampled_jsd,
                      y = mean_downsampled_jsd, x = forcats::fct_reorder(cell_class, mean_downsampled_jsd), group = cell_class, color = cell_class),
                  position = position_dodge2(width = 1), alpha = 0.5, size = 0.5, stroke = 0.25) +
  scale_x_discrete(limits = cell_class_plot_order, labels = cell_class_plot_labels) +
  scale_y_continuous(name = "Cell distance (Jensen-Shannon Distance)") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1")) +
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

pdf(paste0(dir, "Plots/cell_figure_plots/jsd_downsample_20_cell_class_pointrange.pdf"), height = 5, width = 7.5)
print(jsd_cell_class_pointrange)
dev.off()

##########################
# Progenitor point range #
##########################

# <- Light Dark ->
#C6DBEFFF, #9ECAE1FF, #6BAED6FF, #4292C6FF, #2171B5FF, #08519CFF, #08306BFF

div_stage_plot_order <- c("15", "28", "50", "100", "200", "350", "600")
div_stage_plot_labels <- c("15 cell", "28 cell", "50 cell", "100 cell", "200 cell", "350 cell", "600+ cell")

jsd_div_stage_pointrange <- cell_data_downsample %>%
  filter(! cell_type %in% "ABarppxa") %>%
  arrange(mean_downsampled_jsd) %>%
  ggplot() + 
  geom_boxplot(aes(x = div_stage, y = mean_downsampled_jsd), outlier.shape = NA) +
  geom_pointrange(aes(ymin = ci_lower_downsampled_jsd, ymax = ci_upper_downsampled_jsd,
                      y = mean_downsampled_jsd, x = div_stage, group = div_stage, color = div_stage,
                      alpha = div_stage),
                  position = position_dodge2(width = 1), size = 0.5, stroke = 0.25) +
  scale_y_continuous(name = "Cell distance (Jensen-Shannon Distance)") +
  scale_x_discrete(name = "Division stage", limits = div_stage_plot_order, labels = div_stage_plot_labels) +
  scale_color_manual(values = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF")) +
  scale_alpha_manual(values = c("15" = 0.9, "28" = 0.65, "50" = 0.5, "100" = 0.5, "200" = 0.5, "350" = 0.5, "600" = 0.5)) +
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
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pdf(paste0(dir, "Plots/cell_figure_plots/jsd_downsample_20_div_stage_pointrange.pdf"), height = 5, width = 7.5)
print(jsd_div_stage_pointrange)
dev.off()

pdf(paste0(dir, "Plots/cell_figure_plots/jsd_downsample_20_cell_class_div_stage.pdf"), height = 5, width = 12.5)
plot_grid(jsd_cell_class_pointrange, jsd_div_stage_pointrange,
          nrow = 1, align = "hv")
dev.off()
