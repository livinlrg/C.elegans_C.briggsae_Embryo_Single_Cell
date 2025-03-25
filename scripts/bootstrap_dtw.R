

library(dtw)
library(parallel)
library(dplyr)
library(purrr)
library(ggplot2)

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

cel_cell_gene_list_term <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cel_cell_gene_list.rds"))
cbr_cell_gene_list_term <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cbr_cell_gene_list.rds"))

cel_cell_gene_list_pro <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cel_cell_gene_list_pro.rds"))
cbr_cell_gene_list_pro <- readRDS(paste0(dir, "Objects/BootObjects_CellCorrection_out/cbr_cell_gene_list_pro.rds"))

tpm_matrix <- readRDS(paste0(dir, "Objects/tpm_matrix_time_list_filt_cell_bg.rds"))

shared_genes <- rownames(tpm_matrix[["C.elegans"]])[rownames(tpm_matrix[["C.elegans"]]) %in% rownames(tpm_matrix[["C.briggsae"]])]

lineage_branches <- readRDS(paste0(dir, "Objects/lineage_branches.rds"))
CellTableTimeBins <- readRDS(paste0(dir, "Objects/CellTableTimeBins.rds"))

lineage_genes <- list()
for(cur_lineage in names(lineage_branches)) {
  cel_temp <- tpm_matrix[["C.elegans"]][shared_genes, colnames(tpm_matrix[["C.elegans"]]) %in% lineage_branches[[cur_lineage]]$Lineage]
  cbr_temp <- tpm_matrix[["C.briggsae"]][shared_genes, colnames(tpm_matrix[["C.briggsae"]]) %in% lineage_branches[[cur_lineage]]$Lineage]
  
  lineage_genes[[cur_lineage]] <- shared_genes[apply(cel_temp, 1, max) > 80 & apply(cbr_temp, 1, max) > 80]
}

lineage_genes_to_compute <- unique(unlist(lineage_genes))

gene_lineages <- list()
for(this.gene in lineage_genes_to_compute) {
  gene_lineages[[this.gene]] <- names(lineage_genes)[sapply(lineage_genes, function(x) this.gene %in% x)]
}

lineage_branches_subset <- list()
for(cur_lineage in names(lineage_branches)) {
  lineage_branches_subset[[cur_lineage]] <- lineage_branches[[cur_lineage]][lineage_branches[[cur_lineage]]$Lineage %in% colnames(tpm_matrix[[1]]),]
}

compute_dtw_metrics <- function(series1, series2) {
  # Max normalization
  max_series1 <- series1 / max(series1, na.rm = TRUE)
  max_series2 <- series2 / max(series2, na.rm = TRUE)
  
  # Compute raw distance metrics
  raw_dist_maxnorm <- sum(abs(max_series1 - max_series2), na.rm = TRUE)
  
  if (anyNA(max_series1) || anyNA(max_series2)) {
    message("NA in series")
    return(list(raw_dist_maxnorm = NA, dtw_dist_maxnorm = NA))
  }
  
  # Compute DTW distances
  dtw_dist_maxnorm <- dtw(max_series1, max_series2, keep = FALSE, distance.only = TRUE)$distance
  
  return(
    list(raw_dist_maxnorm = raw_dist_maxnorm,
         dtw_dist_maxnorm = dtw_dist_maxnorm
    )
  )
}

CellTableTimeBinsFilt <- CellTableTimeBins[! is.na(CellTableTimeBins$TPMName),]

distance_list <- mclapply(seq_along(lineage_genes_to_compute), function(idx) {
  tryCatch({
    this.gene <- batch_genes[idx]
    
    if (idx %% 100 == 0) {
      message(paste0(this.gene, ": ", idx, " of ", length(batch_genes)))
    }
    
    cel_tpm_matrix <- cbind(cel_cell_gene_list_pro[[this.gene]], cel_cell_gene_list_term[[this.gene]])
    cbr_tpm_matrix <- cbind(cbr_cell_gene_list_pro[[this.gene]], cbr_cell_gene_list_term[[this.gene]])
    
    cel_temp_gene_tpm <- data.frame(t(cel_tpm_matrix))
    cel_temp_gene_tpm$cell <- rownames(cel_temp_gene_tpm)
    
    cbr_temp_gene_tpm <- data.frame(t(cbr_tpm_matrix))
    cbr_temp_gene_tpm$cell <- rownames(cbr_temp_gene_tpm)
    
    cel_tpm_matrix <- left_join(CellTableTimeBinsFilt[,c("Lineage", "TPMName")], cel_temp_gene_tpm, by = c("TPMName" = "cell"))
    rownames(cel_tpm_matrix) <- cel_tpm_matrix$Lineage
    
    cbr_tpm_matrix <- left_join(CellTableTimeBinsFilt[,c("Lineage", "TPMName")], cbr_temp_gene_tpm, by = c("TPMName" = "cell"))
    rownames(cbr_tpm_matrix) <- cbr_tpm_matrix$Lineage
    
    tpm_matrix_list <- list()
    tpm_matrix_list[["C.elegans"]] <- t(cel_tpm_matrix[, -c(1, 2)])
    tpm_matrix_list[["C.briggsae"]] <- t(cbr_tpm_matrix[, -c(1, 2)])
    
    raw_dist <- matrix(nrow = 100, ncol = length(lineage_branches))
    dtw_dist <-  matrix(nrow = 100, ncol = length(lineage_branches))
    
    rownames(raw_dist) <- seq_len(1000)
    rownames(dtw_dist) <- seq_len(1000)
    
    colnames(raw_dist) <- names(lineage_branches)
    colnames(dtw_dist) <- names(lineage_branches)
    
    for(i in seq_len(100)) {
      p = tpm_matrix_list[["C.elegans"]][i,]
      q = tpm_matrix_list[["C.briggsae"]][i,]
      
      for(cur_lineage in gene_lineages[[this.gene]]) {
        cur_cell_types <- lineage_branches_subset[[cur_lineage]]$Lineage
        
        temp <- compute_dtw_metrics(p[cur_cell_types], q[cur_cell_types])
        
        raw_dist[i, cur_lineage] <- temp[[1]]
        dtw_dist[i, cur_lineage] <- temp[[2]]
      }
    }
    return(list(raw_dist, dtw_dist))
  }, error = function(e) {
    message("Encountered error: ", e$message)
    # Return something that signals an error for this gene
    return(NULL)
  })
}, mc.cores = 54)

saveRDS(distance_list, paste0(dir, "Objects/distance_list_dtw_100.rds"))

names(distance_list) <- lineage_genes_to_compute

# raw_dist
# i x lineage
# dtw_dist
# i x lineage

# generate CI's per gene per lineage
distance_ratio_quant_list <- lapply(distance_list, function(cur_gene_matrix) {
  apply(cur_gene_matrix[[2]]/cur_gene_matrix[[1]], 2, function(x) quantile(x, c(0.05, 0.95), na.rm = TRUE))
})

distance_ratio_median_list <- lapply(distance_list, function(cur_gene_matrix) {
  apply(cur_gene_matrix[[2]]/cur_gene_matrix[[1]], 2, function(x) median(x,na.rm = TRUE))
})

distance_diff_quant_list <- lapply(distance_list, function(cur_gene_matrix) {
  apply(abs(cur_gene_matrix[[1]] - cur_gene_matrix[[2]]), 2, function(x) quantile(x, c(0.05, 0.95), na.rm = TRUE))
})

distance_diff_median_list <- lapply(distance_list, function(cur_gene_matrix) {
  apply(abs(cur_gene_matrix[[1]] - cur_gene_matrix[[2]]), 2, function(x) median(x,na.rm = TRUE))
})

filtered_genes <- lapply(distance_ratio_quant_list, function(x) {
  colnames(x)[which(x[2,] < 0.5)]
})

filtered_genes_diff <- lapply(distance_diff_quant_list, function(x) {
  colnames(x)[which(x[1,] > 0.5)]
})

saveRDS(distance_ratio_quant_list, paste0(dir, "Objects/distance_ratio_quant_list.rds"))
saveRDS(distance_ratio_median_list, paste0(dir, "Objects/distance_ratio_median_list.rds"))
saveRDS(distance_diff_quant_list, paste0(dir, "Objects/distance_diff_quant_list.rds"))
saveRDS(distance_diff_median_list, paste0(dir, "Objects/distance_diff_median_list.rds"))
saveRDS(filtered_genes, paste0(dir, "Objects/filtered_genes_dtw.rds"))

distance_ratio_quant_list <- readRDS(paste0(dir, "Objects/distance_ratio_quant_list.rds"))
distance_ratio_median_list <- readRDS(paste0(dir, "Objects/distance_ratio_median_list.rds"))
distance_diff_quant_list <- readRDS(paste0(dir, "Objects/distance_diff_quant_list.rds"))
distance_diff_median_list <- readRDS(paste0(dir, "Objects/distance_diff_median_list.rds"))
filtered_genes <- readRDS(paste0(dir, "Objects/filtered_genes_dtw.rds"))

distance_ratio_range_matrix <- do.call(rbind, lapply(distance_ratio_quant_list, function(x) {
  abs(x[1,] - x[2,])
}))

distance_diff_matrix <- do.call(rbind, distance_diff_median_list)
summary(apply(distance_diff_matrix, 2, function(x) quantile(x, c(0.95), na.rm = TRUE)))

distance_ratio_matrix <- do.call(rbind, distance_ratio_median_list)
summary(apply(distance_ratio_matrix, 2, function(x) quantile(x, c(0.05), na.rm = TRUE)))

summary(apply(distance_diff_matrix, 2, function(x) ecdf(x)(0.5)))
summary(apply(distance_ratio_matrix, 2, function(x) ecdf(x)(0.67)))

##############
# Time plots #
##############

early_cell_tpm_list <- readRDS(paste0(dir, "Objects/early_cell_tpm_list.rds"))
tpm_matrix_lower_ci <- readRDS(paste0(dir, "Objects/tpm_matrix_lower_ci_time_list_filt_cell_bg.rds"))
tpm_matrix_upper_ci <- readRDS(paste0(dir, "Objects/tpm_matrix_upper_ci_time_list_filt_cell_bg.rds"))

tpm_matrix[["C.elegans"]] <- cbind(tpm_matrix[["C.elegans"]], early_cell_tpm_list[["C.elegans"]])
tpm_matrix[["C.briggsae"]] <- cbind(tpm_matrix[["C.briggsae"]], early_cell_tpm_list[["C.briggsae"]])

lineage_branches_subset <- list()
for(cur_lineage in names(lineage_branches)) {
  lineage_branches_subset[[cur_lineage]] <- lineage_branches[[cur_lineage]][lineage_branches[[cur_lineage]]$Lineage %in% colnames(tpm_matrix[[1]]),]
}

rownames(CellTableTimeBins) <- CellTableTimeBins$Lineage
lineage_data <- data.frame(lineage_name = names(lineage_branches_subset),
                           terminal_name = CellTableTimeBins[names(lineage_branches_subset),]$TerminalDatasetName,
                           tpm_name = CellTableTimeBins[names(lineage_branches_subset),]$TPMName)


plot_lineage_tpm <- function(cur_gene, cur_lineage_name, tpm_matrix, tpm_matrix_lower_ci, tpm_matrix_upper_ci, lineage_branches) {
  cells_to_plot <- lineage_branches[[cur_lineage_name]]$Lineage
  cells_to_plot <- cells_to_plot[cells_to_plot != "28_cell_or_earlier"]
  
  df <- data.frame(tpm = c(tpm_matrix[["C.elegans"]][cur_gene, cells_to_plot],
                          tpm_matrix[["C.briggsae"]][cur_gene, cells_to_plot]),
                  tpm_lower = c(tpm_matrix_lower_ci[["C.elegans"]][cur_gene, cells_to_plot],
                                tpm_matrix_lower_ci[["C.briggsae"]][cur_gene, cells_to_plot]),
                  tpm_upper = c(tpm_matrix_upper_ci[["C.elegans"]][cur_gene, cells_to_plot],
                                tpm_matrix_upper_ci[["C.briggsae"]][cur_gene, cells_to_plot]),
                  species = c(rep("C. elegans", length(cells_to_plot)),
                              rep("C. briggsae", length(cells_to_plot))),
                  time_point = factor(c(cells_to_plot, cells_to_plot), levels = rev(cells_to_plot)))

  p1 <- df %>%
    ggplot(aes(x = time_point, y = tpm, group = species, color = species)) +
    geom_pointrange(aes(ymin = tpm_lower, ymax = tpm_upper)) +
    geom_line() +
    scale_color_manual(limits = c("C. elegans", "C. briggsae"),
                       values = c("#009E73", "#56B4E9")) +
    scale_x_discrete() +
    scale_y_continuous(name = "TPM") +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color = NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid.major.y = element_line(color="grey80", size=0.25),
          panel.grid.minor.y = element_line(color="grey80", size=0.05),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.position = "top")
  
  df <- df %>%
    group_by(species) %>%
    mutate(tpm_lower = tpm_lower/max(tpm),
           tpm_upper = tpm_upper/max(tpm)) %>%
    mutate(tpm = tpm/max(tpm))
  
  p2 <- df %>%
    ggplot(aes(x = time_point, y = tpm, group = species, color = species)) +
    geom_pointrange(aes(ymin = tpm_lower, ymax = tpm_upper)) +
    geom_line() +
    scale_color_manual(limits = c("C. elegans", "C. briggsae"),
                       values = c("#009E73", "#56B4E9")) +
    scale_x_discrete() +
    scale_y_continuous(name = "Normalized TPM") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          axis.title.x = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color = NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid.major.y = element_line(color="grey80", size=0.25),
          panel.grid.minor.y = element_line(color="grey80", size=0.05),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.position = "none")
  
  p_out <- plot_grid(p1, p2, nrow = 2, rel_heights = c(0.7, 1))
  
  return(p_out)
}

for(cur_gene in names(filtered_genes)) {
  if(length(filtered_genes[[cur_gene]]) > 0) {
    for(cur_lineage in filtered_genes[[cur_gene]]) {
      
      main_title <- ggdraw() +
        draw_label(
          paste0(cur_gene, ": ", unique(lineage_data[which(lineage_data$lineage_name == cur_lineage), "tpm_name"])),
          fontface = 'bold', x = 0, hjust = 0) +
        theme(plot.margin = margin(5, 5, 5, 7))
      ratio <- ggdraw() +
        draw_label(
          paste0("Ratio: ", round(distance_ratio_median_list[[cur_gene]][cur_lineage], 2),
                 " Lower 95%: ", round(distance_ratio_quant_list[[cur_gene]][1,cur_lineage], 2),
                 " Upper 95: ", round(distance_ratio_quant_list[[cur_gene]][2,cur_lineage], 2)),
          x = 0, hjust = 0) +
        theme(plot.margin = margin(5, 5, 5, 7))
      diff <- ggdraw() +
        draw_label(
          paste0("Diff.: ", round(distance_diff_median_list[[cur_gene]][cur_lineage], 2),
                 " Lower 95%: ", round(distance_diff_quant_list[[cur_gene]][1,cur_lineage], 2),
                 " Upper 95: ", round(distance_diff_quant_list[[cur_gene]][2,cur_lineage], 2)),
          x = 0, hjust = 0) +
        theme(plot.margin = margin(5, 5, 5, 7))
      
      main_plot <- plot_lineage_tpm(cur_gene, cur_lineage,
                                    tpm_matrix, tpm_matrix_lower_ci, tpm_matrix_upper_ci,
                                    lineage_branches_subset)
      
      pdf(paste0(dir, "Plots/time_plots/", cur_gene, "_", gsub(":", ".", cur_lineage), ".pdf"), width = 5, height = 5)
      print(plot_grid(main_title,
                      ratio,
                      diff,
                      main_plot,
                      nrow = 4,
                      rel_heights = c(0.05, 0.05, 0.05, 1)))
      dev.off()
    }
  }
}

distance_ratio_quant_list <- readRDS(paste0(dir, "Objects/distance_ratio_quant_list.rds"))
distance_ratio_median_list <- readRDS(paste0(dir, "Objects/distance_ratio_median_list.rds"))
distance_diff_quant_list <- readRDS(paste0(dir, "Objects/distance_diff_quant_list.rds"))
distance_diff_median_list <- readRDS(paste0(dir, "Objects/distance_diff_median_list.rds"))
filtered_genes <- readRDS(paste0(dir, "Objects/filtered_genes_dtw.rds"))

distance_diff_matrix <- do.call(rbind, distance_diff_median_list)
distance_ratio_matrix <- do.call(rbind, distance_ratio_median_list)

temp_distance_df <- lapply(colnames(distance_diff_matrix), function(cur_lineage) {
  temp_df = data.frame(gene = rownames(distance_diff_matrix),
                       lineage = rep(cur_lineage, nrow(distance_diff_matrix)),
                       ratio_median = distance_ratio_matrix[,cur_lineage],
                       ratio_upper_ci = unlist(lapply(distance_ratio_quant_list, function(x) x[2, cur_lineage])),
                       ratio_lower_ci = unlist(lapply(distance_ratio_quant_list, function(x) x[1, cur_lineage])),
                       diff_median = distance_diff_matrix[,cur_lineage],
                       diff_upper_ci = unlist(lapply(distance_diff_quant_list, function(x) x[2, cur_lineage])),
                       diff_lower_ci = unlist(lapply(distance_diff_quant_list, function(x) x[1, cur_lineage])))
  return(temp_df)
})
distance_df <- do.call(rbind, temp_distance_df)
distance_df <- distance_df %>% filter(! is.na(ratio_median))

filtered_genes_high <- lapply(distance_ratio_quant_list, function(x) {
  colnames(x)[which(x[2,] < 0.5)]
})

filtered_genes_medium <- lapply(distance_ratio_quant_list, function(x) {
  colnames(x)[which(x[2,] < 0.67)]
})

distance_df$pass <- "not different"
for(cur_gene in unique(distance_df$gene)) {
  if(length(filtered_genes_medium[[cur_gene]]) > 0) {
    distance_df[which(distance_df$lineage %in% filtered_genes_medium[[cur_gene]] &
                        distance_df$gene == cur_gene),]$pass = "medium"
  }
  
  if(length(filtered_genes_high[[cur_gene]]) > 0) {
    distance_df[which(distance_df$lineage %in% filtered_genes_high[[cur_gene]] &
                        distance_df$gene == cur_gene),]$pass = "high"
  }
}

pdf(paste0(dir, "Plots/test.pdf"), width = 15, height = 15)
ggMarginal(
distance_df %>%
  filter(terminal_name == "ADF") %>%
  ggplot(aes(x = ratio_median, y = diff_median,
             ymin = diff_lower_ci, ymax = diff_upper_ci,
             xmax = ratio_upper_ci, xmin = ratio_lower_ci,
             color = pass,
             alpha = ifelse(pass %in% c("high", "medium"), 1, 0.25))) +
  # geom_point(size = 0.1, alpha = 0.5)) +
  geom_pointrange(size = 0.1, linewidth = 0.3) +
  geom_errorbarh(size = 0.3) +
  scale_x_continuous(name = "Time warped distance / raw distance") +
  scale_y_continuous(name = "| Raw distance - Time warped distance |") +
  scale_color_manual(name = "Confidence",
                     breaks = c("high", "medium", "not different"),
                     values = c("#D55E00", "#CC79A7", "black"),
                     labels = c("High (>0.5)", "Medium (>0.67)", "Not different")) +
  guides(color = "legend", alpha = "none") +
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color = NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major.y = element_line(color="grey80", size=0.25),
        panel.grid.minor.y = element_line(color="grey80", size=0.05),
        panel.grid.major.x = element_line(color="grey80", size=0.25),
        panel.grid.minor.x = element_line(color="grey80", size=0.05),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.position = "right")
)
dev.off()

distance_df <- left_join(distance_df, lineage_data, by = c("lineage" = "lineage_name"))
saveRDS(distance_df, paste0(dir, "Objects/dtw_distance_df.rds"))

pdf(paste0(dir, "Plots/test.pdf"), width = 5, height = 2.5)
ggplot() +
  geom_density(data = distance_df,
               aes(x = ratio_median), color = "black", alpha = 0.4) +
  geom_density(data = distance_df[distance_df$pass %in% c("medium"),],
               aes(x = ratio_median), color = "orange", alpha = 0.6) +
  geom_density(data = distance_df[distance_df$pass %in% c("high"),],
               aes(x = ratio_median), color = "red", alpha = 0.6) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 0.67, linetype = "dashed", color = "orange") +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color = NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.position = "none")
dev.off()

ecdf(distance_df$ratio_median)(0.5) * 100
ecdf(distance_df$ratio_median)(0.67) * 100

quantile(distance_df$ratio_median, c(0.95, 0.05), na.rm = TRUE)