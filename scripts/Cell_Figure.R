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

gene_data <- readRDS(paste0(dir, "Objects/gene_data_joint_cell_bg.rds"))
cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))
cell_data_bins <- readRDS(paste0(dir, "Objects/cell_data_joint_cell_bg.rds"))

gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))

TPMListBootstrapMean_term <- readRDS(paste0(dir, "Objects/TPMListBootstrapMean_CellCorrection.rds"))
BinarizeBootstrapList_term <- readRDS(paste0(dir, "Objects/BinarizeBootstrapList_CellCorrection.rds"))

TPMListBootstrapMean_term_subset <- list()
TPMListBootstrapMean_term_subset[["C.elegans"]] <- TPMListBootstrapMean_term[["C.elegans"]][gene_data$gene,]
TPMListBootstrapMean_term_subset[["C.briggsae"]] <- TPMListBootstrapMean_term[["C.briggsae"]][gene_data$gene,]

TPMListBootstrap_term <- readRDS(paste0(dir, "Objects/TPMTimeBinsListBootstrap_CellCorrection.rds"))

TPMListBootstrap_term_subset <- list()
TPMListBootstrap_term_subset[["C.elegans"]] <- TPMListBootstrap_term[["C.elegans"]][gene_data$gene,]
TPMListBootstrap_term_subset[["C.briggsae"]] <- TPMListBootstrap_term[["C.briggsae"]][gene_data$gene,]

TPMListBootstrap_pro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))
BinarizeBootstrapListPro <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListPro_CellCorrection.rds"))

TPMListBootstrap_pro_subset <- list()
TPMListBootstrap_pro_subset[["C.elegans"]] <- TPMListBootstrap_pro[["C.elegans"]][gene_data$gene,]
TPMListBootstrap_pro_subset[["C.briggsae"]] <- TPMListBootstrap_pro[["C.briggsae"]][gene_data$gene,]

term_markers <- readRDS(paste0(dir, "Objects/TermMarker_list_filt_metadata_cell_bg.rds"))
pro_markers <- readRDS(paste0(dir, "Objects/ProMarker_list_filt_metadata_cell_bg.rds"))
joint_markers <- readRDS(paste0(dir, "Objects/joint_markers_cell_bg.rds"))

deg <- readRDS(paste0(dir, "Objects/DEG_df_filt_metadata_cell_bg.rds"))

#############
# TPM Plots #
#############

comp_TPM <- log2((TPMListBootstrapMean_term_subset[["C.elegans"]] + 1) / (TPMListBootstrapMean_term_subset[["C.briggsae"]] + 1))
mean_TPM <- (TPMListBootstrapMean_term_subset[["C.elegans"]] + TPMListBootstrapMean_term_subset[["C.briggsae"]])/2

TPMPlot_List <- list()
for(cell_type in c("BWM_middle", "ASG")) {
  print(cell_type)
  
  cell_plot_df <- data.frame(gene = rownames(TPMListBootstrapMean_term_subset[["C.elegans"]]))
  rownames(cell_plot_df) <- cell_plot_df$gene
  
  cell_plot_df$cel_marker <- "Not a marker"
  cell_plot_df[cell_plot_df$gene %in% term_markers[["C.elegans"]][term_markers[["C.elegans"]]$cell_type == cell_type,]$gene,]$cel_marker <- "C. elegans marker"
  cell_plot_df$cbr_marker <- "Not a marker"
  cell_plot_df[cell_plot_df$gene %in% term_markers[["C.briggsae"]][term_markers[["C.briggsae"]]$cell_type == cell_type,]$gene,]$cbr_marker <- "C. briggsae marker"
  
  cell_plot_df$comb_marker <- cell_plot_df$cel_marker
  cell_plot_df$comb_marker <- ifelse(cell_plot_df$cel_marker == "C. elegans marker" & cell_plot_df$cbr_marker == "C. briggsae marker",
                                     "Joint marker", ifelse(cell_plot_df$cbr_marker == "C. briggsae marker", "C. briggsae marker", cell_plot_df$comb_marker))
  
  TPMPlot_List[[cell_type]] <- ggplot() +
    geom_point(aes_(x = mean_TPM[cell_plot_df$comb_marker == "Not a marker", cell_type] + 1,
                    y = comp_TPM[cell_plot_df$comb_marker == "Not a marker", cell_type],
                    color = cell_plot_df[cell_plot_df$comb_marker == "Not a marker",]$comb_marker),
               alpha = 0.5, size = 1, stroke = 0.3) +
    geom_point(aes_(x = mean_TPM[cell_plot_df$comb_marker != "Not a marker", cell_type] + 1,
                    y = comp_TPM[cell_plot_df$comb_marker != "Not a marker", cell_type],
                    color = cell_plot_df[cell_plot_df$comb_marker != "Not a marker",]$comb_marker),
               alpha = 0.5, size = 1, stroke = 0.3) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    scale_color_manual(name = "Markers", breaks = c("C. elegans marker", "C. briggsae marker", "Joint marker", "Not a marker"),
                       labels = c("C. elegans marker", "C. briggsae marker", "Joint marker", "Neither"),
                       values = c("#009E73", "#56B4E9", "#CC79A7", "black")) +
    scale_x_continuous(name = paste0("Mean ", cell_type," TPM"), limits = c(1, 46000), trans = "log2") +
    scale_y_continuous(name = paste0("Log2 ", cell_type, " TPM"), limits = c(-11.5, 11.5)) +
    theme(legend.title= element_blank(),
          legend.position = "none",
          rect = element_rect(fill = "transparent"),
          # axis.title.x = element_text(color = "#009E73", size = 14),
          # axis.title.y = element_text(color = "#56B4E9", size = 14),
          axis.text = element_text(size = 12),
          axis.line = element_line(color="grey80", size=1),
          panel.grid.major = element_line(color="grey80", size=0.25),
          panel.grid.minor = element_line(color="grey80", size=0.05),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          plot.margin = unit(c(-10, -10, 10, 10), "pt"))
  
  comp_TPM_hist <- ggplot() +
    geom_histogram(aes_(y = comp_TPM[mean_TPM[, cell_type] > 100, cell_type]), bins = 80) + 
    scale_y_continuous(limits = c(-11.5, 11.5)) + 
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_rect(fill='transparent'),
          plot.background = element_rect(fill='transparent', color=NA),
          panel.spacing = unit(0, "cm"),
          plot.margin = unit(c(10, 10, 10, 10), "pt"))
  
  TPMPlot_List[[cell_type]] <- plot_grid(TPMPlot_List[[cell_type]],
                                         comp_TPM_hist,
                                         ncol = 2, nrow = 1, align = "h", axis = "rlbt",
                                         rel_heights = c(0.25, 1), rel_widths = c(1, 0.25))
}

pdf(paste0(dir, "Plots/cell_figure_plots/TPM_ASG_BWM_middle.pdf"), width = 4, height = 8)
plot_grid(TPMPlot_List[["BWM_middle"]],
          TPMPlot_List[["ASG"]], align = "hv", nrow = 2)
dev.off()

##############
# jsd matrix #
##############

## Matrix version
pseudocount = 1/dim(TPMListBootstrapMean_term_subset[["C.elegans"]])[1]
jsd_matrix <- matrix(nrow = dim(TPMListBootstrapMean_term_subset[["C.elegans"]])[2],
                     ncol = dim(TPMListBootstrapMean_term_subset[["C.briggsae"]])[2])
rownames(jsd_matrix) <- colnames(TPMListBootstrapMean_term_subset[["C.elegans"]])
colnames(jsd_matrix) <- colnames(TPMListBootstrapMean_term_subset[["C.briggsae"]])
for(ele_cell in rownames(jsd_matrix)) {
  for(bri_cell in colnames(jsd_matrix)) {
    p = TPMListBootstrapMean_term_subset[["C.elegans"]][,ele_cell] + pseudocount
    q = TPMListBootstrapMean_term_subset[["C.briggsae"]][,bri_cell] + pseudocount
    jsd_matrix[ele_cell, bri_cell] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
}

## Fix the cell classes for plotting
ListOfCellTypes <- list()
for(cur_cell_class in unique(cell_data$cell_class)[unique(cell_data$cell_class) != "progenitor"]) {
  print(cur_cell_class)
  ListOfCellTypes[[cur_cell_class]] <- cell_data[cell_data$cell_class == cur_cell_class,]$cell_type
  
  if(cur_cell_class != "Germline") {
    ListOfCellTypes[[cur_cell_class]] <- ListOfCellTypes[[cur_cell_class]][hclust(as.dist(jsd_matrix[ListOfCellTypes[[cur_cell_class]], ListOfCellTypes[[cur_cell_class]]]))$order]
  }
}

saveRDS(ListOfCellTypes, paste0(dir, "Objects/ListOfCellTypes_Mean.rds"))

## Want to order the terminal cell types
OrderForPlot <- c(ListOfCellTypes[["Germline"]],
                  ListOfCellTypes[["Intestine"]],
                  c("BWM_headrow1_in", "BWM_headrow1_out", "BWM_headrow2_in", "BWM_headrow2_out", "BWM_anterior", "BWM_middle", "BWM_posterior", "BWM_far_posterior", "mu_int_mu_anal", "mu_sph"), # Muscle
                  c(ListOfCellTypes[["Mesoderm"]]), # Mesoderm
                  c(ListOfCellTypes[["Glia and excretory"]]), # Glia and Excretory
                  c("Rectal_gland", "Anterior_arcade_cell", "Posterior_arcade_cell", "MC", "B_F_K_Kp_U_Y", "hyp1_hyp2", "hyp1V", "g1A", "g1P", "g2", "mc1", "mc2", "mc3", "pm1_pm2", "pm3_pm4_pm5", "pm6", "pm7", "pm8", "Pharyngeal_intestinal_valve"), # Pharynx and rectal
                  c("P_cells","hyp4_hyp5_hyp6", "hyp7_AB_lineage", "hyp7_C_lineage", "Seam_cells", "T", "Tail_hypodermis"), # Hypodermis and seam
                  c(ListOfCellTypes[["Ciliated neurons"]]), # Ciliated neurons
                  c(ListOfCellTypes[["Non-ciliated neurons"]])) # Non-ciliated neurons

saveRDS(OrderForPlot, paste0(dir, "Objects/OrderForPlot.rds"))
OrderForPlot <- readRDS(paste0(dir, "Objects/OrderForPlot.rds"))

# Plot order for matrices
plot_df <- data.frame(OrderForPlot)
rownames(plot_df) <- plot_df$OrderForPlot

plot_df$cell_class <- NA
for(tissue in names(ListOfCellTypes)) {
  plot_df[ListOfCellTypes[[tissue]],]$cell_class <- tissue
}

plot_df$cell_type <- rownames(plot_df)
plot_df$cell_class <- factor(plot_df$cell_class)
plot_df <- plot_df[plot_df$cell_type %in% OrderForPlot,]

cols = list(cell_class = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                          'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                          'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1"))

pdf(paste0(dir, "Plots/cell_figure_plots/jsd_matrix_ordered_rownames.pdf"), width = 20, height = 20)
pheatmap(jsd_matrix[OrderForPlot, OrderForPlot],
         annotation_row = plot_df[,c(2), drop = FALSE],
         annotation_colors = cols,
         color = rev(inferno(1000)),
         breaks = seq(0.325, 0.85, 0.525/1000),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         dendrogram = "none",
         scale = "none",
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_names_row = FALSE,
         labels_row_side = "right",
         labels_col_side = "right",
         border_color = "black", 
         border_width = 1,
         legend = TRUE,
         legend_title = "Jensen-Shannon Distance")
dev.off()

pdf(paste0(dir, "Plots/cell_figure_plots/jsd_matrix_ordered.pdf"), width = 7, height = 5)
pheatmap(jsd_matrix[OrderForPlot, OrderForPlot],
         annotation_row = plot_df[,c(2), drop = FALSE],
         annotation_colors = cols,
         color = rev(inferno(1000)),
         breaks = seq(0.325, 0.85, 0.525/1000),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         dendrogram = "none",
         scale = "none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_names_row = FALSE,
         labels_row_side = "right",
         labels_col_side = "right",
         border_color = "black", 
         border_width = 1,
         legend = FALSE,
         legend_title = "Jensen-Shannon Distance")
dev.off()


########################
# jsd matrix time bins #
########################

## Matrix version
pseudocount = 1/dim(TPMListBootstrap_term_subset[["C.elegans"]])[1]
jsd_matrix_bins <- matrix(nrow = dim(TPMListBootstrap_term_subset[["C.elegans"]])[2],
                     ncol = dim(TPMListBootstrap_term_subset[["C.briggsae"]])[2])
rownames(jsd_matrix_bins) <- colnames(TPMListBootstrap_term_subset[["C.elegans"]])
colnames(jsd_matrix_bins) <- colnames(TPMListBootstrap_term_subset[["C.briggsae"]])
for(ele_cell in rownames(jsd_matrix_bins)) {
  for(bri_cell in colnames(jsd_matrix_bins)) {
    p = TPMListBootstrap_term_subset[["C.elegans"]][,ele_cell] + pseudocount
    q = TPMListBootstrap_term_subset[["C.briggsae"]][,bri_cell] + pseudocount
    jsd_matrix_bins[ele_cell, bri_cell] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
}

# Plot order for matrices
plot_df <- data.frame(unlist(cell_data[OrderForPlot,]$cell_type_bin))
rownames(plot_df) <- unlist(cell_data[OrderForPlot,]$cell_type_bin)

plot_df$cell_class <- cell_data_bins[rownames(plot_df),]$cell_class

plot_df$cell_type <- rownames(plot_df)
plot_df$cell_class <- factor(plot_df$cell_class)
# plot_df <- plot_df[plot_df$cell_type %in% unlist(cell_data[OrderForPlot,]$cell_type_bins),]

cols = list(cell_class = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                           'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                           'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1"))

pdf(paste0(dir, "Plots/cell_figure_plots/jsd_matrix_bins_ordered_rownames.pdf"), width = 40, height = 40)
pheatmap(jsd_matrix_bins[unlist(cell_data[OrderForPlot,]$cell_type_bin), unlist(cell_data[OrderForPlot,]$cell_type_bin)],
         annotation_row = plot_df[,c(2), drop = FALSE],
         annotation_colors = cols,
         color = rev(inferno(1000)),
         breaks = seq(0.325, 0.85, 0.525/1000),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         dendrogram = "none",
         scale = "none",
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_names_row = FALSE,
         labels_row_side = "right",
         labels_col_side = "right",
         border_color = "black", 
         border_width = 1,
         legend = TRUE,
         legend_title = "Jensen-Shannon Distance")
dev.off()

######################################
# C. elegans terminal matrix version #
######################################
pseudocount = 1/dim(TPMListBootstrapMean_term[["C.elegans"]])[1]
jsd_matrix_cel <- matrix(nrow = dim(TPMListBootstrapMean_term[["C.elegans"]])[2],
                     ncol = dim(TPMListBootstrapMean_term[["C.elegans"]])[2])
rownames(jsd_matrix_cel) <- colnames(TPMListBootstrapMean_term[["C.elegans"]])
colnames(jsd_matrix_cel) <- colnames(TPMListBootstrapMean_term[["C.elegans"]])
for(cel_cell in rownames(jsd_matrix_cel)) {
  for(cbr_cell in colnames(jsd_matrix_cel)) {
    p = TPMListBootstrapMean_term[["C.elegans"]][,cel_cell] + pseudocount
    q = TPMListBootstrapMean_term[["C.elegans"]][,cbr_cell] + pseudocount
    jsd_matrix_cel[cel_cell, cbr_cell] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
}

pdf(paste0(dir, "Plots/supplemental_plots/jsd_matrix_ordered_elegans.pdf"), width = 7, height = 5)
pheatmap(jsd_matrix_cel[OrderForPlot, OrderForPlot],
         annotation_row = plot_df[,c(2), drop = FALSE],
         annotation_colors = cols,
         color = rev(inferno(1000)),
         breaks = seq(0.15, 0.85, 0.7/1000),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         dendrogram = "none",
         scale = "none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_names_row = FALSE,
         labels_row_side = "right",
         labels_col_side = "right",
         border_color = "black", 
         border_width = 1,
         legend = FALSE,
         legend_title = "Jensen-Shannon Distance")
dev.off()

#######################################
# C. briggsae terminal matrix version #
#######################################
pseudocount = 1/dim(TPMListBootstrapMean_term[["C.briggsae"]])[1]
jsd_matrix_cbr <- matrix(nrow = dim(TPMListBootstrapMean_term[["C.briggsae"]])[2],
                         ncol = dim(TPMListBootstrapMean_term[["C.briggsae"]])[2])
rownames(jsd_matrix_cbr) <- colnames(TPMListBootstrapMean_term[["C.briggsae"]])
colnames(jsd_matrix_cbr) <- colnames(TPMListBootstrapMean_term[["C.briggsae"]])
for(cel_cell in rownames(jsd_matrix_cbr)) {
  for(cbr_cell in colnames(jsd_matrix_cbr)) {
    p = TPMListBootstrapMean_term[["C.briggsae"]][,cel_cell] + pseudocount
    q = TPMListBootstrapMean_term[["C.briggsae"]][,cbr_cell] + pseudocount
    jsd_matrix_cbr[cel_cell, cbr_cell] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
}

pdf(paste0(dir, "Plots/supplemental_plots/jsd_matrix_ordered_briggsae.pdf"), width = 7, height = 5)
pheatmap(jsd_matrix_cbr[OrderForPlot, OrderForPlot],
         annotation_row = plot_df[,c(2), drop = FALSE],
         annotation_colors = cols,
         color = rev(inferno(1000)),
         breaks = seq(0.15, 0.85, 0.7/1000),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         dendrogram = "none",
         scale = "none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_names_row = FALSE,
         labels_row_side = "right",
         labels_col_side = "right",
         border_color = "black",
         border_width = 1,
         legend = FALSE,
         legend_title = "Jensen-Shannon Distance")
dev.off()

#############################
# Pro + Term matrix version #
#############################
TPMJoint <- list()
TPMJoint[["C.elegans"]] <- cbind(TPMListBootstrap_pro_subset[["C.elegans"]][,sort(colnames(TPMListBootstrap_pro_subset[["C.elegans"]]))], TPMListBootstrapMean_term_subset[["C.elegans"]][,OrderForPlot])
TPMJoint[["C.briggsae"]] <- cbind(TPMListBootstrap_pro_subset[["C.briggsae"]][,sort(colnames(TPMListBootstrap_pro_subset[["C.elegans"]]))], TPMListBootstrapMean_term_subset[["C.briggsae"]][,OrderForPlot])

pseudocount = 1/dim(TPMJoint[["C.elegans"]])[2]
jsd_matrix_joint <- matrix(nrow = dim(TPMJoint[["C.elegans"]])[2],
                           ncol = dim(TPMJoint[["C.briggsae"]])[2])
rownames(jsd_matrix_joint) <- colnames(TPMJoint[["C.elegans"]])
colnames(jsd_matrix_joint) <- colnames(TPMJoint[["C.briggsae"]])
for(cel_cell in rownames(jsd_matrix_joint)) {
  for(cbr_cell in colnames(jsd_matrix_joint)) {
    p = TPMJoint[["C.elegans"]][,cel_cell] + pseudocount
    q = TPMJoint[["C.briggsae"]][,cbr_cell] + pseudocount
    jsd_matrix_joint[cel_cell, cbr_cell] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
}

# check the matching rate
summary(apply(jsd_matrix_joint, 1, min) == diag(jsd_matrix_joint))

# index of the match
min.match.1 = apply(jsd_matrix_joint, 1, which.min)
min.match.2 = apply(jsd_matrix_joint, 2, which.min)

min_match_df <- data.frame(original_cell_type = colnames(jsd_matrix_joint),
           min_match_1 = colnames(jsd_matrix_joint)[min.match.1],
           min_match_1_jsd = apply(jsd_matrix_joint, 1, min),
           min_match_2 = colnames(jsd_matrix_joint)[min.match.2],
           min_match_2_jsd = apply(jsd_matrix_joint, 2, min))

min_match_df$original_cell_class <- cell_data[match(min_match_df$original_cell_type, cell_data$cell_type),]$cell_class
min_match_df$min_match_1_cell_class <- cell_data[match(min_match_df$min_match_1, cell_data$cell_type),]$cell_class
min_match_df$min_match_2_cell_class <- cell_data[match(min_match_df$min_match_2, cell_data$cell_type),]$cell_class

min_match_df$original_div_stage <- cell_data[match(min_match_df$original_cell_type, cell_data$cell_type),]$div_stage
min_match_df$min_match_1_div_stage <- cell_data[match(min_match_df$min_match_1, cell_data$cell_type),]$div_stage
min_match_df$min_match_2_div_stage <- cell_data[match(min_match_df$min_match_2, cell_data$cell_type),]$div_stage

mismatched_cells <- min_match_df[min_match_df$original_cell_type != min_match_df$min_match_1 &
               min_match_df$original_cell_type != min_match_df$min_match_2,]

mismatched_cells$original_div_stage == mismatched_cells$min_match_1_div_stage | mismatched_cells$original_div_stage == mismatched_cells$min_match_2_div_stage

# 2 out of 125 terminal cell types mismatched both ways
# 12 out of 125 for 1
# 19 out of 125 for 2

# div_stage from cell_data
# then order from OrderForPlot
progenitor_plot_order <- cell_data %>% filter(cell_class == "progenitor") %>% arrange(div_stage, lineage_group) %>% select(cell_type)

cell_data$linaege_group <- factor(cell_data$lineage_group, levels = unique(cell_data$lineage_group))
pro_lineage_groups <- sort(unique(cell_data[cell_data$cell_class == "progenitor",]$lineage_group))

progenitor_plot_order <- c()
for(div_stage in levels(cell_data$div_stage)) {
  for(lineage_group in pro_lineage_groups[pro_lineage_groups %in% cell_data[cell_data$div_stage == div_stage & cell_data$cell_class == "progenitor",]$lineage_group]) {
    cells <- cell_data[cell_data$div_stage == div_stage & cell_data$lineage_group == lineage_group & cell_data$cell_class == "progenitor",]$cell_type
    if(length(cells) > 1) {
      progenitor_plot_order <- c(progenitor_plot_order, cells[hclust(as.dist(jsd_matrix_joint[cells, cells]))$order])
    } else {
      progenitor_plot_order <- c(progenitor_plot_order, cells)
    }
  }
}

joint_plot_order <- c(progenitor_plot_order, OrderForPlot)

saveRDS(joint_plot_order, paste0(dir, "Objects/joint_plot_order.rds"))
joint_plot_order <- readRDS(paste0(dir, "Objects/joint_plot_order.rds"))

# hclust the progenitors within div_stage and lineage_group
cols = list(cell_class = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                           'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                           'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'progenitor' = "grey30"),
            cell_class_lineage_group = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                         'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                         'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
                                         ' ABalpa/ABaraa (Pharynx)' = "#b07aa1", '28 cell or earlier' = "grey30", 'ABpxa/ABarp (Epidermis)' = "#59a14f",
                                         'ABpxp (Neurons)' = "#4e79a7", 'Cx' = "grey30", 'Cxa (Epidermis)' = "#59a14f", 'Cxp/D (Muscle)' = "#76b7b2",
                                         'E (Intestine)' = "#f28e2b", 'MSxa (Pharynx)' = "#b07aa1", 'MSxp (Muscle)' = "#76b7b2", 'Other AB' = "grey30"),
            div_stage = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                          '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF"))

plot_df <- cell_data[joint_plot_order ,c("cell_type", "cell_class", "lineage_group", "div_stage"), drop = FALSE]
plot_df$cell_class_lineage_group <- ifelse(plot_df$cell_class == "progenitor", plot_df$lineage_group, plot_df$cell_class)

for(i in 1:ncol(plot_df)) {
  plot_df[[i]] <- factor(plot_df[[i]], levels = unique(plot_df[[i]]))
}

plot_df$cell_class_lineage_group <- ifelse(plot_df$cell_class == "progenitor", as.character(plot_df$lineage_group), as.character(plot_df$cell_class))
for(i in 1:nrow(plot_df)) {
  plot_df[i,]$cell_class_lineage_group <- ifelse(plot_df[i,]$cell_class == "progenitor", as.character(plot_df[i,]$lineage_group), as.character(plot_df[i,]$cell_class))
}

pdf(paste0(dir, "Plots/supplemental_plots/jsd_matrix_joint.pdf"), width = 10, height = 10)
pheatmap(jsd_matrix_joint[joint_plot_order, joint_plot_order],
         annotation_row = plot_df[joint_plot_order, c("div_stage", "cell_class_lineage_group"), drop = FALSE],
         annotation_colors = cols,
         color = rev(inferno(1000)),
         breaks = seq(0.35, 0.80, 0.45/1000),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         dendrogram = "none",
         scale = "none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_names_row = FALSE,
         labels_row_side = "right",
         labels_col_side = "right",
         border_color = "black", 
         border_width = 1,
         annotation_legend = FALSE,
         legend = FALSE,
         legend_title = "Jensen-Shannon Distance",
         legend_position = "left")
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/jsd_matrix_joint_legend.pdf"), width = 10, height = 10)
pheatmap(jsd_matrix_joint[joint_plot_order, joint_plot_order],
         annotation_row = plot_df[joint_plot_order, c("div_stage", "cell_class_lineage_group"), drop = FALSE],
         annotation_colors = cols,
         color = rev(inferno(1000)),
         breaks = seq(0.35, 0.80, 0.45/1000),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         dendrogram = "none",
         scale = "none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_names_row = FALSE,
         labels_row_side = "right",
         labels_col_side = "right",
         border_color = "black", 
         border_width = 1,
         annotation_legend = TRUE,
         legend = FALSE,
         legend_title = "Jensen-Shannon Distance",
         legend_position = "left")
dev.off()


############################################
# Within species Pro + Term matrix version #
############################################

TPMJoint <- list()
TPMJoint[["C.elegans"]] <- cbind(TPMListBootstrap_pro[["C.elegans"]][,sort(colnames(TPMListBootstrap_pro[["C.elegans"]]))], TPMListBootstrapMean_term[["C.elegans"]][,OrderForPlot])
TPMJoint[["C.briggsae"]] <- cbind(TPMListBootstrap_pro[["C.briggsae"]][,sort(colnames(TPMListBootstrap_pro[["C.elegans"]]))], TPMListBootstrapMean_term[["C.briggsae"]][,OrderForPlot])

pseudocount = 1/dim(TPMJoint[["C.elegans"]])[2]
jsd_matrix_cel <- matrix(nrow = dim(TPMJoint[["C.elegans"]])[2],
                           ncol = dim(TPMJoint[["C.elegans"]])[2])
rownames(jsd_matrix_cel) <- colnames(TPMJoint[["C.elegans"]])
colnames(jsd_matrix_cel) <- colnames(TPMJoint[["C.elegans"]])
for(cel_cell in rownames(jsd_matrix_cel)) {
  for(cbr_cell in colnames(jsd_matrix_cel)) {
    p = TPMJoint[["C.elegans"]][,cel_cell] + pseudocount
    q = TPMJoint[["C.elegans"]][,cbr_cell] + pseudocount
    jsd_matrix_cel[cel_cell, cbr_cell] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
}

pseudocount = 1/dim(TPMJoint[["C.briggsae"]])[2]
jsd_matrix_cbr <- matrix(nrow = dim(TPMJoint[["C.briggsae"]])[2],
                         ncol = dim(TPMJoint[["C.briggsae"]])[2])
rownames(jsd_matrix_cbr) <- colnames(TPMJoint[["C.briggsae"]])
colnames(jsd_matrix_cbr) <- colnames(TPMJoint[["C.briggsae"]])
for(cel_cell in rownames(jsd_matrix_cbr)) {
  for(cbr_cell in colnames(jsd_matrix_cbr)) {
    p = TPMJoint[["C.briggsae"]][,cel_cell] + pseudocount
    q = TPMJoint[["C.briggsae"]][,cbr_cell] + pseudocount
    jsd_matrix_cbr[cel_cell, cbr_cell] <- sqrt(js_divg(p = p/sum(p), q = q/sum(q)))
  }
}


pdf(paste0(dir, "Plots/supplemental_plots/jsd_matrix_cel.pdf"), width = 10, height = 10)
pheatmap(jsd_matrix_cel[joint_plot_order, joint_plot_order],
         annotation_row = plot_df[joint_plot_order, c("div_stage", "cell_class_lineage_group"), drop = FALSE],
         annotation_colors = cols,
         color = rev(inferno(1000)),
         breaks = seq(0.2, 0.825, 0.625/1000),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         dendrogram = "none",
         scale = "none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_names_row = FALSE,
         labels_row_side = "right",
         labels_col_side = "right",
         border_color = "black", 
         border_width = 1,
         annotation_legend = FALSE,
         legend = FALSE,
         legend_title = "Jensen-Shannon Distance",
         legend_position = "left")
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/jsd_matrix_cbr.pdf"), width = 10, height = 10)
pheatmap(jsd_matrix_cbr[joint_plot_order, joint_plot_order],
         annotation_row = plot_df[joint_plot_order, c("div_stage", "cell_class_lineage_group"), drop = FALSE],
         annotation_colors = cols,
         color = rev(inferno(1000)),
         breaks = seq(0.2, 0.825, 0.625/1000),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         dendrogram = "none",
         scale = "none",
         show_rownames = FALSE,
         show_colnames = FALSE,
         annotation_names_row = FALSE,
         labels_row_side = "right",
         labels_col_side = "right",
         border_color = "black", 
         border_width = 1,
         annotation_legend = FALSE,
         legend = FALSE,
         legend_title = "Jensen-Shannon Distance",
         legend_position = "left")
dev.off()

###########################
# Off-diagonal comparison #
###########################

jsd_df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(jsd_df) <- c("cell_type_1", "cell_type_2", "jsd")
for(i in 1:nrow(jsd_matrix_joint)) {
  for(j in 1:ncol(jsd_matrix_joint)) {
    jsd_df <- rbind(jsd_df, data.frame(cell_type_1 = joint_plot_order[i],
                                       cell_type_2 = joint_plot_order[j],
                                       jsd = jsd_matrix_joint[i,j]))
  }
}

jsd_df <- left_join(jsd_df, cell_data[,c("cell_type", "cell_class", "div_stage", "lineage_group")], by = c("cell_type_1" = "cell_type"))
jsd_df <- left_join(jsd_df, cell_data[,c("cell_type", "cell_class", "div_stage", "lineage_group")], by = c("cell_type_2" = "cell_type"), suffix = c("_1", "_2"))

jsd_df$comparison_type <- ifelse(jsd_df$cell_type_1 == jsd_df$cell_type_2, "within", "between")

jsd_df_plot <- jsd_df %>% group_by(div_stage_1, div_stage_2, comparison_type) %>%
  summarise(jsd_mean = mean(jsd), jsd_sd = sd(jsd), jsd_median = median(jsd), jsd_iqr = IQR(jsd))

## Within group cis versus trans
pdf(paste0(dir, "Plots/cell_figure_plots/jsd_matrix_cis_trans_comparison.pdf"), height = 5, width = 5)
ggplot(jsd_df_plot[jsd_df_plot$div_stage_1 == jsd_df_plot$div_stage_2,], aes(x = div_stage_1, y = jsd_mean, ymin = jsd_mean - jsd_sd, ymax = jsd_mean + jsd_sd,
                        color = comparison_type, label = paste0(div_stage_1, "_", div_stage_2))) +
  geom_pointrange(position = position_jitterdodge(dodge.width = 0.25)) +
  # geom_text_repel(max.overlaps = Inf) +
  scale_x_discrete(name = "Division stage", labels = c("15 cell", "28 cell", "50 cell", "100 cell", "200 cell", "350 cell", "600+ cell")) +
  scale_y_continuous(name = "Cell distance (Jensen-Shannon Distance)") +
  scale_color_manual(values = c('within' = "#0d663e", 'between' = "#c75f2b"), limits = c('within', "between"), labels = c("Same cell type", "Diff. cell type")) +
  theme(legend.title = element_blank(),
        legend.position = "top",
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
dev.off()


pdf(paste0(dir, "Plots/cell_figure_plots/jsd_matrix_div_stage_comparison.pdf"), height = 5, width = 5)
ggplot(jsd_df_plot, aes(x = div_stage_1, y = jsd_mean, ymin = jsd_mean - jsd_sd, ymax = jsd_mean + jsd_sd,
                        color = div_stage_2, label = paste0(div_stage_1, "_", div_stage_2))) +
  geom_pointrange(position = position_jitterdodge()) +
  geom_text_repel(max.overlaps = Inf) +
  scale_x_discrete(name = "Division stage", labels = c("15 cell", "28 cell", "50 cell", "100 cell", "200 cell", "350 cell", "600+ cell")) +
  scale_y_continuous(name = "Cell distance (Jensen-Shannon Distance") +
  scale_color_manual(values = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF")) +
  # scale_alpha_manual(values = c("15" = 0.9, "28" = 0.65, "50" = 0.5, "100" = 0.5, "200" = 0.5, "350" = 0.5, "600" = 0.5)) +
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
dev.off()

########################
# JSD versus DEG plots #
########################

pdf(paste0(dir, "Plots/cell_figure_plots/deg_jsd.pdf"), height = 5, width = 5)
ggplot(data = cell_data, 
       aes(x = deg,
           y = jsd_median,
           label = cell_type,
           color = cell_class)) + 
  geom_point(alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  annotate("text", x = c(-Inf), y = c(Inf), hjust = 1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  scale_x_continuous(name = "Number of DEG between species") +
  scale_y_continuous(name = "Cell distance (Jensen-Shannon Distance)") +
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

temp_cell_data <- cell_data
temp_cell_data$mean_expressed <- apply(temp_cell_data[,c("genes_detected_bootstrap_cel", "genes_detected_bootstrap_cbr")], 1, mean)
temp_cell_data$deg_over_expressed_boot <- temp_cell_data$deg / temp_cell_data$mean_expressed

reg <- lm(formula = jsd_median ~ deg_over_expressed_boot,
        data = temp_cell_data)

#get intercept and slope value
coeff <- coefficients(reg)          
intercept <-coeff[1]
slope <- coeff[2]

pdf(paste0(dir, "Plots/cell_figure_plots/deg_over_expressed_boot_jsd.pdf"), height = 5, width = 5)
ggplot(data = temp_cell_data, 
       aes(x = deg_over_expressed_boot,
           y = jsd_median,
           label = cell_type,
           color = cell_class)) + 
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype="dashed") +
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  annotate("text", x = c(-Inf), y = c(Inf), hjust = -0.1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  scale_x_continuous(name = expression(frac("Number of DEG between species", "Number of genes detected (Lower TPM CI)"))) +
  scale_y_continuous(name = "Cell distance(Jensen-Shannon Distance)") +
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

##############
# Gini plots #
##############

reg <- lm(formula = jsd_median ~ cel_gini_median,
          data = cell_data_bins[! cell_data_bins$cell_type %in% c("ABarppxa", "hyp3"),])                      

#get intercept and slope value
coeff <- coefficients(reg)          
intercept <-coeff[1]
slope <- coeff[2]

pdf(paste0(dir, "Plots/cell_figure_plots/cel_gini_jsd_cell.pdf"), height = 8.5, width = 7.5)
ggplot(data = cell_data_bins[! cell_data_bins$cell_type %in% c("ABarppxa", "hyp3"),], aes(x = cel_gini_median,
                             y = jsd_median,
                             label = cell_type, color = cell_class)) + 
  geom_point() +
  geom_text_repel(max.overlaps = 10, show.legend = FALSE) +
  annotate("text", x = c(-Inf), y = c(Inf), hjust = -0.1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype="dashed") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'Progenitor' = "grey20")) +
  scale_x_continuous(name = "Gini coefficient (C. elegans)") +
  scale_y_continuous(name = "Jensen-Shannon Distance") +
  guides(color = guide_legend(ncol = 4), label = "none") +
  theme(legend.title= element_blank(),
        legend.position = "top",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color = NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

reg <- lm(formula = jsd_median ~ cbr_gini_median,
          data = cell_data_bins[! cell_data_bins$cell_type %in% c("ABarppxa", "hyp3"),])

#get intercept and slope value
coeff <- coefficients(reg)          
intercept <-coeff[1]
slope <- coeff[2]

pdf(paste0(dir, "Plots/cell_figure_plots/cbr_gini_jsd_cell.pdf"), height = 8.5, width = 7.5)
ggplot(data = cell_data_bins[! cell_data_bins$cell_type %in% c("ABarppxa", "hyp3"),], aes(x = cbr_gini_median,
                             y = jsd_median,
                             label = cell_type,
                             color = cell_class)) + 
  geom_point() +
  geom_text_repel(max.overlaps = 10, show.legend = FALSE) +
  annotate("text", x = c(-Inf), y = c(Inf), hjust = -0.1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype="dashed") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'Progenitor' = "grey20")) +
  scale_x_continuous(name = "Gini coefficient (C. briggsae)") +
  scale_y_continuous(name = "Jensen-Shannon Distance") +
  guides(color = guide_legend(ncol = 4), label = "none") +
  theme(legend.title= element_blank(),
        legend.position = "top",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color = NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

cell_data$both_gini_median <- ((cell_data$cel_gini_median + cell_data$cbr_gini_median)/2)

reg <- lm(formula = jsd_median ~ both_gini_median,
          data = cell_data)                      

#get intercept and slope value
coeff <- coefficients(reg)          
intercept <-coeff[1]
slope <- coeff[2]

pdf(paste0(dir, "Plots/cell_figure_plots/both_gini_jsd_cell.pdf"), height = 8.5, width = 7.5)
ggplot(data = cell_data, aes(x = both_gini_median,
                             y = jsd_median,
                             label = cell_type,
                             color = cell_class)) + 
  geom_point() +
  geom_text_repel(max.overlaps = 5, show.legend = FALSE) +
  annotate("text", x = c(-Inf), y = c(Inf), hjust = -0.1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype="dashed") +
  scale_x_continuous(name = "Gini coefficient (Both)") +
  scale_y_continuous(name = "Cell Distance (Jensen-Shannon Distance)") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'Progenitor' = "grey20")) +
  guides(color = guide_legend(ncol = 4), label = "none") +
  theme(legend.title= element_blank(),
        legend.position = "top",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color = NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

########################
# Terminal point range #
########################

cell_class_plot_order <- c("Germline", "Muscle", "Intestine", "Hypodermis and seam", "Mesoderm", "Pharynx and rectal", "Ciliated neurons", "Glia and excretory", "Non-ciliated neurons", "progenitor")
cell_class_plot_labels <- c("Germline", "Muscle", "Intestine", "Hypodermis and seam", "Mesoderm", "Pharynx and rectal", "Ciliated neurons", "Glia and excretory", "Non-ciliated neurons", "Progenitors")

jsd_cell_class_pointrange <- cell_data %>%
  filter(! cell_type %in% "ABarppxa") %>%
  arrange(jsd_median) %>%
  ggplot() + 
  geom_boxplot(aes(x = forcats::fct_reorder(cell_class, jsd_median), y = jsd_median), outlier.shape=NA) +
  geom_pointrange(aes(ymin = jsd_lower, ymax = jsd_upper,
                      y = jsd_median, x = forcats::fct_reorder(cell_class, jsd_median), group = cell_class, color = cell_class),
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

pdf(paste0(dir, "Plots/cell_figure_plots/jsd_cell_class_pointrange.pdf"), height = 5, width = 7.5)
print(jsd_cell_class_pointrange)
dev.off()

##########################
# Progenitor point range #
##########################

# <- Light Dark ->
#C6DBEFFF, #9ECAE1FF, #6BAED6FF, #4292C6FF, #2171B5FF, #08519CFF, #08306BFF

div_stage_plot_order <- c("15", "28", "50", "100", "200", "350", "600")
div_stage_plot_labels <- c("15 cell", "28 cell", "50 cell", "100 cell", "200 cell", "350 cell", "600+ cell")

jsd_div_stage_pointrange <- cell_data %>%
  filter(! cell_type %in% "ABarppxa") %>%
  arrange(jsd_median) %>%
  ggplot() +
  geom_boxplot(aes(x = div_stage, y = jsd_median), outlier.shape = NA) +
  geom_pointrange(aes(ymin = jsd_lower, ymax = jsd_upper,
                      y = jsd_median, x = div_stage, group = div_stage, color = div_stage,
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

pdf(paste0(dir, "Plots/cell_figure_plots/jsd_div_stage_pointrange.pdf"), height = 5, width = 7.5)
print(jsd_div_stage_pointrange)
dev.off()

pdf(paste0(dir, "Plots/cell_figure_plots/jsd_cell_class_div_stage.pdf"), height = 5, width = 12.5)
plot_grid(jsd_cell_class_pointrange, jsd_div_stage_pointrange,
          nrow = 1, align = "hv")
dev.off()

pdf(paste0(dir, "Plots/cell_figure_plots/cor_div_stage_pointrange.pdf"), height = 5, width = 7.5)
cell_data %>%
  arrange(cor_median) %>%
  ggplot() + 
  geom_boxplot(aes(x = div_stage, y = cor_median), outlier.shape = NA) +
  geom_pointrange(aes(ymin = cor_lower, ymax = cor_upper,
                      y = cor_median, x = div_stage, group = div_stage, color = div_stage,
                      alpha = div_stage),
                  position = position_dodge2(width = 1), size = 0.5, stroke = 0.25) +
  scale_y_continuous(name = "Cell correlation (Pearson)") +
  scale_x_discrete(name = "Division stage", labels = c("15 cell", "28 cell", "50 cell", "100 cell", "200 cell", "350 cell", "600+ cell")) +
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
dev.off()

pdf(paste0(dir, "Plots/cell_figure_plots/cos_div_stage_pointrange.pdf"), height = 5, width = 7.5)
cell_data %>%
  arrange(cos_median) %>%
  ggplot() + 
  geom_boxplot(aes(x = div_stage, y = cos_median), outlier.shape = NA) +
  geom_pointrange(aes(ymin = cos_lower, ymax = cos_upper,
                      y = cos_median, x = div_stage, group = div_stage, color = div_stage,
                      alpha = div_stage),
                  position = position_dodge2(width = 1), size = 0.5, stroke = 0.25) +
  scale_y_continuous(name = "Cosine angle") +
  scale_x_discrete(name = "Division stage", labels = c("15 cell", "28 cell", "50 cell", "100 cell", "200 cell", "350 cell", "600+ cell")) +
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
dev.off()

######################
# Marker point range #
######################

temp <- cell_data[, c("jsd_median", "cor_median", "cell_type", "cell_class")] %>% filter(cell_type != "ABarppxa")


temp$cor_median <- 1 - temp$cor_median

reg <- lm(formula = jsd_median ~ cor_median,
          data = temp)                  

#get intercept and slope value
coeff <- coefficients(reg)          
intercept <-coeff[1]
slope <- coeff[2]

pdf(paste0(dir, "Plots/cell_figure_plots/jsd_cor.pdf"), height = 5, width = 5)
cell_data %>% filter(! is.na(jsd_median)) %>%
  filter(cell_type != "ABarppxa") %>%
  ggplot(data = ., aes(x = 1 - cor_median,
                       y = jsd_median,
                       label = cell_type,
                       color = plot,
                       fill = plot)) + 
  geom_point() +
  geom_text_repel(max.overlaps = 10, show.legend = FALSE) +
  annotate("text", x = c(-Inf), y = c(Inf), hjust = -0.1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype="dashed") +
  scale_x_continuous(name = "1 - Pearson correlation") +
  scale_y_continuous(name = "Cell distance (Jensen-Shannon Distance)") +
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
  guides(color = guide_legend(ncol = 3)) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

################################
# Shared markers / All markers #
################################

temp <- cell_data[, c("jsd_median", "shared_markers", "either_markers_common", "cell_type", "cell_class")]

temp$shared_over_either <- temp$shared_markers / temp$either_markers_common

reg <- lm(formula = jsd_median ~ shared_over_either,
          data = temp)

#get intercept and slope value
coeff <- coefficients(reg)          
intercept <-coeff[1]
slope <- coeff[2]

pdf(paste0(dir, "Markers/Prop_shared_either_markers_common_jsd.pdf"), height = 5, width = 5)
ggplot(data = temp, aes(x = shared_markers/either_markers_common, y = jsd_median, label = cell_type, color = cell_class)) + 
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  annotate("text", x = c(-Inf), y = c(0.725), hjust = -0.1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype="dashed") +
  scale_x_continuous(name = expression(frac("Shared cell type markers across species",
                                            "Number of unique cell type markers across both species"))) +
  scale_y_continuous(name = "Cell distance (Jensen-Shannon Distance)") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx' = "#b07aa1", 'Rectal cells' = "#ff9da7")) +
  guides(color = guide_legend(ncol = 3), label = "none") +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

#########################
# Aggregate marker info #
#########################

## Number of markers in each cell type
# Number of markers that are private genes to either species
marker_count <- full_join(data.frame(table(joint_markers[["C.elegans"]]$cell_type)),
                          data.frame(table(joint_markers[["C.briggsae"]]$cell_type)),
                          by = "Var1")
colnames(marker_count) <- c("cell_type", "cel_marker_count", "cbr_marker_count")

# add shared markers
marker_count <- left_join(marker_count, 
                          joint_markers[["C.elegans"]] %>%
                            group_by(cell_type) %>%
                            summarise(sum(in_briggsae)),
                          by = c("cell_type" = "cell_type"))

joint_markers_both <- rbind(joint_markers[["C.elegans"]][, c("cell_gene", "cell_type", "gene.type")],
                            joint_markers[["C.briggsae"]][, c("cell_gene", "cell_type", "gene.type")])

# add either markers
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers_both[! duplicated(joint_markers_both$cell_gene),]$cell_type)),
                          by = c("cell_type" = "Var1"))

## Also need # of markers that are from 1:1 versus private
# cel_markers_common
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.elegans"]][joint_markers[["C.elegans"]]$gene.type == "common",]$cell_type)),
                          by = c("cell_type" = "Var1"))

# cbr_markers_common
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.briggsae"]][joint_markers[["C.briggsae"]]$gene.type == "common",]$cell_type)),
                          by = c("cell_type" = "Var1"))

# cel_markers_private
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.elegans"]][joint_markers[["C.elegans"]]$gene.type != "common",]$cell_type)),
                          by = c("cell_type" = "Var1"))

# cbr_markers_private
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.briggsae"]][joint_markers[["C.briggsae"]]$gene.type != "common",]$cell_type)),
                          by = c("cell_type" = "Var1"))

# either_markers_common
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers_both[! duplicated(joint_markers_both$cell_gene) & joint_markers_both$gene.type == "common",]$cell_type)),
                          by = c("cell_type" = "Var1"))

# 1:2
# 1:many
# many:many

# cel_markers_private_category 1:2
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.elegans"]][#joint_markers[["C.elegans"]]$gene.type != "common" &
                                                                          joint_markers[["C.elegans"]]$cel_OG_count == 1 &
                                                                          joint_markers[["C.elegans"]]$cbr_OG_count == 2,]$cell_type)),
                          by = c("cell_type" = "Var1"))

# cel_markers_private_category 1:many
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.elegans"]][#joint_markers[["C.elegans"]]$gene.type != "common" &
                                                                          joint_markers[["C.elegans"]]$cel_OG_count == 1 &
                                                                          joint_markers[["C.elegans"]]$cbr_OG_count > 2,]$cell_type)),
                          by = c("cell_type" = "Var1"))

# cel_markers_private_category many:many
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.elegans"]][#joint_markers[["C.elegans"]]$gene.type != "common" &
                                                                          joint_markers[["C.elegans"]]$cel_OG_count > 2 &
                                                                          joint_markers[["C.elegans"]]$cbr_OG_count > 2,]$cell_type)),
                          by = c("cell_type" = "Var1"))

# cel_markers_private_category any:0
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.elegans"]][#joint_markers[["C.elegans"]]$gene.type != "common" &
                                                                          joint_markers[["C.elegans"]]$cel_OG_count > 0 &
                                                                          joint_markers[["C.elegans"]]$cbr_OG_count == 0,]$cell_type)),
                          by = c("cell_type" = "Var1"))

# cbr_markers_private_category 1:2
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.briggsae"]][#joint_markers[["C.briggsae"]]$gene.type != "common" &
                                                                          joint_markers[["C.briggsae"]]$cel_OG_count == 2 &
                                                                          joint_markers[["C.briggsae"]]$cbr_OG_count == 1,]$cell_type)),
                          by = c("cell_type" = "Var1"))

# cbr_markers_private_category 1:many
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.briggsae"]][#joint_markers[["C.briggsae"]]$gene.type != "common" &
                                                                           joint_markers[["C.briggsae"]]$cel_OG_count > 2 &
                                                                           joint_markers[["C.briggsae"]]$cbr_OG_count == 1,]$cell_type)),
                          by = c("cell_type" = "Var1"))

# cbr_markers_private_category many:many
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.briggsae"]][#joint_markers[["C.briggsae"]]$gene.type != "common" &
                                                                           joint_markers[["C.briggsae"]]$cel_OG_count > 2 &
                                                                           joint_markers[["C.briggsae"]]$cbr_OG_count > 2,]$cell_type)),
                          by = c("cell_type" = "Var1"))

# cbr_markers_private_category many:0
marker_count <- left_join(marker_count,
                          data.frame(table(joint_markers[["C.briggsae"]][#joint_markers[["C.briggsae"]]$gene.type != "common" &
                                                                           joint_markers[["C.briggsae"]]$cel_OG_count == 0 &
                                                                           joint_markers[["C.briggsae"]]$cbr_OG_count > 0,]$cell_type)),
                          by = c("cell_type" = "Var1"))

marker_count <- left_join(marker_count, cell_data[,c("cell_type", "cell_class", "div_stage", "lineage_group", "embryo_time")])

colnames(marker_count) <- c("cell_type", "cel_markers", "cbr_markers", "shared_markers", "either_markers_with_species_specific",
                            "cel_markers_common", "cbr_markers_common", "cel_markers_private", "cbr_markers_private",
                            "either_markers_common", "cel_markers_1_2", "cel_markers_1_many", "cel_markers_many_many", "cel_markers_many_0",
                            "cbr_markers_1_2", "cbr_markers_1_many", "cbr_markers_many_many", "cbr_markers_many_0",
                            "cell_class", "div_stage", "lineage_group", "embryo_time")
marker_count["germline",]$embryo_time <- 800

rownames(marker_count) <- marker_count$cell_type

saveRDS(marker_count, paste0(dir, "Objects/marker_count.rds"))
marker_count <- readRDS(paste0(dir, "Objects/marker_count.rds"))

####################
# Plot marker info #
####################

marker_count$plot <- marker_count$cell_class
marker_count[which(marker_count$cell_class == "progenitor"),]$plot <- as.character(marker_count[which(marker_count$cell_class == "progenitor"),]$div_stage)

reg <- lm(formula = cbr_markers ~ cel_markers,
          data = marker_count)

#get intercept and slope value
coeff <- coefficients(reg)          
intercept <-coeff[1]
slope <- coeff[2]

pdf(paste0(dir, "Plots/cell_figure_plots/num_markers.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = cel_markers,
                                y = cbr_markers,
                                label = cell_type,
                                color = plot,
                                fill = plot)) + 
  geom_point(pch = 21, alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  annotate("text", x = c(-Inf), y = c(Inf), hjust = -0.1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype="dashed") +
  scale_x_continuous(name = expression(paste(italic("C. elegans"), " marker count")), limits = c(min(min(marker_count$cel_markers), min(marker_count$cbr_markers)),
                                                                                                 max(max(marker_count$cel_markers), max(marker_count$cbr_markers)))) +
  scale_y_continuous(name = expression(paste(italic("C. briggsae"), " marker count")), limits = c(min(min(marker_count$cel_markers), min(marker_count$cbr_markers)),
                                                                                                  max(max(marker_count$cel_markers), max(marker_count$cbr_markers)))) +
  guides(color = guide_legend(ncol = 3), label = "none") +
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

marker_count$mean_markers <- rowMeans(marker_count[,c("cel_markers", "cbr_markers")])

pdf(paste0(dir, "Plots/cell_figure_plots/num_markers_embryo_time.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = embryo_time,
                                y = mean_markers,
                                label = cell_type,
                                color = div_stage)) +
  geom_jitter(alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = "Median embryo time") +
  scale_y_continuous(name = "Mean marker count") +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF")) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

pdf(paste0(dir, "Plots/cell_figure_plots/num_markers_embryo_time_cell_class.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = embryo_time,
                                y = mean_markers,
                                label = cell_type,
                                color = ifelse(marker_count$cell_class == "progenitor", unfactor(marker_count$div_stage), marker_count$cell_class))) +
  geom_jitter(alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = "Median embryo time") +
  scale_y_continuous(name = "Marker count") +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
                                '15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF",
                                "600" = "#08306BFF")) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()


marker_count$mean_markers_private <- rowMeans(marker_count[,c("cel_markers_private", "cbr_markers_private")])

pdf(paste0(dir, "Plots/cell_figure_plots/num_private_markers_embryo_time.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = embryo_time,
                                y = mean_markers_private,
                                label = cell_type,
                                color = div_stage)) +
  geom_jitter(alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = "Median embryo time") +
  scale_y_continuous(name = "Mean non-1:1 marker count") +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF")) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()


pdf(paste0(dir, "Plots/cell_figure_plots/num_private_markers_embryo_time_cell_class.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = embryo_time,
                                y = mean_markers_private,
                                label = cell_type,
                                color = ifelse(marker_count$cell_class == "progenitor", marker_count$lineage_group, marker_count$cell_class))) +
  geom_jitter(alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = "Median embryo time") +
  scale_y_continuous(name = "Mean non-1:1 marker count") +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
                                ' ABalpa/ABaraa (Pharynx)' = "#b07aa1", '28 cell or earlier' = "grey30", 'ABpxa/ABarp (Epidermis)' = "#59a14f",
                                'ABpxp (Neurons)' = "#e15759", 'Cx' = "grey30", 'Cxa (Epidermis)' = "#59a14f", 'Cxp/D (Muscle)' = "#76b7b2",
                                'E (Intestine)' = "#f28e2b", 'MSxa (Pharynx)' = "#b07aa1", 'MSxp (Muscle)' = "#76b7b2", 'Other AB' = "grey30")) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()


pdf(paste0(dir, "Plots/cell_figure_plots/num_private_markers_embryo_time_cell_class.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = embryo_time,
                                y = mean_markers_private,
                                label = cell_type,
                                color = ifelse(marker_count$cell_class == "progenitor", unfactor(marker_count$div_stage), marker_count$cell_class))) +
  geom_jitter(alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = "Median embryo time") +
  scale_y_continuous(name = "Non-1:1 marker count") +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
                                '15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF",
                                "600" = "#08306BFF")) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

marker_count$mean_markers_common <- rowMeans(marker_count[,c("cel_markers_common", "cbr_markers_common")])

pdf(paste0(dir, "Plots/cell_figure_plots/num_common_markers_embryo_time_cell_class.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = embryo_time,
                                y = mean_markers_common,
                                label = cell_type,
                                color = ifelse(marker_count$cell_class == "progenitor", unfactor(marker_count$div_stage), marker_count$cell_class))) +
  geom_jitter(alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = "Median embryo time") +
  scale_y_continuous(name = "1:1 marker count") +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
                                '15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF",
                                "600" = "#08306BFF")) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()


marker_count$cel_prop_private <- marker_count$cel_markers_private / marker_count$cel_markers
marker_count$cbr_prop_private <- marker_count$cbr_markers_private / marker_count$cbr_markers

rownames(marker_count) <- marker_count$cell_type

pdf(paste0(dir, "Plots/cell_figure_plots/proportion_private_markers.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = cel_prop_private * 100,
                                y = cbr_prop_private * 100,
                                label = cell_type,
                                color = cell_class)) + 
  geom_point(alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = expression(paste(italic("C. elegans"), " percent non 1:1 markers"))) +
  scale_y_continuous(name = expression(paste(italic("C. briggsae"), " percent non 1:1 markers"))) +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'progenitor' = "grey30")) +
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

marker_count$mean_prop_private <- rowMeans(marker_count[,c("cel_prop_private", "cbr_prop_private")])

pdf(paste0(dir, "Plots/cell_figure_plots/proportion_private_markers_embryo_time.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = embryo_time,
                                y = mean_prop_private * 100,
                                label = cell_type,
                                color = div_stage)) +
  geom_jitter(alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = "Median embryo time") +
  scale_y_continuous(name = "Percent non-1:1 markers") +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF")) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text = element_text(size = 14),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

pdf(paste0(dir, "Plots/cell_figure_plots/proportion_private_markers_div_stage.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = cel_prop_private * 100,
                                y = cbr_prop_private * 100,
                                label = cell_type,
                                color = div_stage)) + 
  geom_point(alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = expression(paste(italic("C. elegans"), " percent non 1:1 markers"))) +
  scale_y_continuous(name = expression(paste(italic("C. briggsae"), " percent non 1:1 markers"))) +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF")) +
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

pdf(paste0(dir, "Plots/cell_figure_plots/proportion_private_markers_1_2.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = cel_markers_1_2/cel_markers * 100,
                                y = cbr_markers_1_2/cbr_markers * 100,
                                label = cell_type, color = cell_class)) + 
  geom_point(alpha = 0.8) +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = expression(paste(italic("C. elegans"), " percent 1:2 markers"))) +
  scale_y_continuous(name = expression(paste(italic("C. briggsae"), " percent 1:2 markers"))) +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'progenitor' = "grey30")) +
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

pdf(paste0(dir, "Plots/cell_figure_plots/proportion_private_markers_1_many.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = cel_markers_1_many/cel_markers * 100,
                                y = cbr_markers_1_many/cbr_markers * 100,
                                label = cell_type, color = cell_class)) + 
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = expression(paste(italic("C. elegans"), " percent 1:many markers"))) +
  scale_y_continuous(name = expression(paste(italic("C. briggsae"), " percent 1:many markers"))) +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'progenitor' = "grey30")) +
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

pdf(paste0(dir, "Plots/cell_figure_plots/proportion_private_markers_many_many.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = cel_markers_many_many/cel_markers * 100,
                                y = cbr_markers_many_many/cbr_markers * 100,
                                label = cell_type, color = cell_class)) + 
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = expression(paste(italic("C. elegans"), " percent many:many markers"))) +
  scale_y_continuous(name = expression(paste(italic("C. briggsae"), " percent many:many markers"))) +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'progenitor' = "grey30")) +
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

pdf(paste0(dir, "Plots/cell_figure_plots/proportion_private_markers_1_0.pdf"), height = 5, width = 5)
ggplot(data = marker_count, aes(x = cel_markers_many_0/cel_markers * 100,
                                y = cbr_markers_many_0/cbr_markers * 100,
                                label = cell_type, color = cell_class)) + 
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  scale_x_continuous(name = expression(paste(italic("C. elegans"), " percent 1:0 markers"))) +
  scale_y_continuous(name = expression(paste(italic("C. briggsae"), " percent 1:0 markers"))) +
  guides(color = guide_legend(ncol = 3), label = "none") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'progenitor' = "grey30")) +
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

######################
# Marker enrichments #
######################

term_markers <- readRDS(paste0(dir, "Objects/TermMarker_list_metadata_cell_bg.rds"))
pro_markers <- readRDS(paste0(dir, "Objects/ProMarker_list_metadata_cell_bg.rds"))

# Pull enrichment functions from 

############################
# Maker sharing versus jsd #
############################

term_markers <- readRDS(paste0(dir, "Objects/TermMarker_list_metadata_cell_bg.rds"))
pro_markers <- readRDS(paste0(dir, "Objects/ProMarker_list_metadata_cell_bg.rds"))

## Markers
term_markers[["C.elegans"]]$cel_tpm_log2fc_just_pro <- NA
term_markers[["C.elegans"]]$cbr_tpm_log2fc_just_pro <- NA

term_markers[["C.briggsae"]]$cel_tpm_log2fc_just_pro <- NA
term_markers[["C.briggsae"]]$cbr_tpm_log2fc_just_pro <- NA

joint_markers <- list()
joint_markers[["C.elegans"]] <- rbind(term_markers[["C.elegans"]][,colnames(pro_markers[["C.elegans"]])], pro_markers[["C.elegans"]])
joint_markers[["C.briggsae"]] <- rbind(term_markers[["C.briggsae"]][,colnames(pro_markers[["C.briggsae"]])], pro_markers[["C.briggsae"]])


cel_markers <- joint_markers[["C.elegans"]][joint_markers[["C.elegans"]]$p_val_adj.cel < 0.05 &
                                              joint_markers[["C.elegans"]]$avg_log2FC.cel >= 1 &
                                              joint_markers[["C.elegans"]]$cel_tpm_log2fc >= 1 &
                                              joint_markers[["C.elegans"]]$cel_tpm >= 80,]
cel_markers_1to1 <- cel_markers[!is.na(cel_markers$cbr_tpm),]
cel_markers_non_1to1 <- cel_markers[is.na(cel_markers$cbr_tpm) & !is.na(cel_markers$OG) & cel_markers$cbr_OG_count != 0,]
cel_markers_private <-  cel_markers[is.na(cel_markers$cbr_tpm) & (is.na(cel_markers$OG) | cel_markers$cbr_OG_count == 0) ,]
cel_markers_almost <- joint_markers[["C.elegans"]][joint_markers[["C.elegans"]]$p_val_adj.cel < 0.05 &
                                                     joint_markers[["C.elegans"]]$avg_log2FC.cel >= .5 &
                                                     joint_markers[["C.elegans"]]$cel_tpm_log2fc >= .5 &
                                                     joint_markers[["C.elegans"]]$cel_tpm >= 40,]

cbr_markers <- joint_markers[["C.briggsae"]][joint_markers[["C.briggsae"]]$p_val_adj.cbr < 0.05 &
                                               joint_markers[["C.briggsae"]]$avg_log2FC.cbr >= 1 &
                                               joint_markers[["C.briggsae"]]$cbr_tpm_log2fc >= 1 &
                                               joint_markers[["C.briggsae"]]$cbr_tpm >= 80,]
cbr_markers_1to1 <- cbr_markers[!is.na(cbr_markers$cel_tpm),] ##defined usine LaDeanna list)
cbr_markers_non_1to1 <- cbr_markers[is.na(cbr_markers$cel_tpm) & !is.na(cbr_markers$OG) & cbr_markers$cel_OG_count != 0,]
cbr_markers_private <-  cbr_markers[is.na(cbr_markers$cel_tpm) & (is.na(cbr_markers$OG) | cbr_markers$cel_OG_count == 0) ,]
cbr_markers_almost <- joint_markers[["C.briggsae"]][joint_markers[["C.briggsae"]]$p_val_adj.cbr < 0.05 &
                                                      joint_markers[["C.briggsae"]]$avg_log2FC.cbr >= .5 &
                                                      joint_markers[["C.briggsae"]]$cbr_tpm_log2fc >= .5 &
                                                      joint_markers[["C.briggsae"]]$cbr_tpm >= 40,]

shared_markers <- intersect(rownames(cel_markers), rownames(cbr_markers))        
cel_shared_markers <- cel_markers_1to1[shared_markers,]
cbr_shared_markers <- cbr_markers_1to1[shared_markers,]

cel_almost_shared_markers <- cel_markers_1to1[intersect(rownames(cel_markers), rownames(cbr_markers_almost)),]
cbr_almost_shared_markers <- cbr_markers_1to1[intersect(rownames(cbr_markers), rownames(cel_markers_almost)),]

cel_divergent_markers <- cel_markers_1to1[cel_markers_1to1$cbr_tpm < 80 & cel_markers_1to1$cbr_tpm_log2fc < 0 & cel_markers_1to1$cel_tpm/(cel_markers_1to1$cbr_tpm + 1) > 5,]
cbr_divergent_markers <- cbr_markers_1to1[cbr_markers_1to1$cel_tpm < 80 & cbr_markers_1to1$cel_tpm_log2fc < 0 & cbr_markers_1to1$cbr_tpm/(cbr_markers_1to1$cel_tpm + 1) > 5,]

cel_almost_shared_markers_counts <- data.frame(cel_almost_shared_markers %>% group_by(cell_type) %>% summarize(n = n()) %>% arrange(desc(n)))
cel_almost_shared_markers_counts <- cel_almost_shared_markers_counts[! is.na(cel_almost_shared_markers_counts$cell_type),]
rownames(cel_almost_shared_markers_counts) <- cel_almost_shared_markers_counts$cell_type

cbr_almost_shared_markers_counts <- data.frame(cbr_almost_shared_markers %>% group_by(cell_type) %>% summarize(n = n()) %>% arrange(desc(n)))
cbr_almost_shared_markers_counts <- cbr_almost_shared_markers_counts[! is.na(cbr_almost_shared_markers_counts$cell_type),]
rownames(cbr_almost_shared_markers_counts) <- cbr_almost_shared_markers_counts$cell_type
# 
# cell_data$almost_shared_markers <- NA
# cell_data[rownames(almost_shared_markers_counts),]$almost_shared_markers <- almost_shared_markers_counts$n

cel_markers_1to1_count <- data.frame(cel_markers_1to1 %>% group_by(cell_type) %>% summarize(cel_markers_1to1_count = n()))
cel_markers_1to1_count <- cel_markers_1to1_count[! is.na(cel_markers_1to1_count$cell_type),]
rownames(cel_markers_1to1_count) <- cel_markers_1to1_count$cell_type

cbr_markers_1to1_count <- data.frame(cbr_markers_1to1 %>% group_by(cell_type) %>% summarize(cbr_markers_1to1_count = n()))
cbr_markers_1to1_count <- cbr_markers_1to1_count[! is.na(cbr_markers_1to1_count$cell_type),]
rownames(cbr_markers_1to1_count) <- cbr_markers_1to1_count$cell_type

temp <- data.frame(cell_type = cel_almost_shared_markers_counts$cell_type,
           cel_almost_shared_markers_count = cel_almost_shared_markers_counts$n,
           cbr_almost_shared_markers_count = cbr_almost_shared_markers_counts[cel_almost_shared_markers_counts$cell_type,]$n,
           cel_markers_1to1_count = cel_markers_1to1_count[cel_almost_shared_markers_counts$cell_type,]$cel_markers_1to1_count,
           cbr_markers_1to1_count = cbr_markers_1to1_count[cel_almost_shared_markers_counts$cell_type,]$cbr_markers_1to1_count)
rownames(temp) <- temp$cell_type

temp$fraction_almost_shared <- ((temp$cel_almost_shared_markers_count / temp$cel_markers_1to1_count) + (temp$cbr_almost_shared_markers_count / temp$cbr_markers_1to1_count))/2

cell_data$fraction_almost_shared <- NA
cell_data[temp$cell_type,]$fraction_almost_shared <- temp$fraction_almost_shared

reg <- lm(formula = jsd_median ~ fraction_almost_shared,
          data = cell_data[cell_data$cell_class != "progenitor",])

#get intercept and slope value
coeff <- coefficients(reg)
intercept <- coeff[1]
slope <- coeff[2]

pdf(paste0(dir, "Plots/cell_figure_plots/fraction_almost_shared_markers_jsd.pdf"), height = 5, width = 5)
cell_data %>% filter(cell_type != "ABarppxa") %>%
  filter(cell_class != "progenitor") %>%
ggplot(aes(x = fraction_almost_shared, y = jsd_median, label = cell_type, color = cell_class)) +
  geom_point() +
  geom_text_repel(show.legend = FALSE) +
  annotate("text", x = c(0.65), y = c(0.6), hjust = -0.1, vjust = 1, color = "black",
           label = bquote("Adj. " ~ R ^ 2 ~ " = " ~ .(round(summary(reg)$adj.r.squared, 2)))) +
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype = "dashed") +
  scale_y_continuous(name = "Cell distance (Jensen-Shannon Distance)") +
  scale_x_continuous(name = "Fraction of cell type markers shared") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1", 'progenitor' = "grey30")) +
  guides(color = guide_legend(ncol = 3), label = "none") +
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

## Mean fraction markers divergent

cel_divergent_markers_count <- data.frame(cel_divergent_markers %>% group_by(cell_type) %>% summarize(n = n()) %>% arrange(desc(n)))
cel_divergent_markers_count <- cel_divergent_markers_count[! is.na(cel_divergent_markers_count$cell_type),]
rownames(cel_divergent_markers_count) <- cel_divergent_markers_count$cell_type

cbr_divergent_markers_count <- data.frame(cbr_divergent_markers %>% group_by(cell_type) %>% summarize(n = n()) %>% arrange(desc(n)))
cbr_divergent_markers_count <- cbr_divergent_markers_count[! is.na(cbr_divergent_markers_count$cell_type),]
rownames(cbr_divergent_markers_count) <- cbr_divergent_markers_count$cell_type

temp <- data.frame(cell_type = cel_divergent_markers_count$cell_type,
                   cel_divergent_markers_count = cel_divergent_markers_count$n,
                   cbr_divergent_markers_count = cbr_divergent_markers_count[cel_divergent_markers_count$cell_type,]$n,
                   cel_markers_1to1_count = cel_markers_1to1_count[cel_divergent_markers_count$cell_type,]$cel_markers_1to1_count,
                   cbr_markers_1to1_count = cbr_markers_1to1_count[cel_divergent_markers_count$cell_type,]$cbr_markers_1to1_count)
rownames(temp) <- temp$cell_type

temp$fraction_divergent <- ((temp$cel_divergent_markers_count / temp$cel_markers_1to1_count) + (temp$cbr_divergent_markers_count / temp$cbr_markers_1to1_count))/2

cell_data$fraction_divergent <- NA
cell_data[temp$cell_type,]$fraction_divergent <- temp$fraction_divergent

cell_class_plot_order <- c("Germline", "Hypodermis and seam", "Intestine", "Muscle",  "Ciliated neurons", "Pharynx and rectal", "Non-ciliated neurons", "Glia and excretory", "Mesoderm", "progenitor")
cell_class_plot_labels <- c("Germline", "Hypodermis and seam", "Intestine", "Muscle",  "Ciliated neurons", "Pharynx and rectal", "Non-ciliated neurons", "Glia and excretory", "Mesoderm", "Progenitors")

pdf(paste0(dir, "Plots/cell_figure_plots/fraction_divergent_markers.pdf"), height = 4, width = 5)
cell_data %>%
  filter(! cell_type %in% "ABarppxa") %>%
  arrange(fraction_divergent) %>%
  ggplot() + 
  geom_boxplot(aes(x = forcats::fct_reorder(cell_class, fraction_divergent), y = fraction_divergent), outlier.shape=NA) +
  geom_point(aes(y = fraction_divergent, x = forcats::fct_reorder(cell_class, fraction_divergent), group = cell_class, color = plot),
                  position = position_dodge2(width = 1), alpha = 0.5, size = 1.5, stroke = 0.5) +
  scale_x_discrete(limits = cell_class_plot_order, labels = cell_class_plot_labels) +
  scale_y_continuous(name = "Fraction of divergent markers") +
  scale_color_manual(values = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF",
                                'Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
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
dev.off()

#################
# Model fitting #
#################
cell_data$gini_both <- (cell_data$cel_gini_median + cell_data$cbr_gini_median) / 2
cell_data$umi_both <- (cell_data$cel_median_umi + cell_data$cbr_median_umi) / 2
cell_data$both_ci_detected <- (cell_data$genes_detected_bootstrap_cel + cell_data$genes_detected_bootstrap_cbr) / 2

olsrr::ols_step_best_subset(lm(jsd_median ~ ., cell_data[,c("jsd_median", "cell_class", "umi_both", "gini_both", "min_cell_count", "div_stage", "embryo_time", "both_ci_detected")]))

car::Anova(lm(jsd_median ~ cell_class + both_ci_detected + min_cell_count + umi_both + gini_both, cell_data[cell_data$cell_class != "progenitor",]))
car::Anova(lm(jsd_median ~ both_ci_detected + min_cell_count + umi_both + gini_both + div_stage, cell_data[cell_data$cell_class == "progenitor",]))

############################
# Shared:Divergent markers #
############################
marker_types <- readRDS("~jmurr/scRNAseq/briggsae/WS290/out/CellClass_Marker_Category_Counts.rds")

background_gene_set <- unique(unlist(marker_types))

marker_types_list <- list()
for(cat in names(marker_types)) {
  print(cat)
  marker_types_list[[cat]] <- list()
  
  for(cell_class in names(marker_types[[cat]])) {
    print(cell_class)
    
    rgs <- marker_types[[cat]][[cell_class]]
    
    marker_types_list[[cat]][[cell_class]] <- wormcat_enrichment(rgs, background_gene_set, gff_list_mrna[["elegans"]])
    
    for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
      marker_types_list[[cat]][[cell_class]][[tier]]$type <- cat
      marker_types_list[[cat]][[cell_class]][[tier]]$cell_class <- cell_class
      marker_types_list[[cat]][[cell_class]][[tier]]$tier <- tier
    }
  }
}

marker_types_df <- do.call(rbind, do.call(rbind, do.call(rbind, marker_types_list))) %>%
  complete(Category, type, cell_class)
marker_types_df[is.na(marker_types_df$Fold),]$Bonferroni <- 1
marker_types_df[is.na(marker_types_df$Fold),]$Fold <- 0

marker_types_df$signif <- ""
marker_types_df$signif <- ifelse(marker_types_df$Bonferroni < 0.05, "*", marker_types_df$signif)
marker_types_df$signif <- ifelse(marker_types_df$Bonferroni < 0.005, "**", marker_types_df$signif)
marker_types_df$signif <- ifelse(marker_types_df$Bonferroni < 0.0005, "***", marker_types_df$signif)

col_order <- unique(marker_types_df$cell_class)

pdf(paste0(dir, "Plots/big_boy_marker_enrich_all_class.pdf"), height = 10, width = 80)
marker_types_df %>%
  ggplot(aes(x = Category, y = cell_class, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  facet_wrap(~dataset, scale="free_y", nrow=4)+
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = col_order) +
  scale_fill_gradient(name = "Fold\nEnrich.", low = "white", high="blue",
                      breaks = c(0, 1, 2, 3, 4, 5),
                      labels = c("0", "1", "2", "3", "4", ">5"),
                      limits = c(0, 5),
                      oob = scales::squish) +
  # scale_fill_viridis() +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid = element_line(color="grey80", size=0.25),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.position = "none") +
  facet_wrap(~type, nrow = 5)
dev.off()


cats_to_plot <- unique(marker_types_df[marker_types_df$tier %in% "WormCat.1",]$Category)

pdf(paste0(dir, "Plots/big_boy_marker_enrich_all_class_tier_1.pdf"), height = 10, width = 7)
marker_types_df %>%
  filter(Category %in% cats_to_plot) %>%
  ggplot(aes(x = Category, y = cell_class, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  facet_wrap(~dataset, scale="free_y", nrow=4)+
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = col_order) +
  scale_fill_gradient(name = "Fold\nEnrich.", low = "white", high="blue",
                      breaks = c(0, 1, 2, 3, 4, 5),
                      labels = c("0", "1", "2", "3", "4", ">5"),
                      limits = c(0, 5),
                      oob = scales::squish) +
  # scale_fill_viridis() +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid = element_line(color="grey80", size=0.25),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.position = "none") +
  facet_wrap(~type, nrow = 5)
dev.off()

col_order <- c("Ciliated neurons", "Non-ciliated neurons", "Glia and excretory", "Hypodermis and seam", "Pharynx and rectal", "Muscle", "Mesoderm", "Intestine",  "Germline", "progenitor")
col_labels <- c("Ciliated neurons", "Non-ciliated neurons", "Glia and excretory", "Hypodermis and seam", "Pharynx and rectal", "Muscle", "Mesoderm", "Intestine",  "Germline", "Progenitors")
cat_order <- c("Cilia", "Neuronal function", "Extracellular material", "Ribosome", "Muscle function",  "Stress response", "Peroxisome", "Cell cycle")

pdf(paste0(dir, "Plots/big_boy_marker_enrich_all_class_tier_1_filtered.pdf"), height = 6, width = 3.5)
marker_types_df %>%
  filter(Category %in% c("Muscle function", "Neuronal function",
                         "Ribosome", "Cilia", "Cell cycle",
                         "Extracellular material", "Stress response", "Peroxisome")) %>%
  filter(type %in% c("Marker, 1:1 shared", "Marker, 1:1 divergent")) %>%
  ggplot(aes(x = Category, y = cell_class, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  facet_wrap(~dataset, scale="free_y", nrow=4)+
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45), limits = cat_order) +
  scale_y_discrete(limits = col_order, labels = col_labels) +
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
        legend.position = "none") +
  facet_wrap(~type, nrow = 2)
dev.off()

pdf(paste0(dir, "Plots/big_boy_marker_enrich_all_class_tier_1_filtered_horizontal.pdf"), height = 3.5, width = 6)
marker_types_df %>%
  filter(Category %in% c("Muscle function", "Neuronal function",
                         "Ribosome", "Cilia", "Cell cycle",
                         "Extracellular material", "Stress response", "Peroxisome")) %>%
  filter(type %in% c("Marker, 1:1 shared", "Marker, 1:1 divergent")) %>%
  mutate(type = factor(type, levels = c("Marker, 1:1 shared", "Marker, 1:1 divergent"))) %>%
  ggplot(aes(x = Category, y = cell_class, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  facet_wrap(~dataset, scale="free_y", nrow=4)+
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45), limits = cat_order) +
  scale_y_discrete(limits = col_order, labels = col_labels) +
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
        legend.position = "none") +
  facet_wrap(~type, ncol = 2)
dev.off()

saveRDS(marker_types_df, paste0(dir, "Objects/marker_types_df.rds"))

###################
# Non-1:1 markers #
###################

background_gene_set <- unique(unlist(cel_markers$gene))

cel_markers <- left_join(cel_markers, cell_data[,c("cell_type", "div_stage", "cell_class")], by = "cell_type")
cel_markers <- cel_markers[which(cel_markers$cell_type != "hyp3"),]

cbr_markers <- left_join(cbr_markers, cell_data[,c("cell_type", "div_stage", "cell_class")], by = "cell_type")
cbr_markers <- cbr_markers[which(cbr_markers$cell_type != "hyp3"),]

rownames(gff_list_mrna[["elegans"]]) <- ifelse(! is.na(gff_list_mrna[["elegans"]]$cds_gene_name), gff_list_mrna[["elegans"]]$cds_gene_name, gff_list_mrna[["elegans"]]$gene_name)
data.frame(table(gff_list_mrna[["elegans"]][background_gene_set, "WormCat.1"]))

marker_div_stage_list <- list()
for(cur_div_stage in unique(cel_markers$div_stage)) {
  print(cur_div_stage)
  
  rgs <- unique(cel_markers[cel_markers$div_stage %in% cur_div_stage,]$gene)
  marker_div_stage_list[[cur_div_stage]] <- wormcat_enrichment(rgs, background_gene_set, gff_list_mrna[["elegans"]])
  
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    marker_div_stage_list[[cur_div_stage]][[tier]]$div_stage <- cur_div_stage
    marker_div_stage_list[[cur_div_stage]][[tier]]$tier <- tier
  }
}

marker_div_stage_list <- transpose(marker_div_stage_list)

marker_div_stage_list_filt <- list()
for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
  marker_div_stage_list_filt[[tier]] <- do.call(rbind, marker_div_stage_list[[tier]])
  temp_filt <- data.frame(marker_div_stage_list_filt[[tier]] %>%
                            filter(Bonferroni < 0.05) %>%
                            select(Category) %>%
                            unique())[,1]
  
  marker_div_stage_list_filt[[tier]] <- marker_div_stage_list_filt[[tier]] %>%
    filter(Category %in% temp_filt) %>%
    complete(Category, div_stage)
  
  marker_div_stage_list_filt[[tier]][is.na(marker_div_stage_list_filt[[tier]]$Fold),]$Fold <- 0
  marker_div_stage_list_filt[[tier]][is.na(marker_div_stage_list_filt[[tier]]$Bonferroni),]$Bonferroni <- 1
  
  marker_div_stage_list_filt[[tier]]$signif <- ""
  marker_div_stage_list_filt[[tier]]$signif <- ifelse(marker_div_stage_list_filt[[tier]]$Bonferroni < 0.05,
                                                      "*",
                                                      marker_div_stage_list_filt[[tier]]$signif)
  marker_div_stage_list_filt[[tier]]$signif <- ifelse(marker_div_stage_list_filt[[tier]]$Bonferroni < 0.005,
                                                       "**",
                                                       marker_div_stage_list_filt[[tier]]$signif)
  marker_div_stage_list_filt[[tier]]$signif <- ifelse(marker_div_stage_list_filt[[tier]]$Bonferroni < 0.0005,
                                                       "***",
                                                       marker_div_stage_list_filt[[tier]]$signif)
}

pdf(paste0(dir, "Plots/cell_figure_plots/marker_div_stage.pdf"), height = 4, width = 2)
do.call(rbind, marker_div_stage_list_filt) %>%
  filter(Category %in% c("Development", "Cell cycle", "DNA: replication", "Proteolysis proteasome: E3: F box",
                         "Transcription: chromatin structure: histone",
                         "Transcription: general machinery Transcription: general machinery: RNA Pol II")) %>%
  ggplot(aes(x = Category, y = div_stage, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = c("15", "28", "50", "100", "200", "350", "600")) +
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
        legend.position = "none")
dev.off()

background_gene_set <- unique(unlist(cel_markers$gene))

rownames(gff_list_mrna[["elegans"]]) <- ifelse(! is.na(gff_list_mrna[["elegans"]]$cds_gene_name), gff_list_mrna[["elegans"]]$cds_gene_name, gff_list_mrna[["elegans"]]$gene_name)
data.frame(table(gff_list_mrna[["elegans"]][background_gene_set, "WormCat.1"]))

marker_div_stage_list <- list()
for(cur_div_stage in unique(cel_markers$div_stage)) {
  print(cur_div_stage)
  
  rgs <- unique(cel_markers[cel_markers$div_stage %in% cur_div_stage,]$gene)
  marker_div_stage_list[[cur_div_stage]] <- wormcat_enrichment(rgs, background_gene_set, gff_list_mrna[["elegans"]])
  
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    marker_div_stage_list[[cur_div_stage]][[tier]]$div_stage <- cur_div_stage
    marker_div_stage_list[[cur_div_stage]][[tier]]$tier <- tier
  }
}

marker_div_stage_list <- transpose(marker_div_stage_list)

marker_div_stage_list_filt <- list()
for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
  marker_div_stage_list_filt[[tier]] <- do.call(rbind, marker_div_stage_list[[tier]])
  temp_filt <- data.frame(marker_div_stage_list_filt[[tier]] %>%
                            filter(Bonferroni < 0.05) %>%
                            select(Category) %>%
                            unique())[,1]
  
  marker_div_stage_list_filt[[tier]] <- marker_div_stage_list_filt[[tier]] %>%
    filter(Category %in% temp_filt) %>%
    complete(Category, div_stage)
  
  marker_div_stage_list_filt[[tier]][is.na(marker_div_stage_list_filt[[tier]]$Fold),]$Fold <- 0
  marker_div_stage_list_filt[[tier]][is.na(marker_div_stage_list_filt[[tier]]$Bonferroni),]$Bonferroni <- 1
  
  marker_div_stage_list_filt[[tier]]$signif <- ""
  marker_div_stage_list_filt[[tier]]$signif <- ifelse(marker_div_stage_list_filt[[tier]]$Bonferroni < 0.05,
                                                      "*",
                                                      marker_div_stage_list_filt[[tier]]$signif)
  marker_div_stage_list_filt[[tier]]$signif <- ifelse(marker_div_stage_list_filt[[tier]]$Bonferroni < 0.005,
                                                      "**",
                                                      marker_div_stage_list_filt[[tier]]$signif)
  marker_div_stage_list_filt[[tier]]$signif <- ifelse(marker_div_stage_list_filt[[tier]]$Bonferroni < 0.0005,
                                                      "***",
                                                      marker_div_stage_list_filt[[tier]]$signif)
}

pdf(paste0(dir, "Plots/cell_figure_plots/marker_div_stage.pdf"), height = 4, width = 2)
do.call(rbind, marker_div_stage_list_filt) %>%
  filter(Category %in% c("Development", "Cell cycle", "DNA: replication", "Proteolysis proteasome: E3: F box",
                         "Transcription: chromatin structure: histone",
                         "Transcription: general machinery Transcription: general machinery: RNA Pol II")) %>%
  ggplot(aes(x = Category, y = div_stage, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = c("15", "28", "50", "100", "200", "350", "600")) +
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
        legend.position = "none")
dev.off()

background_gene_set <- unique(unlist(cel_markers_non_1to1$gene))

cel_markers_non_1to1 <- left_join(cel_markers_non_1to1, cell_data[,c("cell_type", "div_stage")], by = "cell_type")
cel_markers_non_1to1 <- cel_markers_non_1to1[which(cel_markers_non_1to1$cell_type != "hyp3"),]

marker_non_1to1_div_stage_list <- list()
for(cur_div_stage in unique(cel_markers_non_1to1$div_stage)) {
  print(cur_div_stage)
  
  rgs <- unique(cel_markers_non_1to1[cel_markers_non_1to1$div_stage %in% cur_div_stage,]$gene)
  marker_non_1to1_div_stage_list[[cur_div_stage]] <- wormcat_enrichment(rgs, background_gene_set, gff_list_mrna[["elegans"]])
  
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    marker_non_1to1_div_stage_list[[cur_div_stage]][[tier]]$div_stage <- cur_div_stage
    marker_non_1to1_div_stage_list[[cur_div_stage]][[tier]]$tier <- tier
  }
}

marker_non_1to1_div_stage_list <- transpose(marker_non_1to1_div_stage_list)

marker_non_1to1__div_stage_list_filt <- list()
for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
  marker_non_1to1__div_stage_list_filt[[tier]] <- do.call(rbind, marker_non_1to1_div_stage_list[[tier]])
  temp_filt <- data.frame(marker_non_1to1__div_stage_list_filt[[tier]] %>%
                            filter(Bonferroni < 0.05) %>%
                            select(Category) %>%
                            unique())[,1]
  
  marker_non_1to1__div_stage_list_filt[[tier]] <- marker_non_1to1__div_stage_list_filt[[tier]] %>%
    filter(Category %in% temp_filt) %>%
    complete(Category, div_stage)
  
  marker_non_1to1__div_stage_list_filt[[tier]][is.na(marker_non_1to1__div_stage_list_filt[[tier]]$Fold),]$Fold <- 0
  marker_non_1to1__div_stage_list_filt[[tier]][is.na(marker_non_1to1__div_stage_list_filt[[tier]]$Bonferroni),]$Bonferroni <- 1
  
  marker_non_1to1__div_stage_list_filt[[tier]]$signif <- ""
  marker_non_1to1__div_stage_list_filt[[tier]]$signif <- ifelse(marker_non_1to1__div_stage_list_filt[[tier]]$Bonferroni < 0.05,
                                                      "*",
                                                      marker_non_1to1__div_stage_list_filt[[tier]]$signif)
  marker_non_1to1__div_stage_list_filt[[tier]]$signif <- ifelse(marker_non_1to1__div_stage_list_filt[[tier]]$Bonferroni < 0.005,
                                                      "**",
                                                      marker_non_1to1__div_stage_list_filt[[tier]]$signif)
  marker_non_1to1__div_stage_list_filt[[tier]]$signif <- ifelse(marker_non_1to1__div_stage_list_filt[[tier]]$Bonferroni < 0.0005,
                                                      "***",
                                                      marker_non_1to1__div_stage_list_filt[[tier]]$signif)
}

pdf(paste0(dir, "Plots/cell_figure_plots/marker_non_1to1_div_stage.pdf"), height = 4, width = 1)
do.call(rbind, marker_non_1to1__div_stage_list_filt) %>%
  filter(tier == "WormCat.3") %>%
  ggplot(aes(x = Category, y = div_stage, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = c("15", "28", "50", "100", "200", "350", "600")) +
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
        legend.position = "none")
dev.off()


background_gene_set <- unique(unlist(cel_markers_private$gene))

cel_markers_private <- left_join(cel_markers_private, cell_data[,c("cell_type", "div_stage")], by = "cell_type")
cel_markers_private <- cel_markers_private[which(cel_markers_private$cell_type != "hyp3"),]

marker_private_div_stage_list <- list()
for(cur_div_stage in unique(cel_markers_private$div_stage)) {
  print(cur_div_stage)
  
  rgs <- unique(cel_markers_private[cel_markers_private$div_stage %in% cur_div_stage,]$gene)
  marker_private_div_stage_list[[cur_div_stage]] <- wormcat_enrichment(rgs, background_gene_set, gff_list_mrna[["elegans"]])
  
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    marker_private_div_stage_list[[cur_div_stage]][[tier]]$div_stage <- cur_div_stage
    marker_private_div_stage_list[[cur_div_stage]][[tier]]$tier <- tier
  }
}

marker_private_div_stage_list <- transpose(marker_private_div_stage_list)

marker_private_div_stage_list_filt <- list()
for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
  marker_private_div_stage_list_filt[[tier]] <- do.call(rbind, marker_private_div_stage_list[[tier]])
  temp_filt <- data.frame(marker_private_div_stage_list_filt[[tier]] %>%
                            filter(Bonferroni < 0.05) %>%
                            select(Category) %>%
                            unique())[,1]
  
  marker_private_div_stage_list_filt[[tier]] <- marker_private_div_stage_list_filt[[tier]] %>%
    filter(Category %in% temp_filt) %>%
    complete(Category, div_stage)
  
  marker_private_div_stage_list_filt[[tier]][is.na(marker_private_div_stage_list_filt[[tier]]$Fold),]$Fold <- 0
  marker_private_div_stage_list_filt[[tier]][is.na(marker_private_div_stage_list_filt[[tier]]$Bonferroni),]$Bonferroni <- 1
  
  marker_private_div_stage_list_filt[[tier]]$signif <- ""
  marker_private_div_stage_list_filt[[tier]]$signif <- ifelse(marker_private_div_stage_list_filt[[tier]]$Bonferroni < 0.05,
                                                                "*",
                                                              marker_private_div_stage_list_filt[[tier]]$signif)
  marker_private_div_stage_list_filt[[tier]]$signif <- ifelse(marker_private_div_stage_list_filt[[tier]]$Bonferroni < 0.005,
                                                                "**",
                                                              marker_private_div_stage_list_filt[[tier]]$signif)
  marker_private_div_stage_list_filt[[tier]]$signif <- ifelse(marker_private_div_stage_list_filt[[tier]]$Bonferroni < 0.0005,
                                                                "***",
                                                              marker_private_div_stage_list_filt[[tier]]$signif)
}

pdf(paste0(dir, "Plots/cell_figure_plots/marker_private_div_stage.pdf"), height = 4, width = 7)
do.call(rbind, marker_private_div_stage_list_filt) %>%
  ggplot(aes(x = Category, y = div_stage, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = c("15", "28", "50", "100", "200", "350", "600")) +
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
        legend.position = "none")
dev.off()

### joint non 1:1 and private
comb_private_markers <- rbind(cel_markers_non_1to1, cel_markers_private)
background_gene_set <- unique(unlist(comb_private_markers$gene))

comb_private_markers <- comb_private_markers[which(comb_private_markers$cell_type != "hyp3"),]

marker_comb_private_div_stage_list <- list()
for(cur_div_stage in unique(comb_private_markers$div_stage)) {
  print(cur_div_stage)
  
  rgs <- unique(comb_private_markers[comb_private_markers$div_stage %in% cur_div_stage,]$gene)
  marker_comb_private_div_stage_list[[cur_div_stage]] <- wormcat_enrichment(rgs, background_gene_set, gff_list_mrna[["elegans"]])
  
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    marker_comb_private_div_stage_list[[cur_div_stage]][[tier]]$div_stage <- cur_div_stage
    marker_comb_private_div_stage_list[[cur_div_stage]][[tier]]$tier <- tier
  }
}

marker_comb_private_div_stage_list <- transpose(marker_comb_private_div_stage_list)

marker_comb_private_div_stage_list_filt <- list()
for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
  marker_comb_private_div_stage_list_filt[[tier]] <- do.call(rbind, marker_comb_private_div_stage_list[[tier]])
  temp_filt <- data.frame(marker_comb_private_div_stage_list_filt[[tier]] %>%
                            filter(Bonferroni < 0.05) %>%
                            select(Category) %>%
                            unique())[,1]
  
  marker_comb_private_div_stage_list_filt[[tier]] <- marker_comb_private_div_stage_list_filt[[tier]] %>%
    filter(Category %in% temp_filt) %>%
    complete(Category, div_stage)
  
  marker_comb_private_div_stage_list_filt[[tier]][is.na(marker_comb_private_div_stage_list_filt[[tier]]$Fold),]$Fold <- 0
  marker_comb_private_div_stage_list_filt[[tier]][is.na(marker_comb_private_div_stage_list_filt[[tier]]$Bonferroni),]$Bonferroni <- 1
  
  marker_comb_private_div_stage_list_filt[[tier]]$signif <- ""
  marker_comb_private_div_stage_list_filt[[tier]]$signif <- ifelse(marker_comb_private_div_stage_list_filt[[tier]]$Bonferroni < 0.05,
                                                              "*",
                                                              marker_comb_private_div_stage_list_filt[[tier]]$signif)
  marker_comb_private_div_stage_list_filt[[tier]]$signif <- ifelse(marker_comb_private_div_stage_list_filt[[tier]]$Bonferroni < 0.005,
                                                              "**",
                                                              marker_comb_private_div_stage_list_filt[[tier]]$signif)
  marker_comb_private_div_stage_list_filt[[tier]]$signif <- ifelse(marker_comb_private_div_stage_list_filt[[tier]]$Bonferroni < 0.0005,
                                                              "***",
                                                              marker_comb_private_div_stage_list_filt[[tier]]$signif)
}

pdf(paste0(dir, "Plots/cell_figure_plots/marker_comb_private_div_stage.pdf"), height = 4, width = 7)
do.call(rbind, marker_comb_private_div_stage_list_filt) %>%
  ggplot(aes(x = Category, y = div_stage, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = c("15", "28", "50", "100", "200", "350", "600")) +
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
        legend.position = "none")
dev.off()


### omega
cel_markers <- left_join(cel_markers, gff_list_mrna[["elegans"]][,c("cds_gene_name", "omega", "Cutter_Ka")], by = c("gene" = "cds_gene_name"))
cbr_markers <- left_join(cbr_markers, gff_list_mrna[["briggsae"]][,c("cds_gene_name", "omega", "Cutter_Ka")],
                         by = c("gene" = "cds_gene_name"),
                         na_matches = "never")

joint_omega_markers <- rbind(cel_markers, cbr_markers)
joint_omega_markers <- joint_omega_markers[! duplicated(joint_omega_markers$cell_gene),]

cell_data$plot <- cell_data$cell_class
cell_data[which(cell_data$cell_class == "progenitor"),]$plot <- as.character(cell_data[which(cell_data$cell_class == "progenitor"),]$div_stage)

cell_class_marker_omega_plot <- left_join(cell_data, data.frame(joint_omega_markers %>%
                                  group_by(cell_type) %>%
                                  summarize(median_omega = median(omega, na.rm = TRUE),
                                            median_Cutter_Ka = median(Cutter_Ka, na.rm = TRUE)))) %>%
  filter(! is.na(cell_type)) %>%
  filter(! cell_type %in% "ABarppxa") %>%
  arrange(median_omega) %>%
  ggplot() + 
  geom_boxplot(aes(x = forcats::fct_reorder(cell_class, median_omega), y = median_omega), outlier.shape=NA) +
  geom_point(aes(y = median_omega, x = forcats::fct_reorder(cell_class, median_omega), group = cell_class, color = plot),
             position = position_dodge2(width = 1), alpha = 0.5, size = 1.5, stroke = 0.5) +
  # scale_x_discrete(limits = cell_class_plot_order, labels = cell_class_plot_labels) +
  scale_y_continuous(name = "Cell type marker dN/dS") +
  scale_color_manual(values = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                               '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF",
                               'Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
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

pdf(paste0(dir, "Plots/cell_figure_plots/marker_omega.pdf"), height = 4, width = 5)
print(cell_class_marker_omega_plot)
dev.off()

div_stage_marker_omega_plot <- left_join(cell_data, data.frame(joint_omega_markers %>%
                                  group_by(cell_type) %>%
                                  summarize(median_omega = median(omega, na.rm = TRUE),
                                            median_Cutter_Ka = median(Cutter_Ka, na.rm = TRUE)))) %>%
  filter(! is.na(cell_type)) %>%
  filter(! cell_type %in% "ABarppxa") %>%
  arrange(median_omega) %>%
  ggplot() + 
  geom_boxplot(aes(x = forcats::fct_reorder(div_stage, median_omega), y = median_omega), outlier.shape=NA) +
  geom_point(aes(y = median_omega, x = forcats::fct_reorder(div_stage, median_omega), group = cell_class,
                 color = ifelse(cell_class == "progenitor", unfactor(div_stage), cell_class)),
             position = position_dodge2(width = 1), alpha = 0.5, size = 1.5, stroke = 0.5) +
  scale_x_discrete(limits = c("15", "28", "50", "100", "200", "350", "600")) +
  scale_y_continuous(name = "Cell type marker dN/dS") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
                                '15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF",
                                "600" = "#08306BFF")) +
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

pdf(paste0(dir, "Plots/cell_figure_plots/marker_omega_div_stage.pdf"), height = 5, width = 5)
print(div_stage_marker_omega_plot)
dev.off()

pdf(paste0(dir, "Plots/cell_figure_plots/marker_omega_cell_class_div_stage.pdf"), height = 5, width = 12.5)
plot_grid(cell_class_marker_omega_plot, div_stage_marker_omega_plot,
          nrow = 1, align = "hv")
dev.off()


## marker KA

cell_class_marker_ka_plot <- left_join(cell_data, data.frame(joint_omega_markers %>%
                                                                  group_by(cell_type) %>%
                                                                  summarize(median_omega = median(omega, na.rm = TRUE),
                                                                            median_Cutter_Ka = median(Cutter_Ka, na.rm = TRUE)))) %>%
  filter(! is.na(cell_type)) %>%
  filter(! cell_type %in% "ABarppxa") %>%
  arrange(median_Cutter_Ka) %>%
  ggplot() + 
  geom_boxplot(aes(x = forcats::fct_reorder(cell_class, median_Cutter_Ka), y = median_Cutter_Ka), outlier.shape=NA) +
  geom_point(aes(y = median_Cutter_Ka, x = forcats::fct_reorder(cell_class, median_Cutter_Ka), group = cell_class, color = plot),
             position = position_dodge2(width = 1), alpha = 0.5, size = 1.5, stroke = 0.5) +
  scale_y_continuous(name = "Cell type marker Ka") +
  scale_color_manual(values = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                               '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF",
                               'Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
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

pdf(paste0(dir, "Plots/cell_figure_plots/marker_ka.pdf"), height = 5, width = 5)
print(cell_class_marker_ka_plot)
dev.off()

div_stage_marker_ka_plot <- left_join(cell_data, data.frame(joint_omega_markers %>%
                                                                 group_by(cell_type) %>%
                                                                 summarize(median_omega = median(omega, na.rm = TRUE),
                                                                           median_Cutter_Ka = median(Cutter_Ka, na.rm = TRUE)))) %>%
  filter(! is.na(cell_type)) %>%
  filter(! cell_type %in% "ABarppxa") %>%
  arrange(median_Cutter_Ka) %>%
  ggplot() + 
  geom_boxplot(aes(x = forcats::fct_reorder(div_stage, median_Cutter_Ka), y = median_Cutter_Ka), outlier.shape=NA) +
  geom_point(aes(y = median_Cutter_Ka, x = forcats::fct_reorder(div_stage, median_Cutter_Ka), group = cell_class,
                 color = ifelse(cell_class == "progenitor", unfactor(div_stage), cell_class)),
             position = position_dodge2(width = 1), alpha = 0.5, size = 1.5, stroke = 0.5) +
  scale_x_discrete(limits = c("15", "28", "50", "100", "200", "350", "600")) +
  scale_y_continuous(name = "Cell type marker Ka") +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
                                '15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF",
                                "600" = "#08306BFF")) +
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

pdf(paste0(dir, "Plots/cell_figure_plots/marker_ka_div_stage.pdf"), height = 5, width = 5)
print(div_stage_marker_ka_plot)
dev.off()

pdf(paste0(dir, "Plots/cell_figure_plots/marker_ka.pdf"), height = 5, width = 12.5)
plot_grid(cell_class_marker_ka_plot, div_stage_marker_ka_plot,
          nrow = 1, align = "hv")
dev.off()

pdf(paste0(dir, "Plots/cell_figure_plots/marker_omega_ka.pdf"), height = 10, width = 12.5)
plot_grid(cell_class_marker_omega_plot, div_stage_marker_omega_plot,
          cell_class_marker_ka_plot, div_stage_marker_ka_plot,
          nrow = 2, align = "hv")
dev.off()

##########################
# Marker type enrichment #
##########################

cel_markers$cell_class_div_stage <- cel_markers$cell_class
cel_markers[which(cel_markers$cell_class == "progenitor"),]$cell_class_div_stage <- as.character(cel_markers[which(cel_markers$cell_class == "progenitor"),]$div_stage)

cbr_markers$cell_class_div_stage <- cbr_markers$cell_class
cbr_markers[which(cbr_markers$cell_class == "progenitor"),]$cell_class_div_stage <- as.character(cbr_markers[which(cbr_markers$cell_class == "progenitor"),]$div_stage)

rownames(cel_markers) <- cel_markers$cell_gene
rownames(cbr_markers) <- cbr_markers$cell_gene

shared_markers <- intersect(rownames(cel_markers), rownames(cbr_markers))        

marker_cats_list <- list()
for(cell_class in sort(unique(cel_markers$cell_class_div_stage))) {
  marker_cats_list[[cell_class]] <- list()
  
  if(cell_class == "600") {
    temp_markers <- cel_markers[which(cel_markers$div_stage == cell_class),]
  } else {
    temp_markers <- cel_markers[which(cel_markers$cell_class_div_stage == cell_class),]
  }
  
  marker_cats_list[[cell_class]][["1:1"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1,]$gene)
  marker_cats_list[[cell_class]][["1:1 shared"]] <- unique(temp_markers[temp_markers$cell_gene %in% shared_markers & temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1,]$gene)
  marker_cats_list[[cell_class]][["1:1 divergent"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1 &
                                                                      temp_markers$cbr_tpm < 80 &
                                                                      temp_markers$cbr_tpm_log2fc < 0 &
                                                                      temp_markers$cel_tpm/(temp_markers$cbr_tpm + 1) > 5,]$gene)
  marker_cats_list[[cell_class]][["1:2"]] <- unique(temp_markers[(temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 2) | (temp_markers$cel_OG_count == 2 & temp_markers$cbr_OG_count == 1),]$gene)
  marker_cats_list[[cell_class]][["1:many"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count > 2,]$gene)
  marker_cats_list[[cell_class]][["many:many"]] <- unique(temp_markers[temp_markers$cel_OG_count >= 2 & temp_markers$cbr_OG_count >= 2,]$gene)
  marker_cats_list[[cell_class]][["0:n"]] <- unique(temp_markers[temp_markers$cel_OG_count > 0 & temp_markers$cbr_OG_count == 0,]$gene)
}

# Do one just for progenitors
marker_cats_list[["Progenitors"]] <- list()
temp_markers <- cel_markers[which(cel_markers$cell_class == "progenitor"),]

marker_cats_list[["Progenitors"]][["1:1"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1,]$gene)
marker_cats_list[["Progenitors"]][["1:1 shared"]] <- unique(temp_markers[temp_markers$cell_gene %in% shared_markers & temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1,]$gene)
marker_cats_list[["Progenitors"]][["1:1 divergent"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1 &
                                                                           temp_markers$cbr_tpm < 80 &
                                                                           temp_markers$cbr_tpm_log2fc < 0 &
                                                                           temp_markers$cel_tpm/(temp_markers$cbr_tpm + 1) > 5,]$gene)
marker_cats_list[["Progenitors"]][["1:2"]] <- unique(temp_markers[(temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 2) | (temp_markers$cel_OG_count == 2 & temp_markers$cbr_OG_count == 1),]$gene)
marker_cats_list[["Progenitors"]][["1:many"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count > 2,]$gene)
marker_cats_list[["Progenitors"]][["many:many"]] <- unique(temp_markers[temp_markers$cel_OG_count >= 2 & temp_markers$cbr_OG_count >= 2,]$gene)
marker_cats_list[["Progenitors"]][["0:n"]] <- unique(temp_markers[temp_markers$cel_OG_count > 0 & temp_markers$cbr_OG_count == 0,]$gene)

# Do one just for amphid neurons
amphid_ciliated <- c("ADL", "AFD", "ADF", "ASI", "ASJ", "ASK", "AWA", "AWB", "AWC", "AWD", "ASE", "ASG", "ASH")

marker_cats_list[["Amphid ciliated neurons"]] <- list()
temp_markers <- cel_markers[which(cel_markers$cell_type %in% amphid_ciliated),]

marker_cats_list[["Amphid ciliated neurons"]][["1:1"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1,]$gene)
marker_cats_list[["Amphid ciliated neurons"]][["1:1 shared"]] <- unique(temp_markers[temp_markers$cell_gene %in% shared_markers & temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1,]$gene)
marker_cats_list[["Amphid ciliated neurons"]][["1:1 divergent"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1 &
                                                                              temp_markers$cbr_tpm < 80 &
                                                                              temp_markers$cbr_tpm_log2fc < 0 &
                                                                              temp_markers$cel_tpm/(temp_markers$cbr_tpm + 1) > 5,]$gene)
marker_cats_list[["Amphid ciliated neurons"]][["1:2"]] <- unique(temp_markers[(temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 2) | (temp_markers$cel_OG_count == 2 & temp_markers$cbr_OG_count == 1),]$gene)
marker_cats_list[["Amphid ciliated neurons"]][["1:many"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count > 2,]$gene)
marker_cats_list[["Amphid ciliated neurons"]][["many:many"]] <- unique(temp_markers[temp_markers$cel_OG_count >= 2 & temp_markers$cbr_OG_count >= 2,]$gene)
marker_cats_list[["Amphid ciliated neurons"]][["0:n"]] <- unique(temp_markers[temp_markers$cel_OG_count > 0 & temp_markers$cbr_OG_count == 0,]$gene)

marker_cats_list[["Non-amphid ciliated neurons"]] <- list()
temp_markers <- cel_markers[which(cel_markers$cell_class %in% "Ciliated neurons" & ! cel_markers$cell_type %in% amphid_ciliated),]

marker_cats_list[["Non-amphid ciliated neurons"]][["1:1"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1,]$gene)
marker_cats_list[["Non-amphid ciliated neurons"]][["1:1 shared"]] <- unique(temp_markers[temp_markers$cell_gene %in% shared_markers & temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1,]$gene)
marker_cats_list[["Non-amphid ciliated neurons"]][["1:1 divergent"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 1 &
                                                                                          temp_markers$cbr_tpm < 80 &
                                                                                          temp_markers$cbr_tpm_log2fc < 0 &
                                                                                          temp_markers$cel_tpm/(temp_markers$cbr_tpm + 1) > 5,]$gene)
marker_cats_list[["Non-amphid ciliated neurons"]][["1:2"]] <- unique(temp_markers[(temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count == 2) | (temp_markers$cel_OG_count == 2 & temp_markers$cbr_OG_count == 1),]$gene)
marker_cats_list[["Non-amphid ciliated neurons"]][["1:many"]] <- unique(temp_markers[temp_markers$cel_OG_count == 1 & temp_markers$cbr_OG_count > 2,]$gene)
marker_cats_list[["Non-amphid ciliated neurons"]][["many:many"]] <- unique(temp_markers[temp_markers$cel_OG_count >= 2 & temp_markers$cbr_OG_count >= 2,]$gene)
marker_cats_list[["Non-amphid ciliated neurons"]][["0:n"]] <- unique(temp_markers[temp_markers$cel_OG_count > 0 & temp_markers$cbr_OG_count == 0,]$gene)


marker_cats_list <- transpose(marker_cats_list)
background_gene_set <- unique(unlist(marker_cats_list))

marker_types_list <- list()
for(cat in names(marker_cats_list)) {
  print(cat)
  marker_types_list[[cat]] <- list()
  
  for(cell_class in names(marker_cats_list[[cat]])) {
    print(cell_class)
    
    rgs <- marker_cats_list[[cat]][[cell_class]]
    
    if(length(rgs) < 1) {
      for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
        temp <- data.frame(matrix(nrow = 0, ncol = 9))
        colnames(temp) <- c("Category", "RGS", "AC", "Fold", "PValue", "Bonferroni", "type", "cell_class", "tier")
        marker_types_list[[cat]][[cell_class]][[tier]] <- temp
      }
    } else {
      marker_types_list[[cat]][[cell_class]] <- wormcat_enrichment(rgs, background_gene_set, gff_list_mrna[["elegans"]])
      
      for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
        marker_types_list[[cat]][[cell_class]][[tier]]$type <- cat
        marker_types_list[[cat]][[cell_class]][[tier]]$cell_class <- cell_class
        marker_types_list[[cat]][[cell_class]][[tier]]$tier <- tier
      }
    }
  }
}

marker_types_df <- do.call(rbind, do.call(rbind, do.call(rbind, marker_types_list))) %>%
  complete(Category, type, cell_class)
marker_types_df[is.na(marker_types_df$Fold),]$Bonferroni <- 1
marker_types_df[is.na(marker_types_df$Fold),]$Fold <- 0

marker_types_df$signif <- ""
marker_types_df$signif <- ifelse(marker_types_df$Bonferroni < 0.05, "*", marker_types_df$signif)
marker_types_df$signif <- ifelse(marker_types_df$Bonferroni < 0.005, "**", marker_types_df$signif)
marker_types_df$signif <- ifelse(marker_types_df$Bonferroni < 0.0005, "***", marker_types_df$signif)

col_order <- c("Ciliated neurons", "Amphid ciliated neurons", "Non-amphid ciliated neurons", "Non-ciliated neurons",
               "Glia and excretory", "Hypodermis and seam",
               "Pharynx and rectal", "Muscle",
               "Mesoderm", "Intestine",
               "Germline", "600",
               "350", "200",
               "100", "50",
               "28", "15", "Progenitors")

pdf(paste0(dir, "Plots/marker_all_class_div_stage_cat.pdf"), height = 22.5, width = 80)
marker_types_df %>%
  mutate(type = factor(type, levels = c("1:1", "1:1 shared", "1:1 divergent", "1:2", "1:many", "many:many", "0:n"))) %>%
  ggplot(aes(x = Category, y = cell_class, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  facet_wrap(~dataset, scale="free_y", nrow=4)+
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = col_order) +
  scale_fill_gradient(name = "Fold\nEnrich.", low = "white", high="blue",
                      breaks = c(0, 1, 2, 3, 4, 5),
                      labels = c("0", "1", "2", "3", "4", ">5"),
                      limits = c(0, 5),
                      oob = scales::squish) +
  # scale_fill_viridis() +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid = element_line(color="grey80", size=0.25),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.position = "none") +
  facet_wrap(~type, nrow = 7)
dev.off()

tier_one <- unique(marker_types_df[marker_types_df$tier %in% "WormCat.1",]$Category)

pdf(paste0(dir, "Plots/marker_tier_1_class_div_stage_cat.pdf"), height = 20, width = 6.5)
marker_types_df %>%
  mutate(type = factor(type, levels = c("1:1", "1:1 shared", "1:1 divergent", "1:2", "1:many", "many:many", "0:n"))) %>%
  filter(Category %in% tier_one) %>%
  ggplot(aes(x = Category, y = cell_class, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  facet_wrap(~dataset, scale="free_y", nrow=4)+
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = col_order) +
  scale_fill_gradient(name = "Fold\nEnrich.", low = "white", high="blue",
                      breaks = c(0, 1, 2, 3, 4, 5),
                      labels = c("0", "1", "2", "3", "4", ">5"),
                      limits = c(0, 5),
                      oob = scales::squish) +
  # scale_fill_viridis() +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid = element_line(color="grey80", size=0.25),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.position = "none") +
  facet_wrap(~type, nrow = 7)
dev.off()

cats <- c("Development", "Cell cycle", "DNA: replication", "Proteolysis proteasome: E3: F box",
          "Transcription: chromatin structure: histone",
          "Transcription: general machinery Transcription: general machinery: RNA Pol II")

col_order <- c("600",
               "350", "200",
               "100", "50",
               "28", "15")

pdf(paste0(dir, "Plots/marker_subset_class_div_stage_cat.pdf"), height = 3, width = 6)
marker_types_df %>%
  mutate(type = factor(type, levels = c("1:1", "1:1 shared", "1:1 divergent", "1:2", "many:many", "0:n"))) %>%
  filter(type %in% c("1:1", "1:2", "many:many", "0:n")) %>%
  filter(Category %in% cats) %>%
  filter(cell_class %in% col_order) %>%
  ggplot(aes(x = Category, y = cell_class, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  facet_wrap(~dataset, scale="free_y", nrow=4)+
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45),
                   breaks = c("Cell cycle", "Development", "DNA: replication",
                              "Proteolysis proteasome: E3: F box",
                              "Transcription: chromatin structure: histone"),
                   labels = c("Cell cycle", "Development", "DNA replication",
                              "F box",
                              "Histone")) +
  scale_y_discrete(limits = rev(col_order)) +
  scale_fill_gradient(name = "Fold\nEnrich.", low = "white", high="blue",
                      breaks = c(0, 1, 2, 3, 4, 5),
                      labels = c("0", "1", "2", "3", "4", ">5"),
                      limits = c(0, 5),
                      oob = scales::squish) +
  # scale_fill_viridis() +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid = element_line(color="grey80", size=0.25),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.position = "none") +
  facet_wrap(~type, ncol = 4)
dev.off()

col_order <- c("Amphid ciliated neurons", "Non-amphid ciliated neurons", "Non-ciliated neurons",
               "Glia and excretory", "Hypodermis and seam",
               "Pharynx and rectal", "Muscle",
               "Mesoderm", "Intestine",
               "Germline", "Progenitors")

cats <- c("Development", "mRNA functions: processing",
          "Transcription: general machinery: RNA Pol II",
          "Transcription: chromatin structure: histone",
          "Transcription factor: NHR", "Stress response",
          "Transmembrane protein: seven transmembrane receptor",
          "Proteolysis proteasome: E3: F box", "Unassigned")

pdf(paste0(dir, "Plots/marker_subset_progenitors.pdf"), height = 3, width = 6.5)
marker_types_df %>%
  mutate(type = factor(type, levels = c("1:1", "1:1 shared", "1:1 divergent", "1:2", "many:many", "0:n"))) %>%
  filter(type %in% c("1:1", "many:many", "0:n")) %>%
  filter(Category %in% cats) %>%
  filter(cell_class %in% col_order) %>%
  mutate(Category = factor(Category, levels = cats)) %>%
  ggplot(aes(x = Category, y = cell_class, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  facet_wrap(~dataset, scale="free_y", nrow=4)+
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45),
                   breaks = c("Transcription: general machinery: RNA Pol II", "Development", "mRNA functions: processing",
                              "Proteolysis proteasome: E3: F box",
                              "Transcription: chromatin structure: histone",
                              "Transcription factor: NHR", "Stress response",
                              "Transmembrane protein: seven transmembrane receptor", "Unassigned"),
                   labels = c("RNA Pol II", "Development", "mRNA processing",
                              "F box",
                              "Histone", "NHR",
                              "Stress response", "GPCR", "Unassigned")) +
  scale_y_discrete(limits = rev(col_order)) +
  scale_fill_gradient(name = "Fold\nEnrich.", low = "white", high="blue",
                      breaks = c(0, 1, 2, 3, 4, 5),
                      labels = c("0", "1", "2", "3", "4", ">5"),
                      limits = c(0, 5),
                      oob = scales::squish) +
  # scale_fill_viridis() +
  theme(axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid = element_line(color="grey80", size=0.25),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.position = "none") +
  facet_wrap(~type, ncol = 4)
dev.off()

marker_bar_plot <- data.frame(do.call(rbind, transpose(lapply(marker_cats_list, function(x) lapply(x, length)))))
cell_class_rownames <- rownames(marker_bar_plot)
marker_bar_plot <- data.frame(apply(marker_bar_plot, 2, function(x) as.numeric(unlist(x))))

colnames(marker_bar_plot) <- c("one_one", "one_one_shared", "one_one_divergent", "one_two", "one_many", "many_many", "zero_n")

marker_bar_plot$cell_class <- cell_class_rownames

marker_bar_plot_melt <- melt(marker_bar_plot, id.vars = "cell_class")
marker_bar_plot_melt <- marker_bar_plot_melt[! marker_bar_plot_melt$variable %in% c("one_one_shared", "one_one_divergent"),]

#, limits = c(0, 40), oob = scales::squish
#"one_one",  "1:1", 
pdf(paste0(dir, "Plots/cell_figure_plots/percent_marker_types.pdf"), height = 3, width = 6)
marker_bar_plot_melt %>% group_by(cell_class) %>%
  mutate(value = value / sum(value) * 100) %>%
  filter(! cell_class %in% c("Ciliated neurons", "Progenitors")) %>%
  filter(! variable %in% "one_one") %>%
  ggplot(aes(x = cell_class, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_discrete(limits =  c(rev(c("600",
                                     "350", "200",
                                     "100", "50",
                                     "28", "15")), "Germline", "Amphid ciliated neurons", "Non-amphid ciliated neurons", "Non-ciliated neurons",
                                           "Glia and excretory", "Hypodermis and seam",
                                           "Pharynx and rectal", "Muscle",
                                           "Mesoderm", "Intestine"),
                   guide = guide_axis(angle = 45)) +
  scale_y_continuous(name = "Percent of marker types") +
  scale_fill_discrete(name = "Marker type",
                      breaks = c("one_two", "one_many", "many_many", "zero_n"),
                      labels = c("1:2", "1:many", "many:many", "n:0")) + 
  theme(axis.title.x = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid = element_blank(),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        legend.position = "right")
dev.off()

###########################
# All by all markers plot #
###########################

cols = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
  'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
  'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
  '15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
  '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF",
  "600" = "#08306BFF")[cell_data$plot]

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

pdf(paste0(dir, "Plots/cell_figure_plots/all_by_all_markers.pdf"), height = 10, width = 10)
left_join(cell_data, data.frame(joint_omega_markers %>%
                                  group_by(cell_type) %>%
                                  summarize(median_omega = median(omega, na.rm = TRUE),
                                            median_Cutter_Ka = median(Cutter_Ka, na.rm = TRUE)))) %>%
  filter(! is.na(cell_type)) %>%
  filter(! cell_type %in% "ABarppxa") %>%
  filter(cell_class != "progenitor") %>%
  arrange(median_omega) %>%
  select(c("jsd_median", "fraction_almost_shared", "fraction_divergent", "median_omega")) %>%
pairs(~ jsd_median + fraction_almost_shared + fraction_divergent + median_omega, data = .,
      pch = 21, bg = cols, diag.panel = panel.hist)
dev.off()

pdf(paste0(dir, "Plots/cell_figure_plots/all_by_all_markers.pdf"), height = 10, width = 10)
left_join(cell_data, data.frame(joint_omega_markers %>%
                                  group_by(cell_type) %>%
                                  summarize(median_omega = median(omega, na.rm = TRUE),
                                            median_Cutter_Ka = median(Cutter_Ka, na.rm = TRUE)))) %>%
  filter(! is.na(cell_type)) %>%
  filter(! cell_type %in% "ABarppxa") %>%
  filter(cell_class != "progenitor") %>%
  arrange(median_omega) %>%
  select(c("jsd_median", "fraction_almost_shared", "fraction_divergent", "median_omega", "plot")) %>%
ggpairs(., columns = 1:5, aes(color = plot, fill = plot, alpha = 0.5),
        upper = list(continuous = "points"), lower = list(combo = "count")) +
  scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                         'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                         'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
                         '15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                         '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF",
                         "600" = "#08306BFF")) +
  scale_fill_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
                                '15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF",
                                "600" = "#08306BFF"))
dev.off()