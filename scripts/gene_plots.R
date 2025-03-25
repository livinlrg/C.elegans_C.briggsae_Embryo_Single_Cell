library(monocle3)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggtext)
library(Gviz) # for 
library(ggplotify)
library(parallel)
library(dplyr)
library(ggraph)
library(igraph)
library(Seurat)

# Load objects
dir <- "/kimdata/livinlrg/scAnalysis/WS290/"
source(paste0(dir, "../Scripts/WS290/CalcFunctions.R"))

gene_data <- readRDS(paste0(dir, "Objects/gene_data_joint_cell_bg.rds"))
cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))

gene_data$cel_max_tpm <- apply(gene_data[,c("cel_max_tpm_pro", "cel_max_tpm_term")], 1, max)
gene_data$cbr_max_tpm <- apply(gene_data[,c("cbr_max_tpm_pro", "cbr_max_tpm_term")], 1, max)

gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))

TPMListBootstrapMean_term <- readRDS(paste0(dir, "Objects/TPMListBootstrapMean_CellCorrection.rds"))
TPMListBootstrap_pro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))
BinarizeBootstrapListMean_term <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListMean_CellCorrection.rds"))
BinarizeBootstrapListPro <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListPro_CellCorrection.rds"))

TPMJoint <- list()
TPMJoint[["C.elegans"]] <- cbind(TPMListBootstrap_pro[["C.elegans"]][,sort(colnames(TPMListBootstrap_pro[["C.elegans"]]))], TPMListBootstrapMean_term[["C.elegans"]])
TPMJoint[["C.briggsae"]] <- cbind(TPMListBootstrap_pro[["C.briggsae"]][,sort(colnames(TPMListBootstrap_pro[["C.elegans"]]))], TPMListBootstrapMean_term[["C.briggsae"]])

BinarizeList <- list()
BinarizeList[["C.elegans"]] <- cbind(BinarizeBootstrapListPro[["C.elegans"]][,sort(colnames(BinarizeBootstrapListPro[["C.elegans"]]))], BinarizeBootstrapListMean_term[["C.elegans"]])
BinarizeList[["C.briggsae"]] <- cbind(BinarizeBootstrapListPro[["C.briggsae"]][,sort(colnames(BinarizeBootstrapListPro[["C.briggsae"]]))], BinarizeBootstrapListMean_term[["C.briggsae"]])

synteny_filt <- readRDS(file = paste0(dir, "Objects/synteny_filt.rds"))

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240903.rds"))

tpm_matrix <- readRDS(paste0(dir, "Objects/tpm_matrix_list_cell_bg.rds"))

joint_markers <- readRDS(paste0(dir, "Objects/joint_markers_cell_bg.rds"))


######################################
# Generate list of tissue embeddings #
######################################

so_cell_class_list <- readRDS(paste0(dir, "Objects/seurat_objects/so_cell_class_list_20240827.rds"))
so_time_bin_list <- readRDS(paste0(dir, "Objects/seurat_objects/so_time_bin_list_20241010.rds"))
so_rpca <- readRDS(paste0(dir, "Objects/seurat_objects/seurat_rpca_integration_joined_20240823.rds"))

tissue_embeddings <- list()
for(reduction in c(names(so_cell_class_list), names(so_time_bin_list))) {
  if(reduction %in% names(so_cell_class_list)) {
    tissue_embeddings[[reduction]] <- Embeddings(so_cell_class_list[[reduction]], reduction = "umap.rpca")
  } else {
    tissue_embeddings[[reduction]] <- Embeddings(so_time_bin_list[[reduction]], reduction = "umap.rpca")
  }
}

tissue_embeddings[["global"]] <- Embeddings(so_rpca, reduction = "umap.rpca")

rm(so_cell_class_list)
rm(so_time_bin_list)
rm(so_rpca)

tissue_embeddings <- lapply(tissue_embeddings, function(embeddings) {
  return(data.frame(X = embeddings[,1],
                    Y = embeddings[,2],
                    species = colData(cds)[rownames(embeddings),]$species))
})

gene_counts <- readRDS(paste0(dir, "Objects/cello_creation/matrix.rds"))
gene_counts <- gene_counts[,colnames(gene_counts) %in% rownames(tissue_embeddings[["global"]])]

############################
# Plot expression on UMAPs #
############################

pdf(paste0(dir, "Plots/test.pdf"), width = 5, height = 5)
plotExpUmap("daf-19", tissue_embeddings[["global"]][tissue_embeddings[["global"]]$species == "C.elegans",], gene_counts["daf-19",tissue_embeddings[["global"]]$species == "C.elegans"], "C.elegans", "global")
dev.off()

plotExpUmap <- function(gene, tissue_embeddings, gene_counts, cur_species, tissue = "", marker_size = 0.025) {
  print("Plotting UMAP")
  
  ## Plot UMAP
  projection <- data.frame(X = tissue_embeddings[,1],
                           Y = tissue_embeddings[,2],
                           Expression = gene_counts)

  glim <- ceiling(quantile(projection$Expression[projection$Expression != 0], .995))
  if(is.na(glim)) glim = 6
  if(! is.na(glim) && glim < 6) glim = 6
  
  projection <- projection[order(projection$Expression),]
  
  idx_region <- which(projection$Expression > 0)
  na.col = "lightgrey"
  
  p1 <- ggplot(projection, aes(x = X,
                               y = Y)) +
    geom_point(size = marker_size, color = na.col, show.legend = FALSE, stroke=0) +
    geom_point(data = projection[idx_region,], aes(x = X, y = Y, color = Expression + 1),
               size = marker_size, stroke = marker_size / 3,
               alpha = 1) +
    scale_x_continuous(name = "UMAP 1") +
    scale_y_continuous(name = "UMAP 2") +
    scale_color_continuous(type = "viridis", 
                           limits = c(1, glim),
                           oob = scales::squish) +
    labs(title = paste0(ifelse(cur_species == "C.elegans", "C. elegans", "C. briggsae"), " ", tissue)) +
    theme(plot.title = element_text(face = "italic", size = 16),
          rect = element_rect(fill = "transparent"),
          axis.line = element_line(color="grey80", size=1),
          #axis.line.x.bottom = element_line(color="grey80", size=1),
          #axis.line.y.left = element_line(color="grey80", size=1),
          axis.line.y = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
                                                         ends = "last")),
          axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), 
                                                         ends = "last")),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.ticks = element_blank(),
          #text = element_text(family = "sans", color = "black", size=12),
          axis.text = element_blank(),
          legend.title = element_blank(),
          legend.box = "horizontal",
          legend.position = "right",
          #legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          #legend.key = element_rect(colour = "transparent", fill = "transparent")
    )
  return(p1)
}

#########################
# TPM dot plot function #
#########################

pdf(paste0(dir, "Plots/test.pdf"), width = 5, height = 5)
plot_gene_example("daf-19", TPMJoint, BinarizeList, cell_data, "progenitor")
dev.off()

for(cur_gene in c("bath-42", "nphp-4", "chca-1", "rpac-19", "imb-2", "try-1", "acox-1.6")) {
  pdf(paste0(dir, "Plots/gene_figure_plots/gene_plots/", cur_gene,".pdf"), width = 3, height = 3)
  print(plot_gene_example(cur_gene, TPMJoint, BinarizeList, cell_data, "terminal"))
  dev.off()
}

plot_gene_example <- function(gene, TPMList, BinarizeList, cell_data, type) {
  print("Plotting dot plot")
  
  exp_filter <- (BinarizeList[["C.elegans"]][gene,] > 0 | BinarizeList[["C.briggsae"]][gene,] > 0) & (TPMList[["C.elegans"]][gene,] > 40 | TPMList[["C.briggsae"]][gene,] > 40)
  exp_filter <- colnames(exp_filter)[exp_filter]
  
  # Add terminal cell class
  if(type == "terminal") {
    temp <- lapply(TPMList, function(x) {
      df <- data.frame(expression = t(x[gene,] + 1)[,1])
      df$cell_type <- rownames(df)
      df <- left_join(df, cell_data[,c("cell_type", "cell_class")], by = "cell_type") %>% filter(cell_class != "progenitor")
      return(df[df$cell_type %in% exp_filter,])
    })
    
  # Add lineage group
  } else if(type == "progenitor") {
    temp <- lapply(TPMList, function(x) {
      df <- data.frame(expression = t(x[gene,] + 1)[,1])
      df$cell_type <- rownames(df)
      df <- left_join(df, cell_data[,c("cell_type", "cell_class", "div_stage", "lineage_group")], by = "cell_type") %>% filter(cell_class == "progenitor")
      df$cell_class <- df$div_stage
      return(df[df$cell_type %in% exp_filter,])
    })
  }
  
  return(ggplot() +
           geom_point(aes(x = temp[["C.elegans"]]$expression,
                          y = temp[["C.briggsae"]]$expression,
                          color = temp[["C.elegans"]]$cell_class),
                      shape = 21) +
           geom_text_repel(aes(x = temp[["C.elegans"]]$expression,
                               y = temp[["C.briggsae"]]$expression,
                               label = temp[["C.elegans"]]$cell_type,
                               color = temp[["C.elegans"]]$cell_class),
                           show.legend = FALSE ) +
           scale_x_continuous(name = paste0("C. elegans ", gene," TPM"), trans = "log2",
                              limits = c(min(temp[["C.elegans"]]$expression,
                                             temp[["C.briggsae"]]$expression) * 0.9,
                                         max(temp[["C.elegans"]]$expression,
                                             temp[["C.briggsae"]]$expression) * 1.1)) +
           scale_y_continuous(name = paste0("C.briggsae ", gene, " TPM"), trans = "log2",
                              limits = c(min(temp[["C.elegans"]]$expression,
                                             temp[["C.briggsae"]]$expression) * 0.9,
                                         max(temp[["C.elegans"]]$expression,
                                             temp[["C.briggsae"]]$expression) * 1.1)) +
           guides(color = guide_legend(ncol = 3), label = "none") +
           scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                         'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                         'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1",
                                         '15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                                         '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF",
                                         "600" = "#08306BFF")) +
           theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
                 plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                 axis.line = element_line(color="grey80", size=1),
                 axis.title.x = element_text(color = "#009E73", size = 10),
                 axis.title.y = element_text(color = "#56B4E9", size = 10),
                 panel.grid.major = element_line(color="grey80", size=0.25),
                 panel.grid.minor = element_line(color="grey80", size=0.05),
                 legend.title = element_blank(),
                 legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
                 legend.key = element_rect(colour = "transparent", fill = "transparent"),
                 legend.position = "none"))
}

##############################
# Orthology confidence score #
##############################

# ortho scheme
# mutual best blast hit = 5
# synteny = 5
# small OG number = 3
# large OG number = -1

generate_gene_confidence <- function(gene, gene_data, synteny_filt) {
  ## Confidence in orthology scores
  ortho_confidence <- 0
  
  cel_gene_name = gene_data[gene,]$elegans_id
  cbr_gene_name = gene_data[gene,]$briggsae_id
  
  if(gene_data[gene,]$cel_OG_size == 1 & gene_data[gene,]$cel_OG_size == 1) {
    ortho_confidence = ortho_confidence + 3
  } else if (gene_data[gene,]$cel_OG_size > 3 | gene_data[gene,]$cel_OG_size > 3) {
    ortho_confidence = ortho_confidence - 1
  }
  
  if(cel_gene_name %in% synteny_filt$cel_WBG_name & cbr_gene_name %in% synteny_filt$cbr_WBG_name) {
    ortho_confidence = ortho_confidence + 5
    synteny_temp <- synteny_filt[synteny_filt$cel_WBG_name == cel_gene_name & synteny_filt$cbr_WBG_name == cbr_gene_name,]
    
    if(synteny_temp["syntenic"] == "IS_SYNTENIC") {
      ortho_confidence = ortho_confidence + 5
    } else if (synteny_temp["syntenic"] == "NOT_SYNTENIC") {
      ortho_confidence = ortho_confidence - 3
    }
  }
  
  if(is.na(BobOrtho[gene,]$BlastBestElegansGene)) {
    ortho_confidence = ortho_confidence
  } else if(BobOrtho[gene,]$BlastPValueBest == 0 & BobOrtho[gene,]$BlastBestElegansGene == gene) {
    ortho_confidence = ortho_confidence + 5
  } else if(BobOrtho[gene,]$BlastBestElegansGene == gene) {
    ortho_confidence = ortho_confidence + 1
    if(is.na(BobOrtho[gene,]$BlastPValueNextBest)) {
      ortho_confidence = ortho_confidence + 1
    }
    else if((log(BobOrtho[gene,]$BlastPValueNextBest) - log(BobOrtho[gene,]$BlastPValueBest)) > 100 ) {
      ortho_confidence = ortho_confidence + 3
    } else if((log(BobOrtho[gene,]$BlastPValueNextBest) - log(BobOrtho[gene,]$BlastPValueBest)) > 10) {
      ortho_confidence = ortho_confidence + 1
    }
  }
  
  if(is.na(BobOrtho[gene,]$BLATBestElegansGene)) {
    ortho_confidence = ortho_confidence
  } else if(BobOrtho[gene,]$BLATBestElegansGene == gene) {
    ortho_confidence = ortho_confidence + 1
    if(is.na(BobOrtho[gene,]$BLATNextScore)) {
      ortho_confidence = ortho_confidence + 1
    } else if((log2(BobOrtho[gene,]$BLATScore) - log2(BobOrtho[gene,]$BLATNextScore)) > 1.5) {
      ortho_confidence = ortho_confidence + 3
    } else if((log2(BobOrtho[gene,]$BLATScore) - log2(BobOrtho[gene,]$BLATNextScore)) > 1) {
      ortho_confidence = ortho_confidence + 1
    }
  }
  
  if(ortho_confidence < 2) {
    return("low")
  } else if(ortho_confidence < 5) {
    return("medium-low")
  } else if(ortho_confidence < 7) {
    return("medium")
  } else if(ortho_confidence <= 8) {
    return("medium-high")
  } else if(ortho_confidence > 8) {
    return("high")
  }
}

generate_lagene_confidence <- function(gene, gene_data) {
  if(gene_data[gene,"orthology_conf"] == "one:one") {
    return("high")
  } else if (gene_data[gene,"orthology_conf"] == "confident_canonical") {
    return("medium")
  } else if (gene_data[gene,"orthology_conf"] == "canonical") {
    return("low")
  }
}

# Metadata for object
getMode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# List the files in the dir
files <- list.files(paste0(dir, "Plots/gene_plots_fix/"), pattern = ".png")
genes_run <- unlist(lapply(strsplit(files, ".png"), function(x) x[[1]]))
genes_not_run <- gene_data$gene[! gene_data$gene %in% genes_run]

# for(cur_gene in gene_data$gene) {
# mclapply(gene_data$gene, function(cur_gene) {
mclapply(genes_not_run, function(cur_gene) {
    
  print(cur_gene)
  cur_tissue = getMode(c(cell_data[joint_markers[["C.elegans"]][joint_markers[["C.elegans"]]$gene == cur_gene,]$cell_type,]$cell_class,
                         cell_data[joint_markers[["C.briggsae"]][joint_markers[["C.briggsae"]]$gene == cur_gene,]$cell_type,]$cell_class))
  
  if(is.na(cur_tissue)) {
    cur_tissue = "0_300"
  }
  
  if(cur_tissue == "progenitor") {
    mean_embryo_time = mean(c(cell_data[joint_markers[["C.elegans"]][joint_markers[["C.elegans"]]$gene == cur_gene,]$cell_type,]$embryo_time,
                              cell_data[joint_markers[["C.briggsae"]][joint_markers[["C.briggsae"]]$gene == cur_gene,]$cell_type,]$embryo_time),
                            na.rm = TRUE)
    
    if(mean_embryo_time < 150) {
      cur_tissue = "0_150"
    } else if (mean_embryo_time < 200) {
      cur_tissue = "0_200"
    } else if (mean_embryo_time < 250) {
      cur_tissue = "0_250"
    } else {
      cur_tissue = "0_300"
    }
  }
  
  png(paste0(dir, "Plots/gene_plots_fix/", cur_gene, ".png"), width = 1250, height = 1500)
  print(plot_gene_pattern(cur_gene, tissue = cur_tissue, gene_data, tissue_embeddings, gene_counts, TPMJoint, CellTableList, lo))
  dev.off()
}, mc.cores = 60)


############################
# Circle plotting function #
############################

CellTable <- readRDS(paste0(dir, "Objects/CellTable_Names_20240901.rds"))
tpm_matrix_list <- readRDS(paste0(dir, "Objects/tpm_matrix_list_cell_bg.rds"))

CellTableList <- list()

gene_names <- rownames(tpm_matrix_list[["C.elegans"]])
tpm_matrix_list_t <- list()
tpm_matrix_list_t[["C.elegans"]] <- data.frame(t(tpm_matrix_list[["C.elegans"]]))
tpm_matrix_list_t[["C.elegans"]]$Lineage <- rownames(tpm_matrix_list_t[["C.elegans"]])
colnames(tpm_matrix_list_t[["C.elegans"]]) <- c(gene_names, "Lineage")

CellTableList[["C.elegans"]] <- left_join(CellTable, tpm_matrix_list_t[["C.elegans"]], by = c("Lineage" = "Lineage"))
rownames(CellTableList[["C.elegans"]]) <- CellTableList[["C.elegans"]]$Lineage

gene_names <- rownames(tpm_matrix_list[["C.briggsae"]])
tpm_matrix_list_t <- list()
tpm_matrix_list_t[["C.briggsae"]] <- data.frame(t(tpm_matrix_list[["C.briggsae"]]))
tpm_matrix_list_t[["C.briggsae"]]$Lineage <- rownames(tpm_matrix_list_t[["C.briggsae"]])
colnames(tpm_matrix_list_t[["C.briggsae"]]) <- c(gene_names, "Lineage")

CellTableList[["C.briggsae"]] <- left_join(CellTable, tpm_matrix_list_t[["C.briggsae"]], by = c("Lineage" = "Lineage"))
rownames(CellTableList[["C.briggsae"]]) <- CellTableList[["C.briggsae"]]$Lineage

## Add dummy children to Z2/Z3
CellTable_plot <- CellTable[which(CellTable$level > 1),
                            c("Parent", "Lineage", "level")]
CellTable_plot <- rbind(CellTable_plot, data.frame(Parent = c("Z2", "Z2", "Z2", "Z2", "Z2", "Z2", "Z2", "Z2", "Z3", "Z3", "Z3", "Z3", "Z3", "Z3", "Z3", "Z3"),
                                 Lineage = c("Z2.1", "Z2.2", "Z2.3", "Z2.4", "Z2.5", "Z2.6", "Z2.7", "Z2.8", "Z3.1", "Z3.2", "Z3.3", "Z3.4", "Z3.5", "Z3.6", "Z3.7", "Z3.8"),
                                 level = 0))

gr <- graph_from_data_frame(CellTable_plot, directed = TRUE)
cells_ordered <- sort(V(gr)$name)
cells_ordered <- cells_ordered[! cells_ordered %in% c("E", "MS")]
cells_ordered <- c("MS", "E", cells_ordered)
cells_ordered <- cells_ordered[! is.na(cells_ordered)]

gr <- permute(gr, match(V(gr)$name, rev(cells_ordered)))
lo <- create_layout(gr, layout = "partition", circular = TRUE)

plot_circle_pattern <- function(cur_gene, CellTableList, lo) {
  layout_list <- list()
  
  for(cur_species in c("C.elegans", "C.briggsae")) {
    
    layout_list[[cur_species]] <- left_join(lo,
                                            CellTableList[[cur_species]][,c("Lineage", "level", cur_gene)],
                                            by = c("name" = "Lineage"))
    rownames(layout_list[[cur_species]]) <- layout_list[[cur_species]]$name
    colnames(layout_list[[cur_species]])[ncol(layout_list[["C.elegans"]])] <- "tpm"
    
    layout_list[[cur_species]] <- layout_list[[cur_species]][! is.na(layout_list[[cur_species]]$level),]
    layout_list[[cur_species]]$level <- as.numeric(layout_list[[cur_species]]$level)
  }
  
  tpm_max <- max(max(layout_list[["C.elegans"]]$tpm, na.rm = TRUE),
                 max(layout_list[["C.briggsae"]]$tpm, na.rm = TRUE))
  
  cel_plot <- ggraph(layout_list[["C.elegans"]]) +
    geom_node_arc_bar(aes(filter = level > 2, fill = tpm + 1),
                      color = "grey20", size = 0.25) +
    coord_fixed() +
    geom_node_text(aes(filter = level < 4 & level > 2,
                       label = name,
                       angle = node_angle(x, y) + 90), size = 3, colour = 'white') +
    annotate(geom = "text", x = 0.025, y = 0.175, hjust = 0.5, size = 9,
             label = expression(paste(italic("C. elegans"), " TPM"))) +
    scale_fill_viridis(trans = "log2", option = "magma",
                       limits = c(8, tpm_max), oob = scales::squish) +
    guides(alpha = F, size = F, fill = guide_colorbar(title = "")) +
    theme(legend.position = c(0.49,0.425),
          panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          plot.margin = margin(t = 0, r = -100, b = 0, l = -100, unit = "mm"),
          legend.direction = "horizontal",
          legend.justification = "center",
          legend.key.size = unit(9, 'mm'),
          legend.text = element_text(size = 12))
  
  cbr_plot <- ggraph(layout_list[["C.briggsae"]]) +
    geom_node_arc_bar(aes(filter = level > 2, fill = tpm + 1), 
                      color = "grey20", size = 0.25) +
    coord_fixed() +
    geom_node_text(aes(filter = level < 4 & level > 2,
                       label = name,
                       angle = node_angle(x, y) + 90), size = 3, colour = 'white') +
    annotate(geom = "text", x = 0.025, y = 0.175, hjust = 0.5, size = 9,
             label = expression(paste(italic("C. briggsae"), " TPM"))) +
    scale_fill_viridis(trans = "log2", option = "magma",
                       limits = c(8, tpm_max), oob = scales::squish) +
    guides(alpha = F, size = F, fill = guide_colorbar(title = "")) +
    theme(legend.position = c(0.49,0.425),
          panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          plot.margin = margin(t = 0, r = -100, b = 0, l = -100, unit = "mm"),
          legend.direction = "horizontal",
          legend.justification = "center",
          legend.key.size = unit(9, 'mm'),
          legend.text = element_text(size = 12))
  
  # layout_list[["Joint"]] <- layout_list[["C.elegans"]]
  # # layout_list[["Joint"]]$tpm <- log2(layout_list[["C.elegans"]]$tpm + 1) - log2(layout_list[["C.briggsae"]]$tpm + 1)
  # 
  # layout_list[["Joint"]] <- left_join(layout_list[["Joint"]], layout_list[["C.elegans"]][, c("name", "tpm")], suffix = c("", "_cel"), by = c("name" = "name"))
  # layout_list[["Joint"]] <- left_join(layout_list[["Joint"]], layout_list[["C.briggsae"]][, c("name", "tpm")], suffix = c("", "_cbr"), by = c("name" = "name"))
  # 
  # layout_list[["Joint"]]$max_tpm <- apply(layout_list[["Joint"]], 1, function(x) {
  #   if(! is.na(x["tpm_cel"]) & ! is.na(x["tpm_cbr"])) {
  #     return(as.numeric(max(x["tpm_cel"], x["tpm_cbr"], na.rm = TRUE)))
  #   } else {
  #     return(NA)
  #   }
  # })
  # num_cells <- length(layout_list[["Joint"]]$tpm)
  # layout_list[["Joint"]]$tpm <- log2((layout_list[["Joint"]]$tpm_cel + 1) / (layout_list[["Joint"]]$tpm_cbr + 1)) * log2(layout_list[["Joint"]]$max_tpm + 1)
  # 
  # total_tpm_max <- max(layout_list[["Joint"]]$max_tpm, na.rm = TRUE)
  # print(paste("Total max TPM: ", total_tpm_max))
  # display_cutoff <- log2(total_tpm_max) * 5
  # 
  # joint_plot <- ggraph(layout_list[["Joint"]]) +
  #   geom_node_arc_bar(aes(filter = level > 2, fill = tpm),
  #                     color = "grey20", size = 0.25) +
  #   coord_fixed() +
  #   geom_node_text(aes(filter = level < 4 & level > 2,
  #                      label = name,
  #                      angle = node_angle(x, y) + 90), size = 3, colour = 'white') +
  #   annotate(geom = "text", x = 0.025, y = 0.6, hjust = 0.5, size = 8,
  #            label = bquote(italic(.(gene)))) +
  #   annotate(geom = "text", x = 0.025, y = 0.15, hjust = 0.5, size = 5,
  #            label = expression(paste("log2( ", frac(paste(italic("C. elegans")),
  #                                                    paste(italic("C. briggsae"))), " ) * Max"))) +
  #   scale_fill_gradientn(colors = brewer.pal(9, 'RdBu'),
  #                        limits = c(-display_cutoff, display_cutoff),
  #                        oob = scales::squish) +
  #   guides(alpha = F, size = F, fill = guide_colorbar(title = "")) +
  #   theme(legend.position = c(0.49,0.4),
  #         panel.background = element_rect(fill = 'transparent'),
  #         plot.background = element_rect(fill = 'transparent', color = NA),
  #         plot.margin = margin(t = 0, r = -500, b = 0, l = -500, unit = "mm"),
  #         legend.direction = "horizontal",
  #         legend.justification = "center",
  #         legend.key.size = unit(7, 'mm'),
  #         legend.text = element_text(size = 12))
  
  #log2(C. elegans/C. briggsae
  # circle_plots <- plot_grid(cel_plot, joint_plot, cbr_plot, nrow = 3, scale = 1.05)
  circle_plots <- plot_grid(cel_plot, cbr_plot, nrow = 2, scale = 1.05)
}

pdf(paste0(dir, "Plots/supplemental_plots/test.pdf"), height = 15, width = 5)
print(circle_plots)
dev.off()

##############################
# Combined plotting function #
##############################

plot_gene_pattern <- function(gene, tissue = "none", gene_data, embeddings, counts_matrix, TPMList, CellTableList, lo) {
  
  cel_gene_name = gene_data[gene,]$elegans_id
  cbr_gene_name = gene_data[gene,]$briggsae_id
  
  if(! gene %in% gene_data$gene) {
    print("Not a 1:1 ortholog")
    print("Exiting")
    return(NULL)
  }

  plot_tissue = tissue
  
  ## Make global UMAPs
  global_umap_cel <- plotExpUmap(gene, embeddings[["global"]][embeddings[["global"]]$species == "C.elegans",],
                                 counts_matrix[gene, embeddings[["global"]]$species == "C.elegans"],
                                 "C.elegans", marker_size = 0.025)
  global_umap_cbr <- plotExpUmap(gene, embeddings[["global"]][embeddings[["global"]]$species == "C.briggsae",],
                                 counts_matrix[gene, embeddings[["global"]]$species == "C.briggsae"],
                                 "C.briggsae", marker_size = 0.025)
  
  ## Make tissue specific or progenitor UMAPs
  tissue_umap_cel <- plotExpUmap(gene, embeddings[[plot_tissue]][embeddings[[plot_tissue]]$species == "C.elegans",],
                                 counts_matrix[gene, embeddings[["global"]]$species == "C.elegans" & 
                                                 colnames(counts_matrix) %in% rownames(embeddings[[plot_tissue]])], "C.elegans",
                                 tissue, marker_size = 0.9)
  
  tissue_umap_cbr <- plotExpUmap(gene, embeddings[[plot_tissue]][embeddings[[plot_tissue]]$species == "C.briggsae",],
                                 counts_matrix[gene, embeddings[["global"]]$species == "C.briggsae" & 
                                                 colnames(counts_matrix) %in% rownames(embeddings[[plot_tissue]])], "C.briggsae",
                                 tissue, marker_size = 0.9)
  
  print("UMAPs finished")
  
  # TPM expression plot for cell types
  dotPlot_term <- plot_gene_example(gene, TPMList, BinarizeList, cell_data, "terminal")
  # dotPlot_pro <- plot_gene_example(gene, TPMList, BinarizeList, cell_data, "progenitor")
  print("Dot plots finished")
  
  circle_plot <- plot_circle_pattern(gene, CellTableList, lo)
  
  print("Starting histograms")
  gene_data_subset <- gene_data %>% filter(cel_max_tpm > 80 & cbr_max_tpm > 80)
  
  ## Joint
  jsd_plot_joint <- ggplot() +
    geom_density(data = gene_data_subset, aes(x = jsd_median_joint)) +
    geom_rect(data = gene_data[gene,], aes(xmin = jsd_lower_joint,
                                           xmax = jsd_upper_joint, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    geom_vline(aes(xintercept = gene_data[gene, "jsd_median_joint"]), color = "red") +
    scale_x_continuous(name = "Jensen-Shannon distance", limits = c(0, 1)) +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  cel_tau_plot_joint <- ggplot() +
    geom_density(data = gene_data_subset, aes(x = cel_tau_median_joint)) +
    geom_rect(data = gene_data[gene,], aes(xmin = cel_tau_lower_joint,
                                           xmax = cel_tau_upper_joint, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    geom_vline(aes(xintercept = gene_data[gene, "cel_tau_median_joint"]), color = "red") +
    scale_x_continuous(name = expression(paste(italic("C. elegans")," Tau")), limits = c(0, 1)) +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  cbr_tau_plot_joint <- ggplot() +
    geom_density(data = gene_data_subset, aes(x = cbr_tau_median_joint)) +
    geom_rect(data = gene_data[gene,], aes(xmin = cbr_tau_lower_joint,
                                           xmax = cbr_tau_upper_joint, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    geom_vline(aes(xintercept = gene_data[gene, "cbr_tau_median_joint"]), color = "red") +
    scale_x_continuous(name = expression(paste(italic("C. briggsae")," Tau")), limits = c(0, 1)) +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  ## Terminal
  jsd_plot_term <- ggplot() +
    geom_density(data = gene_data_subset, aes(x = jsd_median_term)) +
    geom_rect(data = gene_data[gene,], aes(xmin = jsd_lower_term,
                                           xmax = jsd_upper_term, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    geom_vline(aes(xintercept = gene_data[gene, "jsd_median_term"]), color = "red") +
    scale_x_continuous(name = "Jensen-Shannon distance", limits = c(0, 1)) +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  cel_tau_plot_term <- ggplot() +
    geom_density(data = gene_data_subset, aes(x = cel_tau_median_term)) +
    geom_rect(data = gene_data[gene,], aes(xmin = cel_tau_lower_term,
                                           xmax = cel_tau_upper_term, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    geom_vline(aes(xintercept = gene_data[gene, "cel_tau_median_term"]), color = "red") +
    scale_x_continuous(name = expression(paste(italic("C. elegans")," Tau")), limits = c(0, 1)) +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  cbr_tau_plot_term <- ggplot() +
    geom_density(data = gene_data_subset, aes(x = cbr_tau_median_term)) +
    geom_rect(data = gene_data[gene,], aes(xmin = cbr_tau_lower_term,
                                           xmax = cbr_tau_upper_term, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    geom_vline(aes(xintercept = gene_data[gene, "cbr_tau_median_term"]), color = "red") +
    scale_x_continuous(name = expression(paste(italic("C. briggsae")," Tau")), limits = c(0, 1)) +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  ## Progenitor
  jsd_plot_pro <- ggplot() +
    geom_density(data = gene_data_subset, aes(x = jsd_median_pro)) +
    geom_rect(data = gene_data[gene,], aes(xmin = jsd_lower_pro,
                                           xmax = jsd_upper_pro, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    geom_vline(aes(xintercept = gene_data[gene, "jsd_median_pro"]), color = "red") +
    scale_x_continuous(name = "Jensen-Shannon distance", limits = c(0, 1)) +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  cel_tau_plot_pro <- ggplot() +
    geom_density(data = gene_data_subset, aes(x = cel_tau_median_pro)) +
    geom_rect(data = gene_data[gene,], aes(xmin = cel_tau_lower_pro,
                                           xmax = cel_tau_upper_pro, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    geom_vline(aes(xintercept = gene_data[gene, "cel_tau_median_pro"]), color = "red") +
    scale_x_continuous(name = expression(paste(italic("C. elegans")," Tau")), limits = c(0, 1)) +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  cbr_tau_plot_pro <- ggplot() +
    geom_density(data = gene_data_subset, aes(x = cbr_tau_median_pro)) +
    geom_vline(aes(xintercept = gene_data[gene, "cbr_tau_median_pro"]), color = "red") +
    geom_rect(data = gene_data[gene,], aes(xmin = cbr_tau_lower_pro,
                                           xmax = cbr_tau_upper_pro, ymin = -Inf, ymax = Inf),
              color = "red3", fill = "red", alpha = 0.25) +
    scale_x_continuous(name = expression(paste(italic("C. briggsae")," Tau")), limits = c(0, 1)) +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.title = element_blank(),
          legend.position = "none")
  
  print("Finished histograms")
  
  metric_plots_joint <- plot_grid(jsd_plot_joint,
                                 cel_tau_plot_joint,
                                 cbr_tau_plot_joint,
                                 ncol = 3)
  
  metric_plots_term <- plot_grid(jsd_plot_term,
                                  cel_tau_plot_term,
                                  cbr_tau_plot_term,
                                  ncol = 3)
  
  metric_plots_pro <- plot_grid(jsd_plot_pro,
                                  cel_tau_plot_pro,
                                  cbr_tau_plot_pro,
                                  ncol = 3)
  
  print("Finished metric plots")
  
  # now add the title
  joint_title <- ggdraw() +
    draw_label("Progenitor + Terminal data",
               fontface = 'italic', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  metric_plots_joint <- plot_grid(joint_title, metric_plots_joint, ncol = 1, rel_heights = c(0.05, 1))
  
  term_title <- ggdraw() +
    draw_label("Terminal data",
               fontface = 'italic', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  metric_plots_term <- plot_grid(term_title, metric_plots_term, ncol = 1, rel_heights = c(0.05, 1))
  
  pro_title <- ggdraw() +
    draw_label("Progenitor data",
               fontface = 'italic', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))
  
  metric_plots_pro <- plot_grid(pro_title, metric_plots_pro, ncol = 1, rel_heights = c(0.05, 1))
  
  metric_plots_merge <- plot_grid(metric_plots_joint, metric_plots_term, metric_plots_pro, nrow = 3)
  
  if(cel_gene_name %in% synteny_filt$cel_WBG_name & cbr_gene_name %in% synteny_filt$cbr_WBG_name) {
    synteny_temp <- synteny_filt[synteny_filt$cel_WBG_name == cel_gene_name & synteny_filt$cbr_WBG_name == cbr_gene_name,]
    
    ortho_plot <- ggplot() +
      geom_textbox(aes(x = 0, y = 1,
                       label = paste0("The gene from C. elegans, ", gene, " is ", ifelse(synteny_temp$mutual_best_blat == 1, "mutual best BLAT ", "not mutual best BLAT "),
                                      "with its C. briggsae counterpart and is ", ifelse(synteny_temp$syntenic == "IS_SYNTENIC", "syntenic ", "not syntenic "), 
                                      " with the gene from C. briggsae, ", gene_data[gene,]$briggsae_gene_short_name, ". ",
                                      "The best match from C. elegans to C. briggsae by SW is ",  strsplit(synteny_temp$cel_cbr_sw_best, "#")[[1]][1],
                                      " with the next best match being ", strsplit(synteny_temp$cel_cbr_sw_next_best, "#")[[1]][1], ". ",
                                      "The best match from C. briggsae to C. elegans by SW is ",strsplit(synteny_temp$cbr_cel_sw_best, "#")[[1]][1],
                                      " with the next best match being ", strsplit(synteny_temp$cbr_cel_sw_next_best, "#")[[1]][1], ". ",
                                      "For the best alignment, the percent of the protein aligned for C. elegans is ", round(synteny_temp$cel_percent_align, 2),
                                      " and from C. briggsae is ", round(synteny_temp$cbr_percent_align, 2),
                                      " with ", round(synteny_temp$percent_identity, 2), " percent identity. ",
                                      "There are ", gene_data[gene,]$cel_OG_count, " genes in the orthogroup from C. elegans and ",
                                      gene_data[gene,]$br_OG_count, "genes from C. briggsae. The C. elegans WBgene name is ", cel_gene_name, " and the C. briggsae WBgene is ", cbr_gene_name, ".")),
                   width = unit(0.95, "npc"),
                   height = unit(0.95, "npc"),
                   hjust = 0, vjust = 1) +
      ggtitle(paste0(" There is a ",
                     generate_lagene_confidence(gene, gene_data),
                     " confidence in the orthology assignment")) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_void()
    
  } else {
    print("Huh?")
  }
  
  print("Plotting everything together")
  
  # now add the title
  main_title <- ggdraw() +
    draw_label(
      paste0(gene, ": ", gene_data[gene,]$WormCat.3),
      fontface = 'bold', x = 0, hjust = 0) +
    theme(plot.margin = margin(0, 0, 0, 7))

  umap_plots <- plot_grid(global_umap_cel, global_umap_cbr, tissue_umap_cel, tissue_umap_cbr)
  left_plots <- plot_grid(umap_plots, metric_plots_merge, ncol = 1, rel_heights = c(2, 1))
  right_plots <- plot_grid(dotPlot_term, circle_plot, ncol = 1, rel_heights = c(1, 1.8))
  
  bottom_plots <- plot_grid(left_plots, right_plots, nrow = 1, rel_widths = c(2, 1))

  main_plot <- plot_grid(ortho_plot, bottom_plots, rel_heights = c(0.08, 1.5), ncol = 1, nrow = 2)
  
  return(plot_grid(main_title, main_plot, ncol = 1, rel_heights = c(0.02, 1)))
}

TPMList <- TPMListBootstrapMean_term


png(paste0(dir, "Plots/test.png"), width = 1250, height = 1500)
print(plot_grid(main_title, main_plot, ncol = 1, rel_heights = c(0.02, 1)))
dev.off()

png(paste0(dir, "Plots/gene_figure_plots/hlh-4.png"), width = 1250, height = 1500)
print(plot_gene_pattern("hlh-4", tissue = "Ciliated neurons", gene_data, tissue_embeddings, gene_counts, TPMJoint, CellTableList, lo))
dev.off()

png(paste0(dir, "Plots/gene_figure_plots/lin-59.png"), width = 1250, height = 1500)
print(plot_gene_pattern("lin-59", tissue = "0_300", gene_data, tissue_embeddings, gene_counts, TPMList, CellTableList, lo))
dev.off()

png(paste0(dir, "Plots/gene_figure_plots/rab-7.png"), width = 1250, height = 1500)
print(plot_gene_pattern("rab-7", tissue = "0_300", gene_data, tissue_embeddings, gene_counts, TPMList, CellTableList, lo))
dev.off()

png(paste0(dir, "Plots/gene_figure_plots/atpf-2.png"), width = 1250, height = 1500)
print(plot_gene_pattern("aptf-2", tissue = "0_300", gene_data, tissue_embeddings, gene_counts, TPMList, CellTableList, lo))
dev.off()

png(paste0(dir, "Plots/test.png"), width = 1250, height = 1500)
print(plot_gene_pattern("cept-2", tissue = "0_300", gene_data, tissue_embeddings, gene_counts, TPMList, CellTableList, lo))
dev.off()

