library(monocle3)
library(Seurat)
library(ggplot2)
library(SeuratWrappers)

# Paths
OldPath <- "/kimdata/livinlrg/scAnalysis/BobDataComb/"
ExternalPath <- "/kimdata/livinlrg/scAnalysis/ExternalDataset/"

# Load data
dir <- "/kimdata/livinlrg/scAnalysis/WS290/"
cds <- readRDS(paste0(dir, "Objects/cds_no_bg_20240826.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

CellTable <- readRDS(paste0(dir, "Objects/CellTable_20240826.rds"))
FiltCellTypes <- readRDS(paste0(dir, "Objects/FiltCellTypes_20240820.rds"))
cells_refined_table <- read.delim(paste0(OldPath, "Objects/cells_refined.txt"))

time_bin_vector <- c("lt_100", "100_130", "130_170", "170_210", "210_270", "270_330", "330_390", "390_450", "450_510", "510_580", "580_650", "650_710", "gt_710")

time_bin_df <- data.frame(bins = time_bin_vector,
                                 start = c(0, 100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710),
                                 end = c(100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710, 1000))

# Add each cell type that is a single category to the cells_refined_table
for(cell_type in unique(CellTable$TerminalDatasetName)[! unique(CellTable$TerminalDatasetName) %in%
                                                       c("", "BWM_", "pm3_pm4_pm5a/b/c", "Intestine_", "GLR_1/GLR_2", "mc2a/mc2b")]) {
  if(! cell_type %in% cells_refined_table$cell_type) {
    print(cell_type)
    cells_refined_table <- rbind(cells_refined_table,
                                 data.frame(unique_id = cell_type,
                                            cell_type = cell_type,
                                            time_cutoff = NA,
                                            stringsAsFactors = FALSE))
  }
}

# check on:
names(table(colData(cds)$cell_type))[! names(table(colData(cds)$cell_type)) %in% cells_refined_table$cell_type]
# [1] "DA_DB"                  "hyp1V_and_ant_arc_V"    "MC"                    
# [4] "mu_int_mu_anal_related" "unassigned"             "URB_and_URA"
# joint labels for DA and DB and hyp1v, mu_int_mu_anal_related and ant arc V won't be used

cells_refined_table <- cells_refined_table[cells_refined_table$cell_type != "MC?",]
cells_refined_table <- cells_refined_table[cells_refined_table$cell_type != "URB_and_possibly_URA",]
cells_refined_table <- rbind(cells_refined_table,
      data.frame(unique_id = c("URB_and_URA", "MC"),
           cell_type = c("URB_and_URA", "MC"),
           time_cutoff = NA,
           stringsAsFactors = FALSE))

saveRDS(cells_refined_table, paste0(dir, "Objects/cells_refined_table_20240826.rds"))

## Create time bin tables for the refined cell type labels
# will need to convert smoothed.embryo.time into the time bins
colData(cds)$embryo.time.bin <- NA
for(bin in time_bin_df$bins) {
  colData(cds)[which(colData(cds)$smoothed.embryo.time >= time_bin_df[time_bin_df$bins == bin,]$start &
                 colData(cds)$smoothed.embryo.time < time_bin_df[time_bin_df$bins == bin,]$end),]$embryo.time.bin <- bin
}

# convert to a factor
col_data_df <- colData(cds)[! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"),]
col_data_df$embryo.time.bin <- factor(col_data_df$embryo.time.bin, levels = time_bin_vector)

cel_table_list <- list()
cbr_table_list <- list()

# count the number of cells from each of these time and cell type bins from each species
for(cell_type in unique(cells_refined_table$unique_id)) {
  temp_table <- cells_refined_table[cells_refined_table$unique_id %in% cell_type,,FALSE]
  print(cell_type)
  if(nrow(temp_table) > 1) {
    time_cutoff = max(c(temp_table$time_cutoff, 0), na.rm = TRUE)
    if(time_cutoff > 0) {
      cel_table_list[[cell_type]] <- table(col_data_df[col_data_df$cell_type %in% temp_table$cell_type & col_data_df$smoothed.embryo.time > time_cutoff & col_data_df$species == "C.elegans",]$embryo.time.bin)
      cbr_table_list[[cell_type]] <- table(col_data_df[col_data_df$cell_type %in% temp_table$cell_type & col_data_df$smoothed.embryo.time > time_cutoff & col_data_df$species == "C.briggsae",]$embryo.time.bin)
    } else{ 
      cel_table_list[[cell_type]] <- table(col_data_df[col_data_df$cell_type %in% temp_table$cell_type & col_data_df$species == "C.elegans",]$embryo.time.bin)
      cbr_table_list[[cell_type]] <- table(col_data_df[col_data_df$cell_type %in% temp_table$cell_type & col_data_df$species == "C.briggsae",]$embryo.time.bin)
    }
  } else {
    cel_table_list[[cell_type]] <- table(col_data_df[col_data_df$cell_type %in% temp_table$cell_type & col_data_df$species == "C.elegans",]$embryo.time.bin)
    cbr_table_list[[cell_type]] <- table(col_data_df[col_data_df$cell_type %in% temp_table$cell_type & col_data_df$species == "C.briggsae",]$embryo.time.bin)
  }
}

cel_table <- do.call("rbind", cel_table_list)
cbr_table <- do.call("rbind", cbr_table_list)

write.table(cel_table, sep = "\t", file = paste0(dir, "CellQC/ele_time_bins_table.txt"))
write.table(cbr_table, sep = "\t", file = paste0(dir, "CellQC/bri_time_bins_table.txt"))

base_cutoff = 40
bin_list <- list()
cell_type_bins <- data.frame(cell_type = "",
                             start_bin = "",
                             end_bin = "",
                             elegans_count = 0,
                             briggsae_count = 0,
                             stringsAsFactors = FALSE)

# anything that isn't in cell_refined.txt will be added as a single line
# want to also curate the second column so that the ambiguous cell types aren't worked on
for(cell_type in sort(unique(cells_refined_table$unique_id))) {
  cel_cell_total <- 0
  cbr_cell_total <- 0
  min_cell_count <- min(sum(cel_table[cell_type,]), sum(cbr_table[cell_type,]))
  temp_cell_type_bins <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames <- c("cell_type", "start_bin", "end_bin", "elegans_count", "briggsae_count")
  colnames(temp_cell_type_bins) <- colnames
  
  n_bin = 0
  i = 1
  print(cell_type)
  ## Need to adjust cutoff such that we only end up with a total of three bins
  ## Also need to find a way to modify the structure such that it can accept cells from a linear series of cells > a certain time point
  while(n_bin < 3 & i < 20) {
    # print(paste(i, cutoff))
    cutoff = base_cutoff + (min_cell_count / i)
    temp_cell_type_bins <- data.frame(matrix(nrow = 0, ncol = 5))
    colnames <- c("cell_type", "start_bin", "end_bin", "elegans_count", "briggsae_count")
    cel_cell_total = 0
    cbr_cell_total = 0
    
    start_bin_name = "lt_100"
    for(time_bin in time_bin_vector) {
      # print(paste(cell_type, time_bin, cel_cell_total, cbr_cell_total))
      if(cel_cell_total == 0 & cbr_cell_total == 0) {
        start_bin_name = time_bin
      }
      
      cel_cell_total <- cel_cell_total + cel_table[cell_type, time_bin]
      cbr_cell_total <- cbr_cell_total + cbr_table[cell_type, time_bin]
      
      if(cel_cell_total >= cutoff & cbr_cell_total >= cutoff) {
        temp_cell_type_bins <- rbind(temp_cell_type_bins, data.frame(cell_type = cell_type,
                                                           start_bin = start_bin_name,
                                                           end_bin = time_bin,
                                                           elegans_count = cel_cell_total,
                                                           briggsae_count = cbr_cell_total))
        start_bin_name = time_bin
        cel_cell_total <- 0
        cbr_cell_total <- 0
      }
    }
    n_bin <- nrow(temp_cell_type_bins)
    # print(n_bin)
    i = i + 1
  }
  
  remaining_cells = max(sum(cel_table[cell_type,]) - sum(temp_cell_type_bins$elegans_count),
                        sum(cbr_table[cell_type,]) - sum(temp_cell_type_bins$briggsae_count))
  if(n_bin < 3 & remaining_cells > base_cutoff) {
    cutoff = base_cutoff
    temp_cell_type_bins <- data.frame(matrix(nrow = 0, ncol = 5))
    colnames <- c("cell_type", "start_bin", "end_bin", "elegans_count", "briggsae_count")
    cel_cell_total = 0
    cbr_cell_total = 0
    
    start_bin_name = "lt_100"
    for(time_bin in time_bin_vector) {
      print(paste(cell_type, time_bin, cel_cell_total, cbr_cell_total))
      if(cel_cell_total == 0 & cbr_cell_total == 0) {
        start_bin_name = time_bin
      }
      
      cel_cell_total <- cel_cell_total + cel_table[cell_type, time_bin]
      cbr_cell_total <- cbr_cell_total + cbr_table[cell_type, time_bin]
      
      if(cel_cell_total >= cutoff & cbr_cell_total >= cutoff) {
        temp_cell_type_bins <- rbind(temp_cell_type_bins, data.frame(cell_type = cell_type,
                                                                     start_bin = start_bin_name,
                                                                     end_bin = time_bin,
                                                                     elegans_count = cel_cell_total,
                                                                     briggsae_count = cbr_cell_total))
        start_bin_name = time_bin
        cel_cell_total <- 0
        cbr_cell_total <- 0
      }
    }
  }
  
  # Reworked this section so that it just adjusts the last bin to the maximum time point that there are cells
  if(nrow(temp_cell_type_bins) > 0) {
    bin_number <- which(time_bin_vector == temp_cell_type_bins[nrow(temp_cell_type_bins),]$end_bin)
    if(bin_number != 13) {
      if(sum(cel_table[cell_type, (bin_number + 1):13]) > 0 |
         sum(cbr_table[cell_type, (bin_number + 1):13]) > 0) {
        temp_cell_type_bins[nrow(temp_cell_type_bins),]$end_bin <- time_bin_vector[max(seq(1,13)[cel_table[cell_type,] > 5 | cbr_table[cell_type,] > 5])]
        
        last_bin_names <- time_bin_vector[which(time_bin_vector == temp_cell_type_bins[nrow(temp_cell_type_bins),]$start_bin):which(time_bin_vector == temp_cell_type_bins[nrow(temp_cell_type_bins),]$end_bin)]
        temp_cell_type_bins[nrow(temp_cell_type_bins),]$elegans_count <- sum(cel_table[cell_type, last_bin_names])
        temp_cell_type_bins[nrow(temp_cell_type_bins),]$briggsae_count <- sum(cbr_table[cell_type, last_bin_names])
      }
    }
  }
  cell_type_bins <- rbind(cell_type_bins, temp_cell_type_bins)
}
cell_type_bins <- cell_type_bins[-1,]


# Add cell bins for cell types that don't have a bin at the max and min times
for(cell_type in setdiff(unique(cells_refined_table$unique_id), cell_type_bins$cell_type)) {
  print(cell_type)
  cell_type_bins <- rbind(cell_type_bins, data.frame(cell_type = cell_type,
                                                     start_bin =   time_bin_vector[cel_table[cell_type, time_bin_vector] > 0 | cbr_table[cell_type, time_bin_vector] > 0][1],
                                                     end_bin = tail(time_bin_vector[cel_table[cell_type, time_bin_vector] > 0 | cbr_table[cell_type, time_bin_vector] > 0], n = 1),
                                                     elegans_count = sum(cel_table[cell_type, ]),
                                                     briggsae_count = sum(cbr_table[cell_type, ])))
}

# don't know why these are weird:
# [1] "ADEsh"
# [1] "hyp3"
# [1] "pm8"

start_vector <- c(0, 100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710)
names(start_vector) <- time_bin_vector

end_vector <- c(100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710, 830)
names(end_vector) <- time_bin_vector

for(cell_type in unique(cell_type_bins$cell_type)) {
  pdf(paste0(dir, "CellQC/Embryo_time_plots/", cell_type, "_time_plot.pdf"), width = 5, height = 2.5)
  print(ggplot() +
          geom_histogram(data = data.frame(colData(cds_filt)[colData(cds_filt)$cell_type %in% unique(cells_refined_table[cells_refined_table$unique_id %in% cell_type,]$cell_type),]),
                         aes(x = smoothed.embryo.time, fill = species),
                         alpha = 0.5, position="identity") +
          geom_vline(aes(xintercept = start_vector[cell_type_bins[cell_type_bins$cell_type == cell_type,]$start_bin]), linetype = "dotted") +
          geom_vline(aes(xintercept = end_vector[cell_type_bins[cell_type_bins$cell_type == cell_type,]$end_bin]), linetype = "dotted") +
          scale_x_continuous(name = "Embryo time", limits = c(0, 830)) +
          scale_y_continuous(name = "Cell count") +
          scale_fill_manual(name = "Species", values = c("#009E73", "#56B4E9"), limits = c("C.elegans", "C.briggsae"), labels = c("C. elegans", "C. briggsae")) +
          theme(legend.title= element_blank(),
                rect = element_rect(fill = "transparent"),
                axis.line = element_line(color="grey80", size=1),
                axis.text = element_text(size = 12),
                axis.title = element_text(size = 12),
                panel.grid = element_blank(),
                panel.grid.major.y = element_line(color="grey80", size=0.25),
                panel.grid.minor.y = element_line(color="grey80", size=0.05),
                panel.background = element_rect(fill='transparent'), #transparent panel bg
                plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
                legend.key = element_rect(colour = "transparent", fill = "transparent")))
  dev.off()
}

saveRDS(cell_type_bins, paste0(dir, "Objects/cell_type_bins.rds"))

#pseudotime = >5 and >20
## Add pseudotime bins for germline

# add bins by pseudotime using monocle3
so_cell_class_list <- readRDS(paste0(dir, "Objects/so_cell_class_list_20240827.rds"))

gemline_so <- JoinLayers(so_cell_class_list[["Germline"]])
gemline_so[["umap"]] <- gemline_so[["umap.rpca"]]
germline_cds <- as.cell_data_set(gemline_so)

germline_cds <- cluster_cells(germline_cds)
germline_cds <- learn_graph(germline_cds)

pdf(paste0(dir, "germline_plot.pdf"),  width = 7.5, height = 7.5)
plot_cells(germline_cds,
           color_cells_by = "smoothed.embryo.time",
           label_groups_by_cluster=TRUE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_principal_points = TRUE)
dev.off()

# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, time = 100){
  cell_ids <- which(colData(cds)[, "embryo.time"] < time)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

pn_early <- get_earliest_principal_node(germline_cds, time = 100)

# The auto chosen one wasn't good
germline_cds <- order_cells(germline_cds, root_pr_nodes = "Y_5")

pdf(paste0(dir, "germline_plot.pdf"),  width = 5, height = 5)
plot_cells(germline_cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)
dev.off()

colData(cds)$germline_pseudotime <- NA
colData(cds)[names(pseudotime(germline_cds)),]$germline_pseudotime <- pseudotime(germline_cds)

saveRDS(cds, paste0(dir, "Objects/cds_no_bg_20240828.rds"))

# looks good
pdf(paste0(dir, "germline_plot_pseudotime_packer.pdf"),  width = 5, height = 5)
ggplot(data.frame(colData(cds)[colData(cds)$cell_type == "germline" & colData(cds)$lineage_packer %in% c("Z2/Z3:pseudotime_bin_1",
                                                                                                         "Z2/Z3:pseudotime_bin_2",
                                                                                                         "Z2/Z3:pseudotime_bin_3"),])) +
  geom_boxplot(aes(x = lineage_packer, y = germline_pseudotime), fill = "grey", color = "black") +
  scale_y_continuous(breaks = seq(0, 16, 1))
dev.off()

cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

# Create filters
cell_type_time_bins <- cell_type_bins
cell_type_time_bins <- cell_type_time_bins[cell_type_time_bins$cell_type != "germline",]

cell_type_time_bins <- rbind(cell_type_time_bins, data.frame(cell_type = "germline",
                                                             start_bin = "pseudotime_0",
                                                             end_bin = "pseudotime_10",
                                                             elegans_count = dim(colData(cds_filt)[colData(cds_filt)$species == "C.elegans" &
                                                                                                 colData(cds_filt)$cell_type == "germline" &
                                                                                                 colData(cds_filt)$germline_pseudotime < 5,])[1],
                                                             briggsae_count = dim(colData(cds_filt)[colData(cds_filt)$species == "C.briggsae" &
                                                                                                  colData(cds_filt)$cell_type == "germline" &
                                                                                                  colData(cds_filt)$germline_pseudotime < 5,])[1]))
cell_type_time_bins <- rbind(cell_type_time_bins, data.frame(cell_type = "germline",
                                                             start_bin = "pseudotime_10",
                                                             end_bin = "pseudotime_13",
                                                             elegans_count = dim(colData(cds_filt)[colData(cds_filt)$species == "C.elegans" &
                                                                                                 colData(cds_filt)$cell_type == "germline" &
                                                                                                 colData(cds_filt)$germline_pseudotime >= 10 & colData(cds_filt)$germline_pseudotime < 13.5,])[1],
                                                             briggsae_count = dim(colData(cds_filt)[colData(cds_filt)$species == "C.briggsae" &
                                                                                                  colData(cds_filt)$cell_type == "germline" &
                                                                                                  colData(cds_filt)$germline_pseudotime >= 10 & colData(cds_filt)$germline_pseudotime < 13.5,])[1]))
cell_type_time_bins <- rbind(cell_type_time_bins, data.frame(cell_type = "germline",
                                                             start_bin = "pseudotime_14",
                                                             end_bin = "pseudotime_30",
                                                             elegans_count = dim(colData(cds_filt)[colData(cds_filt)$species == "C.elegans" &
                                                                                                 colData(cds_filt)$cell_type == "germline" &
                                                                                                 colData(cds_filt)$germline_pseudotime >= 13.5 & colData(cds_filt)$germline_pseudotime < 30,])[1],
                                                             briggsae_count = dim(colData(cds_filt)[colData(cds_filt)$species == "C.briggsae" &
                                                                                                  colData(cds_filt)$cell_type == "germline" &
                                                                                                  colData(cds_filt)$germline_pseudotime >= 13.5 & colData(cds_filt)$germline_pseudotime < 30,])[1]))

#Ensure it's all there
# cell_type_time_bins in the cols of the TPMBinsList
# cols of the TPMBinsList in the cell_type_time_bins
summary(paste0(cell_type_time_bins$cell_type, "_", cell_type_time_bins$start_bin, ":", cell_type_time_bins$end_bin) %in% colnames(TPMBinsList[["C.elegans"]])) # missing 0
paste0(cell_type_time_bins$cell_type, "_", cell_type_time_bins$start_bin, ":", cell_type_time_bins$end_bin)[! paste0(cell_type_time_bins$cell_type, "_", cell_type_time_bins$start_bin, ":", cell_type_time_bins$end_bin) %in% colnames(TPMBinsList[["C.elegans"]])]

summary(paste0(cell_type_time_bins$cell_type, "_", cell_type_time_bins$start_bin, ":", cell_type_time_bins$end_bin) %in% colnames(TPMBinsList[["C.briggsae"]])) # missing 0
paste0(cell_type_time_bins$cell_type, "_", cell_type_time_bins$start_bin, ":", cell_type_time_bins$end_bin)[! paste0(cell_type_time_bins$cell_type, "_", cell_type_time_bins$start_bin, ":", cell_type_time_bins$end_bin) %in% colnames(TPMBinsList[["C.briggsae"]])]

summary(colnames(TPMBinsList[["C.elegans"]]) %in% paste0(cell_type_time_bins$cell_type, "_", cell_type_time_bins$start_bin, ":", cell_type_time_bins$end_bin)) # missing 0
summary(colnames(TPMBinsList[["C.briggsae"]]) %in% paste0(cell_type_time_bins$cell_type, "_", cell_type_time_bins$start_bin, ":", cell_type_time_bins$end_bin)) # missing 0

## Generate the name and embryo time bins for the columns
paste0(cell_type_time_bins$cell_type, "_", cell_type_time_bins$start_bin, ":", cell_type_time_bins$end_bin)

saveRDS(cell_type_time_bins, paste0(dir, "Objects/cell_type_time_bins.rds"))

cell_type_time_bins <- readRDS(paste0(dir, "Objects/cell_type_time_bins.rds"))

##################################################################
# Generate a list of barcodes for each bin to be used to bootstrapping
BarcodeBinsList <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  BarcodeBinsList[[cur_species]] <- lapply(unique(cell_type_time_bins[cell_type_time_bins$cell_type != "germline",]$cell_type), function(cell_type) { 
    print(cell_type)
    temp_barcodes_out <- vector(mode = "list", 
                                length = dim(cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,,FALSE])[1])
    for(i in seq(1, dim(cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,,FALSE])[1])) {
      print(paste0(cell_type, ": ", i))
      if(cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,][i,][,ifelse(cur_species == "C.elegans",
                                                                            "elegans_count",
                                                                            "briggsae_count")] > 0) {
        temp_barcodes_out[[i]] <- colnames(cds_filt[rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                                        "ce.unique",
                                                                        "cb.unique") |
                                    rowData(cds_filt)$gene.type == "common",
                                  colData(cds_filt)$species == cur_species &
                                    colData(cds_filt)$cell_type %in% unique(cells_refined_table[cells_refined_table$unique_id %in% cell_type,]$cell_type) &
                                    colData(cds_filt)$embryo.time.bin %in% time_bin_vector[which(time_bin_vector == cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,][i,]$start_bin):
                                                                                                  which(time_bin_vector == cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,][i,]$end_bin)]])
      }
    }
    names(temp_barcodes_out) <- paste0(cell_type, "_", cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$start_bin, ":",
                                       cell_type_time_bins[cell_type_time_bins$cell_type == cell_type,]$end_bin)
    return(temp_barcodes_out)
  })
  names(BarcodeBinsList[[cur_species]]) <- unique(cell_type_time_bins[cell_type_time_bins$cell_type != "germline",]$cell_type)
}

# Add in germline pseudotime bins as well
for(cur_species in c("C.elegans", "C.briggsae")) {
  BarcodeBinsList[[cur_species]][["germline"]] <- list()
  BarcodeBinsList[[cur_species]][["germline"]][[paste0("germline", "_", "pseudotime_0", ":", "pseudotime_10")]] <- colnames(
    cds_filt[rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                   "ce.unique",
                                                   "cb.unique") |
               rowData(cds_filt)$gene.type == "common",
             which(colData(cds_filt)$species == cur_species &
               colData(cds_filt)$cell_type == "germline" &
               colData(cds_filt)$germline_pseudotime < 10)])
  BarcodeBinsList[[cur_species]][["germline"]][[paste0("germline", "_", "pseudotime_10", ":", "pseudotime_13")]] <- colnames(
    cds_filt[rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                   "ce.unique",
                                                   "cb.unique") |
               rowData(cds_filt)$gene.type == "common",
             which(colData(cds_filt)$species == cur_species &
               colData(cds_filt)$cell_type == "germline" &
               colData(cds_filt)$germline_pseudotime >= 10 & colData(cds_filt)$germline_pseudotime < 13.5)])
  BarcodeBinsList[[cur_species]][["germline"]][[paste0("germline", "_", "pseudotime_14", ":", "pseudotime_30")]] <- colnames(
    cds_filt[rowData(cds_filt)$gene.type == ifelse(cur_species == "C.elegans",
                                                   "ce.unique",
                                                   "cb.unique") |
               rowData(cds_filt)$gene.type == "common",
             which(colData(cds_filt)$species == cur_species &
               colData(cds_filt)$cell_type == "germline" &
               colData(cds_filt)$germline_pseudotime >= 13.5 & colData(cds_filt)$germline_pseudotime < 30)])
}

saveRDS(BarcodeBinsList, paste0(dir, "Objects/BarcodeBinsList.rds"))

# Generate the barcode bins list for the progenitors now
ProBarcodeBinsList <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  ProBarcodeBinsList[[cur_species]] <- lapply(unique(colData(cds_filt)$lineage_broad), function(cell_type) { 
    print(cell_type)
    temp_barcodes_out <- list()
    temp_barcodes_out[[cell_type]] <- colnames(cds_filt[, colData(cds_filt)$species == cur_species &
                                             colData(cds_filt)$lineage_broad %in% cell_type])
    return(temp_barcodes_out)
  })
  names(ProBarcodeBinsList[[cur_species]]) <- unique(colData(cds_filt)$lineage_broad)
}

#remove unassigned
ProBarcodeBinsList[["C.elegans"]] <- ProBarcodeBinsList[["C.elegans"]][-which(names(ProBarcodeBinsList[["C.elegans"]]) == "unassigned")]
ProBarcodeBinsList[["C.briggsae"]] <- ProBarcodeBinsList[["C.briggsae"]][-which(names(ProBarcodeBinsList[["C.briggsae"]]) == "unassigned")]

saveRDS(ProBarcodeBinsList, paste0(dir, "Objects/ProBarcodeBinsList.rds"))

JointBarcodeBinsList <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  JointBarcodeBinsList[[cur_species]] <- c(BarcodeBinsList[[cur_species]], ProBarcodeBinsList[[cur_species]])
}

saveRDS(JointBarcodeBinsList, paste0(dir, "Objects/JointBarcodeBinsList.rds"))
