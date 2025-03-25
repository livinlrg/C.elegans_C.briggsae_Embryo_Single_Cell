library(Palo)
library(ggplot2)
library(ggExtra)
library(dplyr)
library(Seurat)
library(monocle3)
library(ggraph)

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240903.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]
rm(cds)

seurat <- readRDS(paste0(dir, "Objects/seurat_objects/seurat_rpca_integration_joined_20240823.rds"))
umap_proj <- data.frame(Embeddings(seurat, reduction = "umap.rpca"))

marker_size = 0.0075
na.col = "lightgrey"

cl <- as.character(colData(cds_filt)$cell_type)

gg_color_hue <- function(n) {
  hues = seq(30, 1000, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(length(unique(cl)))

palopal_terminal <- Palo(position = umap_proj[,c("umaprpca_1", "umaprpca_2")],
                         cluster = cl,
                         palette = pal,
                         rgb_weight = c(3, 4 ,2),
                         init_iter = 1000,
                         refine_iter = 2000,
                         early_stop = 500)


cl <- as.character(colData(cds_filt)$lineage_broad)

gg_color_hue <- function(n) {
  hues = seq(30, 1000, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal <- gg_color_hue(length(unique(cl)))

palopal_lineage <- Palo(position = umap_proj[,c("umaprpca_1", "umaprpca_2")],
                        cluster = cl,
                        palette = pal,
                        rgb_weight = c(3, 4 ,2),
                        init_iter = 1000,
                        refine_iter = 2000,
                        early_stop = 500)

palopal_comb <- c(palopal_terminal, palopal_lineage)

umap_proj$cell_type <- colData(cds_filt)[rownames(umap_proj),]$cell_type
umap_proj$lineage_broad <- colData(cds_filt)[rownames(umap_proj),]$lineage_broad
umap_proj$species <- colData(cds_filt)[rownames(umap_proj),]$species

umap_proj$joint_id <- ifelse(umap_proj$cell_type == "unassigned", umap_proj$lineage_broad, umap_proj$cell_type)

saveRDS(palopal_terminal, paste0(dir, "Objects/palopal_terminal.rds"))
saveRDS(palopal_lineage, paste0(dir, "Objects/palopal_lineage.rds"))
saveRDS(palopal_comb, paste0(dir, "Objects/palopal_comb.rds"))

###########################
#         UMAP Plot       #
###########################
## Plot UMAP

idx_region <- which(umap_proj$joint_id != "unassigned")
marker_size = 0.08
marker_stroke = 0.02
na.col = "lightgrey"

cell_type_umap_elegans <- ggplot(umap_proj[which(umap_proj$species == "C.elegans"),], aes(x = umaprpca_1, y = umaprpca_2)) +
  geom_point(size = marker_size, stroke = marker_stroke, color = na.col, show.legend = FALSE) +
  geom_point(data = umap_proj[intersect(idx_region, which(umap_proj$species == "C.elegans")),],
             aes(x = umaprpca_1, y = umaprpca_2, color = joint_id),
             size = 0.15,
             alpha = 0.8, stroke = 0.05) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = -Inf,
           yend = min(umap_proj$umaprpca_2) + (abs(min(umap_proj$umaprpca_2)) + abs(max(umap_proj$umaprpca_2)))/4, color = "grey80", size = 2) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = min(umap_proj$umaprpca_1) + (abs(min(umap_proj$umaprpca_1)) + abs(max(umap_proj$umaprpca_1)))/4,
           yend = -Inf, color = "grey80", size = 2) +
  scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(umap_proj$umaprpca_1), max(umap_proj$umaprpca_1), by = 0.20),1)) +
  scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(umap_proj$umaprpca_2), max(umap_proj$umaprpca_2), by = 0.20),1)) +
  scale_color_manual(values = palopal_comb) +
  theme(rect = element_rect(fill = "transparent"),
        plot.title = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(hjust=0.1),
        axis.title.y = element_text(hjust=0.1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

cell_type_umap_briggsae <- ggplot(umap_proj[which(umap_proj$species == "C.briggsae"),], aes(x = umaprpca_1, y = umaprpca_2)) +
  geom_point(size = marker_size, stroke = marker_stroke, color = na.col, show.legend=FALSE) +
  geom_point(data = umap_proj[intersect(idx_region, which(umap_proj$species == "C.briggsae")),],
             aes(x = umaprpca_1, y = umaprpca_2, color = joint_id),
             size = 0.15,
             alpha = 0.8, stroke = 0.05) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = -Inf,
           yend = min(umap_proj$umaprpca_2) + (abs(min(umap_proj$umaprpca_2)) + abs(max(umap_proj$umaprpca_2)))/4, color = "grey80", size = 2) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = min(umap_proj$umaprpca_1) + (abs(min(umap_proj$umaprpca_1)) + abs(max(umap_proj$umaprpca_1)))/4,
           yend = -Inf, color = "grey80", size = 2) +
  scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(umap_proj$umaprpca_1), max(umap_proj$umaprpca_1), by = 0.20),1)) +
  scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(umap_proj$umaprpca_2), max(umap_proj$umaprpca_2), by = 0.20),1)) +
  scale_color_manual(values = palopal_comb) +
  theme(rect = element_rect(fill = "transparent"),
        plot.title = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(hjust=0.1),
        axis.title.y = element_text(hjust=0.1),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

pdf(paste0(dir, "Plots/data_figure/cell_type_umap_elegans.pdf"), height = 6, width = 6)
cell_type_umap_elegans
dev.off()

pdf(paste0(dir, "Plots/data_figure/cell_type_umap_briggsae.pdf"), height = 6, width = 6)
cell_type_umap_briggsae
dev.off()

png(paste0(dir, "Plots/data_figure/cell_type_umap_briggsae.png"), width = 1500, height = 750, units = "px", res = 300)
plot_grid(cell_type_umap_elegans, cell_type_umap_briggsae, ncol = 2)
dev.off()

###########
# Cell count
###########

BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))
ProBarcodeBinsList <- readRDS(paste0(dir, "Objects/ProBarcodeBinsList.rds"))
JointBarcodeBinsList <- readRDS(paste0(dir, "Objects/JointBarcodeBinsList.rds"))

cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))

cell_count <- lapply(c("C.elegans", "C.briggsae"), function(cur_species) {
  data.frame(cell_count = unlist(lapply(JointBarcodeBinsList[[cur_species]], function(cell_type) {
    length(unlist(cell_type))
  }
  )),
  cell_type = names(JointBarcodeBinsList[[cur_species]]),
  species = cur_species)
})
names(cell_count) <- c("C.elegans", "C.briggsae")
# cell_count <- do.call(rbind, cell_count)

cell_count[["C.elegans"]]$cell_prop <- (cell_count[["C.elegans"]]$cell_count / sum(cell_count[["C.elegans"]]$cell_count)) * 100
cell_count[["C.briggsae"]]$cell_prop <- (cell_count[["C.briggsae"]]$cell_count / sum(cell_count[["C.briggsae"]]$cell_count)) * 100

CellTypeCount <- left_join(cell_count[["C.elegans"]], cell_count[["C.briggsae"]], by = "cell_type", suffix = c("_cel", "_cbr"))

CellTypeCount <- left_join(CellTypeCount, cell_data[,c("cell_type", "cell_class", "div_stage")], by = "cell_type")
CellTypeCount$plot <- CellTypeCount$cell_class
CellTypeCount[which(CellTypeCount$cell_class == "progenitor"),]$plot <- as.character(CellTypeCount[which(CellTypeCount$cell_class == "progenitor"),]$div_stage)

pdf(paste0(dir, "Plots/data_figure/cell_count_prop.pdf"), height = 5, width = 5)
CellTypeCount %>%
  filter(! cell_type %in% "hyp3") %>%
ggplot(aes(x = cell_prop_cel,
           y = cell_prop_cbr,
           label = cell_type,
           fill = plot,
           color = plot,
           shape = ifelse(cell_class == "progenitor", "pro", "term"))) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, alpha = 0.5, linetype = "dashed") +
  guides(label = "none") + 
  scale_x_continuous(name = "Percent of C. elegans cell types",
                     trans = "log10", limits = c(0.002, 7.5)) +
  scale_y_continuous(name = "Percent of C. briggsae cell types",
                     trans = "log10", limits = c(0.002, 7.5)) +
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
  scale_shape_manual(limits = c("pro", "term"), values = c(22, 21)) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(color = "#009E73"),
        axis.title.y = element_text(color = "#56B4E9"),
        panel.grid.major = element_line(color="grey80", size=0.5),
        panel.grid.minor = element_line(color="grey80", size=0.1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()


# legend
pdf(paste0(dir, "Plots/data_figure/cell_count_prop_legend.pdf"), height = 5, width = 10)
CellTypeCount %>%
  filter(! cell_type %in% "hyp3") %>%
  ggplot(aes(x = cell_prop_cel,
             y = cell_prop_cbr,
             label = cell_type,
             fill = plot,
             color = plot,
             shape = ifelse(cell_class == "progenitor", "pro", "term"))) +
  geom_point(alpha = 0.8) +
  geom_abline(slope = 1, alpha = 0.5, linetype = "dashed") +
  guides(label = "none") + 
  scale_x_continuous(name = "Percent of C. elegans cell types",
                     trans = "log10", limits = c(0.002, 7.5)) +
  scale_y_continuous(name = "Percent of C. briggsae cell types",
                     trans = "log10", limits = c(0.002, 7.5)) +
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
  scale_shape_manual(limits = c("pro", "term"), values = c(22, 21)) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        axis.title.x = element_text(color = "#009E73"),
        axis.title.y = element_text(color = "#56B4E9"),
        panel.grid.major = element_line(color="grey80", size=0.5),
        panel.grid.minor = element_line(color="grey80", size=0.1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()



#############
# Circle plot
#############

gr_list <- list()
layout_list <- list()
CellTableList <- list()

tpm_matrix_list <- readRDS(paste0(dir, "Objects/tpm_matrix_list_filt_cell_bg.rds"))
CellTable <- readRDS(paste0(dir, "Objects/CellTable_Names_20240901.rds"))

cur_speices = "cell_count"
library(ggraph)
library(igraph)

CellTable$dataset <- "No data"
CellTable[CellTable$TerminalDatasetName != "", "dataset"] <- "Annotated"
CellTable[! is.na(CellTable$cel_lineage_specific),"dataset"] <- "C. elegans"
CellTable[! is.na(CellTable$cbr_lineage_specific),"dataset"] <- "C. briggsae"
CellTable[! is.na(CellTable$cel_lineage_specific) & ! is.na(CellTable$cbr_lineage_specific),"dataset"] <- "Annotated"
CellTable[CellTable$Cell == "death", "dataset"] <- "Cell death"

gr <- graph_from_data_frame(CellTable[CellTable$level != "#N/A", c("Parent", "Lineage", "level", "dataset")],
                            directed = TRUE)
cells_ordered <- sort(V(gr)$name)
cells_ordered <- cells_ordered[! cells_ordered %in% c("E", "MS")]
cells_ordered <- c("MS", "E", cells_ordered)
cells_ordered <- cells_ordered[! is.na(cells_ordered)]

gr <- permute(gr, match(V(gr)$name, rev(cells_ordered)))

lo <- create_layout(gr, layout = "partition", circular = TRUE)
lo <- left_join(lo,
                CellTable[,c("Lineage", "dataset", "level")],
                by = c("name" = "Lineage"))
rownames(lo) <- lo$name

size_scale = c(2, 1, 0.5, 0.25, 0)
names(size_scale) = c(1, 2, 3, 4, 5)

lo[lo$name == "P0",]$level <- 1
lo[lo$name == "P0",]$dataset <- "No data"

lo$level <- as.numeric(lo$level)

lo$size_scale <- size_scale[lo$level]

pdf(paste0(dir, "Plots/data_figure/circle_annotations.pdf"), height = 5, width = 7.5)
ggraph(lo) +
  geom_node_arc_bar(aes(#filter = level > 3,
                        fill = dataset), size = 0.25) +
  coord_fixed() +
  geom_node_text(aes(filter = level < 4,
                     label = name,
                     size = size_scale[level],
                     angle = node_angle(x, y) + 90), colour = 'white') +
  annotate(geom = "text", x = 0, y = -1.2, hjust = 0.5, size = 4,
           label = "P0", color = "white") +
  scale_fill_manual(limits = c('Annotated', 'C. elegans', 'C. briggsae', 'No data',  'Cell death'),
                      values = c("#E69F00", "#009E73", "#56B4E9", "grey80",  "grey10")) +
  guides(alpha = F, size = F) +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        plot.margin = margin(t = 0, r = -100, b = 0, l = -100, unit = "mm"),
        legend.text = element_text(size = 12))
dev.off()

## test
gr <- graph_from_data_frame(CellTable[CellTable$level != "#N/A" &
                                        ! CellTable$level %in% c("1", "2"), c("Parent", "Lineage", "level", "dataset")],
                            directed = TRUE)
cells_ordered <- sort(V(gr)$name)
cells_ordered <- cells_ordered[! cells_ordered %in% c("E", "MS")]
cells_ordered <- c("MS", "E", cells_ordered)
cells_ordered <- cells_ordered[! is.na(cells_ordered)]

gr <- permute(gr, match(V(gr)$name, rev(cells_ordered)))

lo <- create_layout(gr, layout = "partition", circular = TRUE)
lo <- left_join(lo,
                CellTable[,c("Lineage", "dataset", "level")],
                by = c("name" = "Lineage"))
rownames(lo) <- lo$name

size_scale = c(2, 1, 0.5, 0.25, 0)
names(size_scale) = c(1, 2, 3, 4, 5)

lo$level <- as.numeric(lo$level)

lo$size_scale <- size_scale[lo$level]

pdf(paste0(dir, "Plots/data_figure/circle_annotations.pdf"), height = 5, width = 7.5)
ggraph(lo) +
  geom_node_arc_bar(aes(filter = level > 2,
    fill = dataset), size = 0.25) +
  coord_fixed() +
  geom_node_text(aes(filter = level < 4,
                     label = name,
                     size = size_scale[level],
                     angle = node_angle(x, y) + 90), colour = 'white') +
  scale_fill_manual(limits = c('Annotated', 'C. elegans', 'C. briggsae', 'No data',  'Cell death'),
                    values = c("#E69F00", "#009E73", "#56B4E9", "grey80",  "grey10")) +
  guides(alpha = F, size = F) +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        plot.margin = margin(t = 0, r = -100, b = 0, l = -100, unit = "mm"),
        legend.text = element_text(size = 12))
dev.off()

## lines

lo_lines <- create_layout(gr, layout = "dendrogram", circular = TRUE)
lo <- left_join(lo,
                CellTable[,c("Lineage", "dataset", "level")],
                by = c("name" = "Lineage"))
rownames(lo) <- lo$name

pdf(paste0(dir, "Plots/data_figure/circle_annotations_lines.pdf"), height = 5, width = 7.5)
ggraph(lo_lines, 'dendogram', circular = FALSE) +
  # geom_node_point() +
  geom_edge_elbow() +
  coord_fixed() +
  # geom_node_text(aes(filter = level < 4,
  #                    label = name,
  #                    size = size_scale[level],
  #                    angle = node_angle(x, y) + 90), colour = 'black') +
  # annotate(geom = "text", x = 0.025, y = 0.175, hjust = 0.5, size = 7,
  #          label = expression(paste(italic("C. elegans"), " TPM"))) +
  # annotate(geom = "text", x = 0, y = -1.2, hjust = 0.5, size = 4,
  #          label = "P0", color = "white") +
  # scale_fill_manual(limits = c('Annotated', 'C. elegans', 'C. briggsae', 'No data',  'Cell death'),
  #                   values = c("#E69F00", "#009E73", "#56B4E9", "grey80",  "grey10")) +
  guides(alpha = F, size = F) +
  theme(legend.position = 'right',
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        plot.margin = margin(t = 0, r = -100, b = 0, l = -100, unit = "mm"),
        legend.text = element_text(size = 12))
dev.off()


# cell table
CellTableList[[cur_species]] <- left_join(CellTable, tpm_matrix_list[[cur_species]][,c("Lineage", gsub("-", ".", gene))], by = c("Lineage" = "Lineage"))
rownames(CellTableList[[cur_species]]) <- CellTableList[[cur_species]]$Lineage
colnames(CellTableList[[cur_species]]) <- c(colnames(CellTableList[[cur_species]])[-length(colnames(CellTableList[[cur_species]]))], "tpm")

gr_list[[cur_species]] <- graph_from_data_frame(CellTableList[[cur_species]][! CellTableList[[cur_species]]$level %in% c("0", "1", "2", "3"),
                                                                             c("Parent", "Lineage", "tpm", "level")],
                                                directed = TRUE)

## Need to rearrange to have P2 be posterior and EMS anterior
cells_ordered <- sort(V(gr_list[[cur_species]])$name)
cells_ordered <- cells_ordered[! cells_ordered %in% c("E", "MS")]
cells_ordered <- c("MS", "E", cells_ordered)
gr_list[[cur_species]] <- permute(gr_list[[cur_species]], match(V(gr_list[[cur_species]])$name, rev(cells_ordered)))

layout_list[[cur_species]] <- create_layout(gr_list[[cur_species]], layout = "partition", circular = TRUE)
layout_list[[cur_species]] <- left_join(layout_list[[cur_species]],
                                        CellTableList[[cur_species]][,c("Lineage", "tpm", "level")],
                                        by = c("name" = "Lineage"))

layout_list[[cur_species]]$depth <- layout_list[[cur_species]]$depth - 1
layout_list[[cur_species]]$level <- as.numeric(layout_list[[cur_species]]$level)
rownames(layout_list[[cur_species]]) <- layout_list[[cur_species]]$name
# layout_list[[cur_species]][layout_list[[cur_species]]$name == "P0",]$level <- 1

tpm_max <- max(max(layout_list[["C.elegans"]]$tpm, na.rm = TRUE),
               max(layout_list[["C.briggsae"]]$tpm, na.rm = TRUE))

layout_list[["C.elegans"]][1,]$name <- ""



ggraph(layout_list[["C.elegans"]]) +
  geom_node_arc_bar(aes(fill = tpm + 1)) +
  coord_fixed() +
  geom_node_text(aes(filter = level < 6,
                     label = name,
                     size = size_scale[level],
                     angle = node_angle(x, y) + 90), colour = 'white') +
  annotate(geom = "text", x = 0.025, y = 0.175, hjust = 0.5, size = 7,
           label = expression(paste(italic("C. elegans"), " TPM"))) +
  annotate(geom = "text", x = 0, y = -1.2, hjust = 0.5, size = 8,
           label = "P0", color = "white") +
  scale_fill_viridis(trans = "log2", option = "magma",
                     limits = c(8, tpm_max), oob = scales::squish) +
  guides(alpha = F, size = F, fill = guide_colorbar(title = "")) +
  theme(legend.position = c(0.49,0.475),
        panel.background = element_rect(fill = 'transparent'),
        plot.background = element_rect(fill = 'transparent', color = NA),
        plot.margin = margin(t = 0, r = -100, b = 0, l = -100, unit = "mm"),
        legend.direction = "horizontal",
        legend.justification = "center",
        legend.key.size = unit(8, 'mm'),
        legend.text = element_text(size = 12))
