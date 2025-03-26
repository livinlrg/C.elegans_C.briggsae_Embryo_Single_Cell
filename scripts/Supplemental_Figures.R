library(ggplot2)
library(ggExtra)
library(ggrepel)
library(dplyr)
library(Seurat)
library(monocle3)
library(ggraph)
library(viridis)
library(VisCello)
library(Palo)
library(ggpp)

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240903.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]
rm(cds)

synteny <- readRDS(paste0(dir, "Objects/synteny_filt.rds"))
gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))
rownames(gff_list_mrna[["elegans"]]) <- ifelse(is.na(gff_list_mrna[["elegans"]]$cds_gene_name), gff_list_mrna[["elegans"]]$gene_name, gff_list_mrna[["elegans"]]$cds_gene_name)
rownames(gff_list_mrna[["briggsae"]]) <- ifelse(is.na(gff_list_mrna[["briggsae"]]$cds_gene_name), gff_list_mrna[["briggsae"]]$gene_name, gff_list_mrna[["briggsae"]]$cds_gene_name)

cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))
CellTable <- readRDS(paste0(dir, "Objects/CellTable_Names_20240901.rds"))

TPMLineageList <- readRDS(paste0(dir, "Objects/TPMLineageList.rds"))

TPMListBootstrapMean_term <- readRDS(paste0(dir, "Objects/TPMListBootstrapMean_CellCorrection.rds"))
TPMListBootstrap_pro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))

TPMJoint <- list()
TPMJoint[["C.elegans"]] <- cbind(TPMListBootstrap_pro[["C.elegans"]][,sort(colnames(TPMListBootstrap_pro[["C.elegans"]]))], TPMListBootstrapMean_term[["C.elegans"]])
TPMJoint[["C.briggsae"]] <- cbind(TPMListBootstrap_pro[["C.briggsae"]][,sort(colnames(TPMListBootstrap_pro[["C.elegans"]]))], TPMListBootstrapMean_term[["C.briggsae"]])

OldPath <- "/kimdata/livinlrg/scAnalysis/BobDataComb/"
og_data <- readRDS(paste0(OldPath, "Objects/WS290/og_data.rds"))

###################################
# Cell bleach time by embryo time #
###################################
CellCountByBin <- data.frame(data.frame(colData(cds_filt)) %>% group_by(species, embryo.time.bin) %>% summarise(count = n()))

##
p1 <- ggplot(CellCountByBin[! is.na(CellCountByBin$embryo.time.bin),]) +
  geom_col(aes(x = factor(embryo.time.bin, levels = c("lt_100", "100_130", "130_170",
                                                      "170_210", "210_270", "270_330",
                                                      "330_390", "390_450", "450_510",
                                                      "510_580", "580_650", "650_710", "gt_710")),
               y = count, fill = species)) +
  facet_grid(cols = vars(species)) +
  #geom_vline(xintercept = Inf, color = "grey80", size = 2) + 
  scale_x_discrete(name = "Embryo Time Bin (Minutes)",
                   labels = c("< 100", "100 to 130", "130 to 170", "170 to 210",
                              "210 to 270", "270 to 330", "330 to 390", "390 to 450",
                              "450 to 510", "510 to 580", "580 to 650",
                              "650 to 710", "> 710")) +
  scale_y_continuous(name = "Cell count",
                     breaks = c(0, 10000, 20000, 30000, 40000, 50000),
                     labels = c("0", "10,000", "20,000", "30,000", "40,000", "50,000")) +
  scale_fill_manual(name = "Species",
                    labels = c("C. briggsae", "C. elegans"),
                    values = c("#56B4E9", "#009E73")) + 
  theme(legend.position = "none",
        legend.text = (element_text(face = "italic")),
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_line(color = "grey80"),
        axis.ticks.x = element_blank(),
        strip.text = element_text(face = "italic"),
        strip.background = element_rect(color = "grey80", fill = "white", size=1.5, linetype = "solid"),
        plot.margin = margin(5.5, 5.5, 0, 5.5, "pt"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color="grey80", size = 0.1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))

pdf(paste0(dir, "Plots/supplemental_plots/EmbryoTimeDensityPlot_Bins.pdf"), height = 2.5, width = 7.5)
print(p1)
dev.off()

# Bleach timing
colData(cds_filt)$bleach_time <- NA

colData(cds_filt)[colData(cds_filt)$dataset3 == "batch_120",]$bleach_time <- "410 + 120 min." # 120 + 410
colData(cds_filt)[colData(cds_filt)$dataset3 == "batch_180",]$bleach_time <- "180 + 0 min." # + 180
colData(cds_filt)[colData(cds_filt)$dataset3 == "batch_240",]$bleach_time <- "410 + 240 min." # 240 + 410
colData(cds_filt)[colData(cds_filt)$dataset3 == "batch_300",]$bleach_time <- "0 + 300 min." # good
colData(cds_filt)[colData(cds_filt)$dataset3 == "batch_360",]$bleach_time <- "180 + 180 min." # 180 + 180
colData(cds_filt)[colData(cds_filt)$dataset3 == "batch_400",]$bleach_time <- "0 + 400 min." # good
colData(cds_filt)[colData(cds_filt)$dataset3 == "batch_500",]$bleach_time <- "0 + 500 min." # good
colData(cds_filt)[colData(cds_filt)$dataset3 == "Ce_ceh9_300_minutes",]$bleach_time <- "0 + 300 min."
colData(cds_filt)[colData(cds_filt)$dataset3 == "Ce_ceh9_600_minutes",]$bleach_time <- "0 + 600 min."
colData(cds_filt)[colData(cds_filt)$dataset3 == "Ce_M03D44_300_minutes",]$bleach_time <- "0 + 300 min."
colData(cds_filt)[colData(cds_filt)$dataset3 == "Ce_M03D44_500_minutes",]$bleach_time <- "0 + 500 min."
colData(cds_filt)[colData(cds_filt)$dataset3 == "Ce_mec3_300_minutes",]$bleach_time <- "0 + 300 min."
colData(cds_filt)[colData(cds_filt)$dataset3 == "Ce_mec3_600_minutes",]$bleach_time <- "0 + 600 min."
colData(cds_filt)[colData(cds_filt)$dataset3 == "Ce_wt_600_minutes",]$bleach_time <- "0 + 600 min."
colData(cds_filt)[colData(cds_filt)$dataset3 == "Waterston_300_minutes",]$bleach_time <- "0 + 300 min."
colData(cds_filt)[colData(cds_filt)$dataset3 == "Waterston_400_minutes",]$bleach_time <- "0 + 400 min."
colData(cds_filt)[colData(cds_filt)$dataset3 == "Waterston_500_1_minutes",]$bleach_time <- "0 + 500 min."
colData(cds_filt)[colData(cds_filt)$dataset3 == "Waterston_500_2_minutes",]$bleach_time <- "0 + 500 min."

colData(cds_filt)[colData(cds_filt)$dataset3 == "Murray_b01",]$bleach_time <- "Asynchronous"
colData(cds_filt)[colData(cds_filt)$dataset3 == "Murray_b02",]$bleach_time <- "Asynchronous"
colData(cds_filt)[colData(cds_filt)$dataset3 == "Murray_r17",]$bleach_time <- "Asynchronous"

BleachTime <- data.frame(colData(cds_filt)) %>%
  group_by(bleach_time, embryo.time.bin, species) %>%
  summarise(n = n()) %>%
  BleachTime %>% group_by(species) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup() %>%
  tidyr::complete(bleach_time, embryo.time.bin, species)

BleachTime$species <- factor(BleachTime$species, levels = c("C.briggsae", "C.elegans"),
                             labels = c("C. briggsae", "C. elegans"))

p2 <- ggplot() +
  geom_tile(data = data.frame(BleachTime)[which(BleachTime$bleach_time != "NA" &
                                                  BleachTime$embryo.time.bin != "NA"),],
            aes(y = factor(bleach_time, levels = rev(c("Asynchronous", "180 + 0 min.", "410 + 120 min.",  
                                                       "180 + 180 min." ,"410 + 240 min.", "0 + 300 min.", 
                                                       "0 + 400 min.", "0 + 500 min.", "0 + 600 min."))),
                x = factor(embryo.time.bin, levels = c("lt_100", "100_130", "130_170",
                                                       "170_210", "210_270", "270_330",
                                                       "330_390", "390_450", "450_510",
                                                       "510_580", "580_650", "650_710", "gt_710")),
                fill = freq), color = "grey80") +
  facet_grid(cols = vars(species)) +
  #geom_vline(xintercept = Inf, color = "grey80", size = 2) + 
  scale_y_discrete(name = "Collection time\n(Worm + Embryo ageing)") +
  scale_x_discrete(name = "Embryo Time Bin (Minutes)",
                   labels = c("< 100", "100 to 130", "130 to 170", "170 to 210",
                              "210 to 270", "270 to 330", "330 to 390", "390 to 450",
                              "450 to 510", "510 to 580", "580 to 650",
                              "650 to 710", "> 710")) +
  scale_fill_gradientn(name = "Percent",
                       na.value = "#F0F0F0",
                       colors = viridis(4),
                       #values = scales::rescale(c(0.05, 0.1, 0.2, 0.3, 0.4)), 
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
        plot.margin = margin(0, 5.5, 5.5, 5.5, "pt"),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.position = "bottom",
        legend.text = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))

pdf(paste0(dir, "Plots/supplemental_plots/EmbryoTimeTile_Bleach.pdf"), height = 2.5, width = 6)
print(p2)
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/EmbryoTimeBins_Facet.pdf"), height = 5, width = 7.5)
plot_grid(p1, p2, ncol = 1,  align = 'v', axis = 'lr', rel_heights = c(1, 1.5)) 
dev.off()


################################
# Plot orthology relationships #
################################

# first collect all of the 1:1's
# then label them according to the OG information about their status
ortho_summary <- rbind(data.frame(species = "C.elegans",
                                 id = rowData(cds_filt)[rowData(cds_filt)$gene.type %in% c("common", "ce.unique"),c("elegans_id")],
                                 gene_short_name = rowData(cds_filt)[rowData(cds_filt)$gene.type %in% c("common", "ce.unique"),c("elegans_gene_short_name")],
                                 gene.type = rowData(cds_filt)[rowData(cds_filt)$gene.type %in% c("common", "ce.unique"),c("gene.type")],
                                 orthology_conf = rowData(cds_filt)[rowData(cds_filt)$gene.type %in% c("common", "ce.unique"),c("orthology_conf")]),
                      data.frame(species = "C.briggsae",
                                 id = rowData(cds_filt)[rowData(cds_filt)$gene.type %in% c("common", "cb.unique"),c("briggsae_id")],
                                 gene_short_name = rowData(cds_filt)[rowData(cds_filt)$gene.type %in% c("common", "cb.unique"),c("briggsae_gene_short_name")],
                                 gene.type = rowData(cds_filt)[rowData(cds_filt)$gene.type %in% c("common", "cb.unique"),c("gene.type")],
                                 orthology_conf = rowData(cds_filt)[rowData(cds_filt)$gene.type %in% c("common", "cb.unique"),c("orthology_conf")]))

ortho_summary$cel_OG_count <- NA
ortho_summary$cbr_OG_count <- NA

for(i in 1:nrow(ortho_summary)) {
  cel_cur_gene = ortho_summary[i,"id"]
  cbr_cur_gene = ortho_summary[i,"id"]
  
  if(cel_cur_gene %in% gff_list_mrna[["elegans"]]$gene_name) {
    ortho_summary[i, "cel_OG_count"] <- gff_list_mrna[["elegans"]][which(gff_list_mrna[["elegans"]]$gene_name == cel_cur_gene), "cel_OG_count"]
    ortho_summary[i, "cbr_OG_count"] <- gff_list_mrna[["elegans"]][which(gff_list_mrna[["elegans"]]$gene_name == cel_cur_gene), "cbr_OG_count"]
  } else if (cbr_cur_gene %in% gff_list_mrna[["briggsae"]]$gene_name) {
    ortho_summary[i, "cel_OG_count"] <- gff_list_mrna[["briggsae"]][which(gff_list_mrna[["briggsae"]]$gene_name == cbr_cur_gene), "cel_OG_count"]
    ortho_summary[i, "cbr_OG_count"] <- gff_list_mrna[["briggsae"]][which(gff_list_mrna[["briggsae"]]$gene_name == cbr_cur_gene), "cbr_OG_count"]
  } else {
    print("Not in")
  }
}

ortho_summary$orthology <- ifelse(ortho_summary$cel_OG_count == 1 & ortho_summary$cbr_OG_count == 1, "1:1", NA)
ortho_summary$orthology <- ifelse(ortho_summary$cel_OG_count == 1 & ortho_summary$cbr_OG_count == 2, "1:2", ortho_summary$orthology)
ortho_summary$orthology <- ifelse(ortho_summary$cel_OG_count == 2 & ortho_summary$cbr_OG_count == 1, "2:1", ortho_summary$orthology)
ortho_summary$orthology <- ifelse(ortho_summary$cel_OG_count > 2  & ortho_summary$cbr_OG_count == 1, "many:1", ortho_summary$orthology)
ortho_summary$orthology <- ifelse(ortho_summary$cel_OG_count == 1 & ortho_summary$cbr_OG_count > 2, "1:many", ortho_summary$orthology)
ortho_summary$orthology <- ifelse(ortho_summary$cel_OG_count >= 2 & ortho_summary$cbr_OG_count >= 2, "many:many", ortho_summary$orthology)
ortho_summary$orthology <- ifelse(ortho_summary$cel_OG_count == 0 & ortho_summary$cbr_OG_count > 0, "private", ortho_summary$orthology)
ortho_summary$orthology <- ifelse(ortho_summary$cel_OG_count > 0 & ortho_summary$cbr_OG_count == 0, "private", ortho_summary$orthology)

ortho_summary[is.na(ortho_summary$orthology),]$orthology_conf <- "one:one"
ortho_summary[is.na(ortho_summary$orthology),]$orthology <- "1:1"

ortho_summary_df <- data.frame(ortho_summary %>% group_by(orthology, species, gene.type) %>% summarise(count = n()) %>% ungroup() %>% arrange(orthology, species, gene.type))

for(ortho in unique(ortho_summary_df$orthology)) {
  for(species in unique(ortho_summary_df$species)) {
    if(ortho != "1:1") {
      ortho_summary_df[ortho_summary_df$species == species &
                         ortho_summary_df$orthology == ortho &
                         ortho_summary_df$gene.type != "common",]$count <- sum(ortho_summary_df[ortho_summary_df$species == species &
                                                                                              ortho_summary_df$orthology == ortho,]$count)
      ortho_summary_df[ortho_summary_df$species == species &
                         ortho_summary_df$orthology == ortho &
                         ortho_summary_df$gene.type == "common",]$count <- ortho_summary_df[ortho_summary_df$species == species &
                                                                                          ortho_summary_df$orthology == ortho &
                                                                                          ortho_summary_df$gene.type != "common",]$count - ortho_summary_df[ortho_summary_df$species == species &
                                                                                                                                                          ortho_summary_df$orthology == ortho &
                                                                                                                                                          ortho_summary_df$gene.type == "common",]$count
    } else {
      remove_this_amount = ortho_summary_df[ortho_summary_df$species == species &
                                              ortho_summary_df$orthology == ortho &
                                              ortho_summary_df$gene.type != "common",]$count
      
      ortho_summary_df[ortho_summary_df$species == species &
                         ortho_summary_df$orthology == ortho &
                         ortho_summary_df$gene.type != "common",]$count <- sum(ortho_summary_df[ortho_summary_df$species == species &
                                                                                                  ortho_summary_df$orthology == ortho,]$count)
      ortho_summary_df[ortho_summary_df$species == species &
                         ortho_summary_df$orthology == ortho &
                         ortho_summary_df$gene.type == "common",]$count <- ortho_summary_df[ortho_summary_df$species == species &
                                                                                              ortho_summary_df$orthology == ortho &
                                                                                              ortho_summary_df$gene.type != "common",]$count - remove_this_amount
    }
  }
}

OrthoPlot <- ggplot() +
  geom_bar(data = ortho_summary_df[ortho_summary_df$gene.type != "common",],
           aes(y = count, x = factor(orthology), group = species),
           color = "black", fill = "orange", stat = "identity", position = position_dodge(), size = 0.35) +
  geom_bar(data = ortho_summary_df[ortho_summary_df$gene.type == "common",],
           aes(y = count, x = factor(orthology), fill = species, group = species),
           color = "black", stat = "identity", position = position_dodge(), size = 0.35) +
  scale_x_discrete(name = "Orthology",
                   limits = c("1:1", "1:2", "2:1", "1:many", "many:1", "many:many", "private"),
                   labels = c("1 to 1", "1 to 2", "2 to 1", "1 to many", "many to 1", "many to many", "Private")) +
  scale_y_continuous(name = "Number of Genes") +
  scale_fill_manual(name = "Species", breaks = c("C.briggsae", "C.elegans"),
                    labels = c("C. briggsae", "C. elegans"),
                    values = c("#56B4E9", "#009E73")) +
  theme(legend.position = "right",
        legend.text = (element_text(face = "italic")),
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(color="grey80", size = 0.1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pdf(paste0(dir, "Plots/supplemental_plots/Gene_orthology.pdf"), height = 3, width = 5)
print(OrthoPlot)
dev.off()

##################
# UMAP Time Plot #
##################

embeddings <- readRDS(paste0(dir, "Objects/UMAP_rpca_coordinates_global.rds"))

## Plot UMAP
projection <- data.frame(X = embeddings[,1],
                         Y = embeddings[,2],
                         species = colData(cds_filt)[rownames(embeddings),]$species,
                         embryo_time_bin = colData(cds_filt)[rownames(embeddings),]$embryo.time.bin,
                         batch = colData(cds_filt)[rownames(embeddings),]$dataset3,
                         bleach_time = colData(cds_filt)[rownames(embeddings),]$bleach_time)

projection <- projection[! is.na(projection$embryo_time_bin),]

projection$dataset <- NA
projection[which(projection$batch == "Ce_ceh9_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Ce_ceh9_600_minutes"),]$dataset <- "600 min."
projection[which(projection$batch == "Ce_M03D44_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Ce_M03D44_500_minutes"),]$dataset <- "500 min."
projection[which(projection$batch == "Ce_mec3_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Ce_mec3_600_minutes"),]$dataset <- "600 min."
projection[which(projection$batch == "Ce_wt_600_minutes"),]$dataset <- "600 min."
projection[which(projection$batch == "Waterston_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Waterston_400_minutes"),]$dataset <- "400 min."
projection[which(projection$batch == "Waterston_500_1_minutes"),]$dataset <- "500 min."
projection[which(projection$batch == "Waterston_500_2_minutes"),]$dataset <- "500 min."
projection[which(projection$batch == "Murray_b01"),]$dataset <- "0 min."
projection[which(projection$batch == "Murray_b02"),]$dataset <- "0 min."
projection[which(projection$batch == "Murray_r17"),]$dataset <- "0 min."

embryo_time_plot <- ggplot(projection, aes(x = X,
                       y = Y)) +
  geom_point(data = projection, aes(x = X,
                                    y = Y,
                                    color = factor(embryo_time_bin, levels = c("lt_100", "100_130", "130_170",
                                                                               "170_210", "210_270", "270_330",
                                                                               "330_390", "390_450", "450_510",
                                                                               "510_580", "580_650", "650_710", "gt_710"))),
             size = 0.12,
             stroke = 0.04,
             alpha = 0.8) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = -Inf,
           yend = min(projection$Y) + (abs(min(projection$Y)) + abs(max(projection$Y)))/4,
           color = "grey80", size = 2) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = min(projection$X) + (abs(min(projection$X)) + abs(max(projection$X)))/4,
           yend = -Inf,
           color = "grey80", size = 2) +
  scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(projection$X), max(projection$X), by = 0.20),1)) +
  scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(projection$Y), max(projection$Y), by = 0.20),1)) +
  scale_colour_viridis_d(name = "Estimated\nembryo time\n(minutes)",
                         labels = c("< 100", "100 to 130", "130 to 170", "170 to 210",
                                    "210 to 270", "270 to 330", "330 to 390", "390 to 450",
                                    "450 to 510", "510 to 580", "580 to 650",
                                    "650 to 710", "> 710")) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme(rect = element_rect(fill = "transparent"),
        axis.line = element_blank(),
        # axis.line = element_line(color="grey80", size=1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = "right",
        axis.title.x = element_text(hjust=0.1),
        axis.title.y = element_text(hjust=0.1),
        legend.background = element_blank(),
        legend.key = element_blank())

pdf(paste0(dir, "Plots/supplemental_plots/SmoothedEmbryoTime_UMAP.pdf"),  width = 10, height = 10)
print(embryo_time_plot)
dev.off()

#####################
# UMAP Species Plot #
#####################
species_plot <- ggplot(projection, aes(x = X,
                       y = Y,
                       color = species)) +
  geom_point(size = 0.12,
             alpha = 0.8, stroke = 0.04) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = -Inf,
           yend = min(projection$X) + (abs(min(projection$Y)) + abs(max(projection$Y)))/4, color = "grey80", size = 2) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = min(projection$X) + (abs(min(projection$X)) + abs(max(projection$X)))/4,
           yend = -Inf, color = "grey80", size = 2) +
  scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(projection$X), max(projection$X), by = 0.20),1)) +
  scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(projection$Y), max(projection$Y), by = 0.20),1)) +
  # facet_grid(cols = vars(species), scales = "free") +
  scale_color_manual(name = "Species", values = c("#56B4E9", "#009E73")) + 
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme(rect = element_rect(fill = "transparent"),
        plot.title = element_text(face = "italic", size = 32),
        axis.line = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(hjust=0.1),
        axis.title.y = element_text(hjust=0.1),
        legend.position = "right",
        # strip.background = element_blank(),
        # strip.text.x = element_text(angle = 0, hjust = 0, size = 20, face = "italic"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # text = element_text(family = "sans", color = "black", size=12),
        # panel.grid.major = element_line("grey80", size = 0.25, linetype = "solid"),
        # panel.grid.minor = element_line("grey80", size = 0.05, linetype = "solid"),
        # axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        legend.background = element_blank(),
        legend.key = element_blank())

pdf(paste0(dir, "Plots/supplemental_plots/SpeciesPlot.pdf"), width = 10, height = 10)
print(species_plot)
dev.off()

###################
# UMAP Batch Plot #
###################

plan("multicore", workers = 6)
options(future.globals.maxSize = 30000 * 1024^2)

cel_so <- CreateSeuratObject(counts = counts(cds)[rowData(cds)[rowData(cds)$gene.type %in% c("common", "ce.unique"),]$gene_short_name,
                                                     ! colData(cds)$filter %in% c("doublet", "mutant", "damaged", "stressed") &
                                                    colData(cds)$species == "C.elegans"],
                                project = "Nov_2024",
                                min.features = 0,
                                min.cells = 0)

cel_so <- AddMetaData(cel_so, data.frame(colData(cds[rowData(cds)[rowData(cds)$gene.type %in% c("common", "ce.unique"),]$gene_short_name,
                                                     ! colData(cds)$filter %in% c("doublet", "mutant", "damaged", "stressed") &
                                                       colData(cds)$species == "C.elegans"])))
cel_so.list <- split(cel_so, f = cel_so$dataset3)

cel_so.list <- NormalizeData(cel_so.list)
cel_so.list <- FindVariableFeatures(cel_so.list)
cel_so.list <- ScaleData(cel_so.list)
cel_so.list <- RunPCA(cel_so.list, npcs = 300)

saveRDS(cel_so.list, paste0(dir, "Objects/seurat_objects/cel_seurat_pre_integration_20241117.rds"))
cel_so.list <- readRDS(paste0(dir, "Objects/seurat_objects/cel_seurat_pre_integration_20241117.rds"))

# mainstream integration
cel_so.list <- IntegrateLayers(object = cel_so.list,
                               method = RPCAIntegration,
                               orig.reduction = "pca",
                               new.reduction = "integrated.rpca",
                               dims = 1:300,
                               verbose = TRUE)

cel_so.list <- IntegrateLayers(object = cel_so.list,
                               method = CCAIntegration,
                               orig.reduction = "pca",
                               new.reduction = "integrated.cca",
                               dims = 1:300,
                               verbose = TRUE)

cel_so.list <- FindNeighbors(cel_so.list, reduction = "integrated.cca", dims = 1:300)
cel_so.list <- FindClusters(cel_so.list, resolution = 2, cluster.name = "cca_clusters")
cel_so.list <- RunUMAP(cel_so.list, reduction = "integrated.cca", dims = 1:300, reduction.name = "umap.cca", seed.use = 20240826)
cel_so.list <- RunUMAP(cel_so.list, dims = 1:300, reduction.name = "umap", seed.use = 20240826)

cel_so.rpca <- JoinLayers(cel_so.list)

embeddings <- Embeddings(cel_so.rpca, reduction = "umap")

projection <- data.frame(X = embeddings[,1],
                         Y = embeddings[,2],
                         species = colData(cds_filt)[rownames(embeddings),]$species,
                         embryo_time_bin = colData(cds_filt)[rownames(embeddings),]$embryo.time.bin,
                         batch = colData(cds_filt)[rownames(embeddings),]$dataset3,
                         bleach_time = colData(cds_filt)[rownames(embeddings),]$bleach_time)

projection$dataset <- NA
projection[which(projection$batch == "Ce_ceh9_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Ce_ceh9_600_minutes"),]$dataset <- "600 min."
projection[which(projection$batch == "Ce_M03D44_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Ce_M03D44_500_minutes"),]$dataset <- "500 min."
projection[which(projection$batch == "Ce_mec3_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Ce_mec3_600_minutes"),]$dataset <- "600 min."
projection[which(projection$batch == "Ce_wt_600_minutes"),]$dataset <- "600 min."
projection[which(projection$batch == "Waterston_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Waterston_400_minutes"),]$dataset <- "400 min."
projection[which(projection$batch == "Waterston_500_1_minutes"),]$dataset <- "500 min."
projection[which(projection$batch == "Waterston_500_2_minutes"),]$dataset <- "500 min."
projection[which(projection$batch == "Murray_b01"),]$dataset <- "0 min."
projection[which(projection$batch == "Murray_b02"),]$dataset <- "0 min."
projection[which(projection$batch == "Murray_r17"),]$dataset <- "0 min."

cel_batch_plot <- ggplot(projection[projection$species == "C.elegans",], aes(x = X,
                                       y = Y,
                                       color = batch)) +
  geom_point(size = 0.12,
             alpha = 0.6, stroke = 0.04) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = -Inf,
           yend = min(projection$X) + (abs(min(projection$Y)) + abs(max(projection$Y)))/4, color = "grey80", size = 2) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = min(projection$X) + (abs(min(projection$X)) + abs(max(projection$X)))/4,
           yend = -Inf, color = "grey80", size = 2) +
  scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(projection$X), max(projection$X), by = 0.20),1)) +
  scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(projection$Y), max(projection$Y), by = 0.20),1)) +
  # scale_color_manual(name = "Collection Time\n(Embryo ageing)",
  #                      breaks = c("0 min.", "300 min.", "400 min.", "500 min.", "600 min."),
  #                      values = rev(viridis(30)[seq(1, 25, by = 6)])) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme(rect = element_rect(fill = "transparent"),
        plot.title = element_text(face = "italic", size = 32),
        axis.line = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(hjust=0.1),
        axis.title.y = element_text(hjust=0.1),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank())

pdf(paste0(dir, "Plots/supplemental_plots/CelBatchPlot.pdf"), width = 10, height = 10)
print(cel_batch_plot)
dev.off()

cbr_so <- CreateSeuratObject(counts = counts(cds)[rowData(cds)[rowData(cds)$gene.type %in% c("common", "ce.unique"),]$gene_short_name,
                                                  ! colData(cds)$filter %in% c("doublet", "mutant", "damaged", "stressed") &
                                                    colData(cds)$species == "C.elegans"],
                             project = "Nov_2024",
                             min.features = 0,
                             min.cells = 0)

cbr_so <- AddMetaData(cbr_so, data.frame(colData(cds[rowData(cds)[rowData(cds)$gene.type %in% c("common", "ce.unique"),]$gene_short_name,
                                                     ! colData(cds)$filter %in% c("doublet", "mutant", "damaged", "stressed") &
                                                       colData(cds)$species == "C.elegans"])))
cbr_so.list <- split(cbr_so, f = cbr_so$dataset3)

cbr_so.list <- NormalizeData(cbr_so.list)
cbr_so.list <- FindVariableFeatures(cbr_so.list)
cbr_so.list <- ScaleData(cbr_so.list)
cbr_so.list <- RunPCA(cbr_so.list, npcs = 300)

saveRDS(cbr_so.list, paste0(dir, "Objects/seurat_objects/cbr_seurat_pre_integration_20241117.rds"))
cbr_so.list <- readRDS(paste0(dir, "Objects/seurat_objects/cbr_seurat_pre_integration_20241117.rds"))

# mainstream integration
cbr_so.list <- IntegrateLayers(object = cbr_so.list,
                               method = RPCAIntegration,
                               orig.reduction = "pca",
                               new.reduction = "integrated.rpca",
                               dims = 1:300,
                               verbose = TRUE)

cbr_so.list <- IntegrateLayers(object = cbr_so.list,
                               method = CCAIntegration,
                               orig.reduction = "pca",
                               new.reduction = "integrated.cca",
                               dims = 1:300,
                               verbose = TRUE)

cbr_so.list <- FindNeighbors(cbr_so.list, reduction = "integrated.cca", dims = 1:300)
cbr_so.list <- FindClusters(cbr_so.list, resolution = 2, cluster.name = "cca_clusters")
cbr_so.list <- RunUMAP(cbr_so.list, reduction = "integrated.cca", dims = 1:300, reduction.name = "umap.cca", seed.use = 20240826)
cbr_so.list <- RunUMAP(cbr_so.list, dims = 1:300, reduction.name = "umap", seed.use = 20240826)

cbr_so.list <- JoinLayers(cbr_so.list)

embeddings <- Embeddings(cbr_so.list, reduction = "umap")

projection <- data.frame(X = embeddings[,1],
                         Y = embeddings[,2],
                         species = colData(cds_filt)[rownames(embeddings),]$species,
                         embryo_time_bin = colData(cds_filt)[rownames(embeddings),]$embryo.time.bin,
                         batch = colData(cds_filt)[rownames(embeddings),]$dataset3,
                         bleach_time = colData(cds_filt)[rownames(embeddings),]$bleach_time)

projection$dataset <- NA
projection[which(projection$batch == "Ce_ceh9_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Ce_ceh9_600_minutes"),]$dataset <- "600 min."
projection[which(projection$batch == "Ce_M03D44_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Ce_M03D44_500_minutes"),]$dataset <- "500 min."
projection[which(projection$batch == "Ce_mec3_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Ce_mec3_600_minutes"),]$dataset <- "600 min."
projection[which(projection$batch == "Ce_wt_600_minutes"),]$dataset <- "600 min."
projection[which(projection$batch == "Waterston_300_minutes"),]$dataset <- "300 min."
projection[which(projection$batch == "Waterston_400_minutes"),]$dataset <- "400 min."
projection[which(projection$batch == "Waterston_500_1_minutes"),]$dataset <- "500 min."
projection[which(projection$batch == "Waterston_500_2_minutes"),]$dataset <- "500 min."
projection[which(projection$batch == "Murray_b01"),]$dataset <- "0 min."
projection[which(projection$batch == "Murray_b02"),]$dataset <- "0 min."
projection[which(projection$batch == "Murray_r17"),]$dataset <- "0 min."

cbr_batch_plot <- ggplot(projection[projection$species == "C.elegans",], aes(x = X,
                                                                             y = Y,
                                                                             color = batch)) +
  geom_point(size = 0.12,
             alpha = 0.6, stroke = 0.04) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = -Inf,
           yend = min(projection$X) + (abs(min(projection$Y)) + abs(max(projection$Y)))/4, color = "grey80", size = 2) +
  annotate(geom = "segment", x = -Inf, y = -Inf,
           xend = min(projection$X) + (abs(min(projection$X)) + abs(max(projection$X)))/4,
           yend = -Inf, color = "grey80", size = 2) +
  scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(projection$X), max(projection$X), by = 0.20),1)) +
  scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(projection$Y), max(projection$Y), by = 0.20),1)) +
  # scale_color_manual(name = "Collection Time\n(Embryo ageing)",
  #                      breaks = c("0 min.", "300 min.", "400 min.", "500 min.", "600 min."),
  #                      values = rev(viridis(30)[seq(1, 25, by = 6)])) +
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme(rect = element_rect(fill = "transparent"),
        plot.title = element_text(face = "italic", size = 32),
        axis.line = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(hjust=0.1),
        axis.title.y = element_text(hjust=0.1),
        legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank())

pdf(paste0(dir, "Plots/supplemental_plots/CbrBatchPlot.pdf"), width = 10, height = 10)
print(cbr_batch_plot)
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/UMAP_batch_plots.pdf"), width = 10, height = 5)
plot_grid(cel_batch_plot, cbr_batch_plot, ncol = 2, labels = c("C. elegans", "C. briggsae"))
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/UMAP_plots.pdf"), width = 10, height = 10)
plot_grid(embryo_time_plot, species_plot,
          cel_batch_plot, cbr_batch_plot, ncol = 2)
dev.off()

#####################
# UMAP Tissue Plots #
#####################

## Read in the coordinates from the clist and the cell identities in the eset then
## plot the cell identities using the centroids of the values from the clist to label with geom_text_repel
eset <- readRDS(paste0(dir, "Objects/cello_creation/cello_internal/eset.rds"))
clist <- readRDS(paste0(dir, "Objects/cello_creation/cello_internal/clist.rds"))

plot_list <- list()

terminal_umap_names <- c("Muscle", "Hypodermis and seam", 
                         "Ciliated neurons", "Non-ciliated neurons", "Pharynx and rectal",
                         "Glia and excretory", "Intestine", "Mesoderm", "Germline")

names(terminal_umap_names) <- c("muscle", "hypodermis_and_seam",
                                "ciliated_neurons", "nonciliated_neurons",
                                "pharynx_rectal", "glia_excretory",
                                "intestine", "mesoderm", "germline")

cell_type_list <- list()

cell_type_list[["Ciliated neurons"]] <- c("ADL", "AWA", "ADF", "AFD", "ADE", "IL1",
                     "IL2", "OLL", "OLQ", "BAG", "URX", "CEP",
                     "ASJ", "ASI", "ASG", "ASK", "ASE", "ASH", "AWB", "AWC")
cell_type_list[["Non-ciliated neurons"]] <- c("AIY", "AIA", "AIB", "AIZ", "AIM", "AIN",
                        "ALA", "ALM_BDU", "ALM_PLM", "AUA", "AVA", "AVB",
                        "AVD", "AVD", "AVE", "AVF", "AVG",
                        "AVG", "AVH", "AVJ", "AVK", "AVL",
                        "DA", "DA_DB", "DB", "DD", "DVA",
                        "DVC", "PVP", "PVQ_and_possibly_PVC",
                        "PVR", "PVT",  "RIA", "RIB",
                        "RIC", "RID", "RIM", "RIS",
                        "RIV", "RMD", "RME", "SIA", "SIB", "SMB", "SMD",
                        "URB_and_URA", "CAN",
                        "FLP", "RIH", "RIG",  "RMG", "I5")
cell_type_list[["Muscle"]] <- c("BWM_anterior1", "BWM_anterior2", "BWM_far_posterior", "BWM_headrow1", "BWM_headrow1_2", 
            "BWM_headrow1_in", "BWM_headrow1_out", "BWM_headrow2", "BWM_headrow2_in", "BWM_headrow2_out", 
            "BWM_middle", "BWM_posterior_late", "BWM_posterior1", "BWM_posterior2", "mu_int_mu_anal",
            "mu_int_mu_anal_related", "mu_sph", "BWM_anterior_early", "BWM_middle_early", "BWM_posterior_early")
cell_type_list[["Mesoderm"]] <- c("Coelomocytes", "hmc", "hmc_and_hmc_homolog", "Z1_Z4",
              "hmc_homolog", "GLR_1", "GLR_2", "M_cell")
cell_type_list[["Intestine"]] <- c("Intestine_anterior",
               "Intestine_far_posterior", "Intestine_middle", "Intestine_middle_and_posterior",
               "Intestine_early", "Intestine_early_far_posterior")
cell_type_list[["Hypodermis and seam"]] <- c("T", "hyp3", 
                       "hyp4_hyp5_hyp6", "hyp7_AB_lineage", "hyp7_C_lineage",
                       "Hypodermis", "Tail_hypodermis", "Seam_cells", "P_cells",
                       "Seam_hyp_early_cells", "Seam_cells_early")
cell_type_list[["Pharynx and rectal"]] <- c("B_F_K_Kp_U_Y", "Rectal_gland",
                                        "hyp1_hyp2", "hyp1V", "hyp1V_and_ant_arc_V", "Pharyngeal_intestinal_valve",
                                        "pm1_pm2", "pm3_pm4_pm5a", "pm3_pm4_pm5b", "pm3_pm4_pm5c", "pm6",
                                        "pm7", "pm8", "mc1",
                                        "mc2a", "mc2b", "mc3", "Anterior_arcade_cell", 
                                        "Arcade_cell", "Posterior_arcade_cell", "g1A", "g1P", "g2", "MC",
                                        "early_hyp1_hyp2_cell", "early_arcade_cell")
cell_type_list[["Germline"]] <- c("germline")
cell_type_list[["Glia and excretory"]] <- c("Excretory_cell", "Excretory_duct_and_pore", "Excretory_gland", "AMsh", "AMso", "CEPsh", "ADEsh", "ILsh_OLLsh_OLQsh", "ILso", "CEPso")

Biobase::pData(eset)$barcodes <- rownames(Biobase::pData(eset))
pdata <- Biobase::pData(eset)
for(terminal_umap in terminal_umap_names) {
  print(terminal_umap)
  
  # First find the centroid of the cell types coordinates
  cell_types <- unique(pdata[rownames(data.frame(clist[[terminal_umap]]@proj)),]$cell_type)
  cell_types <- cell_types[cell_types != "unannotated" &
                             cell_types %in% cell_type_list[[terminal_umap]]]
  
  centroid_df <- data.frame(cell_type = cell_types, umap1 = NA, umap2 = NA)
  rownames(centroid_df) <- centroid_df$cell_type
  
  temp_umap <- data.frame(clist[[terminal_umap]]@proj[[1]])
  temp_umap$barcodes <- rownames(temp_umap)
  temp_umap <- left_join(temp_umap, pdata[temp_umap$barcodes,], by = c("barcodes" = "barcodes"))
  
  temp_umap[! temp_umap$cell_type %in% cell_types, "cell_type"] <- "unannotated"
  
  for(cur_cell_type in cell_types) {
    temp_pData <- pdata[pdata$cell_type %in% cur_cell_type,]
    # Find centroid for this cell type
    centroid_df[cur_cell_type,]$umap1 <- mean(temp_umap[temp_umap$barcodes %in% rownames(temp_pData),][,"UMAP_1"])
    centroid_df[cur_cell_type,]$umap2 <- mean(temp_umap[temp_umap$barcodes %in% rownames(temp_pData),][,"UMAP_2"])
  }
  
  species_plot <- temp_umap %>% arrange(sample(rownames(temp_umap))) %>%
    ggplot(data = ., aes(x = UMAP_1,
                         y = UMAP_2,
                         color = species)) +
    geom_point(size = 0.3,
               alpha = 0.8,
               stroke = 0.1) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = -Inf,
             yend = min(temp_umap$UMAP_2) + (abs(min(temp_umap$UMAP_2)) + abs(max(temp_umap$UMAP_2)))/4, color = "grey80", size = 2) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = min(temp_umap$UMAP_1) + (abs(min(temp_umap$UMAP_1)) + abs(max(temp_umap$UMAP_1)))/4,
             yend = -Inf, color = "grey80", size = 2) +
    scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(temp_umap$UMAP_1), max(temp_umap$UMAP_1), by = 0.20),1)) +
    scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(temp_umap$UMAP_2), max(temp_umap$UMAP_2), by = 0.20),1)) +
    scale_color_manual(name = "Species", values = c("#56B4E9", "#009E73")) + 
    guides(colour = guide_legend(override.aes = list(size=4))) +
    theme(rect = element_rect(fill = "transparent"),
          plot.title = element_text(face = "italic", size = 32),
          axis.line = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title.x = element_text(hjust=0.1),
          axis.title.y = element_text(hjust=0.1),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank())
  
  time_plot <- temp_umap %>% filter(! is.na(embryo.time.bin)) %>%
    ggplot(aes(x = UMAP_1,
               y = UMAP_2, 
               color = factor(embryo.time.bin, levels = c("lt_100", "100_130", "130_170",
                                                          "170_210", "210_270", "270_330",
                                                          "330_390", "390_450", "450_510",
                                                          "510_580", "580_650", "650_710", "gt_710")))) +
    geom_point(size = 0.3,
               stroke = 0.1,
               alpha = 0.6) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = -Inf,
             yend = min(temp_umap$UMAP_2) + (abs(min(temp_umap$UMAP_2)) + abs(max(temp_umap$UMAP_2)))/4,
             color = "grey80", size = 2) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = min(temp_umap$UMAP_1) + (abs(min(temp_umap$UMAP_1)) + abs(max(temp_umap$UMAP_1)))/4,
             yend = -Inf,
             color = "grey80", size = 2) +
    scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(temp_umap$UMAP_1), max(temp_umap$UMAP_1), by = 0.20),1)) +
    scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(temp_umap$UMAP_2), max(temp_umap$UMAP_2), by = 0.20),1)) +
    scale_colour_viridis_d(name = "Estimated\nembryo time\n(minutes)",
                           breaks = c("lt_100", "100_130", "130_170",
                                      "170_210", "210_270", "270_330",
                                      "330_390", "390_450", "450_510",
                                      "510_580", "580_650", "650_710", "gt_710"),
                           labels = c("< 100", "100 to 130", "130 to 170", "170 to 210",
                                      "210 to 270", "270 to 330", "330 to 390", "390 to 450",
                                      "450 to 510", "510 to 580", "580 to 650",
                                      "650 to 710", "> 710")) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    theme(rect = element_rect(fill = "transparent"),
          axis.line = element_blank(),
          # axis.line = element_line(color="grey80", size=1),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "right",
          axis.title.x = element_text(hjust=0.1),
          axis.title.y = element_text(hjust=0.1),
          # legend.key.size = unit(10, 'cm'), #change legend key size
          # legend.key.height = unit(8, 'cm'), #change legend key height
          # legend.key.width = unit(6, 'cm'), #change legend key width
          # legend.title = element_text(size=72), #change legend title font size
          # legend.text = element_text(size=60),
          legend.background = element_blank(),
          legend.key = element_blank())
  
  ## Add cell type colors
  cl <- temp_umap$cell_type
  
  gg_color_hue <- function(n) {
    hues = seq(30, 1000, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  pal <- gg_color_hue(length(unique(cl)))
  
  palopal_temp <- Palo(position = temp_umap[,c("UMAP_1", "UMAP_2")],
                       cluster = cl,
                       palette = pal,
                       rgb_weight = c(3, 4 ,2),
                       init_iter = 1000,
                       refine_iter = 2000,
                       early_stop = 500)
  
  palopal_temp["unannotated"] = "lightgrey"
  
  centroid_df$is_centroid = TRUE
  centroid_df <- rbind(centroid_df, data.frame(cell_type = temp_umap$cell_type,
                                               umap1 = temp_umap$UMAP_1,
                                               umap2 = temp_umap$UMAP_2,
                                               is_centroid = FALSE))
  centroid_df$label <- ifelse(centroid_df$is_centroid, centroid_df$cell_type, "")
  
  cell_type_plot <- centroid_df %>% 
    mutate(cell_type = ifelse(cell_type == "unannotated", NA, cell_type)) %>%
    arrange(cell_type) %>%
    ggplot(aes(x = umap1,
               y = umap2,
               color = cell_type,
               label = label)) +
    geom_point(data = centroid_df[! centroid_df$is_centroid,],
               size = 0.3,
               stroke = 0.1) +
    geom_text_repel(max.overlaps = Inf) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = -Inf,
             yend = min(temp_umap$UMAP_2) + (abs(min(temp_umap$UMAP_2)) + abs(max(temp_umap$UMAP_2)))/4, color = "grey80", size = 2) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = min(temp_umap$UMAP_1) + (abs(min(temp_umap$UMAP_1)) + abs(max(temp_umap$UMAP_1)))/4,
             yend = -Inf, color = "grey80", size = 2) +
    scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(temp_umap$UMAP_1), max(temp_umap$UMAP_1), by = 0.20),1)) +
    scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(temp_umap$UMAP_2), max(temp_umap$UMAP_2), by = 0.20),1)) +
    scale_color_manual(values = palopal_temp, na.value = "grey") +
    ggtitle(paste0(terminal_umap, " UMAP [2D]")) +
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
  
  plot_list[[terminal_umap]] <- plot_grid(cell_type_plot,
                                          plot_grid(time_plot,
                                                    species_plot, nrow = 2),
                                          ncol = 2, rel_widths = c(1, 0.5))
}

for(terminal_umap in names(terminal_umap_names)) {
  print(terminal_umap)
  pdf(paste0(dir, "Plots/supplemental_plots/tissue_umaps/", terminal_umap, ".pdf"), width = 12, height = 9)
  print(plot_list[[terminal_umap_names[terminal_umap]]])
  dev.off()
}


progenitor_umap_names <- c("0 to 150 embryo time",
                           "0 to 200 embryo time",
                           "0 to 250 embryo time",
                           "0 to 300 embryo time")

names(progenitor_umap_names) <- c("0_150",
                                  "0_200",
                                  "0_250",
                                  "0_300")

for(progenitor_umap in progenitor_umap_names) {
  print(progenitor_umap)
  
  # First find the centroid of the cell types coordinates
  cell_types <- unique(pdata[rownames(data.frame(clist[[progenitor_umap]]@proj)),]$lineage_broad)
  cell_types <- cell_types[cell_types != "unannotated"]
  
  centroid_df <- data.frame(cell_type = cell_types, umap1 = NA, umap2 = NA)
  rownames(centroid_df) <- centroid_df$cell_type
  
  temp_umap <- data.frame(clist[[progenitor_umap]]@proj[[1]])
  temp_umap$barcodes <- rownames(temp_umap)
  temp_umap <- left_join(temp_umap, pdata[temp_umap$barcodes,], by = c("barcodes" = "barcodes"))
  
  for(cur_cell_type in cell_types) {
    temp_pData <- pdata[pdata$lineage_broad %in% cur_cell_type,]
    # Find centroid for this cell type
    centroid_df[cur_cell_type,]$umap1 <- mean(temp_umap[temp_umap$barcodes %in% rownames(temp_pData),][,"UMAP_1"])
    centroid_df[cur_cell_type,]$umap2 <- mean(temp_umap[temp_umap$barcodes %in% rownames(temp_pData),][,"UMAP_2"])
  }
  
  species_plot <- temp_umap %>% arrange(sample(rownames(temp_umap))) %>%
    ggplot(data = ., aes(x = UMAP_1,
                         y = UMAP_2,
                         color = species)) +
    geom_point(size = 0.3,
               stroke = 0.1,
               alpha = 0.6) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = -Inf,
             yend = min(temp_umap$UMAP_2) + (abs(min(temp_umap$UMAP_2)) + abs(max(temp_umap$UMAP_2)))/4, color = "grey80", size = 6) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = min(temp_umap$UMAP_1) + (abs(min(temp_umap$UMAP_1)) + abs(max(temp_umap$UMAP_1)))/4,
             yend = -Inf, color = "grey80", size = 6) +
    scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(temp_umap$UMAP_1), max(temp_umap$UMAP_1), by = 0.20),1)) +
    scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(temp_umap$UMAP_2), max(temp_umap$UMAP_2), by = 0.20),1)) +
    scale_color_manual(name = "Species", values = c("#56B4E9", "#009E73")) + 
    guides(colour = guide_legend(override.aes = list(size=8))) +
    theme(rect = element_rect(fill = "transparent"),
          plot.title = element_text(face = "italic", size = 32),
          axis.line = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title.x = element_text(hjust=0.1, size = 32),
          axis.title.y = element_text(hjust=0.1, size = 32),
          legend.title = element_text(size = 32),
          legend.text = element_text(size = 32),
          legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_blank(),
          legend.key = element_blank())
  
  time_plot <- temp_umap %>% filter(! is.na(embryo.time.bin)) %>%
    ggplot(aes(x = UMAP_1,
               y = UMAP_2, 
               color = factor(embryo.time.bin, levels = c("lt_100", "100_130", "130_170",
                                                          "170_210", "210_270", "270_330",
                                                          "330_390", "390_450", "450_510",
                                                          "510_580", "580_650", "650_710", "gt_710")))) +
    geom_point(size = 0.3,
               stroke = 0.1,
               alpha = 0.6) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = -Inf,
             yend = min(temp_umap$UMAP_2) + (abs(min(temp_umap$UMAP_2)) + abs(max(temp_umap$UMAP_2)))/4,
             color = "grey80", size = 6) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = min(temp_umap$UMAP_1) + (abs(min(temp_umap$UMAP_1)) + abs(max(temp_umap$UMAP_1)))/4,
             yend = -Inf,
             color = "grey80", size = 6) +
    scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(temp_umap$UMAP_1), max(temp_umap$UMAP_1), by = 0.20),1)) +
    scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(temp_umap$UMAP_2), max(temp_umap$UMAP_2), by = 0.20),1)) +
    scale_colour_viridis_d(name = "Estimated\nembryo time\n(minutes)",
                           breaks = c("lt_100", "100_130", "130_170",
                                      "170_210", "210_270", "270_330",
                                      "330_390", "390_450", "450_510",
                                      "510_580", "580_650", "650_710", "gt_710"),
                           labels = c("< 100", "100 to 130", "130 to 170", "170 to 210",
                                      "210 to 270", "270 to 330", "330 to 390", "390 to 450",
                                      "450 to 510", "510 to 580", "580 to 650",
                                      "650 to 710", "> 710")) +
    guides(colour = guide_legend(override.aes = list(size = 8))) +
    theme(rect = element_rect(fill = "transparent"),
          axis.line = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "right",
          axis.title.x = element_text(hjust=0.1, size = 32),
          axis.title.y = element_text(hjust=0.1, size = 32),
          legend.title = element_text(size = 32),
          legend.text = element_text(size = 32),
          legend.background = element_blank(),
          legend.key = element_blank())
  
  ## Add cell type colors
  cl <- temp_umap[temp_umap$lineage_broad != "unannotated",]$lineage_broad
  
  cl_count <- table(cl)
  
  cl <- ifelse(cl %in% names(cl_count)[cl_count > 4], cl, "unannotated")
  
  gg_color_hue <- function(n) {
    hues = seq(30, 1000, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  pal <- gg_color_hue(length(unique(cl)))

  palopal_temp <- Palo(position = temp_umap[,c("UMAP_1", "UMAP_2")],
                       cluster = cl,
                       palette = pal,
                       rgb_weight = c(3, 4 ,2),
                       init_iter = 1000,
                       refine_iter = 2000,
                       early_stop = 500)

  palopal_temp[palopal_temp %in% c("unassigned", "unannotated")] <- "lightgrey"
  
  centroid_df$is_centroid = TRUE
  centroid_df <- rbind(centroid_df, data.frame(cell_type = temp_umap$lineage_broad,
                                               umap1 = temp_umap$UMAP_1,
                                               umap2 = temp_umap$UMAP_2,
                                               is_centroid = FALSE))
  centroid_df$label <- ifelse(centroid_df$is_centroid, centroid_df$cell_type, "")
  
  cell_type_plot <- centroid_df %>% 
    mutate(cell_type = ifelse(cell_type %in% c("unannotated", "unassigned"), NA, cell_type)) %>%
    mutate(cell_type = ifelse(cell_type %in% cl, cell_type, NA)) %>%
    mutate(label = ifelse(label %in% c("unannotated", "unassigned"), NA, label)) %>%
    mutate(label = ifelse(label %in% cl, label, NA)) %>%
    arrange(cell_type) %>%
    ggplot(aes(x = umap1,
               y = umap2,
               color = cell_type,
               label = label)) +
    geom_point(data = centroid_df[! centroid_df$is_centroid,],
               size = 0.3,
               stroke = 0.1) +
    geom_text_repel(max.overlaps = Inf,
                    force = 10,
                    max.time = 50, max.iter = 1000000) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = -Inf,
             yend = min(temp_umap$UMAP_2) + (abs(min(temp_umap$UMAP_2)) + abs(max(temp_umap$UMAP_2)))/4, color = "grey80", size = 6) +
    annotate(geom = "segment", x = -Inf, y = -Inf,
             xend = min(temp_umap$UMAP_1) + (abs(min(temp_umap$UMAP_1)) + abs(max(temp_umap$UMAP_1)))/4,
             yend = -Inf, color = "grey80", size = 6) +
    scale_x_continuous(name = "UMAP 1", breaks = round(seq(min(temp_umap$UMAP_1), max(temp_umap$UMAP_1), by = 0.20),1),
                       limits = c(min(temp_umap$UMAP_1) * 1.1, max(temp_umap$UMAP_1) * 1.1)) +
    scale_y_continuous(name = "UMAP 2", breaks = round(seq(min(temp_umap$UMAP_2), max(temp_umap$UMAP_2), by = 0.20),1),
                       limits = c(min(temp_umap$UMAP_2) * 1.1, max(temp_umap$UMAP_2) * 1.1)) +
    scale_color_manual(values = palopal_temp, na.value = "lightgrey") +
    ggtitle(paste0(progenitor_umap, " UMAP [2D]")) +
    theme(rect = element_rect(fill = "transparent"),
          plot.title = element_blank(),
          axis.line = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title.x = element_text(hjust=0.1, size = 32),
          axis.title.y = element_text(hjust=0.1, size = 32),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  plot_list[[progenitor_umap]] <- plot_grid(cell_type_plot,
                                            plot_grid(time_plot,
                                                      species_plot, nrow = 2),
                                            ncol = 2, rel_widths = c(1, 0.5))
}

# work on resizing stuffs
pdf(paste0(dir, "Plots/supplemental_plots/time_umaps/progenitor_0_150.pdf"), width = 16, height = 12)
plot_list[["0 to 150 embryo time"]]
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/time_umaps/progenitor_0_200.pdf"), width = 24, height = 18)
plot_list[["0 to 200 embryo time"]]
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/time_umaps/progenitor_0_250.pdf"), width = 24, height = 18)
plot_list[["0 to 250 embryo time"]]
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/time_umaps/progenitor_0_300.pdf"), width = 24, height = 18)
plot_list[["0 to 300 embryo time"]]
dev.off()



# The y-axis will be the mean embryo time and the x-axis will be the birth and death times

plot_df <- left_join(CellTable[! is.na(CellTable$MergedDatasetName),], cell_data, by = join_by(MergedDatasetName == cell_type)) %>%
  mutate(d_time = replace(d_time, is.na(d_time), 720)) %>%
  filter(! is.na(embryo_time))

plot_df$d_time <- as.numeric(plot_df$d_time)
plot_df$br_time <- as.numeric(plot_df$br_time)

plot_df$mean_time <- (plot_df$br_time + plot_df$d_time) / 2

reg <- lm(formula = mean_time ~ embryo_time,
          data = plot_df)

#get intercept and slope value
coeff <- coefficients(reg)          
intercept <-coeff[1]
slope <- coeff[2]

pdf(paste0(dir, "Plots/supplemental_plots/embryo_time_division_time.pdf"), width = 5.5, height = 5)
plot_df %>%
ggplot() +
  geom_segment(aes(x = as.numeric(br_time),
                   y = as.numeric(embryo_time),
                   xend = as.numeric(d_time),
                   yend = as.numeric(embryo_time), color = div_stage),
               lineend = "butt") +
  geom_point(aes(x = as.numeric(br_time),
                 y = as.numeric(embryo_time),
                 color = div_stage),
             size = 1) +
  geom_point(aes(x = as.numeric(d_time),
                 y = as.numeric(embryo_time),
                 color = div_stage),
             size = 1) +
  geom_abline(intercept = intercept,
              slope = slope,
              color="black", linetype="dashed") +
  scale_x_continuous(name = "Cell type birth and death time") +
  scale_y_continuous(name = "Median embryo time") +
  scale_color_manual(values = c('15' = "#C6DBEFFF", '28' = "#9ECAE1FF", '50' = "#6BAED6FF",
                       '100' = "#4292C6FF", '200' = "#2171B5FF", '350' = "#08519CFF", "600" = "#08306BFF")) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        rect = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()


######################################
# fkh-8 and daf-19 expression timing #
######################################

time_bin_vector <- c("lt_100", "100_130", "130_170", "170_210", "210_270", "270_330", "330_390", "390_450", "450_510", "510_580", "580_650", "650_710", "gt_710")

time_bin_df <- data.frame(bins = time_bin_vector,
                          start = c(0, 100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710),
                          end = c(100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710, 1000))

core_cilia <- c("arl-13", "bbs-1", "bbs-2", "bbs-4", "bbs-5",
                "bbs-8", "bbs-9", "ccep-290", "cfap-36", "che-10",
                "che-11", "che-12", "che-13", "che-2", "che-3",
                "clic-1", "daf-10", "dli-1", "dyci-1", "dyf-1",
                "dyf-11", "dyf-13", "dyf-19", "dyf-2", "dyf-3",
                "dyf-5", "dyf-6", "dyf-7", "dylt-2", "peli-1",
                "hyls-1", "ift-139", "ift-20", "ift-74", "ift-81",
                "ifta-1", "ifta-2", "jbts-14", "kap-1", "klp-11",
                "klp-20", "mks-1", "mks-2", "mks-3", "mks-5",
                "mks-6", "mksr-1", "mksr-2", "nphp-1", "nphp-2",
                "nphp-4", "nubp-1", "rbg-3", "osm-1", "osm-12",
                "osm-3", "osm-5", "osm-6", "arl-3", "arl-6",
                "rab-8", "rpi-2", "osta-1", "tba-5", "tba-6",
                "tbb-4", "tmem-107", "tmem-138", "tmem-17", "tmem-218",
                "tmem-231", "tub-1", "xbx-1")

dtw_cilia <- c("C54C6.6", "kap-1", "nubp-1", "che-10",
               "dyf-6", "nphp-1", "dyf-13", "osm-12",
               "osm-6", "dyci-1", "ifta-2", "mks-6",
               "osm-5", "R02E12.4", "ccep-290", "che-13",
               "che-2", "dyf-19", "dyf-2", "dyf-3",
               "mks-2", "tmem-218")

core_cilia_cel_exp <- TPMLineageList[["C.elegans"]][core_cilia[core_cilia %in% rownames(TPMLineageList[["C.briggsae"]])], grepl("Ciliated neurons", colnames(TPMLineageList[["C.elegans"]]))]
core_cilia_cel_exp$gene <- rownames(core_cilia_cel_exp)
core_cilia_cbr_exp <- TPMLineageList[["C.briggsae"]][core_cilia[core_cilia %in% rownames(TPMLineageList[["C.briggsae"]])], grepl("Ciliated neurons", colnames(TPMLineageList[["C.briggsae"]]))]
core_cilia_cbr_exp$gene <- rownames(core_cilia_cbr_exp)

core_cilia_cel_exp$species <- "C.elegans"
core_cilia_cbr_exp$species <- "C.briggsae"
core_cilia_exp <- rbind(core_cilia_cel_exp[,c(paste0("Ciliated neurons_", time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$bins), "species", "gene")],
                       core_cilia_cbr_exp[,c(paste0("Ciliated neurons_", time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$bins), "species", "gene")])

core_cilia_exp_melt <- reshape2::melt(core_cilia_exp, id.vars = c("species", "gene"))

head(core_cilia_exp_melt)
core_cilia_exp_melt$line_type = "solid"
core_cilia_exp_melt[core_cilia_exp_melt$gene %in% dtw_cilia,]$line_type <- "dashed"

pdf(paste0(dir, "Plots/supplemental_plots/cilia_expression.pdf"), width = 5, height = 5.5)
core_cilia_exp_melt %>%
  group_by(species, gene) %>%
  mutate(norm_expression = value / max(value, na.rm = TRUE)) %>%
  mutate(norm_expression_dir = ifelse(species == "C.briggsae", -norm_expression, norm_expression)) %>%
  ggplot(aes(x = variable, y = norm_expression_dir,
             color = line_type,
             group = paste0(species, gene),
             linetype = line_type)) +
  geom_line(size = 1, alpha = 0.5) +
  scale_x_discrete(name = "Embryo time bin", limits = paste0("Ciliated neurons_", time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$bins),
                   labels = paste0(time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$start, " to ",
                                   time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$end)) +
  scale_y_continuous(name = "Max normalized expression") +
  scale_linetype_manual(values = c('solid' = "solid",
                                   'dashed' = "dashed")) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        rect = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        # panel.grid.major.x = element_blank(),
        # panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/cilia_expression_exp_mean.pdf"), width = 5, height = 5.5)
core_cilia_exp_melt %>%
  group_by(species, variable) %>%
  summarise(mean_exp = mean(value, na.rm = TRUE)) %>%
  mutate(norm_expression = mean_exp / max(mean_exp, na.rm = TRUE)) %>%
  # mutate(norm_expression_dir = ifelse(species == "C.briggsae", -norm_expression, norm_expression)) %>%
  ggplot(aes(x = variable, y = norm_expression,
             color = species,
             group = species)) +
  geom_line(size = 1, alpha = 0.5) +
  scale_x_discrete(name = "Embryo time bin", limits = paste0("Ciliated neurons_", time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$bins),
                   labels = paste0(time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$start, " to ",
                                   time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$end)) +
  scale_y_continuous(name = "Max normalized expression") +
  # scale_linetype_manual(values = c('solid' = "solid",
  #                                  'dashed' = "dashed")) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        rect = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        # panel.grid.major.x = element_blank(),
        # panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

dtw_cilia_cel_exp <- TPMLineageList[["C.elegans"]][dtw_cilia[dtw_cilia %in% rownames(TPMLineageList[["C.briggsae"]])], grepl("Ciliated neurons", colnames(TPMLineageList[["C.elegans"]]))]
dtw_cilia_cbr_exp <- TPMLineageList[["C.briggsae"]][dtw_cilia[dtw_cilia %in% rownames(TPMLineageList[["C.briggsae"]])], grepl("Ciliated neurons", colnames(TPMLineageList[["C.briggsae"]]))]

timing_gene_df <- rbind(rbind(data.frame(gene = "daf-19",
                                               time = c("100_130", "130_170",
                                                        "170_210", "210_270", "270_330",
                                                        "330_390", "390_450", "450_510",
                                                        "510_580", "580_650", "650_710"),
                                               expression = unlist(TPMLineageList[["C.elegans"]][c("daf-19"), grepl("Ciliated neurons", colnames(TPMLineageList[["C.elegans"]]))]),
                                               species = "C.elegans"),
                                    data.frame(gene = "fkh-8",
                                               time = c("100_130", "130_170",
                                                        "170_210", "210_270", "270_330",
                                                        "330_390", "390_450", "450_510",
                                                        "510_580", "580_650", "650_710"),
                                               expression = unlist(TPMLineageList[["C.elegans"]][c("fkh-8"), grepl("Ciliated neurons", colnames(TPMLineageList[["C.elegans"]]))]),
                                               species = "C.elegans")),
                              rbind(data.frame(gene = "daf-19",
                                               time = c("lt_100", "100_130", "130_170",
                                                        "170_210", "210_270", "270_330",
                                                        "330_390", "390_450", "450_510",
                                                        "510_580", "580_650", "650_710", "gt_710"),
                                               expression = unlist(TPMLineageList[["C.briggsae"]][c("daf-19"), grepl("Ciliated neurons", colnames(TPMLineageList[["C.briggsae"]]))]),
                                               species = "C.briggsae"),
                                    data.frame(gene = "fkh-8",
                                               time = c("lt_100", "100_130", "130_170",
                                                        "170_210", "210_270", "270_330",
                                                        "330_390", "390_450", "450_510",
                                                        "510_580", "580_650", "650_710", "gt_710"),
                                               expression = unlist(TPMLineageList[["C.briggsae"]][c("fkh-8"), grepl("Ciliated neurons", colnames(TPMLineageList[["C.briggsae"]]))]),
                                               species = "C.briggsae"))) %>%
  group_by(species, gene) %>%
  mutate(norm_expression = expression / max(expression)) %>% 
  rbind(., rbind(data.frame(gene = "core_cilia",
                            time = c("lt_100", "100_130", "130_170",
                                     "170_210", "210_270", "270_330",
                                     "330_390", "390_450", "450_510",
                                     "510_580", "580_650", "650_710", "gt_710"),
                            expression = colMeans(dtw_cilia_cbr_exp),
                            norm_expression = rowMeans(apply(dtw_cilia_cbr_exp, 1, function(x) x/max(x)))/max(rowMeans(apply(dtw_cilia_cbr_exp, 1, function(x) x/max(x)))),
                            species = "C.briggsae"),
                 data.frame(gene = "core_cilia",
                            time = c("100_130", "130_170",
                                     "170_210", "210_270", "270_330",
                                     "330_390", "390_450", "450_510",
                                     "510_580", "580_650", "650_710"),
                            expression = colMeans(dtw_cilia_cel_exp),
                            norm_expression = rowMeans(apply(dtw_cilia_cel_exp, 1, function(x) x/max(x)))/max(rowMeans(apply(dtw_cilia_cel_exp, 1, function(x) x/max(x)))),
                            species = "C.elegans"))) %>%
  filter(time != "lt_100" & time != "gt_710")

pdf(paste0(dir, "Plots/supplemental_plots/daf_19_fkh_8_expression.pdf"), width = 5, height = 5.5)
timing_gene_df %>%
  ggplot(aes(x = time, y = norm_expression,
             color = paste0(species, "_", gene),
             group = paste0(species, gene),
             linetype = gene)) +
  geom_line(size = 1) +
  scale_color_manual(name = "Species", values = c('C.briggsae_daf-19' = "#56B4E9",
                                                  'C.briggsae_fkh-8' = "#56B4E9",
                                                  'C.elegans_daf-19' = "#009E73",
                                                  'C.elegans_fkh-8' = "#009E73",
                                                  'C.elegans_core_cilia' = "orange",
                                                  'C.briggsae_core_cilia' = "red")) +
  scale_linetype_manual(name = "Gene", values = c('daf-19' = "solid",
                                                  'fkh-8' = "dashed",
                                                  'core_cilia' = "solid")) +
  scale_x_discrete(name = "Embryo time bin", limits = time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$bins,
                   labels = paste0(time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$start, " to ",
                                   time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$end)) +
  scale_y_continuous(name = "Max normalized expression") +
  theme(legend.title = element_blank(),
        legend.position = "top",
        rect = element_rect(fill = "transparent"),
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

pdf(paste0(dir, "Plots/supplemental_plots/daf_19_fkh_8_expression_grouped.pdf"), width = 5, height = 5.5)
timing_gene_df %>%
  ggplot(aes(x = time, y = norm_expression,
             color = species,
             group = paste0(species, "_", gene))) +
  facet_wrap(~gene, nrow = 3) + 
  geom_line(size = 1.5) +
  scale_color_manual(name = "Species", values = c('C.briggsae' = "#56B4E9",
                                                  'C.elegans' = "#009E73")) +
  scale_x_discrete(name = "Embryo time bin", limits = time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$bins,
                   labels = paste0(time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$start, " to ",
                                   time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$end)) +
  scale_y_continuous(name = "Max normalized expression") +
  theme(legend.title = element_blank(),
        legend.position = "top",
        rect = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.line = element_line(color="grey80", size=1),
        strip.background = element_rect(fill = "transparent"),
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


CN_TPMLineageList <- list()
CN_TPMLineageList[["C.elegans"]] <- TPMLineageList[["C.elegans"]][,grepl("Ciliated neurons", colnames(TPMLineageList[["C.elegans"]]))]
CN_TPMLineageList[["C.briggsae"]] <- TPMLineageList[["C.briggsae"]][,grepl("Ciliated neurons", colnames(TPMLineageList[["C.briggsae"]]))]

timing_sets_df <- rbind(
  rbind(
    rbind(data.frame(gene = "daf-19",
                     time = c("100_130", "130_170",
                              "170_210", "210_270", "270_330",
                              "330_390", "390_450", "450_510",
                              "510_580", "580_650", "650_710"),
                     expression = unlist(CN_TPMLineageList[["C.elegans"]][c("daf-19"),]),
                     norm_expression = unlist(CN_TPMLineageList[["C.elegans"]][c("daf-19"),]) / max(CN_TPMLineageList[["C.elegans"]]["daf-19",], CN_TPMLineageList[["C.briggsae"]]["daf-19",]),
                     species = "C.elegans"),
          data.frame(gene = "daf-19",
                     time = c("lt_100", "100_130", "130_170",
                              "170_210", "210_270", "270_330",
                              "330_390", "390_450", "450_510",
                              "510_580", "580_650", "650_710", "gt_710"),
                     expression = unlist(CN_TPMLineageList[["C.briggsae"]][c("daf-19"),]),
                     norm_expression = unlist(CN_TPMLineageList[["C.briggsae"]][c("daf-19"),]) / max(CN_TPMLineageList[["C.elegans"]]["daf-19",], CN_TPMLineageList[["C.briggsae"]]["daf-19",]),
                     species = "C.briggsae")),
    rbind(data.frame(gene = "fkh-8",
                     time = c("100_130", "130_170",
                              "170_210", "210_270", "270_330",
                              "330_390", "390_450", "450_510",
                              "510_580", "580_650", "650_710"),
                     expression = unlist(CN_TPMLineageList[["C.elegans"]][c("fkh-8"),]),
                     norm_expression = unlist(CN_TPMLineageList[["C.elegans"]][c("fkh-8"),]) / max(CN_TPMLineageList[["C.elegans"]]["fkh-8",], CN_TPMLineageList[["C.briggsae"]]["fkh-8",]),
                     species = "C.elegans"),
          data.frame(gene = "fkh-8",
                     time = c("lt_100", "100_130", "130_170",
                              "170_210", "210_270", "270_330",
                              "330_390", "390_450", "450_510",
                              "510_580", "580_650", "650_710", "gt_710"),
                     expression = unlist(CN_TPMLineageList[["C.briggsae"]][c("fkh-8"),]),
                     norm_expression = unlist(CN_TPMLineageList[["C.briggsae"]][c("fkh-8"),]) / max(CN_TPMLineageList[["C.elegans"]]["fkh-8",], CN_TPMLineageList[["C.briggsae"]]["fkh-8",]),
                     species = "C.briggsae"))),
  rbind(data.frame(gene = "core_cilia",
                   time = c("100_130", "130_170",
                            "170_210", "210_270", "270_330",
                            "330_390", "390_450", "450_510",
                            "510_580", "580_650", "650_710"),
                   expression = colMeans(core_cilia_cel_exp[,grepl("Ciliated", colnames(core_cilia_cel_exp))]),
                   norm_expression = colMeans(core_cilia_cel_exp[,grepl("Ciliated", colnames(core_cilia_cel_exp))]) / max(colMeans(core_cilia_cel_exp[,grepl("Ciliated", colnames(core_cilia_cel_exp))]), colMeans(core_cilia_cbr_exp[,grepl("Ciliated", colnames(core_cilia_cbr_exp))])),
                   species = "C.elegans"),
        data.frame(gene = "core_cilia",
                   time = c("lt_100", "100_130", "130_170",
                            "170_210", "210_270", "270_330",
                            "330_390", "390_450", "450_510",
                            "510_580", "580_650", "650_710", "gt_710"),
                   expression = colMeans(core_cilia_cbr_exp[,grepl("Ciliated", colnames(core_cilia_cbr_exp))]),
                   norm_expression = colMeans(core_cilia_cbr_exp[,grepl("Ciliated", colnames(core_cilia_cbr_exp))]) / max(colMeans(core_cilia_cel_exp[,grepl("Ciliated", colnames(core_cilia_cel_exp))]), colMeans(core_cilia_cbr_exp[,grepl("Ciliated", colnames(core_cilia_cbr_exp))])),
                   species = "C.briggsae")))


pdf(paste0(dir, "Plots/supplemental_plots/daf_19_fkh_8_expression.pdf"), width = 5, height = 5.5)
timing_sets_df %>%
  ggplot(aes(x = time, y = norm_expression,
             color = gene,
             group = paste0(species, gene),
             linetype = species)) +
  geom_line(size = 0.5) +
  scale_color_manual(name = "Gene set", values = c("daf-19" = "#CC79A7", "fkh-8" = "#0072B2", "core_cilia" = "#D55E00")) +
  scale_linetype_manual(name = "Species", values = c('C.elegans' = "solid",
                                                  'C.briggsae' = "dashed")) +
  scale_x_discrete(name = "Embryo time bin", limits = time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$bins,
                   labels = paste0(time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$start, " to ",
                                   time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$end)) +
  scale_y_continuous(name = "Max normalized expression") +
  guides(color = guide_legend(title = "Gene set"),
         linetype = guide_legend(title = "Species", nrow = 2)) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        rect = element_rect(fill = "transparent"),
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

## cep-1/gld-1
timing_sets_df <- rbind(
  rbind(data.frame(gene = "cep-1",
                   time = c("100_130", "130_170",
                            "170_210", "210_270", "270_330",
                            "330_390", "390_450", "450_510",
                            "510_580", "580_650", "650_710"),
                   expression = unlist(TPMLineageList[["C.elegans"]][c("cep-1"),]),
                   norm_expression = unlist(TPMLineageList[["C.elegans"]][c("cep-1"),]) / max(TPMLineageList[["C.elegans"]]["cep-1",], TPMLineageList[["C.briggsae"]]["cep-1",]),
                   species = "C.elegans"),
        data.frame(gene = "cep-1",
                   time = c("lt_100", "100_130", "130_170",
                            "170_210", "210_270", "270_330",
                            "330_390", "390_450", "450_510",
                            "510_580", "580_650", "650_710", "gt_710"),
                   expression = unlist(TPMLineageList[["C.briggsae"]][c("cep-1"),]),
                   norm_expression = unlist(TPMLineageList[["C.briggsae"]][c("cep-1"),]) / max(TPMLineageList[["C.elegans"]]["cep-1",], TPMLineageList[["C.briggsae"]]["cep-1",]),
                   species = "C.briggsae")),
  rbind(data.frame(gene = "gld-1",
                   time = c("100_130", "130_170",
                            "170_210", "210_270", "270_330",
                            "330_390", "390_450", "450_510",
                            "510_580", "580_650", "650_710"),
                   expression = unlist(TPMLineageList[["C.elegans"]][c("gld-1"),]),
                   norm_expression = unlist(TPMLineageList[["C.elegans"]][c("gld-1"),]) / max(TPMLineageList[["C.elegans"]]["gld-1",], TPMLineageList[["C.briggsae"]]["gld-1",]),
                   species = "C.elegans"),
        data.frame(gene = "gld-1",
                   time = c("lt_100", "100_130", "130_170",
                            "170_210", "210_270", "270_330",
                            "330_390", "390_450", "450_510",
                            "510_580", "580_650", "650_710", "gt_710"),
                   expression = unlist(TPMLineageList[["C.briggsae"]][c("gld-1"),]),
                   norm_expression = unlist(TPMLineageList[["C.briggsae"]][c("gld-1"),]) / max(TPMLineageList[["C.elegans"]]["gld-1",], TPMLineageList[["C.briggsae"]]["gld-1",]),
                   species = "C.briggsae")))


pdf(paste0(dir, "Plots/supplemental_plots/daf_19_fkh_8_expression.pdf"), width = 5, height = 5.5)
timing_sets_df %>%
  ggplot(aes(x = time, y = norm_expression,
             color = gene,
             group = paste0(species, gene),
             linetype = species)) +
  geom_line(size = 0.5) +
  scale_color_manual(name = "Gene set", values = c("daf-19" = "#CC79A7", "fkh-8" = "#0072B2", "core_cilia" = "#D55E00")) +
  scale_linetype_manual(name = "Species", values = c('C.elegans' = "solid",
                                                     'C.briggsae' = "dashed")) +
  scale_x_discrete(name = "Embryo time bin", limits = time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$bins,
                   labels = paste0(time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$start, " to ",
                                   time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$end)) +
  scale_y_continuous(name = "Max normalized expression") +
  guides(color = guide_legend(title = "Gene set"),
         linetype = guide_legend(title = "Species", nrow = 2)) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        rect = element_rect(fill = "transparent"),
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

###################
# log2(Tau) x JSD #
###################

pdf(paste0(dir, "Plots/supplemental_plots/tau_log2_jsd.pdf"), width = 6, height = 5)
gene_data %>%
  filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>%
  ggplot(aes(x = jsd_median_joint, y = log2(cel_tau_median_joint/cbr_tau_median_joint), color = mean_tau_joint)) +
  geom_point() +
  scale_color_viridis(name = "Mean tau") +
  scale_x_continuous(name = "Cell distance (JSD)") +
  scale_y_continuous(name = "log2(C. elegans Tau / C. briggsae Tau)") +
  theme(legend.position = "right",
        rect = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/tau_log2_tau"), width = 6, height = 5)
gene_data %>%
  filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>%
  ggplot(aes(x = mean_tau_joint, y = log2(cel_tau_median_joint/cbr_tau_median_joint), color = jsd_median_joint)) +
  geom_point() +
  scale_color_viridis(name = "Cell distance") +
  scale_x_continuous(name = "Mean tau") +
  scale_y_continuous(name = "log2(C. elegans Tau / C. briggsae Tau)") +
  theme(legend.position = "right",
        rect = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major = element_line(color="grey80", size=0.25),
        panel.grid.minor = element_line(color="grey80", size=0.05),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"))
dev.off()

############################################################
# Synteny + Phylostrata + Chromosome + Recombination region #
############################################################

PS.values <- unlist((gene_data %>% group_by(PS.Name) %>% summarise(PS.number = mean(PS.Value)) %>% filter(! is.na(PS.number)) %>% arrange(PS.number))[,1])

lm(jsd_median_joint ~ PS.Name, gene_data %>% filter(cel_max_tpm > 80 & cbr_max_tpm > 80))

ps_plot <- gene_data %>%
  filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>%
  filter(! is.na(PS.Name)) %>%
  ggplot(aes(x = PS.Name, y = jsd_median_joint)) +
  geom_boxplot() +
  scale_y_continuous(name = "Gene distance (Jensen-Shannon distance)") +
  scale_x_discrete(name = "Phylostratigraphic age", limits = PS.values) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
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

pdf(paste0(dir, "Plots/supplemental_plots/ps_value.pdf"), width = 8, height = 5)
print(ps_plot)
dev.off()

syteny_plot <- gene_data %>%
  filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>%
  filter(! is.na(syntenic)) %>%
  ggplot(aes(x = syntenic, y = jsd_median_joint)) +
  geom_boxplot() +
  scale_y_continuous(name = "Gene distance (Jensen-Shannon distance)", limits = c(0, 1)) +
  scale_x_discrete(name = "Gene synteny", limits = c("IS_SYNTENIC", "NOT_SYNTENIC"), labels = c("Syntenic", "Not syntenic")) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
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

pdf(paste0(dir, "Plots/supplemental_plots/syteny.pdf"), width = 2.5, height = 5)
print(syteny_plot)
dev.off()

gene_data_temp <- gene_data
gene_data_temp$cel_RecombRegion <- NA
gene_data_temp$cbr_RecombRegion <- NA
for(cur_gene in gene_data_temp$gene) {
  gene_data_temp[cur_gene,]$cel_RecombRegion <- gff_list_mrna[["elegans"]][which(gff_list_mrna[["elegans"]]$cds_gene_name == cur_gene),]$RecombRegion
  gene_data_temp[cur_gene,]$cbr_RecombRegion <- gff_list_mrna[["briggsae"]][which(gff_list_mrna[["briggsae"]]$cds_gene_name == cur_gene),]$RecombRegion
}

recomb_regions <- unique(gene_data$RecombRegion)
recomb_regions <- recomb_regions[!is.na(recomb_regions)]

recomb_simple_list <- list("Tip", "Tip", "Center", "Arm", "Arm")
names(recomb_simple_list) <- c("LeftTip", "RightTip", "Center", "LeftArm", "RightArm")

gene_data_temp$cel_simple_recomb_region <- NA
gene_data_temp$cbr_simple_recomb_region <- NA
for(recomb in recomb_regions) {
  recomb_simple <- unique(unlist(lapply(strsplit(recomb, "_"), function(x) x[2])))
  gene_data_temp[which(gene_data_temp$cel_RecombRegion == recomb),]$cel_simple_recomb_region <- unlist(recomb_simple_list[recomb_simple])
  gene_data_temp[which(gene_data_temp$cbr_RecombRegion == recomb),]$cbr_simple_recomb_region <- unlist(recomb_simple_list[recomb_simple])
}

# check if they maintain the same recombination region or not
gene_data_temp$recomb_region_match <- gene_data_temp$cel_RecombRegion == gene_data_temp$cbr_RecombRegion
gene_data_temp$cel_simple_recomb_region_match <- ifelse(gene_data_temp$recomb_region_match, gene_data_temp$cel_simple_recomb_region, "Translocation")

recomb_plot <- gene_data_temp %>%
  filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>%
  filter(! is.na(cel_simple_recomb_region_match)) %>%
  ggplot(aes(x = cel_simple_recomb_region_match, y = jsd_median_joint)) +
  geom_boxplot() +
  stat_summary(fun.data = n_fun, geom = "text", angle = 90) +
  scale_y_continuous(name = "Gene distance (Jensen-Shannon distance)", limits = c(0, 1)) +
  scale_x_discrete(name = "Gene recombination region", limits = c("Tip", "Arm", "Center", "Translocation")) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
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


pdf(paste0(dir, "Plots/supplemental_plots/same_recomb.pdf"), width = 5, height = 5)
print(recomb_plot)
dev.off()


maternal_plot <- gene_data %>%
  filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>%
  filter(! is.na(maternal)) %>%
  ggplot(aes(x = maternal, y = jsd_median_joint)) +
  geom_boxplot() +
  scale_y_continuous(name = "Gene distance (Jensen-Shannon distance)", limits = c(0, 1)) +
  scale_x_discrete(name = "Maternal inheritance", limits = c(TRUE, FALSE), labels = c("Maternally inherited", "Not maternally inherited")) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
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

pdf(paste0(dir, "Plots/supplemental_plots/jsd_maternal.pdf"), width = 2.5, height = 5)
print(maternal_plot)
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/random_associations.pdf"), width = 10, height = 10)
plot_grid(plot_grid(syteny_plot, maternal_plot, recomb_plot, ncol = 3, align = "hv", rel_widths = c(0.5, 0.5, 1), axis = "tblr"),
          ps_plot, ncol = 1, align = "v")
dev.off()


pairwise.wilcox.test(gene_data[which(gene_data_$cel_max_tpm > 80 & gene_data$cbr_max_tpm > 80),]$jsd_median_joint,
                     gene_data[which(gene_data_$cel_max_tpm > 80 & gene_data$cbr_max_tpm > 80),]$syntenic, p.adj = "bonferroni")


pairwise.wilcox.test(gene_data[which(gene_data$cel_max_tpm > 80 & gene_data$cbr_max_tpm > 80),]$jsd_median_joint,
                     gene_data[which(gene_data$cel_max_tpm > 80 & gene_data$cbr_max_tpm > 80),]$maternal, p.adj = "bonferroni")

pairwise.wilcox.test(gene_data_temp[which(gene_data_temp$cel_max_tpm > 80 & gene_data_temp$cbr_max_tpm > 80),]$jsd_median_joint,
                     gene_data_temp[which(gene_data_temp$cel_max_tpm > 80 & gene_data_temp$cbr_max_tpm > 80),]$cel_simple_recomb_region_match, p.adj = "bonferroni")

cor(gene_data[which(gene_data$cel_max_tpm > 80 & gene_data$cbr_max_tpm > 80),]$jsd_median_joint,
    gene_data[which(gene_data$cel_max_tpm > 80 & gene_data$cbr_max_tpm > 80),]$PS.Value, use = "na.or.complete")

# pearson r of 0.42

temp <- left_join(data.frame(table((gene_data %>%
  filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>%
  filter(! is.na(PS.Name)) %>%
  filter(mean_tau_joint < 0.4))$PS.Name)),
  data.frame(table((gene_data %>%
           filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>%
           filter(! is.na(PS.Name)))$PS.Name)), by = "Var1")

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
    annotate(geom = "text", x = 0.025, y = 0.175, hjust = 0.5, size = 7,
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
          legend.key.size = unit(7, 'mm'),
          legend.text = element_text(size = 12))
  
  cbr_plot <- ggraph(layout_list[["C.briggsae"]]) +
    geom_node_arc_bar(aes(filter = level > 2, fill = tpm + 1), 
                      color = "grey20", size = 0.25) +
    coord_fixed() +
    geom_node_text(aes(filter = level < 4 & level > 2,
                       label = name,
                       angle = node_angle(x, y) + 90), size = 3, colour = 'white') +
    annotate(geom = "text", x = 0.025, y = 0.175, hjust = 0.5, size = 7,
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
          legend.key.size = unit(7, 'mm'),
          legend.text = element_text(size = 12))
  
  layout_list[["Joint"]] <- layout_list[["C.elegans"]]
  # layout_list[["Joint"]]$tpm <- log2(layout_list[["C.elegans"]]$tpm + 1) - log2(layout_list[["C.briggsae"]]$tpm + 1)
  
  layout_list[["Joint"]] <- left_join(layout_list[["Joint"]], layout_list[["C.elegans"]][, c("name", "tpm")], suffix = c("", "_cel"), by = c("name" = "name"))
  layout_list[["Joint"]] <- left_join(layout_list[["Joint"]], layout_list[["C.briggsae"]][, c("name", "tpm")], suffix = c("", "_cbr"), by = c("name" = "name"))
  
  layout_list[["Joint"]]$max_tpm <- apply(layout_list[["Joint"]], 1, function(x) {
    if(! is.na(x["tpm_cel"]) & ! is.na(x["tpm_cbr"])) {
      return(as.numeric(max(x["tpm_cel"], x["tpm_cbr"], na.rm = TRUE)))
    } else {
      return(NA)
    }
  })
  num_cells <- length(layout_list[["Joint"]]$tpm)
  layout_list[["Joint"]]$tpm <- log2((layout_list[["Joint"]]$tpm_cel + 1) / (layout_list[["Joint"]]$tpm_cbr + 1)) * log2(layout_list[["Joint"]]$max_tpm + 1)
  
  total_tpm_max <- max(layout_list[["Joint"]]$max_tpm, na.rm = TRUE)
  print(paste("Total max TPM: ", total_tpm_max))
  display_cutoff <- log2(total_tpm_max) * 5
  
  joint_plot <- ggraph(layout_list[["Joint"]]) +
    geom_node_arc_bar(aes(filter = level > 2, fill = tpm),
                      color = "grey20", size = 0.25) +
    coord_fixed() +
    geom_node_text(aes(filter = level < 4 & level > 2,
                       label = name,
                       angle = node_angle(x, y) + 90), size = 3, colour = 'white') +
    annotate(geom = "text", x = 0.025, y = 0.6, hjust = 0.5, size = 8,
             label = bquote(italic(.(gene)))) +
    annotate(geom = "text", x = 0.025, y = 0.15, hjust = 0.5, size = 5,
             label = expression(paste("log2( ", frac(paste(italic("C. elegans")),
                                                     paste(italic("C. briggsae"))), " ) * Max"))) +
    scale_fill_gradientn(colors = brewer.pal(9, 'RdBu'),
                         limits = c(-display_cutoff, display_cutoff),
                         oob = scales::squish) +
    guides(alpha = F, size = F, fill = guide_colorbar(title = "")) +
    theme(legend.position = c(0.49,0.4),
          panel.background = element_rect(fill = 'transparent'),
          plot.background = element_rect(fill = 'transparent', color = NA),
          plot.margin = margin(t = 0, r = -500, b = 0, l = -500, unit = "mm"),
          legend.direction = "horizontal",
          legend.justification = "center",
          legend.key.size = unit(7, 'mm'),
          legend.text = element_text(size = 12))
  
  #log2(C. elegans/C. briggsae
  circle_plots <- plot_grid(cel_plot, joint_plot, cbr_plot, ncol = 3, scale = 1.05)
}

pdf(paste0(dir, "Plots/supplemental_plots/cep-1.pdf"), height = 5, width = 15)
print(plot_circle_pattern("cep-1", CellTableList, lo))
dev.off()

pdf(paste0(dir, "Plots/supplemental_plots/gld-1.pdf"), height = 5, width = 15)
print(plot_circle_pattern("gld-1", CellTableList, lo))
dev.off()


####
# Gene duplications
####

load(paste0(dir, "Objects/BobDownload_20241120/ele_files_for_fig_5.RData"))
# ele_2_1_better_gene
# ele_2_1_poorer_gene
# genes_1_1_ele
# genes_m_m_ele
# tau_ele

load(paste0(dir, "Objects/BobDownload_20241120/bri_files_for_fig_5.RData"))
# bri_1_2_better_gene
# bri_1_2_poorer_gene
# genes_1_1_bri
# genes_m_m_bri
# tau_bri

cel_dup = tau_ele  %>% mutate(og_class2 = ifelse(gene %in%
                                                   genes_1_1_ele, "1_1",
                                                 ifelse(gene %in% ele_2_1_better_gene, "2:1 better",
                                                        ifelse(gene %in%
                                                                 ele_2_1_poorer_gene, "2:1 poorer",
                                                               ifelse( gene %in% genes_m_m_ele, "m:m",
                                                                       "other")))))

cbr_dup = tau_bri  %>% mutate(og_class2 = ifelse(gene_name %in%
                                                   genes_1_1_bri, "1_1",
                                                 ifelse(gene_name %in% bri_1_2_better_gene, "2:1 better",
                                                        ifelse(gene_name %in%
                                                                 bri_1_2_poorer_gene, "2:1 poorer",
                                                               ifelse(
                                                                 gene_name %in% genes_m_m_bri, "m:m",
                                                                 "other")))))

cel_dup$species = "C. elegans"
cbr_dup$species = "C. briggsae"

min_exp <- min(c(cel_dup$sum, cbr_dup$sum))
max_exp <- max(c(cel_dup$sum, cbr_dup$sum))

rownames(cel_dup) <- cel_dup$gene
rownames(cbr_dup) <- cbr_dup$joint_name

cel_dup$cel_max_tpm <- apply(TPMJoint[["C.elegans"]][rownames(cel_dup),], 1, max)
cbr_dup$cbr_max_tpm <- apply(TPMJoint[["C.briggsae"]][rownames(cbr_dup),], 1, max)

cel_dup$og <- gff_list_mrna[["elegans"]][cel_dup$gene,]$OG
cbr_dup$og <- gff_list_mrna[["briggsae"]][cbr_dup$joint_name,]$OG

cel_dup$og_min_exp <- NA
cbr_dup$og_min_exp <- NA

cel_dup$og_min_tau <- NA
cbr_dup$og_min_tau <- NA

all_og <- unique(c(cel_dup$og, cbr_dup$og))
for(cur_og in all_og[! is.na(all_og)]) {
  og_min_exp <- min(cel_dup[which(cel_dup$og == cur_og), "cel_max_tpm"], cbr_dup[which(cbr_dup$og == cur_og), "cbr_max_tpm"], na.rm = TRUE)
  og_min_tau <- min(cel_dup[which(cel_dup$og == cur_og), "tau_ele"], cbr_dup[which(cbr_dup$og == cur_og), "tau_bri"], na.rm = TRUE)
  
  if(cur_og %in% cel_dup$og) {
    cel_dup[which(cel_dup$og == cur_og),]$og_min_exp <- og_min_exp
    cel_dup[which(cel_dup$og == cur_og),]$og_min_tau <- og_min_tau
  }
  
  if(cur_og %in% cbr_dup$og) {
    cbr_dup[which(cbr_dup$og == cur_og),]$og_min_exp <- og_min_exp
    cbr_dup[which(cbr_dup$og == cur_og),]$og_min_tau <- og_min_tau
  }
}

cel_dup$og_max_sum <- NA
cbr_dup$og_max_sum <- NA
for(cur_og in all_og[! is.na(all_og)]) {
  if(cur_og %in% cel_dup$og) {
    cel_dup[which(cel_dup$og == cur_og),]$og_max_sum <- max(cel_dup[which(cel_dup$og == cur_og),]$sum)
  }
  
  if(cur_og %in% cbr_dup$og) {
    cbr_dup[which(cbr_dup$og == cur_og),]$og_max_sum <- max(cbr_dup[which(cbr_dup$og == cur_og),]$sum)
  }
}

cel_dup_filt <- cel_dup[which(cel_dup$og_min_exp > 80),]
cbr_dup_filt <- cbr_dup[which(cbr_dup$og_min_exp > 80),]

subset_og <- unique(c(cel_dup_filt$og, cbr_dup_filt$og))
for(cur_og in subset_og[! is.na(subset_og)]) {
  if(sum(cel_dup_filt[which(cel_dup_filt$og == cur_og),]$og_class2 %in% "2:1 better") > 0) {
    print(cel_dup_filt[which(cel_dup_filt$og == cur_og),]$gene)
    print(cel_dup_filt[which(cel_dup_filt$og == cur_og),]$og_class2)
    print(cel_dup_filt[which(cel_dup_filt$og == cur_og),]$tau_ele)
  }
}

# 2:1 in general are expressed at lower levels (Fine)
wilcox.test(cel_dup[which(cel_dup$og_class2 %in% c("2:1 better", "2:1 poorer")),]$sum, cel_dup[which(cel_dup$og_class2 == "1_1"),]$sum)
wilcox.test(cbr_dup[which(cbr_dup$og_class2 %in% c("2:1 better", "2:1 poorer")),]$sum, cbr_dup[which(cbr_dup$og_class2 == "1_1"),]$sum)

wilcox.test(cel_dup[which(cel_dup$og_class2 %in% c("2:1 better")),]$og_max_sum, cel_dup[which(cel_dup$og_class2 == "1_1"),]$sum)
wilcox.test(cbr_dup[which(cbr_dup$og_class2 %in% c("2:1 better")),]$og_max_sum, cbr_dup[which(cbr_dup$og_class2 == "1_1"),]$sum)

# The better match is expressed at higher levels (Fine)
wilcox.test(cel_dup[which(cel_dup$og_class2 == "2:1 better"),]$sum, cel_dup[which(cel_dup$og_class2 == "2:1 poorer"),]$sum)
wilcox.test(cbr_dup[which(cbr_dup$og_class2 == "2:1 better"),]$sum, cbr_dup[which(cbr_dup$og_class2 == "2:1 poorer"),]$sum)

# how many reduced (need to fix)
# summary(cel_dup[which(cel_dup$og_class2 == "2:1 better"),]$sum / cel_dup[which(cel_dup$og_class2 == "2:1 poorer"),]$sum > 1.5)
# summary(cel_dup[which(cel_dup$og_class2 == "2:1 poorer"),]$sum / cel_dup[which(cel_dup$og_class2 == "2:1 better"),]$sum > 1.5)
# 
# summary(cbr_dup[which(cbr_dup$og_class2 == "2:1 better"),]$sum / cbr_dup[which(cbr_dup$og_class2 == "2:1 poorer"),]$sum > 1.5)
# summary(cbr_dup[which(cbr_dup$og_class2 == "2:1 poorer"),]$sum / cbr_dup[which(cbr_dup$og_class2 == "2:1 better"),]$sum > 1.5)

# Now filtered for only those sets of duplications where all of the genes in the og are expressed at > 80
wilcox.test(cel_dup_filt[which(cel_dup_filt$og_class2 %in% c("2:1 better", "2:1 poorer")),]$tau_ele, cel_dup_filt[which(cel_dup_filt$og_class2 == "1_1"),]$tau_ele)
wilcox.test(cbr_dup_filt[which(cbr_dup_filt$og_class2 %in% c("2:1 better", "2:1 poorer")),]$tau_bri, cbr_dup_filt[which(cbr_dup_filt$og_class2 == "1_1"),]$tau_bri)

wilcox.test(cel_dup_filt[which(cel_dup_filt$og_class2 %in% c("2:1 better")),]$og_min_tau, cel_dup_filt[which(cel_dup_filt$og_class2 == "1_1"),]$tau_ele)
wilcox.test(cbr_dup_filt[which(cbr_dup_filt$og_class2 %in% c("2:1 better")),]$og_min_tau, cbr_dup_filt[which(cbr_dup_filt$og_class2 == "1_1"),]$tau_bri)

wilcox.test(cel_dup_filt[which(cel_dup_filt$og_class2 == "2:1 better"),]$tau_ele, cel_dup_filt[which(cel_dup_filt$og_class2 == "2:1 poorer"),]$tau_ele)
wilcox.test(cbr_dup_filt[which(cbr_dup_filt$og_class2 == "2:1 better"),]$tau_bri, cbr_dup_filt[which(cbr_dup_filt$og_class2 == "2:1 poorer"),]$tau_bri)

cel_exp_plot <- cel_dup %>%
  ggplot(aes(x = sum, color = og_class2)) +
  geom_density() +
  geom_vline(xintercept = cel_dup["sod-2", "sum"], color = "red") +
  geom_vline(xintercept = cel_dup["sod-3", "sum"], color = "red") +
  scale_x_continuous(name = "Total expression (TPM)", trans = "log10",
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     limits = c(1, max_exp),
                     oob = scales::squish) +
  scale_y_continuous(name = "C. elegans") +
  scale_color_manual(values = c('#CC6677','#88CCEE','#DDCC77','#117733', '#332288'),
                     labels = c("1:1", "2:1 better match", "2:1 poorer match", "many:many", "Uncategorized")) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major.x = element_line(color="grey80", size=0.25),
        panel.grid.minor.x = element_line(color="grey80", size=0.05),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        plot.margin = margin(t = 0, r = 2.5, b = 0, l = 2.5, unit = "mm"))

cbr_exp_plot <- cbr_dup %>%
  ggplot(aes(x = sum, color = og_class2)) +
  geom_density() +
  geom_vline(xintercept = cbr_dup["sod-2","sum"], color = "red") +
  scale_x_continuous(name = "Total expression (TPM)", trans = "log10",
                     labels = scales::trans_format("log10", scales::math_format(10^.x)),
                     limits = c(1, max_exp),
                     oob = scales::squish) +
  scale_y_continuous(name = "C. briggsae") +
  scale_color_manual(values = c('#CC6677','#88CCEE','#DDCC77','#117733', '#332288'),
                     labels = c("1:1", "2:1 better match", "2:1 poorer match", "many:many", "Uncategorized")) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major.x = element_line(color="grey80", size=0.25),
        panel.grid.minor.x = element_line(color="grey80", size=0.05),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        plot.margin = margin(t = 0, r = 2.5, b = 0, l = 2.5, unit = "mm"))

min_tau <- min(c(cel_dup_filt$tau_ele, cbr_dup_filt$tau_bri))

cel_tau_plot <- cel_dup_filt %>%
  ggplot(aes(x = tau_ele, color = og_class2)) +
  geom_density() +
  geom_vline(xintercept = cel_dup_filt["sod-2","tau_ele"], color = "red") +
  geom_vline(xintercept = cel_dup_filt["sod-3","tau_ele"], color = "red") +
  scale_x_continuous(name = "Tau", limits = c(min_tau, 1)) +
  scale_y_continuous(name = "C. elegans") +
  scale_color_manual(values = c('#CC6677','#88CCEE','#DDCC77','#117733', '#332288'),
                     labels = c("1:1", "2:1 better match", "2:1 poorer match", "many:many", "Uncategorized")) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major.x = element_line(color="grey80", size=0.25),
        panel.grid.minor.x = element_line(color="grey80", size=0.05),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        plot.margin = margin(t = 0, r = 2.5, b = 0, l = 2.5, unit = "mm"))

cbr_tau_plot <- cbr_dup_filt %>%
  ggplot(aes(x = tau_bri, color = og_class2)) +
  geom_density() +
  geom_vline(xintercept = cbr_dup_filt["sod-2","tau_bri"], color = "red") +
  scale_x_continuous(name = "Tau", limits = c(min_tau, 1)) +
  scale_y_continuous(name = "C. briggsae") +
  scale_color_manual(values = c('#CC6677','#88CCEE','#DDCC77','#117733', '#332288'),
                     labels = c("1:1", "2:1 better match", "2:1 poorer match", "many:many", "Uncategorized")) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major.x = element_line(color="grey80", size=0.25),
        panel.grid.minor.x = element_line(color="grey80", size=0.05),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        plot.margin = margin(t = 0, r = 2.5, b = 0, l = 2.5, unit = "mm"))

pdf(paste0(dir, "Plots/supplemental_plots/2_1_expression_tau_patterns_fixed.pdf"), height = 2.5, width = 7)
plot_grid(plot_grid(cel_exp_plot, cbr_exp_plot, nrow = 2, rel_heights = c(0.8, 1)),
          plot_grid(cel_tau_plot, cbr_tau_plot, nrow = 2, rel_heights = c(0.8, 1)), nrow = 1)
dev.off()


################
# ehn-3 figure #
################

OrderForPlot <- readRDS(paste0(dir, "Objects/OrderForPlot.rds"))

# ehn-3, ztf-16, sdz-12
cel_tpm <- data.frame(t(TPMListBootstrapMean_term[["C.elegans"]][c("ehn-3", "ztf-16", "sdz-12"),]))
cel_tpm$cell_type <- rownames(cel_tpm)
cel_tpm$species <- "C. elegans"
cbr_tpm <- data.frame(t(TPMListBootstrapMean_term[["C.briggsae"]][c("ehn-3", "ztf-16", "sdz-12"),]))
cbr_tpm$cell_type <- rownames(cel_tpm)
cbr_tpm$species <- "C. briggsae"

tpm <- rbind(cel_tpm, cbr_tpm)
plot_tpm <- tpm %>% melt(value.name = "tpm", id = c("cell_type", "species"))

"CBG19502"

cel_bin <- data.frame(t(BinarizeBootstrapListMean_term[["C.elegans"]][c("ehn-3", "ztf-16", "sdz-12"),]))
cel_bin$cell_type <- rownames(cel_bin)
cel_bin$species <- "C. elegans"

cbr_bin <- data.frame(t(BinarizeBootstrapListMean_term[["C.briggsae"]][c("ehn-3", "ztf-16", "sdz-12"),]))
cbr_bin$cell_type <- rownames(cbr_bin)
cbr_bin$species <- "C. briggsae"

bin <- rbind(cel_bin, cbr_bin)
plot_bin <- bin %>% melt(value.name = "bin", id = c("cell_type", "species"))

plot_order <- (left_join(plot_tpm, plot_bin) %>%
  filter(bin > 10) %>%
  select(cell_type) %>%
  unique() %>% arrange(match(cell_type, OrderForPlot)))[,1]

# plot_tpm <- rbind(plot_tpm, data.frame(cell_type = (plot_tpm %>% group_by(cell_type, species) %>% summarise(sum_tpm = sum(tpm)))[, "cell_type", drop = TRUE],
#                            species = (plot_tpm %>% group_by(cell_type, species) %>% summarise(sum_tpm = sum(tpm)))[, "species", drop = TRUE],
#                            variable = "Sum",
#                            tpm = (plot_tpm %>% group_by(cell_type, species) %>% summarise(sum_tpm = sum(tpm)))[, "sum_tpm", drop = TRUE]))

pdf(paste0(dir, "Plots/supplemental_plots/ehn_3_orthologs_tpm.pdf"), height = 3, width = 10)
plot_tpm %>%
  filter(cell_type %in% plot_order) %>%
  arrange(match(cell_type, plot_order)) %>%
ggplot(aes(x = cell_type, y = paste0(species, "_", variable), fill = tpm + 1)) +
  geom_tile(color = "grey80") + 
  scale_fill_viridis(name = "TPM", option = "magma", limits = c(1, 9603), trans = "log2") +
  scale_x_discrete(limits = plot_order) +
  scale_y_discrete(limits = rev(c("C. elegans_ehn.3", "C. briggsae_ehn.3",
                              "C. elegans_ztf.16", "C. briggsae_ztf.16",
                              "C. elegans_sdz.12", "C. briggsae_sdz.12")),
                   labels = rev(c("C. elegans ehn-3", "C. briggsae ehn-3",
                              "C. elegans ztf-16", "C. briggsae ztf-16",
                              "C. elegans sdz-12", "C. briggsae sdz-12"))) +
  theme(legend.position = "right",
        rect = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.line = element_line(color="grey80", size=1),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()


#################
# acdh-6 figure #
#################

OrderForPlot <- readRDS(paste0(dir, "Objects/OrderForPlot.rds"))

cel_tpm <- data.frame(t(TPMListBootstrapMean_term[["C.elegans"]][c("acdh-5", "acdh-6"),]))
cel_tpm$cell_type <- rownames(cel_tpm)
cel_tpm$species <- "C. elegans"
cbr_tpm <- data.frame(t(TPMListBootstrapMean_term[["C.briggsae"]][c("acdh-6"),, drop = FALSE]))
cbr_tpm$cell_type <- rownames(cbr_tpm)
cbr_tpm$species <- "C. briggsae"
cbr_tpm$acdh.5 <- NA

tpm <- rbind(cel_tpm, cbr_tpm)
plot_tpm <- tpm %>% melt(value.name = "tpm", id = c("cell_type", "species"))

cel_bin <- data.frame(t(BinarizeBootstrapListMean_term[["C.elegans"]][c("acdh-5", "acdh-6"),]))
cel_bin$cell_type <- rownames(cel_bin)
cel_bin$species <- "C. elegans"

cbr_bin <- data.frame(t(BinarizeBootstrapListMean_term[["C.briggsae"]][c("acdh-6"),, drop = FALSE]))
cbr_bin$cell_type <- rownames(cbr_bin)
cbr_bin$species <- "C. briggsae"
cbr_bin$acdh.5 <- NA

bin <- rbind(cel_bin, cbr_bin)
plot_bin <- bin %>% melt(value.name = "bin", id = c("cell_type", "species"))

plot_order <- (left_join(plot_tpm, plot_bin) %>%
                 filter(bin > 10) %>%
                 select(cell_type) %>%
                 unique() %>% arrange(match(cell_type, OrderForPlot)))[,1]

# plot_tpm <- rbind(plot_tpm, data.frame(cell_type = (plot_tpm %>% group_by(cell_type, species) %>% summarise(sum_tpm = sum(tpm)))[, "cell_type", drop = TRUE],
#                            species = (plot_tpm %>% group_by(cell_type, species) %>% summarise(sum_tpm = sum(tpm)))[, "species", drop = TRUE],
#                            variable = "Sum",
#                            tpm = (plot_tpm %>% group_by(cell_type, species) %>% summarise(sum_tpm = sum(tpm)))[, "sum_tpm", drop = TRUE]))

pdf(paste0(dir, "Plots/supplemental_plots/acdh-5_orthologs_tpm.pdf"), height = 3, width = 10)
plot_tpm %>%
  filter(cell_type %in% plot_order) %>%
  arrange(match(cell_type, plot_order)) %>%
  ggplot(aes(x = cell_type, y = paste0(species, "_", variable), fill = tpm + 1)) +
  geom_tile(color = "grey80") + 
  scale_fill_viridis(name = "TPM", option = "magma", limits = c(1, 1000), trans = "log2") +
  scale_x_discrete(limits = plot_order) +
  scale_y_discrete(limits = rev(c("C. elegans_acdh.6", "C. briggsae_acdh.6",
                                  "C. elegans_acdh.5")),
                   labels = rev(c("C. elegans acdh-6", "C. briggsae acdh-6",
                                  "C. elegans acdh-5"))) +
  theme(legend.position = "right",
        rect = element_rect(fill = "transparent"),
        axis.text = element_text(size = 12),
        axis.line = element_line(color="grey80", size=1),
        axis.title = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()



# Load Orthogroups.txt
gene_duplications_og <- read.table("/kimdata/livinlrg/GenomeAnalysis/orthofinder_sans_uteleia/fastas/OrthoFinder/Results_Apr23/Gene_Duplication_Events/Duplications.tsv",
                                   sep="\t",
                                   header = TRUE)
gene_duplications_og <- gene_duplications_og[gene_duplications_og$Support > 0.5,]

# multiple hit
cur_gene <- "T07F8.1"
# no match

out <- apply(cel_dup[cel_dup$og_class2 %in% c("2:1 better", "2:1 poorer"),], 1, function(cur_row) {
  cur_gene <- unlist(cur_row["gene"])
  temp_1 <- gene_duplications_og[grepl(paste0("CELEG.", gff_list_mrna[["elegans"]][cur_gene, "short_name"], "(,|$)"), gene_duplications_og$Genes.1),]
  temp_2 <- gene_duplications_og[grepl(paste0("CELEG.", gff_list_mrna[["elegans"]][cur_gene, "short_name"], "(,|$)"), gene_duplications_og$Genes.2),]
  
  if(nrow(temp_1) == 1) {
    if(temp_1$Support < 0.5) {
      return("low support")
    } else if(temp_1['Type'] == "Terminal") {
      return("terminal")
    } else if(temp_1['Type'] == "Non-Terminal") {
      return(temp_1['Species.Tree.Node'])
    }
  } else if(nrow(temp_2) == 1) {
    if(temp_2$Support < 0.5) {
      return("low support")
    } else if(temp_2['Type'] == "Terminal") {
      return("terminal")
    } else if(temp_2['Type'] == "Non-Terminal") {
      return(temp_2['Species.Tree.Node'])
    }
  } else if (nrow(temp_1) == 0 & nrow(temp_2) == 0) {
    # print(cur_gene)
    return("No match found")
  } else {
    print(temp_1)
    print(temp_2)
    return("Multiple matches found")
  }
})

cel_dup$age <- NA
cel_dup[names(out),]$age <- unlist(out)

# terminal
cel_dup[c("T07F8.1", "ZK686.6", "K05D4.4", "Y17D7A.4", "cyp-33D3"),]$age <- "terminal"

# n4, n2, n0
table(cel_dup[cel_dup$og_class2 %in% c("2:1 poorer"),]$sum > 80, cel_dup[cel_dup$og_class2 %in% c("2:1 poorer"),]$age)

cel_dup$test_age <- cel_dup$age
cel_dup[which(cel_dup$test_age %in% c("N0", "N2", "N4")),]$test_age <- "old"
cel_dup[which(cel_dup$test_age %in% c("Multiple matches found", "No match found")),]$test_age <- NA

chisq.test(table(cel_dup[cel_dup$og_class2 %in% c("2:1 poorer"),]$sum > 80, cel_dup[cel_dup$og_class2 %in% c("2:1 poorer"),]$test_age))

# they are the same
ogs <- cel_dup[! is.na(cel_dup$age),]$og
ogs <- ogs[! is.na(ogs)]
for(cur_og in ogs) {
  cur_ages <- cel_dup[which(cel_dup$og == cur_og), "age"]
  if(cur_ages[1] != cur_ages[2]) {
    print(cur_og)
  }
}

out <- apply(cbr_dup[cbr_dup$og_class2 %in% c("2:1 better", "2:1 poorer"),], 1, function(cur_row) {
  cur_gene <- unlist(cur_row["joint_name"])
  temp_1 <- gene_duplications_og[grepl(paste0("_", gff_list_mrna[["briggsae"]][cur_gene, "short_name"], "(,|$)"), gene_duplications_og$Genes.1),]
  temp_2 <- gene_duplications_og[grepl(paste0("_", gff_list_mrna[["briggsae"]][cur_gene, "short_name"], "(,|$)"), gene_duplications_og$Genes.2),]
  
  if(nrow(temp_1) == 1) {
    if(temp_1$Support < 0.5) {
      return("low support")
    } else if(temp_1['Type'] == "Terminal") {
      return("terminal")
    } else if(temp_1['Type'] == "Non-Terminal") {
      return(temp_1['Species.Tree.Node'])
    }
  } else if(nrow(temp_2) == 1) {
    if(temp_2$Support < 0.5) {
      return("low support")
    } else if(temp_2['Type'] == "Terminal") {
      return("terminal")
    } else if(temp_2['Type'] == "Non-Terminal") {
      return(temp_2['Species.Tree.Node'])
    }
  } else if (nrow(temp_1) == 0 & nrow(temp_2) == 0) {
    # print(cur_gene)
    return("No match found")
  } else {
    print(temp_1)
    print(temp_2)
    return("Multiple matches found")
  }
})

cbr_dup$age <- NA
cbr_dup[names(out),]$age <- unlist(out)

#terminal
cbr_dup[c("CBG15349", "CBG00609", "Cbr-cwf-19L2",
          "CBG22991", "CBG22992"),]$age <- "terminal"


cbr_dup$test_age <- cbr_dup$age
cbr_dup[which(cbr_dup$test_age %in% c("N0", "N11", "N2", "N5", "N7", "N9")),]$test_age <- "old"
cbr_dup[which(cbr_dup$test_age %in% c("Multiple matches found", "No match found")),]$test_age <- NA

chisq.test(table(cbr_dup[cbr_dup$og_class2 %in% c("2:1 poorer"),]$sum > 80,
                 cbr_dup[cbr_dup$og_class2 %in% c("2:1 poorer"),]$test_age))