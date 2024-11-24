
library(monocle3)
library(Seurat)
library(ggplot2)
library(viridis)
library(reshape2)
library(ggrepel)
library(ggExtra) # For ggmarginal
library(dplyr)
library(cowplot) # for plot_grid
library(ggsankey)
library(tidyr)
library(purrr)

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"
source(paste0(dir, "../Scripts/WS290/CalcFunctions.R"))

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240903.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]
rm(cds)

ListOfCellTypes <- readRDS(paste0(dir, "Objects/ListOfCellTypes_Mean.rds"))

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"
source(paste0(dir, "../Scripts/WS290/CalcFunctions.R"))

gene_data <- readRDS(paste0(dir, "Objects/gene_data_joint_cell_bg.rds"))
cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))

gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))

TPMListBootstrapMean_term <- readRDS(paste0(dir, "Objects/TPMListBootstrapMean_CellCorrection.rds"))
BinarizeBootstrapListMean_term <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListMean_CellCorrection.rds"))
BinarizeBootstrapList_term <- readRDS(paste0(dir, "Objects/BinarizeBootstrapTimeBinsList_CellCorrection.rds"))

TPMListBootstrapMean_term_subset <- list()
TPMListBootstrapMean_term_subset[["C.elegans"]] <- TPMListBootstrapMean_term[["C.elegans"]][gene_data$gene,]
TPMListBootstrapMean_term_subset[["C.briggsae"]] <- TPMListBootstrapMean_term[["C.briggsae"]][gene_data$gene,]

BinarizeBootstrapListMean_term_subset <- list()
BinarizeBootstrapListMean_term_subset[["C.elegans"]] <- BinarizeBootstrapListMean_term[["C.elegans"]][gene_data$gene,]
BinarizeBootstrapListMean_term_subset[["C.briggsae"]] <- BinarizeBootstrapListMean_term[["C.briggsae"]][gene_data$gene,]

TPMListBootstrap_pro <- readRDS(paste0(dir, "Objects/TPMListBootstrapPro_CellCorrection.rds"))
BinarizeBootstrapListPro <- readRDS(paste0(dir, "Objects/BinarizeBootstrapListPro_CellCorrection.rds"))

TPMListBootstrap_pro_subset <- list()
TPMListBootstrap_pro_subset[["C.elegans"]] <- TPMListBootstrap_pro[["C.elegans"]][gene_data$gene,]
TPMListBootstrap_pro_subset[["C.briggsae"]] <- TPMListBootstrap_pro[["C.briggsae"]][gene_data$gene,]

term_markers <- readRDS(paste0(dir, "Objects/TermMarker_list_filt_metadata_cell_bg.rds"))
pro_markers <- readRDS(paste0(dir, "Objects/ProMarker_list_filt_metadata_cell_bg.rds"))
joint_markers <- readRDS(paste0(dir, "Objects/joint_markers_cell_bg.rds"))

deg <- readRDS(paste0(dir, "Objects/DEG_df_filt_metadata_cell_bg.rds"))

# for geom_text_repel
gene_data$plot_name <- NA
gene_data[c("rab-7", "vab-15", "lin-59", "pha-4", "ceh-53", "T04A6.1"),]$plot_name <- c("rab-7", "vab-15", "lin-59", "pha-4", "ceh-53", "T04A6.1")

#############
# Tau plots #
#############

gene_data$max_tpm <- apply(gene_data[,c("max_tpm_term", "max_tpm_pro")], 1, max)
gene_data$cel_max_tpm <- apply(gene_data[,c("cel_max_tpm_pro", "cel_max_tpm_term")], 1, max)
gene_data$cbr_max_tpm <- apply(gene_data[,c("cbr_max_tpm_pro", "cbr_max_tpm_term")], 1, max)

pdf(paste0(dir, "Plots/gene_figure_plots/term_tau.pdf"), width = 7.5, height = 7.5)
ggMarginal(gene_data %>%
             filter(max_tpm_term > 80) %>% 
             arrange(max_tpm_term) %>%
             ggplot(., aes(x = cel_tau_median_term, y = cbr_tau_median_term, color = max_tpm_term, label = plot_name)) +
             geom_abline(intercept = 0.3, slope = 1, linetype = "dashed", size = 1.5) +
             geom_abline(intercept = -0.3, slope = 1, linetype = "dashed", size = 1.5) +
             geom_point(alpha = 0.45, size = 1.5, stroke = 0.5) +
             # geom_text_repel(color = "black", max.overlaps = Inf) +
             scale_x_continuous(name = expression(paste(italic("C. elegans")," Tau")), limits = c(0, 1)) +
             scale_y_continuous(name = expression(paste(italic("C. briggsae"), "Tau")), limits = c(0, 1)) +
             scale_color_viridis(name = "Expression \n (TPM)", trans = "log2") +
             theme(legend.position = c(0.15, 0.8),
                   axis.title.x = element_text(color = "#009E73", size = 24),
                   axis.title.y = element_text(color = "#56B4E9", size = 24),
                   axis.text = element_text(size = 18),
                   legend.text = element_text(size = 12),
                   legend.title = element_text(size = 18),
                   rect = element_rect(fill = "transparent"),
                   axis.line = element_line(color="grey80", size=1),
                   panel.grid.major = element_line(color="grey80", size=0.5),
                   panel.grid.minor = element_line(color="grey80", size=0.1),
                   panel.background = element_rect(fill='transparent'), #transparent panel bg
                   plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
             ),
           type = "histogram", col = "black", alpha = 0.8, bins = 30,
           xparams = list(fill = "#009E73"),
           yparams = list(fill = "#56B4E9"))
dev.off()

pdf(paste0(dir, "Plots/gene_figure_plots/pro_tau.pdf"), width = 7.5, height = 7.5)
ggMarginal(gene_data %>%
             filter(max_tpm_pro > 80) %>% 
             arrange(max_tpm_pro) %>%
             ggplot(., aes(x = cel_tau_median_pro, y = cbr_tau_median_pro, color = max_tpm_pro)) +
             geom_abline(intercept = 0.3, slope = 1, linetype = "dashed", size = 1.5) +
             geom_abline(intercept = -0.3, slope = 1, linetype = "dashed", size = 1.5) +
             geom_point(alpha = 0.45, size = 1.5, stroke = 0.5) +
             scale_x_continuous(name = expression(paste(italic("C. elegans")," Tau")), limits = c(0, 1)) +
             scale_y_continuous(name = expression(paste(italic("C. briggsae"), "Tau")), limits = c(0, 1)) +
             scale_color_viridis(name = "Expression \n (TPM)", trans = "log2") +
             theme(legend.position = c(0.15, 0.8),
                   axis.title.x = element_text(color = "#009E73", size = 24),
                   axis.title.y = element_text(color = "#56B4E9", size = 24),
                   axis.text = element_text(size = 18),
                   legend.text = element_text(size = 12),
                   legend.title = element_text(size = 18),
                   rect = element_rect(fill = "transparent"),
                   axis.line = element_line(color="grey80", size=1),
                   panel.grid.major = element_line(color="grey80", size=0.5),
                   panel.grid.minor = element_line(color="grey80", size=0.1),
                   panel.background = element_rect(fill='transparent'), #transparent panel bg
                   plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
             ),
           type = "histogram", col = "black", alpha = 0.8, bins = 30,
           xparams = list(fill = "#009E73"),
           yparams = list(fill = "#56B4E9"))
dev.off()

gene_data_more_genes <- gene_data
gene_data_more_genes[c("cep-1", "gld-1"),]$plot_name <- c("cep-1", "gld-1")

pdf(paste0(dir, "Plots/gene_figure_plots/joint_tau.pdf"), width = 7.5, height = 7.5)
ggMarginal(gene_data_more_genes %>%
             # filter(max_tpm_term > 80 | max_tpm_pro > 80) %>% 
             filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>% 
             arrange(max_tpm) %>%
             ggplot(., aes(x = cel_tau_median_joint, y = cbr_tau_median_joint, color = max_tpm, label = plot_name), fill = NA) +
             geom_abline(intercept = 0.3, slope = 1, linetype = "dashed", size = 1) +
             geom_abline(intercept = -0.3, slope = 1, linetype = "dashed", size = 1) +
             geom_point(alpha = 0.45, size = 1.5, stroke = 0.5) +
             geom_text_repel(color = "black", size = 6, max.overlaps = Inf) +
             geom_point(data = gene_data_more_genes[! is.na(gene_data_more_genes$plot_name),],
                        aes(x = cel_tau_median_joint, y = cbr_tau_median_joint, fill = max_tpm), pch = 21, color = "black", show.legend = FALSE) +
             scale_x_continuous(name = expression(paste(italic("C. elegans")," Tau")), limits = c(0, 1)) +
             scale_y_continuous(name = expression(paste(italic("C. briggsae"), "Tau")), limits = c(0, 1)) +
             scale_color_viridis(name = "Expression\n(TPM)", trans = "log2", limits = c(80, 65536), oob = scales::squish) +
             scale_fill_viridis(name = "Expression\n(TPM)", trans = "log2", limits = c(80, 65536), oob = scales::squish) +
             theme(legend.position = c(0.15, 0.8),
                   axis.title.x = element_text(color = "#009E73", size = 24),
                   axis.title.y = element_text(color = "#56B4E9", size = 24),
                   axis.text = element_text(size = 18),
                   legend.text = element_text(size = 12),
                   legend.title = element_text(size = 18),
                   rect = element_rect(fill = "transparent"),
                   axis.line = element_line(color="grey80", size=1),
                   panel.grid.major = element_line(color="grey80", size=0.5),
                   panel.grid.minor = element_line(color="grey80", size=0.1),
                   panel.background = element_rect(fill='transparent'), #transparent panel bg
                   plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
             ),
           type = "histogram", col = "black", alpha = 0.8, bins = 30,
           xparams = list(fill = "#009E73"),
           yparams = list(fill = "#56B4E9"))
dev.off()

write.csv(gene_data %>%
            filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>%
            filter(cel_tau_median_joint - cbr_tau_median_joint > 0.3 |
                     cel_tau_median_joint - cbr_tau_median_joint < -0.3),
          paste0(dir, "Tables/tau_difference.csv"), row.names = FALSE, quote = FALSE)

####################
# Tau to jsd plots #
####################

gene_data$mean_tau_joint <- (gene_data$cel_tau_median_joint + gene_data$cbr_tau_median_joint) / 2

pdf(paste0(dir, "Plots/gene_figure_plots/joint_tau_jsd.pdf"), width = 5, height = 5)
gene_data %>%
  filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>% 
  arrange(max_tpm) %>%
ggplot(., aes(x = mean_tau_joint, y = jsd_median_joint, color = max_tpm, label = plot_name)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 0.4, linetype = "dashed", size = 1, color = "red") +
  geom_vline(xintercept = 0.7, linetype = "dashed", size = 1, color = "red") +
  geom_hline(yintercept = 0.45, linetype = "dashed", size = 1, color = "red") +
  geom_hline(yintercept = 0.55, linetype = "dashed", size = 1, color = "red") +
  geom_text_repel(color = "black", size = 6, max.overlaps = Inf) +
  geom_point(data = gene_data[! is.na(gene_data$plot_name),],
             aes(x = mean_tau_joint, y = jsd_median_joint, fill = max_tpm), pch = 21, color = "black", show.legend = FALSE) +
  scale_x_continuous(name = "Tau (Mean between species)", breaks = c(0, 0.2, 0.4, 0.6, 0.80, 1.0)) +
  scale_y_continuous(name = "Joint gene JSD") +
  scale_color_viridis(name = "Expression\n(TPM)", trans = "log2", limits = c(80, 65536), oob = scales::squish) +
  scale_fill_viridis(name = "Expression\n(TPM)", trans = "log2", limits = c(80, 65536), oob = scales::squish) +
  theme(#legend.title= element_blank(),
        legend.position = c(0.15, 0.8),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text = element_text(size = 12),
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major = element_line(color="grey80", size=0.5),
        panel.grid.minor = element_line(color="grey80", size=0.1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  )
dev.off()

gene_data$mean_tau_term <- (gene_data$cel_tau_median_term + gene_data$cbr_tau_median_joint) / 2

pdf(paste0(dir, "Plots/gene_figure_plots/term_tau_jsd.pdf"), width = 5, height = 5)
gene_data %>%
  filter(max_tpm_term > 80) %>% 
  arrange(max_tpm_term) %>%
  ggplot(., aes(x = mean_tau_term, y = jsd_median_term, color = max_tpm_term)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 0.4, linetype = "dashed", size = 1, color = "red") +
  geom_vline(xintercept = 0.7, linetype = "dashed", size = 1, color = "red") +
  geom_hline(yintercept = 0.45, linetype = "dashed", size = 1, color = "red") +
  geom_hline(yintercept = 0.55, linetype = "dashed", size = 1, color = "red") +
  scale_x_continuous(name = "Tau (Mean between species)", breaks = c(0, 0.2, 0.4, 0.6, 0.80, 1.0)) +
  scale_y_continuous(name = "Term. gene JSD") +
  scale_color_viridis(name = "Expression\n(TPM)", trans = "log2", limits = c(80, 65536), oob = scales::squish) +
  # scale_color_manual(name = "Markers", breaks = c("C. elegans marker", "C. briggsae marker", "Conserved marker", "Neither"),
  #                      labels = c("C. elegans marker", "C. briggsae marker", "Conserved marker", "Neither"),
  #                      values = c("#009E73", "#56B4E9", "#CC79A7", "black")) +
  theme(#legend.title= element_blank(),
    legend.position = c(0.15, 0.8),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 12),
    rect = element_rect(fill = "transparent"),
    axis.line = element_line(color="grey80", size=1),
    panel.grid.major = element_line(color="grey80", size=0.5),
    panel.grid.minor = element_line(color="grey80", size=0.1),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  )
dev.off()

gene_data$mean_tau_pro <- (gene_data$cel_tau_median_pro + gene_data$cbr_tau_median_pro) / 2

pdf(paste0(dir, "Plots/gene_figure_plots/pro_tau_jsd.pdf"), width = 5, height = 5)
gene_data %>%
  filter(max_tpm_pro > 80) %>% 
  arrange(max_tpm_pro) %>%
  ggplot(., aes(x = mean_tau_pro, y = jsd_median_pro, color = max_tpm_pro)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = 0.4, linetype = "dashed", size = 1, color = "red") +
  geom_vline(xintercept = 0.7, linetype = "dashed", size = 1, color = "red") +
  geom_hline(yintercept = 0.45, linetype = "dashed", size = 1, color = "red") +
  geom_hline(yintercept = 0.55, linetype = "dashed", size = 1, color = "red") +
  scale_x_continuous(name = "Tau (Mean between species)", breaks = c(0, 0.2, 0.4, 0.6, 0.80, 1.0)) +
  scale_y_continuous(name = "Pro. gene JSD") +
  scale_color_viridis(name = "Expression\n(TPM)", trans = "log2", limits = c(80, 65536), oob = scales::squish) +
  theme(#legend.title= element_blank(),
    legend.position = c(0.15, 0.8),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text = element_text(size = 12),
    rect = element_rect(fill = "transparent"),
    axis.line = element_line(color="grey80", size=1),
    panel.grid.major = element_line(color="grey80", size=0.5),
    panel.grid.minor = element_line(color="grey80", size=0.1),
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  )
dev.off()

############
# jsd hist #
############

pdf(paste0(dir, "Plots/gene_figure_plots/gene_jsd_histogram.pdf"), width = 4.5, height = 2)
gene_data %>%
  filter(max_tpm_term > 80 | max_tpm_pro > 80) %>% 
  arrange(max_tpm) %>%
ggplot() +
  geom_density(aes(x = jsd_median_joint)) +
  geom_rect(data = gene_data[gene_data$gene == "rab-7",],
            aes(xmin = jsd_lower_joint, xmax = jsd_upper_joint, ymin = -Inf, ymax = Inf),
            color = "red3", fill = "red", alpha = 0.25, size = 0.1) +
  geom_rect(data = gene_data[gene_data$gene == "vab-15",],
            aes(xmin = jsd_lower_joint, xmax = jsd_upper_joint, ymin = -Inf, ymax = Inf),
            color = "red3", fill = "red", alpha = 0.25, size = 0.1) +
  geom_rect(data = gene_data[gene_data$gene == "pha-4",],
            aes(xmin = jsd_lower_joint, xmax = jsd_upper_joint, ymin = -Inf, ymax = Inf),
            color = "red3", fill = "red", alpha = 0.25, size = 0.1) +
  geom_rect(data = gene_data[gene_data$gene == "ceh-53",],
            aes(xmin = jsd_lower_joint, xmax = jsd_upper_joint, ymin = -Inf, ymax = Inf),
            color = "red3", fill = "red", alpha = 0.25, size = 0.1) +
  geom_rect(data = gene_data[gene_data$gene == "T04A6.1",],
            aes(xmin = jsd_lower_joint, xmax = jsd_upper_joint, ymin = -Inf, ymax = Inf),
            color = "red3", fill = "red", alpha = 0.25, size = 0.1) +
  geom_vline(aes(xintercept = gene_data[gene_data$gene == "rab-7", "jsd_median_joint"]), color = "red") +
  geom_vline(aes(xintercept = gene_data[gene_data$gene == "vab-15", "jsd_median_joint"]), color = "red") +
  geom_vline(aes(xintercept = gene_data[gene_data$gene == "pha-4", "jsd_median_joint"]), color = "red") +
  geom_vline(aes(xintercept = gene_data[gene_data$gene == "ceh-53", "jsd_median_joint"]), color = "red") +
  geom_vline(aes(xintercept = gene_data[gene_data$gene == "T04A6.1", "jsd_median_joint"]), color = "red") +
  scale_x_continuous(name = "Jensen-Shannon Distance", limits = c(0, 1)) +
  theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        axis.line = element_line(color="grey80", size=1),
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),
        legend.position = "none")
dev.off()


####################
# jsd btwn species #
####################

pdf(paste0(dir, "Plots/gene_figure_plots/term_to_pro_jsd.pdf"), width = 5, height = 5)
ggMarginal(gene_data %>%
             filter(max_tpm_term > 80 & max_tpm_pro > 80) %>% 
             arrange(abs(log2(max_tpm_term/max_tpm_pro))) %>%
             ggplot(aes(x = jsd_median_term, y = jsd_median_pro, color = log2((max_tpm_term/max_tpm_pro)))) +
             geom_point(alpha = 0.65, size = 0.75) +
             scale_x_continuous(name = "Terminal gene JSD", limits = c(0, 1)) +
             scale_y_continuous(name = "Progenitor gene JSD", limits = c(0, 1)) +
             scale_color_gradient2(name = "log2(Term./Pro.)", low = '#000099', mid = '#DCDCDC', high = '#FF0000', midpoint = 0,
                                   limits = c(-5, 5),
                                   breaks = c(-5, -2.5, 0, 2.5, 5),
                                   labels = c("<-5", "-2.5", "0", "2.5", ">5"),
                                   oob = scales::squish) +
             theme(legend.position = c(.65, .12),
                   legend.direction = "horizontal",
                   # legend.title = element_blank(),
                   # axis.title.x = element_text(color = "#009E73", size = 16),
                   # axis.title.y = element_text(color = "#56B4E9", size = 16),
                   axis.text = element_text(size = 12),
                   rect = element_rect(fill = "transparent"),
                   axis.line = element_line(color="grey80", size=1),
                   panel.grid.major = element_line(color="grey80", size=0.5),
                   panel.grid.minor = element_line(color="grey80", size=0.1),
                   panel.background = element_rect(fill='transparent'), #transparent panel bg
                   plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
             ),
           type = "histogram", col = "black", alpha = 0.8,
           xparams = list(bins = 30, fill = "grey80"),
           yparams = list(bins = 30, fill = "grey80"))
dev.off()

####################
# lethal phenotype #
####################

gene_data$lethal <- "None annotated"
gene_data[which(gene_data$RNAi_or_Allele),]$lethal <- "No lethal phenotype"
gene_data[which(gene_data$lethal_direct_inferred),]$lethal <- "Non-specific lethal"
gene_data[which(gene_data$larval_lethal_direct_inferred),]$lethal <- "Larval specific lethal"
gene_data[which(gene_data$embryonic_lethal_direct_inferred),]$lethal <- "Embryonic lethal"

wilcox_list <- list()
for(cat_outer in unique(gene_data$lethal)) {
  for(cat_inner in unique(gene_data$lethal)) {
    print(paste0(cat_outer, "_", cat_inner))
    if(cat_outer == cat_inner) next
    if(paste0(cat_inner, "_", cat_outer) %in% names(wilcox_list)) next
    wilcox_list[[paste0(cat_outer, "_", cat_inner)]] <- wilcox.test(gene_data[which((gene_data$max_tpm_term > 80 & gene_data$max_tpm_pro > 80) & gene_data$lethal == cat_outer),]$jsd_median_joint,
                                                                    gene_data[which((gene_data$max_tpm_term > 80 & gene_data$max_tpm_pro > 80) & gene_data$lethal == cat_inner),]$jsd_median_joint)$p.value
  }
}
p.adjust(unlist(wilcox_list), method = "bonferroni")

pdf(paste0(dir, "Plots/gene_figure_plots/lethal_jsd.pdf"), height = 5, width = 5)
gene_data %>% 
  filter(max_tpm_term > 80 | max_tpm_pro > 80) %>% 
ggplot(aes(x = forcats::fct_reorder(lethal, jsd_median_joint),
           y = jsd_median_joint)) +
  geom_violin() +
  geom_boxplot(width = 0.25) +
  stat_summary(fun.data = n_fun, geom = "text", angle = 0) +
  scale_x_discrete(name = "RNAi + Allele phenotype") +
  scale_y_continuous(name = "Gene distance (jsd)") +
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

##################
# Gene dot plots #
##################
print("Plotting dot plot")

for(gene in c("rab-7", "vab-15", "pha-4", "ceh-53", "T04A6.1")) {
  TPMList <- TPMListBootstrapMean_term_subset
  BinarizeList <- BinarizeBootstrapListMean_term_subset
  
  exp_filter <- (BinarizeList[["C.elegans"]][gene,] > 0 | BinarizeList[["C.briggsae"]][gene,] > 0) & (TPMList[["C.elegans"]][gene,] > 80 | TPMList[["C.briggsae"]][gene,] > 80)

  temp <- lapply(TPMList, function(x) {
    df <- data.frame(expression = t(x[gene,, drop = FALSE] + 1))
    df$cell_type <- rownames(df)
    
    colnames(df) <- c("expression", "cell_type")
    df <- left_join(df, cell_data[,c("cell_type", "cell_class")], by = "cell_type")
    return(df[exp_filter,])
  })
  
  pdf(paste0(dir, "Plots/gene_figure_plots/gene_plots/", gene, "_tpm.pdf"), height = 3.5, width = 3.5)
  print(ggplot() +
          geom_point(aes(x = temp[["C.elegans"]]$expression,
                         y = temp[["C.briggsae"]]$expression,
                         color = temp[["C.elegans"]]$cell_class),
                     shape = 21) +
          geom_text_repel(aes(x = temp[["C.elegans"]]$expression,
                              y = temp[["C.briggsae"]]$expression,
                              label = temp[["C.elegans"]]$cell_type,
                              color = temp[["C.elegans"]]$cell_class),
                          show.legend = FALSE) +
          scale_x_continuous(name = paste("C. elegans ", gene," TPM"),
                             trans = "log2",
                             limits = c(min(temp[["C.elegans"]]$expression,
                                            temp[["C.briggsae"]]$expression) * 0.9,
                                        #oob = scales::squish,
                                        max(temp[["C.elegans"]]$expression,
                                            temp[["C.briggsae"]]$expression) * 1.1)) +
          scale_y_continuous(name = paste("C. briggsae  ", gene, " TPM"),
                             trans = "log2",
                             limits = c(min(temp[["C.elegans"]]$expression,
                                            temp[["C.briggsae"]]$expression) * 0.9,
                                        #oob = scales::squish,
                                        max(temp[["C.elegans"]]$expression,
                                            temp[["C.briggsae"]]$expression) * 1.1)) +
          scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                        'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                        'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1")) +
          theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
                plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
                axis.line = element_line(color="grey80", size=1),
                panel.grid = element_line(color="grey80", size=0.25),
                legend.title = element_blank(),
                legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
                legend.key = element_rect(colour = "transparent", fill = "transparent"),
                legend.position = "none",
                axis.text = element_text(size = 12),
                axis.title.x = element_text(color = "#009E73", size = 16),
                axis.title.y = element_text(color = "#56B4E9", size = 16)))
  dev.off()
  
  pdf(paste0(dir, "Plots/gene_figure_plots/gene_plots/", gene, "_tpm_legend.pdf"), height = 7.5, width = 8)
  print(ggplot() +
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
    scale_color_manual(values = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                                  'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                                  'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1")) +
    guides(color = guide_legend(ncol = 3, byrow=TRUE)) +
    theme(panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          axis.line = element_line(color="grey80", size=1),
          panel.grid = element_line(color="grey80", size=0.25),
          legend.title = element_blank(),
          legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
          legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.position = "top"))
  dev.off()
}

######################
# Supplemental plots #
######################
# supplemental plot comparing gene expression level and jsd_median with colors segmeting into tau categories
pdf(paste0(dir, "Plots/gene_figure_plots/jsd_expression.pdf"), width = 6, height = 5)
ggMarginal(
gene_data %>% filter(max_tpm > 80 & ! is.na(mean_tau_joint)) %>%
  mutate(tau_cat = ifelse(gene_data[gene_data$max_tpm > 80 & ! is.na(gene_data$mean_tau_joint),]$mean_tau_joint <= 0.4,
                          "Broad", ifelse(gene_data[gene_data$max_tpm > 80,]$mean_tau_joint < 0.6, "Specific", "Patterned"))) %>%
  filter(! is.na(tau_cat)) %>%
  ggplot(aes(x = max_tpm, y = jsd_median_joint,
             color = tau_cat)) +
  geom_point(alpha = 0.65, size = 0.25) +
  scale_y_continuous(name = "Jensen-Shannon distance", limits = c(0, 1)) +
  scale_x_continuous(name = "Maximum expression (TPM)", trans = "log2") +
  scale_color_manual(name = "Tau", values = c("#4682B4", "#FFA500", "#B22222"), limits = c("Broad", "Patterned", "Specific")) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
theme(legend.position = "left",
      legend.background = element_blank(),
      legend.box.background = element_blank(),
      legend.key = element_blank(),
      axis.text = element_text(size = 12),
      rect = element_rect(fill = "transparent"),
      axis.line = element_line(color="grey80", size=1),
      panel.grid.major = element_line(color="grey80", size=0.5),
      panel.grid.minor = element_line(color="grey80", size=0.1),
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
), groupColour = TRUE, groupFill = TRUE)
dev.off()

pdf(paste0(dir, "Plots/gene_figure_plots/jsd_synteny.pdf"), width = 6, height = 5)
gene_data %>% filter(cel_max_tpm > 80 & cbr_max_tpm > 80) %>%
  ggplot(aes(x = syntenic, y = jsd_median_joint)) +
  geom_boxplot() +
  scale_y_continuous(name = "Jensen-Shannon distance", limits = c(0, 1)) +
  # scale_x_continuous(name = "Maximum expression (TPM)", trans = "log2") +
  # scale_color_manual(name = "Tau", values = c("#4682B4", "#FFA500", "#B22222"), limits = c("Broad", "Patterned", "Specific")) +
  # guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.position = "right",
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        axis.text = element_text(size = 12),
        rect = element_rect(fill = "transparent"),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major = element_line(color="grey80", size=0.5),
        panel.grid.minor = element_line(color="grey80", size=0.1),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
  )
dev.off()

###############
# Sankey plot #
###############

gene_data_sankey <- data.frame(matrix(nrow = dim(gene_data)[1], ncol = 13))
colnames(gene_data_sankey) <- c("gene",
                                "joint_detected", "progenitor_detected", "terminal_detected",
                                "joint_tau", "joint_conservation",
                                "terminal_tau", "terminal_conservation",
                                "progenitor_tau", "progenitor_conservation",
                                "joint_just_conservation", "terminal_just_conservation", "progenitor_just_conservation")
gene_data_sankey$gene <- gene_data$gene

exp_cutoff = 80
jsd_cutoff = c(0.45, 0.55)
names(jsd_cutoff) <- c("lower", "upper")
tau_cutoff = c(0.4, 0.7)
names(tau_cutoff) <- c("lower", "upper")

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

gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                              (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Broad conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                              (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Patterned conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                              (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Specific conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint <= tau_cutoff["lower"] &
                                                              (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Broad diverged"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint > tau_cutoff["lower"] & gene_data$mean_tau_joint <= tau_cutoff["upper"] &
                                                              (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Patterned diverged"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_joint > tau_cutoff["upper"] &
                                                              (gene_data$cel_max_tpm >= exp_cutoff | gene_data$cbr_max_tpm >= exp_cutoff)),]$gene,]$joint_conservation <- "Specific diverged"

# just conservation
gene_data$joint_just_conservation <- "Neutral"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_joint <= jsd_cutoff["lower"]),]$gene,]$joint_just_conservation <- "Conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_joint >= jsd_cutoff["upper"]),]$gene,]$joint_just_conservation <- "Diverged"


# terminal conservation
gene_data_sankey[gene_data_sankey$terminal_detected %in%
                   c("Both", "C. elegans", "C. briggsae"),]$terminal_conservation <- paste0(gene_data_sankey[gene_data_sankey$terminal_detected %in%
                                                                                                               c("Both", "C. elegans", "C. briggsae"),]$terminal_tau, " ", "neutral")

gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                   (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Broad conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                   (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Patterned conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"] & gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                   (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Specific conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term <= tau_cutoff["lower"] &
                                                                   (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Broad diverged"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term > tau_cutoff["lower"] & gene_data$mean_tau_term <= tau_cutoff["upper"] &
                                                                   (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Patterned diverged"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"] & gene_data$mean_tau_term > tau_cutoff["upper"] &
                                                                   (gene_data$cel_max_tpm_term >= exp_cutoff | gene_data$cbr_max_tpm_term >= exp_cutoff)),]$gene,]$terminal_conservation <- "Specific diverged"

# just conservation
gene_data$terminal_just_conservation <- "Neutral"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term <= jsd_cutoff["lower"]),]$gene,]$terminal_just_conservation <- "Conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_term >= jsd_cutoff["upper"]),]$gene,]$terminal_just_conservation <- "Diverged"

# progenitor conservation
gene_data_sankey[gene_data_sankey$progenitor_detected %in%
                   c("Both", "C. elegans", "C. briggsae"),]$progenitor_conservation <- paste0(gene_data_sankey[gene_data_sankey$progenitor_detected %in%
                                                                                                                 c("Both", "C. elegans", "C. briggsae"),]$progenitor_tau, " ", "neutral")

gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                   (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Broad conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                   (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Patterned conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"] & gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                   (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Specific conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro <= tau_cutoff["lower"] &
                                                                   (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Broad diverged"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro > tau_cutoff["lower"] & gene_data$mean_tau_pro <= tau_cutoff["upper"] &
                                                                   (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Patterned diverged"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"] & gene_data$mean_tau_pro > tau_cutoff["upper"] &
                                                                   (gene_data$cel_max_tpm_pro >= exp_cutoff | gene_data$cbr_max_tpm_pro >= exp_cutoff)),]$gene,]$progenitor_conservation <- "Specific diverged"

# just conservation
gene_data$progenitor_just_conservation <- "Neutral"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro <= jsd_cutoff["lower"]),]$gene,]$progenitor_just_conservation <- "Conserved"
gene_data_sankey[gene_data_sankey$gene %in% gene_data[which(gene_data$jsd_median_pro >= jsd_cutoff["upper"]),]$gene,]$progenitor_just_conservation <- "Diverged"

gene_data_sankey <- left_join(gene_data_sankey, gene_data[,c("gene", "WormCat.1", "WormCat.2", "WormCat.3")])
rownames(gene_data_sankey) <- gene_data_sankey$gene

gene_data_sankey[which(gene_data_sankey$joint_conservation == "NA neutral"),]$joint_conservation <- NA
gene_data_sankey[which(gene_data_sankey$terminal_conservation == "NA neutral"),]$terminal_conservation <- NA
gene_data_sankey[which(gene_data_sankey$progenitor_conservation == "NA neutral"),]$progenitor_conservation <- NA

saveRDS(gene_data_sankey, paste0(dir, "Objects/gene_data_sankey.rds"))
gene_data_sankey <- readRDS(paste0(dir, "Objects/gene_data_sankey.rds"))

##########
# Plot sankey
##########

sankey_df <- gene_data_sankey %>% make_long(progenitor_conservation, progenitor_tau, progenitor_detected, terminal_detected, terminal_tau, terminal_conservation)
sankey_df <- left_join(sankey_df,
                       sankey_df %>%
                         group_by(x, node) %>%
                         tally())

sankey_df$label <- sankey_df$node
sankey_df[sankey_df$node %in% c("Broad conserved", "Specific conserved", "Patterned conserved"),]$label <- "Conserved"
sankey_df[sankey_df$node %in% c("Broad diverged", "Specific diverged", "Patterned diverged"),]$label <- "Diverged"

pdf(paste0(dir, "Plots/gene_figure_plots/Sankey.pdf"), height = 20, width = 30)
sankey_df %>% drop_na(node) %>%
  ggplot(aes(x = x,
             next_x = next_x,
             node = node,
             next_node = next_node,
             fill = factor(node),
             label = paste0(label,"\nn = ", n)),
         width = 5) +
  geom_sankey(flow.alpha = 0.5,
              node.color = 1,
              width = 0.25,
              node.size = 1,
              color = "grey40",
              type = "sankey") +
  geom_sankey_label(size = 10,
                    node.width = 4,
                    node.size = 1) +
  scale_x_discrete(labels = c("Progenitor\nconservation", "Progenitor Tau", "Detected in\nprogenitor cells", "Detected in\nterminal cells", "Terminal Tau", "Terminal\nconservation")) +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 36),
        axis.title = element_blank())
dev.off()

sankey_df <- gene_data_sankey %>% make_long(joint_detected, joint_tau, joint_conservation, terminal_conservation, progenitor_conservation)
sankey_df <- left_join(sankey_df,
                       sankey_df %>%
                         group_by(x, node) %>%
                         tally())

sankey_df$label <- sankey_df$node
# sankey_df[sankey_df$node %in% c("Broad conserved", "Specific conserved", "Patterned conserved"),]$label <- "Conserved"
# sankey_df[sankey_df$node %in% c("Broad diverged", "Specific diverged", "Patterned diverged"),]$label <- "Diverged"

pdf(paste0(dir, "Plots/gene_figure_plots/Sankey.pdf"), height = 20, width = 30)
sankey_df %>% drop_na(node) %>%
  ggplot(aes(x = x,
             next_x = next_x,
             node = node,
             next_node = next_node,
             fill = factor(node),
             label = paste0(label,"\nn = ", n)),
         width = 5) +
  geom_sankey(flow.alpha = 0.5,
              node.color = 1,
              width = 0.25,
              node.size = 1,
              # color = "grey40",
              color = NA,
              type = "sankey") +
  geom_sankey_label(size = 10,
                    node.width = 4,
                    node.size = 1) +
  scale_x_discrete(labels = c("Detected in", "Tau", "Conservation", "Terminal\nconservation", "Progenitor\nconservation")) +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 36),
        axis.title = element_blank())
dev.off()

###########################
# Wormcat Enrichment Plot #
###########################

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


background_gene_set <- gene_data_sankey[gene_data_sankey$terminal_detected %in% c("Both", "C. elegans", "C. briggsae"),]$gene
Cat_WormCat_term <- list()
for(cat in c("Broad conserved", "Patterned conserved", "Specific conserved", "Specific diverged", "Patterned diverged", "Broad diverged")) {
  print(cat)
  rgs <- gene_data_sankey[which(gene_data_sankey$terminal_conservation == cat),]$gene
  
  Cat_WormCat_term[[cat]] <- wormcat_enrichment(rgs, background_gene_set, gff_list_mrna[["elegans"]])
  
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    Cat_WormCat_term[[cat]][[tier]]$type <- cat
  }
}

background_gene_set <- gene_data_sankey[gene_data_sankey$progenitor_detected %in% c("Both", "C. elegans", "C. briggsae"),]$gene
Cat_WormCat_pro <- list()
for(cat in c("Broad conserved", "Patterned conserved", "Specific conserved", "Specific diverged", "Patterned diverged", "Broad diverged")) {
  print(cat)
  rgs <- gene_data_sankey[which(gene_data_sankey$progenitor_conservation == cat),]$gene
  
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

pdf(paste0(dir, "Plots/gene_figure_plots/WormCatEnrichment.pdf"), height = 3.75, width = 4)
plot_grid(pro_wormcat_plot, term_wormcat_plot, align = "v", nrow = 2, rel_heights = c(0.6, 1), axis = 'lr')
dev.off()

n_fun <- function(x) {
  return(data.frame(y = 0, label = length(x)))
}
  
pdf(paste0(dir, "Plots/gene_figure_plots/WormCat_jsd.pdf"), height = 5, width = 8)
gene_data %>%
  filter(max_tpm > 80) %>%
  arrange(jsd_median_joint) %>%
  filter(! WormCat.1 %in% c("NA", "Pseudogene", "Major sperm protein") & ! is.na(WormCat.1)) %>%
  ggplot(aes(x = forcats::fct_reorder(WormCat.1, jsd_median_joint), y = jsd_median_joint)) +
  geom_boxplot() +
  # geom_pointrange(aes(ymin = jsd_lower_joint, ymax = jsd_upper_joint,
  #                     y = jsd_median_joint, x = forcats::fct_reorder(WormCat.1, jsd_median_joint), group = WormCat.1),
  #                 position = position_dodge2(width = 1), alpha = 0.5, size = 0.1, stroke = 0.035) +
  stat_summary(fun.data = n_fun, geom = "text", size = 4, angle = 90, hjust = 0) +
  scale_y_continuous(name = "Gene distance (JSD)") +
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



## put histograms on the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, ...)
}

pdf(paste0(dir, "Plots/gene_figure_plots/Metric_Matrix.pdf"),  width = 10, height = 10)
gene_data %>%
  filter(max_tpm > 80) %>%
  select(c("jsd_median_joint", "pcor_median_joint", "scor_median_joint",
           "cos_median_joint", "cel_tau_median_joint", "cbr_tau_median_joint")) %>%
pairs(.,
      lower.panel = NULL,
      diag.panel = panel.hist,
      cex = 0.1,
      pch = 1,
      col = alpha("black", 0.5),
      labels = c("Jensen-shannon dist.", "Pearson cor.", "Spearman cor.", "Cosine dist.", "C. elegans Tau", "C. briggsae Tau"))
dev.off()

# Modeling
car::Anova(lm(jsd_median_joint ~ log2(max_tpm) * mean_tau_joint + WormCat.1 + syntenic + percent_identity + PS.Value + omega + RecombRegion, gene_data[gene_data$cel_max_tpm > 80 & gene_data$cbr_max_tpm > 80,]))
car::Anova(lm(jsd_median_joint ~ log2(max_tpm) + WormCat.1 + syntenic + percent_identity + PS.Value + omega + RecombRegion, gene_data[gene_data$cel_max_tpm > 80 & gene_data$cbr_max_tpm > 80,]))

lm(jsd_median_joint ~ log2(max_tpm), gene_data[gene_data$cel_max_tpm > 80 & gene_data$cbr_max_tpm > 80,])
lm(jsd_median_joint ~ log2(max_tpm), gene_data[gene_data$max_tpm > 80,])

cor(gene_data[gene_data$max_tpm > 80,]$jsd_median_joint, log2(gene_data[gene_data$max_tpm > 80,]$max_tpm))


lm(jsd_median_joint ~ log2(max_tpm + 1), gene_data)




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
                              limits = c(1,
                                         max(temp[["C.elegans"]]$expression,
                                             temp[["C.briggsae"]]$expression) * 1.1)) +
           scale_y_continuous(name = paste0("C.briggsae ", gene, " TPM"), trans = "log2",
                              limits = c(1, max(temp[["C.elegans"]]$expression,
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



