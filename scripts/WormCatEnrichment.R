
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)

set.seed(2022)

# source 
dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

term_markers <- readRDS(paste0(dir, "Objects/TermMarker_list_filt_metadata_cell_bg.rds"))
pro_markers <- readRDS(paste0(dir, "Objects/ProMarker_list_filt_metadata_cell_bg.rds"))

deg <- readRDS(paste0(dir, "Objects/DEG_df_filt_metadata_cell_bg.rds"))

gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))
rownames(gff_list_mrna[["elegans"]]) <- ifelse(is.na(gff_list_mrna[["elegans"]]$cds_gene_name), gff_list_mrna[["elegans"]]$gene_name, gff_list_mrna[["elegans"]]$cds_gene_name)
rownames(gff_list_mrna[["briggsae"]]) <- ifelse(is.na(gff_list_mrna[["briggsae"]]$cds_gene_name), gff_list_mrna[["briggsae"]]$gene_name, gff_list_mrna[["briggsae"]]$cds_gene_name)

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
  
  print(paste0("Background genes: ", total_background, " RGS genes: ", total_rgs))
  
  out_list <- list()
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    out_list[[tier]] <- .merger_cats(rgs_cat_list[[tier]], background_cat_list[[tier]],
                                     total_background, total_rgs)
    out_list[[tier]] <- .worm_cat_acceptable_pvalues(out_list[[tier]])
  }
  return(out_list)
}


background_gene_set <- unique(deg$gene)
deg_enrichment_list <- list()
for(cur_cell_type in unique(deg$cell_type)) {
  print(cur_cell_type)
  deg_enrichment_list[[cur_cell_type]] <- wormcat_enrichment(unique(deg[deg$cell_type %in% cur_cell_type,]$gene), background_gene_set, gff_list_mrna[["elegans"]])
  
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    deg_enrichment_list[[cur_cell_type]][[tier]]$cell_type <- cur_cell_type
  }
}

deg_enrichment_list_filt <- list()
for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
  deg_enrichment_list_filt[[tier]] <- do.call(rbind, lapply(deg_enrichment_list, function(x) x[[tier]]))
  temp_filt <- data.frame(deg_enrichment_list_filt[[tier]] %>%
    filter(Bonferroni < 0.05) %>%
    select(Category) %>%
    unique())[,1]
  
  deg_enrichment_list_filt[[tier]] <- deg_enrichment_list_filt[[tier]] %>%
    filter(Category %in% temp_filt) %>%
    complete(Category, cell_type)
  
  deg_enrichment_list_filt[[tier]][is.na(deg_enrichment_list_filt[[tier]]$Fold),]$Fold <- 0
  deg_enrichment_list_filt[[tier]][is.na(deg_enrichment_list_filt[[tier]]$Bonferroni),]$Bonferroni <- 1
}

saveRDS(deg_enrichment_list_filt, paste0(dir, "Objects/deg_wormcat_enrichment_cell_bg.R"))

## Markers
term_markers[["C.elegans"]]$cel_tpm_log2fc_just_pro <- NA
term_markers[["C.elegans"]]$cbr_tpm_log2fc_just_pro <- NA

term_markers[["C.briggsae"]]$cel_tpm_log2fc_just_pro <- NA
term_markers[["C.briggsae"]]$cbr_tpm_log2fc_just_pro <- NA

joint_markers <- list()
joint_markers[["C.elegans"]] <- rbind(term_markers[["C.elegans"]][,colnames(pro_markers[["C.elegans"]])], pro_markers[["C.elegans"]])
joint_markers[["C.briggsae"]] <- rbind(term_markers[["C.briggsae"]][,colnames(pro_markers[["C.briggsae"]])], pro_markers[["C.briggsae"]])

# [[cur_species]][[tier]][[cur_cell_type]]

marker_enrichment_list <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  marker_enrichment_list[[cur_species]] <- list()
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    marker_enrichment_list[[cur_species]][[tier]] <- list()
  }
  
  background_gene_set <- unique(joint_markers[[cur_species]]$gene)
  for(cur_cell_type in  unique(joint_markers[[cur_species]]$cell_type)) {
    print(cur_cell_type)
    temp_list <- wormcat_enrichment(unique(joint_markers[[cur_species]][joint_markers[[cur_species]]$cell_type %in% cur_cell_type,]$gene),
                                    background_gene_set,
                                    gff_list_mrna[[ifelse(cur_species == "C.elegans", "elegans", "briggsae")]])
    
    for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
      marker_enrichment_list[[cur_species]][[tier]][[cur_cell_type]] <- temp_list[[tier]]
      if(nrow(marker_enrichment_list[[cur_species]][[tier]][[cur_cell_type]]) > 0) {
        marker_enrichment_list[[cur_species]][[tier]][[cur_cell_type]]$cell_type <- cur_cell_type
      }
    }
  }
}

marker_enrichment_list_filt <- list()
for(cur_species in c("C.elegans", "C.briggsae")) {
  marker_enrichment_list_filt[[cur_species]] <- list()
  
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    marker_enrichment_list_filt[[cur_species]][[tier]] <- do.call(rbind, marker_enrichment_list[[cur_species]][[tier]])
    temp_filt <- data.frame(marker_enrichment_list_filt[[cur_species]][[tier]] %>%
                              filter(Bonferroni < 0.05) %>%
                              select(Category) %>%
                              unique())[,1]
    
    marker_enrichment_list_filt[[cur_species]][[tier]] <- marker_enrichment_list_filt[[cur_species]][[tier]] %>%
      filter(Category %in% temp_filt) %>%
      complete(Category, cell_type)
    
    marker_enrichment_list_filt[[cur_species]][[tier]][is.na(marker_enrichment_list_filt[[cur_species]][[tier]]$Fold),]$Fold <- 0
    marker_enrichment_list_filt[[cur_species]][[tier]][is.na(marker_enrichment_list_filt[[cur_species]][[tier]]$Bonferroni),]$Bonferroni <- 1
    
    marker_enrichment_list_filt[[cur_species]][[tier]]$signif <- ""
    marker_enrichment_list_filt[[cur_species]][[tier]]$signif <- ifelse(marker_enrichment_list_filt[[cur_species]][[tier]]$Bonferroni < 0.05,
                                                                        "*",
                                                                        marker_enrichment_list_filt[[cur_species]][[tier]]$signif)
    marker_enrichment_list_filt[[cur_species]][[tier]]$signif <- ifelse(marker_enrichment_list_filt[[cur_species]][[tier]]$Bonferroni < 0.005,
                                                                        "**",
                                                                        marker_enrichment_list_filt[[cur_species]][[tier]]$signif)
    marker_enrichment_list_filt[[cur_species]][[tier]]$signif <- ifelse(marker_enrichment_list_filt[[cur_species]][[tier]]$Bonferroni < 0.0005,
                                                                        "***",
                                                                        marker_enrichment_list_filt[[cur_species]][[tier]]$signif)
    
  }
}

saveRDS(marker_enrichment_list_filt, paste0(dir, "Objects/marker_wormcat_enrichment_cell_bg.R"))

markers.df.reshape <- dcast(marker_enrichment_list_filt[["C.elegans"]][["WormCat.1"]], Category ~ cell_type, value.var = "Fold", fun.aggregate = mean)
rownames(markers.df.reshape) <- markers.df.reshape$Category

markers.matrix <- as.matrix(markers.df.reshape[-1])

row_order <- rownames(markers.matrix[hclust(dist(markers.matrix, method = "euclidean"), method = "complete")$order,])
col_order <- colnames(markers.matrix[,hclust(dist(t(markers.matrix), method = "euclidean"), method = "complete")$order])

pdf(paste0(dir, "Plots/marker_enrich_all.pdf"), height = 60, width = 7)
marker_enrichment_list_filt[["C.elegans"]][["WormCat.1"]] %>%
  ggplot(aes(x = Category, y = cell_type, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45), limits = row_order) +
  scale_y_discrete(limits = col_order) +
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

joint_plot_order <- readRDS(paste0(dir, "Objects/joint_plot_order.rds"))

plot_df <- cell_data[joint_plot_order ,c("cell_type", "cell_class", "lineage_group", "div_stage"), drop = FALSE]
plot_df$cell_class_lineage_group <- ifelse(plot_df$cell_class == "progenitor", plot_df$lineage_group, plot_df$cell_class)

cats <- (do.call(rbind, marker_enrichment_list_filt[["C.elegans"]]) %>%
           filter(cell_type %in% plot_df[plot_df$cell_class != "progenitor",]$cell_type) %>%
           group_by(Category) %>%
           summarise(signif = sum(Bonferroni < 0.05)) %>%
           filter(signif > 0))[,1, drop = TRUE]

pdf(paste0(dir, "Plots/marker_enrich_terminal.pdf"), height = 17.5, width = 22.5)
do.call(rbind, marker_enrichment_list_filt[["C.elegans"]]) %>%
  filter(cell_type %in% plot_df[plot_df$cell_class != "progenitor",]$cell_type) %>%
  filter(Category %in% cats) %>%
  ggplot(aes(x = Category, y = cell_type, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = rev(plot_df[plot_df$cell_class != "progenitor",]$cell_type)) +
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

cats <- (do.call(rbind, marker_enrichment_list_filt[["C.elegans"]]) %>%
           filter(cell_type %in% plot_df[plot_df$cell_class == "progenitor",]$cell_type) %>%
           group_by(Category) %>%
           summarise(signif = sum(Bonferroni < 0.05)) %>%
           filter(signif > 0))[,1, drop = TRUE]

pdf(paste0(dir, "Plots/marker_enrich_progenitor.pdf"), height = 42.5, width = 20)
do.call(rbind, marker_enrichment_list_filt[["C.elegans"]]) %>%
  filter(cell_type %in% plot_df[plot_df$cell_class == "progenitor",]$cell_type) %>%
  filter(Category %in% cats) %>%
  ggplot(aes(x = Category, y = cell_type, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  scale_y_discrete(limits = rev(plot_df[plot_df$cell_class == "progenitor",]$cell_type)) +
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



