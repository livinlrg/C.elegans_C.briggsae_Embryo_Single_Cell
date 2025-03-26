library(data.table)
library(dplyr)
library(RcisTarget)
library(parallel)
library(ggplot2)
library(tidyr)
library(reshape2)

dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

cel_database_table_motif <- read.table(paste0(dir, "Objects/scenic/cel_motifs.tbl"), header = FALSE, sep = "\t")
cbr_database_table_motif <- read.table(paste0(dir, "Objects/scenic/cbr_motifs.tbl"), header = FALSE, sep = "\t")

cel_database_table_chipseq <- read.table(paste0(dir, "Objects/scenic/cel_chipseq_motif_table.tbl"), header = FALSE, sep = "\t")

chipseq_motifs <- cel_database_table_chipseq[! cel_database_table_chipseq$V1 %in% cel_database_table_motif$V1,]$V1

# Function to read regulatory files
read_in_reg <- function(scanned_reg) {
  return(data.frame(rbindlist(lapply(scanned_reg, function(line) {
    temp_line <- strsplit(gsub(", ", "|", line), ",", fixed = TRUE)
    temp_df <- data.frame(tf = temp_line[[1]][1],
                          motif = temp_line[[1]][2],
                          auc = as.numeric(temp_line[[1]][3]),
                          nes = as.numeric(temp_line[[1]][4]),
                          genes = temp_line[[1]][9],
                          gene_scores = NA,
                          rank_at_max = temp_line[[1]][10])
    temp_df$genes <- gsub("\\", "", temp_df$genes, fixed = TRUE)
    temp_df$genes <- gsub("\"", "", temp_df$genes, fixed = TRUE)
    temp_df$genes <- gsub("[", "", temp_df$genes, fixed = TRUE)
    temp_df$genes <- gsub("]", "", temp_df$genes, fixed = TRUE)
    temp_df$genes <- gsub("'", "", temp_df$genes, fixed = TRUE)
    temp_df$genes <- gsub("(", "", temp_df$genes, fixed = TRUE)
    temp_df$genes <- gsub(")", "", temp_df$genes, fixed = TRUE)
    if(is.na(temp_df$genes)) {
      return(temp_df)
    }
    # print(temp_df$genes)
    temp_genes <- strsplit(temp_df$genes, "|", fixed = TRUE)[[1]]
    # print(temp_genes)
    if(length(temp_genes > 1)) {
      temp_df$gene_scores <- list(temp_genes[seq(from = 2, to = length(temp_genes), by = 2)])
      temp_df$genes <- list(temp_genes[seq(from = 1, to = length(temp_genes), by = 2)])
    }
    return(temp_df)
  }))))
}

# elegans
files_in_dir_cel <- list.files(paste0(dir, "Objects/scenic/cell_based/cel_marathon/ctx/"))

cel_regs <- lapply(files_in_dir_cel, function(file) {
  print(file)
  scanned_reg <- scan(paste0(dir, "Objects/scenic/cell_based/cel_marathon/ctx/", file), what="", sep="\n")
  scanned_reg <- scanned_reg[4:length(scanned_reg)]
  scanned_reg_df <- read_in_reg(scanned_reg)
  return(scanned_reg_df)
})

# briggsae
files_in_dir_cbr <- list.files(paste0(dir, "Objects/scenic/cell_based/cbr_marathon/ctx/"))

cbr_regs <- lapply(files_in_dir_cbr, function(file) {
  print(file)
  scanned_reg <- scan(paste0(dir, "Objects/scenic/cell_based/cbr_marathon/ctx/", file), what="", sep="\n")
  scanned_reg <- scanned_reg[4:length(scanned_reg)]
  scanned_reg_df <- read_in_reg(scanned_reg)
  return(scanned_reg_df)
})

saveRDS(cel_regs, paste0(dir, "Objects/scenic/cell_based/cel_regs.rds"))
saveRDS(cbr_regs, paste0(dir, "Objects/scenic/cell_based/cbr_regs.rds"))

cel_regs <- readRDS(paste0(dir, "Objects/scenic/cell_based/cel_regs.rds"))
cbr_regs <- readRDS(paste0(dir, "Objects/scenic/cell_based/cbr_regs.rds"))

## Break into rows and lists
cel_regs_filt <- lapply(cel_regs, function(reg) {
  temp_df <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(temp_df) <- c("tf", "motif", "genes", "gene_scores")
  for(motif in unique(reg$motif)) {
    temp_df <- rbind(temp_df, data.frame("tf" = unique(reg[reg$motif %in% motif, "tf"]),
                                         "motif" = unique(reg[reg$motif %in% motif, "motif"]),
                                         "genes" = unlist(reg[reg$motif %in% motif, "genes"]),
                                         "gene_scores" = unlist(reg[reg$motif %in% motif, "gene_scores"])))
  }
  return(temp_df)
})

cel_regs_filt <- lapply(seq(1, length(cel_regs_filt)), function(i) {
  cel_regs_filt[[i]]$iteration <- i
  return(cel_regs_filt[[i]])
})

cel_regs_filt_df <- do.call(rbind, cel_regs_filt)

saveRDS(cel_regs_filt_df, paste0(dir, "Objects/scenic/cell_based/cel_regs_filt_df.rds"))
cel_regs_filt_df <- readRDS(paste0(dir, "Objects/scenic/cell_based/cel_regs_filt_df.rds"))

# tf percent in
plot_vector <- sort(rowSums(table(cel_regs_filt_df$tf, cel_regs_filt_df$iteration) > 0))/1000
plot_df = data.frame(tf = names(plot_vector), value = plot_vector)
plot_df$tf <- factor(plot_df$tf, plot_df$tf)

pdf(paste0(dir, "Plots/scenic/cel_tfs_in.pdf"), width = 35, height = 5)
ggplot(plot_df) +
  geom_bar(aes(x = tf,
               y = value), stat = "identity") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "red") +
  scale_y_continuous(name = "Fraction observed", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major.y = element_line(color="grey80", size=0.25),
        panel.grid.minor.y = element_line(color="grey80", size=0.05),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
dev.off()

# Repeat wtih cbr
cbr_regs_filt <- lapply(cbr_regs, function(reg) {
  temp_df <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(temp_df) <- c("tf", "motif", "genes", "gene_scores")
  for(motif in unique(reg$motif)) {
    temp_df <- rbind(temp_df, data.frame("tf" = unique(reg[reg$motif %in% motif, "tf"]),
                                         "motif" = unique(reg[reg$motif %in% motif, "motif"]),
                                         "genes" = unlist(reg[reg$motif %in% motif, "genes"]),
                                         "gene_scores" = unlist(reg[reg$motif %in% motif, "gene_scores"])))
  }
  return(temp_df)
})

cbr_regs_filt <- lapply(seq(1, length(cbr_regs_filt)), function(i) {
  cbr_regs_filt[[i]]$iteration <- i
  return(cbr_regs_filt[[i]])
})

cbr_regs_filt_df <- do.call(rbind, cbr_regs_filt)

saveRDS(cbr_regs_filt_df, paste0(dir, "Objects/scenic/cell_based/cbr_regs_filt_df.rds"))
cbr_regs_filt_df <- readRDS(paste0(dir, "Objects/scenic/cell_based/cbr_regs_filt_df.rds"))

# tf percent in
plot_vector <- sort(rowSums(table(cbr_regs_filt_df$tf, cbr_regs_filt_df$iteration) > 0))/1000
plot_df = data.frame(tf = names(plot_vector), value = plot_vector)
plot_df$tf <- factor(plot_df$tf, plot_df$tf)

pdf(paste0(dir, "Plots/scenic/cbr_tfs_in.pdf"), width = 35, height = 5)
ggplot(plot_df) +
  geom_bar(aes(x = tf,
                 y = value), stat = "identity") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "red") +
  scale_y_continuous(name = "Fraction observed", breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme(legend.title= element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 24),
        axis.line = element_line(color="grey80", size=1),
        panel.grid.major.y = element_line(color="grey80", size=0.25),
        panel.grid.minor.y = element_line(color="grey80", size=0.05),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
        legend.box.background = element_rect(colour = "transparent", fill = "transparent"),
        legend.key = element_rect(colour = "transparent", fill = "transparent"),
        axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
dev.off()



## calculate proportions
cel_regs_df <- cel_regs_filt_df[! duplicated(paste0(cel_regs_filt_df$tf, cel_regs_filt_df$motif, cel_regs_filt_df$genes)),
                                 c("tf", "motif", "genes", "gene_scores")]
cel_regs_df$gene_frac_in <- NA
cel_regs_df$gene_motif_frac_in <- NA
cel_regs_df$motif_frac_in <- NA
cel_regs_df$tf_frac_in <- NA
rownames(cel_regs_df) <- paste(cel_regs_df$tf, cel_regs_df$motif, cel_regs_df$genes, sep = "_")

cel_regs_df$type <- ifelse(cel_regs_df$motif %in% chipseq_motifs, "chipseq", "motif")
cel_regs_filt_df$type <- ifelse(cel_regs_filt_df$motif %in% chipseq_motifs, "chipseq", "motif")

cel_tf_in <- sort(rowSums(table(cel_regs_filt_df[cel_regs_filt_df$type == "motif",]$tf,
                                cel_regs_filt_df[cel_regs_filt_df$type == "motif",]$iteration) > 0)/1000)
cel_motif_in <- sort(rowSums(table(cel_regs_filt_df[cel_regs_filt_df$type == "motif",]$motif,
                                   cel_regs_filt_df[cel_regs_filt_df$type == "motif",]$iteration) > 0)/1000)

for(tf in unique(cel_regs_filt_df[cel_regs_filt_df$type == "motif",]$tf)) {
  cel_regs_df$tf_frac_in[cel_regs_df$tf %in% tf] <- cel_tf_in[tf]
  for(motif in unique(cel_regs_filt_df[cel_regs_filt_df$tf %in% tf & cel_regs_filt_df$type == "motif",]$motif)) {
    print(paste0(tf, "_", motif))
    cel_regs_df$motif_frac_in[cel_regs_df$motif %in% motif & cel_regs_df$tf %in% tf] <- cel_motif_in[motif]
    
    temp_motif_df <- cel_regs_filt_df[cel_regs_filt_df$motif %in% motif &
                                        cel_regs_filt_df$tf %in% tf &
                                        cel_regs_filt_df$type == "motif",]
    temp_gene_in_frac <- sort(rowSums(table(temp_motif_df$genes, temp_motif_df$iteration) > 0)/length(unique(temp_motif_df$iteration)))
    
    cel_regs_df[paste(tf, motif, names(temp_gene_in_frac), sep = "_"),]$gene_motif_frac_in <- temp_gene_in_frac
    
    temp_motif_df <- cel_regs_filt_df[cel_regs_filt_df$tf %in% tf & cel_regs_filt_df$type == "motif",]
    temp_gene_in_frac <- sort(rowSums(table(temp_motif_df$genes, temp_motif_df$iteration) > 0)/length(unique(temp_motif_df$iteration)))
    
    cel_regs_df[paste(tf, motif, names(temp_gene_in_frac), sep = "_"),]$gene_frac_in <- temp_gene_in_frac
  }
}
cel_regs_df <- cel_regs_df[! is.na(cel_regs_df$tf),]
cel_regs_df <- cel_regs_df %>% arrange(cel_regs_df$tf, cel_regs_df$motif, cel_regs_df$genes)
cel_gene_scores_median <- cel_regs_filt_df %>% group_by(tf, motif, genes) %>% summarise(gene_scores = median(as.numeric(gene_scores)))

cel_regs_df <- left_join(cel_regs_df, cel_gene_scores_median, by = c("tf" = "tf", "motif" = "motif", "genes" = "genes"))

## celegans chip-seq
cel_regs_df_chip <- cel_regs_filt_df[! duplicated(paste0(cel_regs_filt_df$tf, cel_regs_filt_df$motif, cel_regs_filt_df$genes)),
                                c("tf", "motif", "genes", "gene_scores")]
cel_regs_df_chip$gene_frac_in <- NA
cel_regs_df_chip$gene_motif_frac_in <- NA
cel_regs_df_chip$motif_frac_in <- NA
cel_regs_df_chip$tf_frac_in <- NA
rownames(cel_regs_df_chip) <- paste(cel_regs_df_chip$tf, cel_regs_df_chip$motif, cel_regs_df_chip$genes, sep = "_")

cel_regs_df_chip$type <- ifelse(cel_regs_df_chip$motif %in% chipseq_motifs, "chipseq", "motif")
cel_tf_in <- sort(rowSums(table(cel_regs_filt_df$tf,
                                cel_regs_filt_df$iteration) > 0)/1000)
cel_motif_in <- sort(rowSums(table(cel_regs_filt_df$motif,
                                   cel_regs_filt_df$iteration) > 0)/1000)
for(tf in unique(cel_regs_filt_df$tf)) {
  cel_regs_df_chip$tf_frac_in[cel_regs_df_chip$tf %in% tf] <- cel_tf_in[tf]
  for(motif in unique(cel_regs_filt_df[cel_regs_filt_df$tf %in% tf,]$motif)) {
    print(paste0(tf, "_", motif))
    cel_regs_df_chip$motif_frac_in[cel_regs_df_chip$motif %in% motif & cel_regs_df_chip$tf %in% tf] <- cel_motif_in[motif]
    
    temp_motif_df <- cel_regs_filt_df[cel_regs_filt_df$motif %in% motif &
                                        cel_regs_filt_df$tf %in% tf,]
    temp_gene_in_frac <- sort(rowSums(table(temp_motif_df$genes, temp_motif_df$iteration) > 0)/length(unique(temp_motif_df$iteration)))
    
    cel_regs_df_chip[paste(tf, motif, names(temp_gene_in_frac), sep = "_"),]$gene_motif_frac_in <- temp_gene_in_frac
    
    temp_motif_df <- cel_regs_filt_df[cel_regs_filt_df$tf %in% tf,]
    temp_gene_in_frac <- sort(rowSums(table(temp_motif_df$genes, temp_motif_df$iteration) > 0)/length(unique(temp_motif_df$iteration)))
    
    cel_regs_df_chip[paste(tf, motif, names(temp_gene_in_frac), sep = "_"),]$gene_frac_in <- temp_gene_in_frac
  }
}
cel_regs_df_chip <- cel_regs_df_chip[! is.na(cel_regs_df_chip$tf),]
cel_regs_df_chip <- cel_regs_df_chip %>% arrange(cel_regs_df_chip$tf, cel_regs_df_chip$motif, cel_regs_df_chip$genes)
cel_gene_scores_median <- cel_regs_filt_df %>% group_by(tf, motif, genes) %>% summarise(gene_scores = median(as.numeric(gene_scores)))

cel_regs_df_chip <- left_join(cel_regs_df_chip, cel_gene_scores_median, by = c("tf" = "tf", "motif" = "motif", "genes" = "genes"))

## now for briggsae
cbr_regs_df <- cbr_regs_filt_df[! duplicated(paste0(cbr_regs_filt_df$tf, cbr_regs_filt_df$motif, cbr_regs_filt_df$genes)),
                                c("tf", "motif", "genes", "gene_scores")]
cbr_regs_df$gene_frac_in <- NA # actually tf_gene_frac_in
cbr_regs_df$gene_motif_frac_in <- NA
cbr_regs_df$motif_frac_in <- NA
cbr_regs_df$tf_frac_in <- NA

rownames(cbr_regs_df) <- paste(cbr_regs_df$tf, cbr_regs_df$motif, cbr_regs_df$genes, sep = "_")

for(tf in unique(cbr_regs_filt_df$tf)) {
  cbr_regs_df$tf_frac_in[cbr_regs_df$tf %in% tf] <- sort(rowSums(table(cbr_regs_filt_df$tf, cbr_regs_filt_df$iteration) > 0)/1000)[tf]
  for(motif in unique(cbr_regs_filt_df[cbr_regs_filt_df$tf %in% tf,]$motif)) {
    print(paste0(tf, "_", motif))
    cbr_regs_df$motif_frac_in[cbr_regs_df$motif %in% motif &
                                cbr_regs_df$tf %in% tf] <- sort(rowSums(table(cbr_regs_filt_df$motif, cbr_regs_filt_df$iteration) > 0)/1000)[motif]
    
    temp_motif_df <- cbr_regs_filt_df[cbr_regs_filt_df$motif %in% motif &
                                        cbr_regs_filt_df$tf %in% tf,]
    temp_gene_in_frac <- sort(rowSums(table(temp_motif_df$genes, temp_motif_df$iteration) > 0)/length(unique(temp_motif_df$iteration)))
    
    cbr_regs_df[paste(tf, motif, names(temp_gene_in_frac), sep = "_"),]$gene_motif_frac_in <- temp_gene_in_frac
    
    temp_motif_df <- cbr_regs_filt_df[cbr_regs_filt_df$tf %in% tf,]
    temp_gene_in_frac <- sort(rowSums(table(temp_motif_df$genes, temp_motif_df$iteration) > 0)/length(unique(temp_motif_df$iteration)))
    
    cbr_regs_df[paste(tf, motif, names(temp_gene_in_frac), sep = "_"),]$gene_frac_in <- temp_gene_in_frac
  }
}
cbr_regs_df <- cbr_regs_df[! is.na(cbr_regs_df$tf),]
cbr_regs_df <- cbr_regs_df %>% arrange(cbr_regs_df$tf, cbr_regs_df$motif, cbr_regs_df$genes)
cbr_gene_scores_median <- cbr_regs_filt_df %>% group_by(tf, motif, genes) %>% summarise(gene_scores = median(as.numeric(gene_scores)))

cbr_regs_df <- left_join(cbr_regs_df, cbr_gene_scores_median, by = c("tf" = "tf", "motif" = "motif", "genes" = "genes"))

## establish the 80% and 40% thresholds for tfs and motifs
# establish the 80% and 20% for genes
cel_regs_df <- cel_regs_df %>% 
  filter(type == "motif") %>%
  mutate(tf_in = ifelse(tf_frac_in < 0.8, ifelse(tf_frac_in >= 0.2, "Medium", "Low"), "High"),
         motif_in = ifelse(motif_frac_in < 0.8, ifelse(motif_frac_in >= 0.2, "Medium", "Low"), "High"),
         gene_in = ifelse(gene_frac_in < 0.8, ifelse(gene_frac_in >= 0.2, "Medium", "Low"), "High"))
cbr_regs_df <- cbr_regs_df %>% mutate(tf_in = ifelse(tf_frac_in < 0.8, ifelse(tf_frac_in >= 0.2, "Medium", "Low"), "High"),
                                      motif_in = ifelse(motif_frac_in < 0.8, ifelse(motif_frac_in >= 0.2, "Medium", "Low"), "High"),
                                      gene_in = ifelse(gene_frac_in < 0.8, ifelse(gene_frac_in >= 0.2, "Medium", "Low"), "High"))
cel_regs_df_chip <- cel_regs_df_chip %>% mutate(tf_in = ifelse(tf_frac_in < 0.8, ifelse(tf_frac_in >= 0.2, "Medium", "Low"), "High"),
                                      motif_in = ifelse(motif_frac_in < 0.8, ifelse(motif_frac_in >= 0.2, "Medium", "Low"), "High"),
                                      gene_in = ifelse(gene_frac_in < 0.8, ifelse(gene_frac_in >= 0.2, "Medium", "Low"), "High"))

data.frame(cel_regs_df %>% filter(tf_in == "High" & genes %in% common_genes) %>% group_by(tf) %>% filter(! duplicated(genes)) %>% summarise(gene_sum = sum(gene_in == "High"), mean_gene_score = mean(gene_scores))) %>% filter(gene_sum > 20)
pat_9_genes <- (cel_regs_df[cel_regs_df$tf %in% "pat-9",] %>% filter(tf_in == "High" & genes %in% common_genes) %>% filter(! duplicated(genes)) %>% filter(gene_in == "High") %>% select(genes))[,1]

cbr_regs_df[cbr_regs_df$tf %in% "pat-9",] %>% filter(tf_in == "High" & genes %in% common_genes & genes %in% pat_9_genes)
data.frame(cbr_regs_df %>% filter(tf_in == "High" & genes %in% common_genes) %>% group_by(tf) %>% filter(! duplicated(genes)) %>% summarise(gene_sum = sum(gene_in == "High"), mean_gene_score = mean(gene_scores))) %>% filter(gene_sum > 20)
data.frame(cel_regs_df_chip %>% filter(tf_in == "High" & genes %in% common_genes) %>% group_by(tf) %>% filter(! duplicated(genes)) %>% summarise(gene_sum = sum(gene_in == "High"), mean_gene_score = mean(gene_scores))) %>% filter(gene_sum > 20)

# cel 24
# cbr 21
# 16 shared
# not in "ceh-36"  "ceh-37"  "fos-1"   "lin-15B" "nhr-91"
unique(cbr_regs_df[cbr_regs_df$tf_in %in% "High",]$tf)[! unique(cbr_regs_df[cbr_regs_df$tf_in %in% "High",]$tf) %in% 
  unique(cel_regs_df[cel_regs_df$tf_in %in% c("High") & cel_regs_df$type == "motif",]$tf)]

# 20 out of 21 conserved
# not in "nhr-91"
unique(cbr_regs_df[cbr_regs_df$tf_in %in% "High",]$tf)[! unique(cbr_regs_df[cbr_regs_df$tf_in %in% "High",]$tf) %in% 
  unique(cel_regs_df[cel_regs_df$tf_in %in% c("High", "Medium") & cel_regs_df$type == "motif",]$tf)]

# 16 out of 24 conserved
# not in "ceh-10"  "lsl-1"   "nhr-109" "nhr-216" "nhr-23"  "nob-1"   "php-3"   "ttx-1"   "unc-30"
unique(cel_regs_df[cel_regs_df$tf_in %in% "High" & cel_regs_df$type == "motif",]$tf)[! unique(cel_regs_df[cel_regs_df$tf_in %in% "High" & cel_regs_df$type == "motif",]$tf) %in% 
  unique(cbr_regs_df[cbr_regs_df$tf_in %in% c("High"),]$tf)]

# 21 out of 24 conserved
# not in "ceh-10"  "nhr-109" "php-3"
unique(cel_regs_df[cel_regs_df$tf_in %in% "High" & cel_regs_df$type == "motif",]$tf)[! unique(cel_regs_df[cel_regs_df$tf_in %in% "High",]$tf) %in% 
  unique(cbr_regs_df[cbr_regs_df$tf_in %in% c("High", "Medium"),]$tf)]

unique(cel_regs_df[cel_regs_df$tf_in %in% c("High", "Medium"),]$tf)[! unique(cel_regs_df[cel_regs_df$tf_in %in% c("High", "Medium"),]$tf) %in% 
                                                         unique(cbr_regs_df[cbr_regs_df$tf_in %in% c("High", "Medium"),]$tf)]

conserved_tfs <- unique(unique(cel_regs_df[cel_regs_df$tf_in %in% "High",]$tf)[unique(cel_regs_df[cel_regs_df$tf_in %in% "High",]$tf) %in% 
                                                                        unique(cbr_regs_df[cbr_regs_df$tf_in %in% c("High", "Medium"),]$tf)],
                        unique(cbr_regs_df[cbr_regs_df$tf_in %in% "High",]$tf)[unique(cbr_regs_df[cbr_regs_df$tf_in %in% "High",]$tf) %in% 
                                                                                 unique(cel_regs_df[cel_regs_df$tf_in %in% c("High", "Medium"),]$tf)])

# load gene data
gene_data <- readRDS(paste0(dir, "Objects/gene_data_joint_cell_bg.rds"))
common_genes <- gene_data$gene

for(tf in conserved_tfs) {
  cel_genes <- cel_regs_df[cel_regs_df$tf %in% tf & cel_regs_df$gene_in %in% c("High"),]$genes
  cel_genes_filt <- cel_genes[cel_genes %in% common_genes]
  cbr_genes <- cbr_regs_df[cbr_regs_df$tf %in% tf & cbr_regs_df$gene_in %in% c("High"),]$genes
  cbr_genes_filt <- cbr_genes[cbr_genes %in% common_genes]
  cbr_genes_unfilt <- cbr_regs_df[cbr_regs_df$tf %in% tf & cbr_regs_df$gene_in %in% c("High", "Medium"),]$genes
  
  print(paste0("TF: ", tf,
               " Cel length: ", length(cel_genes_filt),
               " Cbr filt length: ", length(cbr_genes_filt),
               " Cel in Cbr: ", sum(cel_genes_filt %in% cbr_genes_unfilt) / length(cel_genes_filt),
               " Jaccard: ", jaccard_dist(cel_genes_filt, cbr_genes_filt)))
}

#jaccard distance function
jaccard_dist <- function(x, y) {
  x <- as.character(x)
  y <- as.character(y)
  x_union_y <- length(unique(c(x, y)))
  x_intersect_y <- length(intersect(x, y))
  return(1 - (x_intersect_y / x_union_y))
}

########################
# Check gene threshold #
########################

# Empirically look at the gene_threshold distributions for the in regs
out_genes <- regs_df[regs_df$tf_motif %in% tf_motif & regs_df$gene_frac_in >= gene_threshold & regs_df$genes %in% common_genes,]$genes
out_genes_score <- regs_df[regs_df$tf_motif %in% tf_motif & regs_df$gene_frac_in >= gene_threshold & regs_df$genes %in% common_genes,]$gene_scores

for(regs_df in list(cel_regs_df, cel_regs_df_chip, cbr_regs_df)) {
  temp_regs_df <- regs_df[regs_df$tf_frac_in > 0.8 & regs_df$gene_frac_in > 0.1,]
  
  for(tf in unique(temp_regs_df$tf)) {
    # temp_regs_df <- temp_regs_df[temp_regs_df$tf == tf,]
    print(paste0("TF: ", tf, " Gene frac in: ", (sum(temp_regs_df[temp_regs_df$tf == tf,]$gene_frac_in > 0.8)/length(temp_regs_df[temp_regs_df$tf == tf,]))))
  }
  print(summary(temp_regs_df$gene_frac_in))
}
temp_regs_df <- cbr_regs_df[cbr_regs_df$tf_frac_in > 0.8,]

temp_regs_df <- temp_regs_df %>%
  arrange(tf, gene_frac_in)

# Define the cumulative sum within each tf group
temp <- lapply(unique(temp_regs_df$tf), function(cur_tf) {
  temp_regs <- temp_regs_df[temp_regs_df$tf == cur_tf,]
  temp_regs %>%
    arrange(tf, desc(gene_frac_in)) %>%
    group_by(tf) %>%
    mutate(cumulative_sum = cumsum(gene_frac_in >= gene_frac_in))
})
temp <- data.frame(do.call(rbind, temp))

pdf(paste0(dir, "Plots/scenic/test.pdf"), width = 10, height = 5)
ggplot(temp, aes(x = gene_frac_in, y = cumulative_sum, color = tf)) +
  geom_line() +           # Line plot to connect points
  # geom_point() +          # Points for each threshold
  # scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
  #                    labels = c(1, 0.8, 0.6, 0.4, 0.2, 0)) +
  scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Cumulative Sum of Genes Passing Threshold by TF",
       x = "Threshold",
       y = "Cumulative Sum of Genes") +
  theme_minimal()
dev.off()

###########################
# Create grn output files #
###########################

# original file looks like:
# [1] ",,Enrichment,Enrichment,Enrichment,Enrichment,Enrichment,Enrichment,Enrichment,Enrichment"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
# [2] ",,AUC,NES,MotifSimilarityQvalue,OrthologousIdentity,Annotation,Context,TargetGenes,RankAtMax"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
# [3] "TF,MotifID,,,,,,,,"                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
# [4] "che-1,che-1_M00609_2.00,0.06062529071699116,4.495677307881477,0.0,1.0,gene is directly annotated,\"frozenset({'elegans_WS260_1000_.genes_vs_motifs.rankings', 'activating', 'weight>75.0%'})\",
# \"[('nspa-10', 1.4634721254022351), ('srv-21', 1.2794719822226035), ('ZK354.7', 0.926197359659824), ('M01E10.3', 3.2440600068308614), ('hlh-10', 0.7201306455334311), ('gcy-19', 0.7445089809123647), ('gcy-13', 0.3470247877362786), ('T08G3.6', 0.6539952233373272), ('srh-266', 0.572205180353126), ('F41C3.7', 0.5147284621436841), ('gcy-8', 0.3626958237783603), ('srh-112', 0.3642028323563814), ('T09E11.8', 0.5005412589780202), ('H20J04.1', 0.6057657401427019), ('ZC84.7', 5.064011607137173), ('T02H6.8', 0.4860539885772801), ('kel-10', 0.3759495592600916), ('T22D1.3', 0.3843489573997686), ('nhr-254', 4.460848837788149), ('gcy-5', 0.6376361665637749), ('ceh-36', 1.0742530524559204), ('tag-4', 2.2519306250827453), ('T04B2.7', 0.3140774872573552), ('clec-176', 0.6575581842664913), ('C46F4.3', 0.6122194462034529), ('W10C6.2', 0.6201735935371883), ('ttll-12', 0.8249977364163552), ('odr-10', 4.417513860640181), ('Y116A8A.150', 0.5422255567748224), ('R102.2', 0.3643394451164232), ('mig-39', 0.4181075569604617), ('C29F5.8', 0.4475400239656149), ('C05G5.3', 0.7702318611149106), ('B0285.4', 1.2889239654935565), ('tyr-3', 0.3298755373830985), ('C07B5.3', 0.3127497325265043), ('ZC190.8', 0.4683286535568636), ('W09C3.8', 2.0265668298292705), ('srh-279', 0.5687087549479333), ('lys-3', 0.4259912696133049), ('flp-26', 0.3381053925241072), ('K09C6.3', 0.5602582395046924), ('ugt-56', 4.944743883017949), ('ins-26', 1.037841137506648), ('trx-2', 0.8614433325228168), ('Y69E1A.8', 0.4762393107535524), ('gcy-17', 3.137578234180936), ('che-1', 1.0), ('srz-99', 0.5849726098057242), ('srx-47', 0.4841509281433063), ('srh-111', 0.3049247740931969), ('C32C4.3', 1.783383995308551), ('Y11D7A.8', 0.7911690248380737), ('F58A4.1', 0.4447277897075052), ('T01G6.1', 0.9224694009571932), ('T10B10.9', 0.6461240509011945), ('T10H4.13', 2.8276614229398493), ('str-33', 0.397312104764441), ('twk-4', 0.8939917945000393), ('C38C10.3', 0.3210522204495623), ('ceh-90', 3.594934511956937), ('fbxa-44', 0.776828837715083), ('chil-25', 0.3237656744884648), ('Y53F4B.20', 1.313056693572033), ('F42A6.5', 0.3172347617234859), ('ilys-5', 0.7714801940734264), ('srxa-15', 0.328973749559788), ('F16G10.2', 0.5534783320607537), ('F15D3.4', 0.4926170043422586), ('F58F9.1', 0.5106101769012235), ('srt-11', 0.3079493412718683), ('Y39A1A.9', 0.4806476974253419), ('sri-25', 0.548614875717464), ('C47F8.6', 0.7969247067996484), 
# ('F46F5.8', 0.5241088298829004)]\",3248"
cbr_regs_df$tf_motif <- paste(cbr_regs_df$tf, cbr_regs_df$motif, sep = "_")
cel_regs_df$tf_motif <- paste(cel_regs_df$tf, cel_regs_df$motif, sep = "_")
cel_regs_df_chip$tf_motif <- paste(cel_regs_df_chip$tf, cel_regs_df_chip$motif, sep = "_")

recreate_file <- function(file, regs_df, module_threshold, gene_threshold, common_genes, cur_species) {
  print(file)
  fileConn <- file(file)
  writeLines(",,Enrichment,Enrichment,Enrichment,Enrichment,Enrichment,Enrichment,Enrichment,Enrichment", fileConn)
  cat(",,AUC,NES,MotifSimilarityQvalue,OrthologousIdentity,Annotation,Context,TargetGenes,RankAtMax", file = file, append = TRUE, sep = "\n")
  cat("TF,MotifID,,,,,,,,", file = file, append = TRUE, sep = "\n")
  tfs_to_use <- unique(regs_df[regs_df$tf_frac_in >= module_threshold,]$tf)
  tf_motif_to_use <- unique(paste(regs_df[regs_df$tf %in% tfs_to_use,]$tf, regs_df[regs_df$tf %in% tfs_to_use,]$motif, sep = "_"))
  
  print(tf_motif_to_use)
  for(tf_motif in tf_motif_to_use) {
    out_genes <- regs_df[regs_df$tf_motif %in% tf_motif & regs_df$gene_frac_in >= gene_threshold & regs_df$genes %in% common_genes,]$genes
    out_genes_score <- regs_df[regs_df$tf_motif %in% tf_motif & regs_df$gene_frac_in >= gene_threshold & regs_df$genes %in% common_genes,]$gene_scores
    print(paste0("TF motif: ", tf_motif, " Common gene num.: ", length(unique(out_genes)), " ",
                 "Total gene num.: ", length(unique(regs_df[regs_df$tf_motif %in% tf_motif & regs_df$gene_frac_in >= gene_threshold,]$genes))))
    if(length(unique(out_genes)) < 20) {
      next
    } else {
      # Using mapply to concatenate elements of out_genes and out_genes_score
      gene_print <- mapply(function(gene, score) paste0("('", gene, "', ", score, ")"), out_genes, out_genes_score)
      
      # Concatenating the result into a single string
      gene_print_collapse <- paste0(gene_print, collapse = ", ")
      
      cat(paste(unique(regs_df[regs_df$tf_motif %in% tf_motif,]$tf),
                unique(regs_df[regs_df$tf_motif %in% tf_motif,]$motif),
                "1",
                "5",
                "0.0",
                "1.0",
                "gene is directly annotated",
                paste0("\"frozenset({'", cur_species ,"_WS290_1000_.genes_vs_motifs.rankings', 'activating', 'weight>75.0%'})\""),
                paste0("\"[", gene_print_collapse, "]\""),
                "3000",
                sep = ","), file = file, append = TRUE, sep = "\n")
    }
  }
  close(fileConn)
}

recreate_file(file = paste0(dir, "Objects/scenic/", "elegans", "_motif.aggregate.high.cell_bg.csv"),
              regs_df = cel_regs_df,
              module_threshold = 0.8,
              gene_threshold = 0.6,
              common_genes = common_genes,
              cur_species = "elegans")

recreate_file(file = paste0(dir, "Objects/scenic/", "elegans", "_chipseq_motif.aggregate.high.cell_bg.csv"),
              regs_df = cel_regs_df_chip,
              module_threshold = 0.8,
              gene_threshold = 0.6,
              common_genes = common_genes,
              cur_species = "elegans")

recreate_file(file = paste0(dir, "Objects/scenic/", "briggsae", "_motif.aggregate.high.cell_bg.csv"),
              regs_df = cbr_regs_df,
              module_threshold = 0.8,
              gene_threshold = 0.6,
              common_genes = common_genes,
              cur_species = "briggsae")

# Medium
recreate_file(file = paste0(dir, "Objects/scenic/", "elegans", "_motif.aggregate.medium.cell_bg.csv"),
              regs_df = cel_regs_df,
              module_threshold = 0.4,
              gene_threshold = 0.4,
              common_genes = common_genes,
              cur_species = "elegans")

recreate_file(file = paste0(dir, "Objects/scenic/", "elegans", "_chipseq_motif.aggregate.medium.cell_bg.csv"),
              regs_df = cel_regs_df_chip,
              module_threshold = 0.4,
              gene_threshold = 0.4,
              common_genes = common_genes,
              cur_species = "elegans")

recreate_file(file = paste0(dir, "Objects/scenic/", "briggsae", "_motif.aggregate.medium.cell_bg.csv"),
              regs_df = cbr_regs_df,
              module_threshold = 0.4,
              gene_threshold = 0.4,
              common_genes = common_genes,
              cur_species = "briggsae")

## Output the curated GRN

curate_grn <- function(regs_df, module_threshold, gene_threshold) {
  tfs_to_use <- unique(regs_df[regs_df$tf_frac_in >= module_threshold,]$tf)
  tf_motif_to_use <- unique(paste(regs_df[regs_df$tf %in% tfs_to_use,]$tf, regs_df[regs_df$tf %in% tfs_to_use,]$motif, sep = "_"))
  
  print(tf_motif_to_use)
  for(tf_motif in tf_motif_to_use) {
    out_genes <- regs_df[regs_df$tf_motif %in% tf_motif & regs_df$gene_frac_in >= gene_threshold,]$genes
    
    unique(regs_df[regs_df$tf_motif %in% tf_motif,]$tf)
    unique(regs_df[regs_df$tf_motif %in% tf_motif,]$motif)
  }
}

gene_threshold = 0.6
tf_threshold = 0.8

cel_motif_tfs_to_use <- (cel_regs_df %>%
                              filter(tf_frac_in > tf_threshold & gene_frac_in > gene_threshold) %>%
                              group_by(tf) %>%
                              summarise(gene_count = n_distinct(genes)) %>%
                              filter(gene_count > 20))[,"tf", drop = TRUE]

cbr_motif_tfs_to_use <- (cbr_regs_df %>%
                              filter(tf_frac_in > tf_threshold & gene_frac_in > gene_threshold) %>%
                              group_by(tf) %>%
                              summarise(gene_count = n_distinct(genes)) %>%
                              filter(gene_count > 20))[,"tf", drop = TRUE]

tfs_to_use <- unique(cel_motif_tfs_to_use, cbr_motif_tfs_to_use)

cel_motif_tfs_to_use <- (cel_regs_df %>%
                           filter(tf %in% tfs_to_use & gene_frac_in > gene_threshold) %>%
                           group_by(tf) %>%
                           summarise(gene_count = n_distinct(genes)) %>%
                           filter(gene_count > 20))[,"tf", drop = TRUE]

cbr_motif_tfs_to_use <- (cbr_regs_df %>%
                           filter(tf %in% tfs_to_use & gene_frac_in > gene_threshold) %>%
                           group_by(tf) %>%
                           summarise(gene_count = n_distinct(genes)) %>%
                           filter(gene_count > 20))[,"tf", drop = TRUE]

# cel motif
cel_regs_df_filt <- cel_regs_df[cel_regs_df$tf %in% cel_motif_tfs_to_use & cel_regs_df$gene_frac_in > gene_threshold,]

# cbr motif
cbr_regs_df_filt <- cbr_regs_df[cbr_regs_df$tf %in% cbr_motif_tfs_to_use & cbr_regs_df$gene_frac_in > gene_threshold,]

# cel chip
cel_chip_tfs_to_use <- (cel_regs_df_chip %>%
                          filter(tf_frac_in > tf_threshold & gene_frac_in > gene_threshold) %>%
                          group_by(tf) %>%
                          summarise(gene_count = n_distinct(genes)) %>%
                          filter(gene_count > 20))[,"tf", drop = TRUE]

cel_regs_df_chip_filt <- cel_regs_df_chip[cel_regs_df_chip$tf %in% cel_chip_tfs_to_use & cel_regs_df_chip$gene_frac_in > gene_threshold,]

write.table(cel_regs_df_filt, paste0(dir, "Tables/cel_regs_df_filt.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cbr_regs_df_filt, paste0(dir, "Tables/cbr_regs_df_filt.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(cel_regs_df_chip_filt, paste0(dir, "Tables/cel_regs_df_chip_filt.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(cel_regs_df, paste0(dir, "Objects/scenic/cel_regs_df.rds"))
saveRDS(cbr_regs_df, paste0(dir, "Objects/scenic/cbr_regs_df.rds"))
saveRDS(cel_regs_df_chip, paste0(dir, "Objects/scenic/cel_regs_df_chip.rds"))

saveRDS(cel_regs_df_filt, paste0(dir, "Objects/scenic/cel_regs_df_filt.rds"))
saveRDS(cbr_regs_df_filt, paste0(dir, "Objects/scenic/cbr_regs_df_filt.rds"))
saveRDS(cel_regs_df_chip_filt, paste0(dir, "Objects/scenic/cel_regs_df_chip_filt.rds"))

##############################
# TF expression conservation #
##############################

gene_data <- readRDS(paste0(dir, "Objects/gene_data_joint_cell_bg.rds"))
cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))

cel_tfs <- read.table(paste0(dir, "Objects/scenic/cel_tfs.txt"))[,1]
cbr_tfs <- read.table(paste0(dir, "Objects/scenic/cbr_tfs.txt"))[,1]

summary(gene_data[gene_data$gene %in% cel_tfs,]$max_tpm_term > 80 | gene_data[gene_data$gene %in% cel_tfs,]$max_tpm_pro > 80)
summary(gene_data[gene_data$gene %in% cbr_tfs,]$max_tpm_term > 80 | gene_data[gene_data$gene %in% cbr_tfs,]$max_tpm_pro > 80)


pdf(paste0(dir, "Plots/scenic/tf_expression.pdf"), width = 6, height = 5)
gene_data[gene_data$gene %in% c(cel_tfs, cbr_tfs),] %>%
  ggplot(aes(x = max_tpm_term + 1, y = max_tpm_pro + 1, label = gene, color = cel_OG_count + cbr_OG_count)) +
  geom_point() +
  geom_text_repel() +
  scale_x_continuous(trans = "log2") +
  scale_y_continuous(trans = "log2") +
  scale_color_continuous(name = "OG")
  # scale_color_viridis(trans = "log2")
dev.off()

pdf(paste0(dir, "Plots/scenic/tf_conservation.pdf"), width = 6, height = 5)
gene_data[gene_data$gene %in% c(cel_tfs, cbr_tfs),] %>%
  filter(max_tpm > 80) %>%
  ggplot(aes(x = jsd_median_joint, y = max_tpm, label = gene, color = gene %in% c(cel_regs_df_filt$tf, cbr_regs_df_filt$tf, cel_regs_df_chip_filt$tf))) +
  geom_point() +
  geom_text_repel() +
  scale_x_continuous() +
  scale_y_continuous(trans = "log2") +
  scale_color_discrete(name = "tf")
# scale_color_viridis(trans = "log2")
dev.off()

pdf(paste0(dir, "Plots/scenic/tf_conservation.pdf"), width = 6, height = 5)
gene_data[gene_data$gene %in% c(cel_tfs, cbr_tfs),] %>%
  filter(max_tpm > 80) %>%
  ggplot(aes(x = mean_tau_joint, y = jsd_median_joint, label = gene, color = cel_OG_count + cbr_OG_count)) +
  geom_point() +
  geom_text_repel() +
  scale_x_continuous() +
  scale_y_continuous(trans = "log2") +
  scale_color_continuous(name = "OG")
# scale_color_viridis(trans = "log2")
dev.off()

gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))
rownames(gff_list_mrna[["elegans"]]) <- ifelse(is.na(gff_list_mrna[["elegans"]]$cds_gene_name), gff_list_mrna[["elegans"]]$gene_name, gff_list_mrna[["elegans"]]$cds_gene_name)

og_data <- readRDS("/kimdata/livinlrg/scAnalysis/BobDataComb/Objects/WS290/og_data.rds")

cel_ogs <- gff_list_mrna[["elegans"]][unique(c(cel_tfs)[! c(cel_tfs) %in% gene_data$gene]),]$OG

og_data[cel_ogs,1:9]$BriGenes == 0

cel_tfs <- read.table(paste0(dir, "Objects/scenic/cel_tfs.txt"))[,1]
cbr_tfs <- read.table(paste0(dir, "Objects/scenic/cbr_tfs.txt"))[,1]

# 658 total TF's that are 1:1
# 619 are 'expressed'
gene_data[gene_data$gene %in% cel_tfs,] %>% filter(max_tpm_term > 80 | max_tpm_pro > 80) %>% nrow()

# joint_jsd  term_jsd   pro_jsd
# 0.5084851 0.5431008 0.6143365
gene_data[gene_data$gene %in% cel_tfs,] %>% filter(max_tpm_term > 80 | max_tpm_pro > 80) %>% summarise(joint_jsd = mean(jsd_median_joint),
                                                                                                       term_jsd = mean(jsd_median_term),
                                                                                                       pro_jsd = mean(jsd_median_pro))

# joint_jsd  term_jsd   pro_jsd
# 0.5417206 0.5703673 0.6572025
gene_data[! gene_data$gene %in% cel_tfs,] %>% filter(max_tpm_term > 80 | max_tpm_pro > 80) %>% summarise(joint_jsd = mean(jsd_median_joint),
                                                                                                         term_jsd = mean(jsd_median_term),
                                                                                                         pro_jsd = mean(jsd_median_pro))

gene_data$max_tpm <- max(gene_data$max_tpm_term, gene_data$max_tpm_pro)

###########################
# Load the SCENIC results #
###########################

# load barcode bins from EmbryoTimeBins
BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))
ProBarcodeBinsList <- readRDS(paste0(dir, "Objects/ProBarcodeBinsList.rds"))
JointBarcodeBinsList <- readRDS(paste0(dir, "Objects/JointBarcodeBinsList.rds"))

# load cell tables
CellTable <- readRDS(paste0(dir, "Objects/CellTable_Names_20240901.rds"))
cell_type_time_bins <- readRDS(paste0(dir, "Objects/cell_type_time_bins.rds"))
CellTableTimeBins <- readRDS(paste0(dir, "Objects/CellTableTimeBins.rds"))

# Load the scenic auc files
# Aggregate based on cell_types in the barcode lists
# Create median values and output an auc matrix
# Convert to z-score matrix

z_score <- function(value, dist) {
  return((value - mean(dist))/sd(dist))
}

scenic_auc_list <- list()
scenic_mean_list <- list()
scenic_zscore_list <- list()
scenic_mean_list_out <- list()
for(species in c("elegans", "briggsae")) {
  cur_species = ifelse(species == "elegans", "C.elegans", "C.briggsae")
  for(version in c("cel", "cbr", "cel_chip")) {
    print(paste0(cur_species, " ", version))
    
    scenic_auc_list[[paste0(species, "_", version)]] <- read.csv(paste0(dir, "Objects/scenic/WS290_", species, ".", version, "_grn.pyscenic.high.cell_bg.csv"),
                        header = FALSE, stringsAsFactors = FALSE)
    
    module_names_temp <- unlist(scenic_auc_list[[paste0(species, "_", version)]][1,-1])
    colnames(scenic_auc_list[[paste0(species, "_", version)]]) <- c("Cell", module_names_temp)
    
    # pull barcode names
    barcode_names <- scenic_auc_list[[paste0(species, "_", version)]][, 1]
    
    scenic_auc_list[[paste0(species, "_", version)]] <- scenic_auc_list[[paste0(species, "_", version)]][-1, -1]
    
    # Convert to numeric
    scenic_auc_list[[paste0(species, "_", version)]] <- apply(scenic_auc_list[[paste0(species, "_", version)]], 2, as.numeric)
    rownames(scenic_auc_list[[paste0(species, "_", version)]]) <- barcode_names[-1]
    
    scenic_mean_list[[paste0(species, "_", version)]] <- list()
    # Aggregate based on cell_types in the barcode lists
    for(broad_cell_type in names(JointBarcodeBinsList[[cur_species]])) {
      for(i in 1:length(JointBarcodeBinsList[[cur_species]][[broad_cell_type]])) {
        scenic_mean_list[[paste0(species, "_", version)]][[names(JointBarcodeBinsList[[cur_species]][[broad_cell_type]])[i]]] <- colMeans(scenic_auc_list[[paste0(species, "_", version)]][JointBarcodeBinsList[[cur_species]][[broad_cell_type]][[i]],,drop = FALSE])
      }
    }
    scenic_mean_list_out[[paste0(species, "_", version)]] <- t(do.call(rbind, scenic_mean_list[[paste0(species, "_", version)]]))
    
    temp_zscore_list <- list()
    for(tf in rownames(scenic_mean_list_out[[paste0(species, "_", version)]])) {
      temp_zscore_list[[tf]] <- sapply(scenic_mean_list_out[[paste0(species, "_", version)]][tf,], function(auc_mean) {
        z_score(auc_mean, scenic_mean_list_out[[paste0(species, "_", version)]][tf,])
        })
    }
    scenic_zscore_list[[paste0(species, "_", version)]] <- t(do.call(cbind, temp_zscore_list))
  }
}


col_order <- pheatmap(scenic_zscore_list[["elegans_cel"]])$tree_col$labels[pheatmap(scenic_zscore_list[["elegans_cel"]])$tree_col$order]
row_order <- pheatmap(scenic_zscore_list[["elegans_cel"]])$tree_row$labels[pheatmap(scenic_zscore_list[["elegans_cel"]])$tree_row$order]

pdf(paste0(dir, "Plots/scenic/cbr.cel_grn.heatmap.zscore.pdf"), width = 50, height = 12.5)
pheatmap(scenic_zscore_list[["briggsae_cel"]][row_order,col_order],
         color = rev(inferno(1000)),
         breaks = seq(-2.5, 5, 7.5/1000),
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

## Output the modules as a txt file
cel_regs_df$gene_orthology <- ifelse(cel_regs_df$genes %in% common_genes, "orthologous", "non_orthologous")
cel_regs_df$tf_orthology <- ifelse(cel_regs_df$tf %in% common_genes, "orthologous", "non_orthologous")
cel_regs_df_out <-(cel_regs_df %>% select(c("tf", "motif", "genes", "tf_orthology", "gene_orthology", "tf_frac_in", "gene_motif_frac_in", "gene_frac_in", "motif_frac_in", "gene_scores", "tf_in", "motif_in", "gene_in")))

write.table(cel_regs_df_out, file = paste0(NewBobPath, "SCENIC/cel_regs_df.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

cbr_regs_df$gene_orthology <- ifelse(cbr_regs_df$genes %in% common_genes, "orthologous", "non_orthologous")
cbr_regs_df$tf_orthology <- ifelse(cbr_regs_df$tf %in% common_genes, "orthologous", "non_orthologous")
cbr_regs_df_out <-(cbr_regs_df %>% select(c("tf", "motif", "genes", "tf_orthology", "gene_orthology","tf_frac_in", "motif_frac_in", "gene_frac_in", "gene_motif_frac_in", "gene_scores", "tf_in", "motif_in", "gene_in")))

write.table(cbr_regs_df_out, file = paste0(NewBobPath, "SCENIC/cbr_regs_df.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

######################
# Load into VisCello #
######################
library(VisCello)

pData <- readRDS(file = paste0(dir, "Objects/cello_creation/pdata.rds"))
# fData <- readRDS(file = paste0(dir, "Objects/cello_creation/fdata.rds"))
clist <- readRDS(paste0(dir, "Objects/cello_creation/cello_internal/clist.rds"))

for(version in c("cel", "cbr", "cel_chip")) {
  print(paste0(version))
  
  eset_matrix <- data.matrix(cbind(t(scenic_auc_list[[paste0("elegans", "_", version)]]),
                                   t(scenic_auc_list[[paste0("briggsae", "_", version)]])))
  
  motifs <- rownames(eset_matrix)
  rownames(eset_matrix) <- gsub("(+)", "", motifs, fixed = TRUE)
  
  fData <- data.frame(motif = motifs,
                      motif_short = rownames(eset_matrix))
  rownames(fData) <- fData$motif_short
  
  eset_pData <- pData[colnames(eset_matrix),]
  
  # Cello  object
  eset <- new("ExpressionSet",
              assayData = assayDataNew("environment",
                                       exprs = Matrix(eset_matrix, sparse = TRUE),
                                       norm_exprs = Matrix(eset_matrix, sparse = TRUE)),
              phenoData = new("AnnotatedDataFrame", data = eset_pData),
              featureData = new("AnnotatedDataFrame", data = fData))
  
  clist_out <- list()
  for(i in 1:length(clist)) {
    obj_name = names(clist)[i]
    print(obj_name)
    
    umap_2d <- clist[i][[1]]@proj[1][[1]]
    umap_3d <- clist[i][[1]]@proj[2][[1]]
    
    umap_2d <- umap_2d[rownames(umap_2d) %in% colnames(eset),]
    umap_3d <- umap_3d[rownames(umap_3d) %in% colnames(eset),]
    
    cur_idx = match(rownames(umap_2d), colnames(eset))
    if(sum(is.na(cur_idx)) > 0) {
      print("Fuck")
    }
    cello <- new("Cello", name = obj_name, idx = cur_idx)
    
    cello@proj <- list(assign(paste0(obj_name, " UMAP [2D]"), umap_2d),
                       assign(paste0(obj_name, " UMAP [3D]"), umap_3d))
    
    names(cello@proj) <- c(paste0(obj_name, " UMAP [2D]"), paste0(obj_name, " UMAP [3D]"))
    
    clist_out[[obj_name]] <- cello
  }
  
  dir.create(paste0(dir, "Objects/scenic/VisCello/", version), showWarnings = FALSE)
  saveRDS(eset, file = paste0(dir, "Objects/scenic/VisCello/", version, "/eset.rds"))
  saveRDS(clist_out, file = paste0(dir, "Objects/scenic/VisCello/", version, "/clist.rds"))
}

###############
# Create plot #
###############

joint_plot_order <- readRDS(paste0(dir, "Objects/joint_plot_order.rds"))

CellTable <- readRDS(paste0(dir, "Objects/CellTable_Names_20240901.rds"))
cell_type_time_bins <- readRDS(paste0(dir, "Objects/cell_type_time_bins.rds"))
CellTableTimeBins <- readRDS(paste0(dir, "Objects/CellTableTimeBins.rds"))

# load barcode bins from EmbryoTimeBins
BarcodeBinsList <- readRDS(paste0(dir, "Objects/BarcodeBinsList.rds"))
ProBarcodeBinsList <- readRDS(paste0(dir, "Objects/ProBarcodeBinsList.rds"))
JointBarcodeBinsList <- readRDS(paste0(dir, "Objects/JointBarcodeBinsList.rds"))

cell_data <- readRDS(paste0(dir, "Objects/cell_data_joint_mean_cell_bg.rds"))

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240903.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

time_bin_vector <- c("lt_100", "100_130", "130_170", "170_210", "210_270", "270_330", "330_390", "390_450", "450_510", "510_580", "580_650", "650_710", "gt_710")

time_bin_df <- data.frame(bins = time_bin_vector,
                          start = c(0, 100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710),
                          end = c(100, 130, 170, 210, 270, 330, 390, 450, 510, 580, 650, 710, 1000))

## cell bins
sort(unique(c(cell_type_time_bins$end_bin, cell_type_time_bins$start_bin)))
cell_type_time_bins %>% group_by(cell_type) %>%
  summarise(first = sum(end_bin %in% c("130_170", "170_210", "210_270", "270_330", "330_390")),
            middle = sum(end_bin %in% c("390_450", "450_510", "510_580")), 
            last = sum(end_bin %in% c("580_650", "650_710", "gt_710")))

GetParent <- function(x){
  if(is.na(x) || x=="P") return(NA)
  if(x=="P0" || x=="tossed") return(NA)
  if(is.element(substr(x,nchar(x),nchar(x)),c("a","p","d","v","l","r","x"))) return (substr(x,1,nchar(x)-1))
  if(x=="MSx1" || x=="MSx2") return("MS")
  if(x=="Cx1" || x=="Cx2") return("C")
  if(x=="AB" || x=="P1") return("P0")
  if(x=="EMS" || x=="P2") return("P1")
  if(x=="E" || x=="MS") return("EMS")
  if(x=="C" || x=="P3") return("P2")
  if(x=="D" || x=="P4") return("P3")
  if(x=="Z2" || x=="Z3") return("P4")
  return("P")
}

lineage_list <- list()
for(cur_cell_class in unique(cell_data$cell_class)[unique(cell_data$cell_class) != "progenitor"]) {
  print(cur_cell_class)
  
  for(cur_cell_type in cell_data[cell_data$cell_class %in% cur_cell_class,]$cell_type) {
    cur_lineage_set <- CellTable[which(CellTable$MergedDatasetName == cur_cell_type),]$Lineage
    for(i in 1:length(cur_lineage_set)) {
      cur_lineage <- cur_lineage_set[i]
      lineage_list[[cur_cell_class]] <- c(lineage_list[[cur_cell_class]], cur_lineage)
      while(! is.na(GetParent(cur_lineage))) {
        lineage_list[[cur_cell_class]] <- c(lineage_list[[cur_cell_class]], GetParent(cur_lineage))
        cur_lineage <- GetParent(cur_lineage)
      }
    }
  }
  lineage_list[[cur_cell_class]] <- sort(unique(lineage_list[[cur_cell_class]]))
}

saveRDS(lineage_list, paste0(dir, "Objects/lineage_list.rds"))

cell_class_barcodes <- list()
for(cur_cell_class in unique(cell_data$cell_class[! cell_data$cell_class %in% c("progenitor", "Germline")])) {
  print(cur_cell_class)
  temp_cell_types <- lineage_list[[cur_cell_class]]
  temp_cell_table <- CellTable[temp_cell_types,]
  temp_cell_table <- temp_cell_table[! is.na(temp_cell_table$MergedDatasetName),]
  
  cell_class_barcodes[[cur_cell_class]] <- unlist(c(JointBarcodeBinsList[["C.elegans"]][temp_cell_table$MergedDatasetName],
                                                  JointBarcodeBinsList[["C.briggsae"]][temp_cell_table$MergedDatasetName]))
}

saveRDS(cell_class_barcodes, paste0(dir, "Objects/cell_class_barcodes.rds"))

scenic_zscore_mean_list <- list()
for(species in c("elegans", "briggsae")) {
  cur_species = ifelse(species == "elegans", "C.elegans", "C.briggsae")
  for(version in c("cel", "cbr", "cel_chip")) {
    scenic_zscore_mean_list[[paste0(species, "_", version)]] <- scenic_zscore_list[[paste0(species, "_", version)]][,unlist(cell_data[cell_data$cell_class %in% c("progenitor", "Germline"),]$cell_type_bin)]
    
    for(cur_cell_type in rownames(cell_data[! cell_data$cell_class %in% c("progenitor", "Germline"),])) {
      cur_cell_time_bins <- unlist(cell_data[cur_cell_type,]$cell_type_bin)
      
      temp_zscore <- data.frame(rowMeans(scenic_zscore_list[[paste0(species, "_", version)]][,cur_cell_time_bins, drop = FALSE]))
      colnames(temp_zscore) <- cur_cell_type
      scenic_zscore_mean_list[[paste0(species, "_", version)]] <- cbind(scenic_zscore_mean_list[[paste0(species, "_", version)]],
                                                                        temp_zscore)
    }
  }
}

table(colData(cds_filt[,cell_class_barcodes[["Pharynx and rectal"]]])$embryo.time.bin,
      colData(cds_filt[,cell_class_barcodes[["Pharynx and rectal"]]])$cell_type)
table(colData(cds_filt[,cell_class_barcodes[["Ciliated neurons"]]])$embryo.time.bin)[time_bin_vector]

scenic_auc_df_list <- list()
for(species in c("elegans", "briggsae")) {
  cur_species = ifelse(species == "elegans", "C.elegans", "C.briggsae")
  for(version in c("cel", "cbr", "cel_chip")) {
    dataset = paste0(species, "_", version)
    print(dataset)
    
    time_class_df <- data.frame(matrix(nrow = 0, ncol = 6))
    colnames(time_class_df) <- c("motif", "cell_class", "time_bin", "mean_auc", "median_auc", "dataset")
    
    for(cur_cell_class in unique(cell_data$cell_class[! cell_data$cell_class %in% c("progenitor", "Germline")])) {
      print(cur_cell_class)
      
      temp_cds_pdata <- colData(cds_filt[,cell_class_barcodes[[cur_cell_class]]])
      
      # Find the time bins that both species have
      temp_cds_pdata_cel <- temp_cds_pdata[temp_cds_pdata$species == "C.elegans",]
      temp_cds_pdata_cbr <- temp_cds_pdata[temp_cds_pdata$species == "C.briggsae",]
      
      temp_cds_time_bins_cel <- sort(unique(temp_cds_pdata_cel$embryo.time.bin))
      temp_cds_time_bins_cbr <- sort(unique(temp_cds_pdata_cbr$embryo.time.bin))
      
      temp_cds_time_bins <- intersect(temp_cds_time_bins_cel, temp_cds_time_bins_cbr)
      
      temp_cds_pdata <- temp_cds_pdata[temp_cds_pdata$species == cur_species,]
      
      for(time_bin in time_bin_vector) {
        if(time_bin %in% temp_cds_time_bins) {
          temp_barcodes <- temp_cds_pdata[temp_cds_pdata$embryo.time.bin == time_bin,]$barcode
          
          for(motif in colnames(scenic_auc_list[[dataset]])) {
            temp_auc <- scenic_auc_list[[dataset]][temp_barcodes, motif]
            time_class_df <- rbind(time_class_df, 
                                   data.frame(motif = motif, cell_class = cur_cell_class, time_bin = time_bin,
                                              mean_auc = mean(temp_auc, na.rm = TRUE),
                                              median_auc = median(temp_auc, na.rm = TRUE),
                                              dataset = dataset))
          }
        }
      }
    }
    
    # germline
    time_class_df <- rbind(time_class_df, 
                           data.frame(motif = names(scenic_mean_list[[dataset]]$"germline_pseudotime_0:pseudotime_10"),
                                      cell_class = "Germline", time_bin = "germline_pseudotime_0:pseudotime_10",
                                      mean_auc = unlist(scenic_mean_list[[dataset]]["germline_pseudotime_0:pseudotime_10"]),
                                      median_auc = NA,
                                      dataset = dataset))
    time_class_df <- rbind(time_class_df, 
                           data.frame(motif = names(scenic_mean_list[[dataset]]$"germline_pseudotime_10:pseudotime_13"),
                                      cell_class = "Germline", time_bin = "germline_pseudotime_10:pseudotime_13",
                                      mean_auc = unlist(scenic_mean_list[[dataset]]["germline_pseudotime_10:pseudotime_13"]),
                                      median_auc = NA,
                                      dataset = dataset))
    time_class_df <- rbind(time_class_df, 
                           data.frame(motif = names(scenic_mean_list[[dataset]]$"germline_pseudotime_14:pseudotime_30"),
                                      cell_class = "Germline", time_bin = "germline_pseudotime_14:pseudotime_30",
                                      mean_auc = unlist(scenic_mean_list[[dataset]]["germline_pseudotime_14:pseudotime_30"]),
                                      median_auc = NA,
                                      dataset = dataset))
    scenic_auc_df_list[[dataset]] <- time_class_df
  }
}

#Calculate zscores
for(species in c("elegans", "briggsae")) {
  cur_species = ifelse(species == "elegans", "C.elegans", "C.briggsae")
  for(version in c("cel", "cbr", "cel_chip")) {
    dataset = paste0(species, "_", version)
    print(dataset)
    scenic_auc_df_list[[dataset]]$z_score <- NA
    for(i in 1:nrow(scenic_auc_df_list[[dataset]])) {
      cur_motif = scenic_auc_df_list[[dataset]][i, "motif"]
      scenic_auc_df_list[[dataset]][i, "z_score"] <- z_score(scenic_auc_df_list[[dataset]][i, "mean_auc"],
                                                             scenic_auc_list[[dataset]][,cur_motif])
    }
  }
}

# identify which motifs are present just in the chip set
chip_exclusive <- colnames(scenic_auc_list[["elegans_cel_chip"]])[! colnames(scenic_auc_list[["elegans_cel_chip"]]) %in% 
                                                                    c(colnames(scenic_auc_list[["elegans_cel"]]),
                                                                      colnames(scenic_auc_list[["elegans_cbr"]]))]

# Create a matrix to load to pheatmap
# Cell class
# -> embryo time bin ->
# cel motif auc - cel grn
# cbr motif auc - cel grn
# cel motif auc - cbr grn
# cbr motif auc - cbr grn
# -----------------------
# cel motif auc - cel chip grn
# cbr motif auc - cel chip grn

# each row is a motif, grouped by cell class, ordered by time bin
# cell class, time bin

col_metadata <- data.frame(cell_class = scenic_auc_df_list[["briggsae_cel"]][scenic_auc_df_list[["briggsae_cel"]]$motif == "daf-19(+)",]$cell_class,
                           time_bin = scenic_auc_df_list[["briggsae_cel"]][scenic_auc_df_list[["briggsae_cel"]]$motif == "daf-19(+)",]$time_bin)

out_matrix <- matrix(nrow = 0, ncol = dim(col_metadata)[1])
row_metadata <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(row_metadata) <- c("motif", "species_grn", "species_data", "type")
for(cur_motif in unique(c(scenic_auc_df_list[["elegans_cel"]]$motif,
                        scenic_auc_df_list[["elegans_cbr"]]$motif))) {
  print(cur_motif)
  # cel grn
  # cbr grn
  for(cur_grn in c("cel", "cbr")) {
    if(cur_motif %in% scenic_auc_df_list[[paste0("elegans_", cur_grn)]]$motif) {
      for(cur_species in c("elegans", "briggsae")) {
        out_matrix <- rbind(out_matrix,
                            scenic_auc_df_list[[paste0(cur_species, "_", cur_grn)]][scenic_auc_df_list[[paste0(cur_species, "_", cur_grn)]]$motif == cur_motif,]$z_score)
        row_metadata <- rbind(row_metadata,
                              data.frame(motif = cur_motif,
                                         species_grn = cur_grn,
                                         species_data = cur_species,
                                         type = "motif"))
      }
    }
  }
}

rownames(out_matrix) <- paste0(row_metadata$motif, "_", row_metadata$species_data, "_", row_metadata$species_grn)
colnames(out_matrix) <- paste0(col_metadata$cell_class, "_", col_metadata$time_bin)

rownames(row_metadata) <- paste0(row_metadata$motif, "_", row_metadata$species_data, "_", row_metadata$species_grn)
rownames(col_metadata) <- paste0(col_metadata$cell_class, "_", col_metadata$time_bin)

motif_order <- c("lsl-1(+)",
                 "elt-2(+)", "elt-7(+)", "nhr-109(+)", "tbx-43(+)", # intestine
                 "hlh-1(+)","pat-9(+)", # muscle 
                 "let-381(+)", # mesoderm
                 "ceh-37(+)", # glia and excretory
                 "pha-4(+)", "ceh-22(+)", # pharynx and rectal
                 "elt-1(+)", "elt-3(+)", "nhr-23(+)", # hypodermis and seam
                 "daf-19(+)", "ets-5(+)", "ceh-36(+)", "che-1(+)", # ciliated neurons
                 "unc-3(+)", # non-ciliated neurons
                 "dpl-1(+)", "efl-1(+)", "lin-15B(+)") # misc.
cell_class_order <- c("Germline", "Intestine", "Muscle", "Mesoderm", "Glia and excretory", "Pharynx and rectal", "Hypodermis and seam", "Ciliated neurons", "Non-ciliated neurons")
col_metadata <- col_metadata %>% arrange(match(cell_class,cell_class_order))
row_metadata <- row_metadata %>% arrange(match(motif, motif_order))

cols = list(cell_class = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                           'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                           'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1"),
            species_grn = c('cel' = "#009E73", 'cbr' = "#56B4E9"),
            species_data = c('elegans' = "#009E73", 'briggsae' = "#56B4E9"))

pdf(paste0(dir, "Plots/scenic/motif_z_score_cell_class.pdf"), width = 8, height = 6)
pheatmap(out_matrix[rownames(row_metadata), rownames(col_metadata)],
         annotation_col = col_metadata[,c("cell_class"), drop = FALSE],
         annotation_row = row_metadata[,c("species_data", "species_grn")],
         gaps_col = c(3, 14, 27, 40, 52, 64, 75, 86),
         annotation_colors = cols,
         border_color = "black",
         # border_width = 0.1,
         color = (inferno(1000)),
         breaks = seq(-1, 6, 8/1000),
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         legend = FALSE,
         annotation_legend = FALSE,
         main = "Z-scored AUC")
dev.off()

###
# Now plus the chip-exclusive motifs
col_metadata <- data.frame(cell_class = scenic_auc_df_list[["briggsae_cel"]][scenic_auc_df_list[["briggsae_cel"]]$motif == "daf-19(+)",]$cell_class,
                           time_bin = scenic_auc_df_list[["briggsae_cel"]][scenic_auc_df_list[["briggsae_cel"]]$motif == "daf-19(+)",]$time_bin)

out_matrix <- matrix(nrow = 0, ncol = 98)
row_metadata <- data.frame(matrix(nrow = 0, ncol = 4))
colnames(row_metadata) <- c("motif", "species_grn", "species_data", "type")
for(cur_motif in unique(c(scenic_auc_df_list[["elegans_cel"]]$motif,
                          scenic_auc_df_list[["elegans_cbr"]]$motif))) {
  print(cur_motif)
  # cel grn
  # cbr grn
  for(cur_grn in c("cel", "cbr")) {
    if(cur_motif %in% scenic_auc_df_list[[paste0("elegans_", cur_grn)]]$motif) {
      for(cur_species in c("elegans", "briggsae")) {
        out_matrix <- rbind(out_matrix,
                            scenic_auc_df_list[[paste0(cur_species, "_", cur_grn)]][scenic_auc_df_list[[paste0(cur_species, "_", cur_grn)]]$motif == cur_motif,]$z_score)
        row_metadata <- rbind(row_metadata,
                              data.frame(motif = cur_motif,
                                         species_grn = cur_grn,
                                         species_data = cur_species,
                                         type = "motif"))
      }
    }
  }
}

out_matrix_chip <- out_matrix
row_metadata_chip <- row_metadata

for(cur_motif in chip_exclusive) {
  print(cur_motif)
  for(cur_grn in c("cel_chip")) {
    if(cur_motif %in% scenic_auc_df_list[[paste0("elegans_", cur_grn)]]$motif) {
      for(cur_species in c("elegans", "briggsae")) {
        out_matrix_chip <- rbind(out_matrix_chip,
                            scenic_auc_df_list[[paste0(cur_species, "_", cur_grn)]][scenic_auc_df_list[[paste0(cur_species, "_", cur_grn)]]$motif == cur_motif,]$z_score)
        row_metadata_chip <- rbind(row_metadata_chip,
                              data.frame(motif = cur_motif,
                                         species_grn = cur_grn,
                                         species_data = cur_species,
                                         type = "chip"))
      }
    }
  }
}

rownames(out_matrix_chip) <- paste0(row_metadata_chip$motif, "_", row_metadata_chip$species_data, "_", row_metadata_chip$species_grn)
colnames(out_matrix_chip) <- paste0(col_metadata$cell_class, "_", col_metadata$time_bin)

rownames(row_metadata_chip) <- paste0(row_metadata_chip$motif, "_", row_metadata_chip$species_data, "_", row_metadata_chip$species_grn)
rownames(col_metadata) <- paste0(col_metadata$cell_class, "_", col_metadata$time_bin)

motif_order_chip <- c("lsl-1(+)",
                      "elt-2(+)", "elt-7(+)", "nhr-109(+)", "tbx-43(+)", # intestine
                      "hlh-1(+)","pat-9(+)", # muscle 
                      "let-381(+)", # mesoderm
                      "ceh-37(+)", # glia and excretory
                      "pha-4(+)", "ceh-22(+)", # pharynx and rectal
                      "elt-1(+)", "elt-3(+)", "nhr-23(+)", # hypodermis and seam
                      "daf-19(+)", "ets-5(+)", "ceh-36(+)", "che-1(+)", # ciliated neurons
                      "unc-3(+)", # non-ciliated neurons
                      "dpl-1(+)", "efl-1(+)", "lin-15B(+)", # misc.
                      "madf-6(+)", "madf-8(+)", "ceh-91(+)", "D2030.7(+)",  # germline
                      "nhr-28(+)", "let-607(+)", "nhr-80(+)", "pqm-1(+)",  "crh-2(+)",  # intestine
                      "unc-62(+)", "dsc-1(+)", "hif-1(+)", "unc-120(+)", # muscle
                      "nhr-23(+)", "blmp-1(+)", "nhr-25(+)", "ham-2(+)", # 
                      "fkh-8(+)", "ceh-48(+)",
                      "unc-42(+)") # misc.

cell_class_order <- c("Germline", "Intestine", "Muscle", "Mesoderm", "Glia and excretory", "Pharynx and rectal", "Hypodermis and seam", "Ciliated neurons", "Non-ciliated neurons")
col_metadata <- col_metadata %>% arrange(match(cell_class, cell_class_order))
row_metadata_chip <- row_metadata_chip %>% arrange(match(motif, motif_order_chip))

cols = list(cell_class = c('Ciliated neurons' = "#e15759", 'Germline' = "#bab0ac", 'Glia and excretory' = "#edc948",
                           'Hypodermis and seam' = "#59a14f", 'Intestine' = "#f28e2b", 'Muscle' = "#76b7b2", 'Mesoderm' = "#9c755f",
                           'Non-ciliated neurons' = "#4e79a7", 'Pharynx and rectal' = "#b07aa1"),
            species_grn = c('cel' = "#009E73", 'cbr' = "#56B4E9", 'cel_chip' = "#E69F00"),
            species_data = c('elegans' = "#009E73", 'briggsae' = "#56B4E9"))

# remove hte last two mesoderm and pharynx datasets and the last hyp/seam
col_metadata <- col_metadata[rownames(col_metadata) != paste0("Mesoderm", "_", c("650_710", "gt_710")),]
col_metadata <- col_metadata[rownames(col_metadata) != paste0("Pharynx and rectal", "_", c("650_710", "gt_710")),]
col_metadata <- col_metadata[rownames(col_metadata) != paste0("Hypodermis and seam", "_", c("gt_710")),]

pdf(paste0(dir, "Plots/scenic/motif_z_score_cell_class_chip.pdf"), width = 8, height = 6)
pheatmap(out_matrix_chip[rownames(row_metadata_chip), rownames(col_metadata)],
         border_color = "black",
         annotation_col = col_metadata[,c("cell_class"), drop = FALSE],
         annotation_row = row_metadata_chip[,c("species_data", "species_grn")],
         gaps_col = c(3, 14, 27, 38, 50, 62, 72, 83),
         gaps_row = c(74),
         annotation_colors = cols,
         color = (inferno(100)),
         breaks = seq(-1, 8, 9/100),
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         legend = TRUE,
         fontsize_row = 4,
         annotation_legend = TRUE,
         main = "Z-scored AUC")
dev.off()

# max'd value
pdf(paste0(dir, "Plots/scenic/motif_z_score_cell_class_chip.pdf"), width = 8, height = 6)
pheatmap(t(apply(out_matrix_chip, 1, function(x) x/max(x)))[rownames(row_metadata_chip), rownames(col_metadata)],
         border_color = "black",
         annotation_col = col_metadata[,c("cell_class"), drop = FALSE],
         annotation_row = row_metadata_chip[,c("species_data", "species_grn")],
         gaps_col = c(3, 14, 27, 38, 50, 62, 72, 83),
         gaps_row = c(74),
         annotation_colors = cols,
         color = (inferno(10000)),
         breaks = seq(-0.25, 1, 1.25/10000),
         show_rownames = TRUE,
         show_colnames = FALSE,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         legend = FALSE,
         fontsize_row = 4,
         annotation_legend = TRUE,
         main = "Z-scored AUC")
dev.off()

cor_list <- list()
for(tf in unique(row_metadata_chip$motif)) {
  temp_rowdata <- row_metadata_chip[row_metadata_chip$motif == tf,]
  for(cur_species_grn in unique(temp_rowdata$species_grn)) {
    temp_comparables <- rownames(temp_rowdata[temp_rowdata$species_grn == cur_species_grn,])
    cor_list[[paste0(tf, "_", cur_species_grn)]] <- cor(out_matrix_chip[temp_comparables[1],], out_matrix_chip[temp_comparables[2],])^2
  }
}

###############
# Embryo time #
###############

# grab all of the cells in the ciliated neuron program
# left_join with the cds pdata to grab smoothed.embryo time
# plot

# grab all of the cells in the ciliated neuron program
cell_motif_list_cn <- list()
for(cur_motif in c("daf-19(+)")) {
  print(cur_motif)
  cell_motif_list_cn[[cur_motif]] <- data.frame(cell = cell_class_barcodes[["Ciliated neurons"]],
                                                embryo.time = colData(cds_filt[,cell_class_barcodes[["Ciliated neurons"]]])$smoothed.embryo.time,
                                                embryo.time.bin = colData(cds_filt[,cell_class_barcodes[["Ciliated neurons"]]])$embryo.time.bin,
                                                auc_cel = rbind(scenic_auc_list[["elegans_cel"]], scenic_auc_list[["briggsae_cel"]])[cell_class_barcodes[["Ciliated neurons"]], cur_motif],
                                                auc_cbr = rbind(scenic_auc_list[["elegans_cbr"]], scenic_auc_list[["briggsae_cbr"]])[cell_class_barcodes[["Ciliated neurons"]], cur_motif],
                                                motif = cur_motif,
                                                species = colData(cds_filt[,cell_class_barcodes[["Ciliated neurons"]]])$species)
}

for(cur_motif in c("fkh-8(+)")) {
  print(cur_motif)
  cell_motif_list_cn[[cur_motif]] <- data.frame(cell = cell_class_barcodes[["Ciliated neurons"]],
                                                embryo.time = colData(cds_filt[,cell_class_barcodes[["Ciliated neurons"]]])$smoothed.embryo.time,
                                                embryo.time.bin = colData(cds_filt[,cell_class_barcodes[["Ciliated neurons"]]])$embryo.time.bin,
                                                auc_cel = rbind(scenic_auc_list[["elegans_cel_chip"]], scenic_auc_list[["briggsae_cel_chip"]])[cell_class_barcodes[["Ciliated neurons"]], cur_motif],
                                                auc_cbr = NA,
                                                motif = cur_motif,
                                                species = colData(cds_filt[,cell_class_barcodes[["Ciliated neurons"]]])$species)
}

pdf(paste0(dir, "Plots/scenic/daf_19_fkh_8_time_bins_cel.pdf"), width = 5, height = 5)
do.call(rbind, cell_motif_list_cn[c("fkh-8(+)", "daf-19(+)")]) %>%
  filter(! embryo.time.bin %in% c("gt_710", "lt_100")) %>%
  ggplot(aes(x = embryo.time.bin, y = auc_cel,
             group = paste0(motif, "_", species), color = species,
             linetype = motif, fill = species)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_normal, alpha = 0.5, linetype = "solid") +
  stat_summary(geom = "line", fun.y = mean) +
  stat_summary(geom = "point", fun.y = mean) +
  scale_y_continuous(name = "Regulon activity (AUC)") +
  scale_x_discrete(name = "Embryo time bin", limits = time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$bins,
                   labels = paste0(time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$start, " to ",
                                   time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$end)) +
  scale_color_manual(breaks = c("C.elegans", "C.briggsae"), values = c("#009E73", "#56B4E9")) +
  scale_fill_manual(breaks = c("C.elegans", "C.briggsae"), values = c("#009E73", "#56B4E9")) +
  guides(color = guide_legend(ncol = 1),
         linetype = guide_legend(ncol = 1)) +
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
dev.off()

for(cur_motif in c("elt-7(+)", "nhr-109(+)")) {
  print(cur_motif)
  cell_motif_list_cn[[cur_motif]] <- data.frame(cell = cell_class_barcodes[["Intestine"]],
                                                embryo.time = colData(cds_filt[,cell_class_barcodes[["Intestine"]]])$smoothed.embryo.time,
                                                embryo.time.bin = colData(cds_filt[,cell_class_barcodes[["Intestine"]]])$embryo.time.bin,
                                                auc_cel = rbind(scenic_auc_list[["elegans_cel"]], scenic_auc_list[["briggsae_cel"]])[cell_class_barcodes[["Intestine"]], cur_motif],
                                                auc_cbr = NA, # rbind(scenic_auc_list[["elegans_cbr"]], scenic_auc_list[["briggsae_cbr"]])[cell_class_barcodes[["Ciliated neurons"]], cur_motif],
                                                motif = cur_motif,
                                                species = colData(cds_filt[,cell_class_barcodes[["Intestine"]]])$species)
}

pdf(paste0(dir, "Plots/scenic/elt_7_time_bins_cel.pdf"), width = 5, height = 5)
do.call(rbind, cell_motif_list_cn[c("elt-7(+)")]) %>%
  filter(! embryo.time.bin %in% c("gt_710", "lt_100")) %>%
  ggplot(aes(x = embryo.time.bin, y = auc_cel, group = paste0(species, "_", motif),
             color = paste0(species, "_", motif), fill = paste0(species, "_", motif))) +
  stat_summary(geom="ribbon", fun.data = mean_cl_normal, width=0.1, conf.int=0.95, alpha = 0.25) +
  stat_summary(geom="line", fun.y = mean, linetype = "dashed") +
  stat_summary(geom="point", fun.y = mean) +
  scale_y_continuous(name = "Regulon activity (AUC)") +
  scale_x_discrete(name = "Embryo time bin", limits = time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100", "100_130"),]$bins,
                     labels = paste0(time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$start, " to ",
                                     time_bin_df[! time_bin_df$bins %in% c("gt_710", "lt_100"),]$end)) +
  scale_color_manual(breaks = c("C.elegans_elt-7(+)", "C.briggsae_elt-7(+)"), values = c("#009E73", "#56B4E9")) +
  scale_fill_manual(breaks = c("C.elegans_elt-7(+)", "C.briggsae_elt-7(+)"), values = c("#009E73", "#56B4E9")) +
  theme(legend.title = element_blank(),
        legend.position = "none",
        rect = element_rect(fill = "transparent"),
        # axis.title.x = element_blank(),
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


cel_regs_df_filt
cbr_regs_df_filt
cel_regs_df_chip_filt
TPMLineageList[["C.elegans"]][cel_regs_df_filt[cel_regs_df_filt$tf == "daf-19","genes"], temp_time_class_bin]

temp_time_class_bin <- data.frame(time_class_metadata_df %>%
  filter(cell_class %in% "Ciliated neurons" & species == "C.elegans") %>%
  reframe(class_time = paste0(cell_class, "_", time_bin),
          time_bin = time_bin,
          cell_class = cell_class))

motif_of_interest <- "daf-19"

cel_temp_tpm <- TPMLineageList[["C.elegans"]][cel_regs_df_filt[cel_regs_df_filt$tf == motif_of_interest, "genes"], temp_time_class_bin[,"class_time"]]
colnames(cel_temp_tpm) <- temp_time_class_bin$time_bin
cel_temp_tpm$gene <- rownames(cel_temp_tpm)
cel_temp_tpm$species <- "C.elegans"
cel_temp_tpm$motif <- motif_of_interest

cbr_temp_tpm <- TPMLineageList[["C.briggsae"]][cel_regs_df_filt[cel_regs_df_filt$tf == motif_of_interest,"genes"], temp_time_class_bin[,"class_time"]]
colnames(cbr_temp_tpm) <- temp_time_class_bin$time_bin
cbr_temp_tpm$gene <- rownames(cbr_temp_tpm)
cbr_temp_tpm$species <- "C.briggsae"
cbr_temp_tpm$motif <- motif_of_interest

time_bin_vector <- c("lt_100", "100_130", "130_170", "170_210", "210_270", "270_330", "330_390", "390_450", "450_510", "510_580", "580_650", "650_710", "gt_710")

temp_tpm <- rbind(cel_temp_tpm, cbr_temp_tpm) %>% melt(id.vars = c("gene", "species", "motif")) 

temp_tpm$time <- match(temp_tpm$variable, time_bin_vector)

##############
# Enrichments
##############

gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))
rownames(gff_list_mrna[["elegans"]]) <- ifelse(is.na(gff_list_mrna[["elegans"]]$cds_gene_name), gff_list_mrna[["elegans"]]$gene_name, gff_list_mrna[["elegans"]]$cds_gene_name)
rownames(gff_list_mrna[["briggsae"]]) <- ifelse(is.na(gff_list_mrna[["briggsae"]]$cds_gene_name), gff_list_mrna[["briggsae"]]$gene_name, gff_list_mrna[["briggsae"]]$cds_gene_name)

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

cel_regs_df_filt
cbr_regs_df_filt
cel_regs_df_chip_filt

background_gene_set <- gene_data[gene_data$max_tpm_pro > 80 | gene_data$max_tpm_term > 80,]$gene

enrich_list <- list()
for(tf in unique(c(cel_regs_df_filt$tf, cel_regs_df_chip_filt$tf))) {
  if(tf %in% unique(cel_regs_df_filt$tf)) {
    rgs = cel_regs_df_filt[cel_regs_df_filt$tf == tf,]$genes
  } else {
    rgs = cel_regs_df_chip_filt[cel_regs_df_chip_filt$tf == tf,]$genes
  }
  
  enrich_list[[tf]] <- wormcat_enrichment(rgs, background_gene_set, gff_list_mrna[["elegans"]])
  
  for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
    if(nrow(enrich_list[[tf]][[tier]]) > 0) {
      enrich_list[[tf]][[tier]]$tf <- tf
    }
  }
}


enrich_list_filt <- list()
for(tier in c("WormCat.1", "WormCat.2", "WormCat.3")) {
  enrich_list_filt[[tier]] <- do.call(rbind, lapply(enrich_list, function(x) x[[tier]]))
  temp_filt <- data.frame(enrich_list_filt[[tier]] %>%
                            filter(Bonferroni < 0.05) %>%
                            select(Category) %>%
                            unique())[,1]
  
  enrich_list_filt[[tier]] <- enrich_list_filt[[tier]] %>%
    filter(Category %in% temp_filt) %>%
    complete(Category, tf)
  
  enrich_list_filt[[tier]][is.na(enrich_list_filt[[tier]]$Fold),]$Fold <- 0
  enrich_list_filt[[tier]][is.na(enrich_list_filt[[tier]]$Bonferroni),]$Bonferroni <- 1
  
  enrich_list_filt[[tier]]$signif <- ""
  enrich_list_filt[[tier]]$signif <- ifelse(enrich_list_filt[[tier]]$Bonferroni < 0.05,
                                            "*",
                                            enrich_list_filt[[tier]]$signif)
  enrich_list_filt[[tier]]$signif <- ifelse(enrich_list_filt[[tier]]$Bonferroni < 0.005,
                                            "**",
                                            enrich_list_filt[[tier]]$signif)
  enrich_list_filt[[tier]]$signif <- ifelse(enrich_list_filt[[tier]]$Bonferroni < 0.0005,
                                            "***",
                                            enrich_list_filt[[tier]]$signif)
}

tf_to_use <- unique(enrich_list_filt[["WormCat.1"]][enrich_list_filt[["WormCat.1"]]$signif != "",]$tf)

enrich.df.reshape <- dcast(enrich_list_filt[["WormCat.1"]] %>%
                             filter(tf %in% tf_to_use), Category ~ tf, value.var = "Fold", fun.aggregate = mean)
rownames(enrich.df.reshape) <- enrich.df.reshape$Category

enrich.matrix <- as.matrix(enrich.df.reshape[-1])

row_order <- rownames(enrich.matrix[hclust(dist(enrich.matrix, method = "euclidean"), method = "complete")$order,])
col_order <- colnames(enrich.matrix[,hclust(dist(t(enrich.matrix), method = "euclidean"), method = "complete")$order])

pdf(paste0(dir, "Plots/scenic/enrich_all.pdf"), height = 6, width = 5.25)
enrich_list_filt[["WormCat.1"]] %>%
  filter(tf %in% tf_to_use) %>%
  ggplot(aes(x = Category, y = tf, fill = Fold)) +
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

tf_to_use <- c("daf-19", "fkh-8", "ceh-48", "hlh-1", "pat-9", "efl-1", "dpl-1", "elt-2", "elt-3", "pha-4")
cat_to_use <- c("Cilia", "Neuronal function", "Muscle function", "Cell cycle", "DNA", "Development", "Stress response", "Peroxisome", "Extracellular material", "Transcription factor")
pdf(paste0(dir, "Plots/scenic/enrich_subset.pdf"), height = 4, width = 3)
enrich_list_filt[["WormCat.1"]] %>%
  filter(tf %in% tf_to_use & Category %in% cat_to_use) %>%
  ggplot(aes(x = Category, y = tf, fill = Fold)) +
  geom_tile(color="grey80",  linewidth = 0.1) +
  geom_text(aes(label = signif), vjust = 0.65, hjust = 0.5, size = 2) +
  geom_segment(x = 0, xend = Inf, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = Inf, y = Inf, yend = Inf, color = "grey80", linewidth = 2) +
  geom_segment(x = 0, xend = 0, y = 0, yend = 0, color = "grey80", linewidth = 2) +
  geom_segment(x = Inf, xend = Inf, y = 0, yend = Inf, color = "grey80", linewidth = 2) +
  scale_x_discrete(guide = guide_axis(angle = 45), limits = cat_to_use) +
  scale_y_discrete(limits = tf_to_use) +
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
        legend.position = "top")
dev.off()

