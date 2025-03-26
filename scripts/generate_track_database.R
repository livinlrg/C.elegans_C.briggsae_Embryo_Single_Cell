
# Using chip-seq data from the modencode and modern projects found here:
# https://epic.gs.washington.edu/modERNresource/
# downloaded most recently 2024_09_27
# chip-seq tracks are created using the identified peaks
# then an upstream target bed file is created using the peaks

library(rtracklayer)
library(dplyr)
dir <- "/kimdata/livinlrg/scAnalysis/WS290/"

## for the chip-seq data, the goal will be to seperate the peak file with primary targets into a file for each experiment
peaks <- read.table(paste0(dir, "Objects/scenic/WormTFPeaksPrimaryTargets.tsv"), header = TRUE, sep = "\t")

elegans_gtf <- data.frame(rtracklayer::import("/kimdata/livinlrg/scAnalysis/Reference/elegans/C_elegans.ws290.extendedUTRv2.format.gtf"))
briggsae_gtf <- data.frame(rtracklayer::import("/kimdata/livinlrg/scAnalysis/Reference/briggsae/C_briggsae.ws290.extendedUTRv3.format.gtf"))
gff_list_mrna <- readRDS(paste0(dir, "Objects/gff_list_mrna.rds"))

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_no_bg_20240828.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

rm(cds)

# want to make sure all of the gene names are in my mrna gff file
gene_names <- rowData(cds_filt)[rowData(cds_filt)$gene.type %in% c("common", "ce.unique"),]$gene_short_name
gene_names[! (gene_names %in% gff_list_mrna[["elegans"]]$out_name)]

# first need to generate bedgraph files for each of the chip experiments
# I will likely use the experiment name as the prefix for the bed file
# I will also need to generate a bed file for the primary targets
# if there are no primary targets for that gene, I still likely need to create an entry?

# bedgraph is:
# chrom start end dataValue

dir.create(paste0(dir, "Objects/scenic/bedGraphFiles/"))

for(exp in unique(peaks$experiment)) {
  exp_peaks <- peaks[peaks$experiment == exp,]
  exp_bed <- data.frame(seqnames = exp_peaks$chrom,
                     start = exp_peaks$chromStart,
                     end = exp_peaks$chromEnd,
                     score = exp_peaks$signalValue)
 write.table(exp_bed, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t",
           file = paste0(dir, "Objects/scenic/bedGraphFiles/", exp, ".bedgraph"))
}

peaks_plus <- left_join(peaks, elegans_gtf[,c("gene_id", "protein_id", "transcript_id")],
                        by = c("Transcript" = "transcript_id"),
                        na_matches = "never")
peaks_plus <- peaks_plus[! duplicated(peaks_plus$peakID),]

for(cur_gene in unique(peaks_plus[is.na(peaks_plus$gene_id),]$Gene)) {
  if(cur_gene %in% gff_list_mrna[["elegans"]]$short_name | cur_gene %in% gff_list_mrna[["elegans"]]$out_name) {
    peaks_plus[peaks_plus$Gene %in% cur_gene,]$gene_id <- gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$short_name %in% cur_gene |
                                                                                       gff_list_mrna[["elegans"]]$out_name %in% cur_gene,]$gene_name
  }
  # print(cur_gene %in% gff_list_mrna[["elegans"]]$short_name | cur_gene %in% gff_list_mrna[["elegans"]]$out_name)
}

unique(peaks_plus[is.na(peaks_plus$gene_id),]$Gene)
'
Remaining genes that cannot be matched:
 [1] "K11E4.2"   "C31C9.2"   "C14B1.9"   "ZK112.5"   "K02E10.1"  "T19B4.3"   "R05D3.2"  
 [8] "Y70G10A.2" "bus-8"     "C09B8.5"   "Y71H2AM.2"
'
peaks_plus[peaks_plus$Gene %in% "K11E4.2",]$gene_id <- "WBGene00010774"
peaks_plus[peaks_plus$Gene %in% "C31C9.2",]$gene_id <- "WBGene00007836"
peaks_plus[peaks_plus$Gene %in% "C14B1.9",]$gene_id <- "WBGene00007580"
peaks_plus[peaks_plus$Gene %in% "ZK112.5",]$gene_id <- "WBGene00022662" # was dead and merged into acs-25
peaks_plus[peaks_plus$Gene %in% "K02E10.1",]$gene_id <- "WBGene00019317"
peaks_plus[peaks_plus$Gene %in% "T19B4.3",]$gene_id <- "WBGene00020557"
peaks_plus[peaks_plus$Gene %in% "R05D3.2",]$gene_id <- "WBGene00019877"
peaks_plus[peaks_plus$Gene %in% "Y70G10A.2",]$gene_id <- "WBGene00013498"
peaks_plus[peaks_plus$Gene %in% "bus-8",]$gene_id <- "WBGene00044623" # merged into bus-8b
peaks_plus[peaks_plus$Gene %in% "C09B8.5",]$gene_id <- "WBGene00015624"
peaks_plus[peaks_plus$Gene %in% "Y71H2AM.2",]$gene_id <- "WBGene00022167"

elegans_genome <- read.table(paste0(dir, "Objects/scenic/elegans_genome.bed"),
                             header = TRUE)
rownames(elegans_genome) <- elegans_genome$chr
rownames(gff_list_mrna[["elegans"]]) <- gff_list_mrna[["elegans"]]$gene_name

bed_out <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(bed_out) <- c("chrom", "start", "end", "gene")
for(cur_gene in unique(elegans_gtf$gene_id)) {
  if(cur_gene %in% peaks_plus$gene_id) {
    cur_peak <- peaks_plus[peaks_plus$gene_id %in% cur_gene, , drop = FALSE]
    
    bed_out <- rbind(bed_out, data.frame(chrom = unique(cur_peak$chrom),
                                         start = min(cur_peak$metaPeakStart),
                                         end = max(cur_peak$metaPeakEnd),
                                         gene = gff_list_mrna[["elegans"]][cur_gene,]$out_name))
  } else {
    # print("no peak")
    elegans_gtf_cur_gene <- gff_list_mrna[["elegans"]][cur_gene,]
    if(elegans_gtf_cur_gene$strand == "-") {
      bed_out <- rbind(bed_out, data.frame(chrom = elegans_gtf_cur_gene$seqnames,
                                           start = elegans_gtf_cur_gene$start,
                                           end = elegans_gtf_cur_gene$end + 1000,
                                           gene = elegans_gtf_cur_gene$out_name))
    } else {
      bed_out <- rbind(bed_out, data.frame(chrom = elegans_gtf_cur_gene$seqnames,
                                           start = elegans_gtf_cur_gene$start - 1000,
                                           end = elegans_gtf_cur_gene$end,
                                           gene = elegans_gtf_cur_gene$out_name))
    }
    
    # check for out of bounds issues
    if(bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$start < 0) {
      bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$start = 0
    } else if (bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$end > as.numeric(elegans_genome[as.character(bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$chrom),]$end)) {
      bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$end = as.numeric(elegans_genome[as.character(bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$seqnames),]$end)
    }
  }
}

options("scipen"=100, "digits"=4)
colnames(bed_out) <- c('#chrom', 'chromStart', 'chromEnd', 'name')
bed = bed_out %>%
  mutate(across(everything(), as.character))

write.table(bed, paste0(dir, "Objects/scenic/WS290_track.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## Also try with a bed file that has each meatapeak as a separate line
bed_out <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(bed_out) <- c("chrom", "start", "end", "gene")
for(cur_gene in unique(elegans_gtf$gene_id)) {
  if(cur_gene %in% peaks_plus$gene_id) {
    cur_peak <- peaks_plus[peaks_plus$gene_id %in% cur_gene, , drop = FALSE]
    
    bed_out <- rbind(bed_out, data.frame(chrom = gsub("chr", "", cur_peak$chrom),
                                         start = cur_peak$metaPeakStart,
                                         end = cur_peak$metaPeakEnd,
                                         gene = rep(gff_list_mrna[["elegans"]][cur_gene,]$out_name,
                                                    nrow(cur_peak))))
  } else {
    # print("no peak")
    elegans_gtf_cur_gene <- gff_list_mrna[["elegans"]][cur_gene,]
    if(elegans_gtf_cur_gene$strand == "-") {
      bed_out <- rbind(bed_out, data.frame(chrom = elegans_gtf_cur_gene$seqnames,
                                           start = elegans_gtf_cur_gene$start,
                                           end = elegans_gtf_cur_gene$end + 1000,
                                           gene = elegans_gtf_cur_gene$out_name))
    } else {
      bed_out <- rbind(bed_out, data.frame(chrom = elegans_gtf_cur_gene$seqnames,
                                           start = elegans_gtf_cur_gene$start - 1000,
                                           end = elegans_gtf_cur_gene$end,
                                           gene = elegans_gtf_cur_gene$out_name))
    }
    
    # check for out of bounds issues
    if(bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$start < 0) {
      bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$start = 0
    } else if (bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$end > as.numeric(elegans_genome[as.character(bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$chrom),]$end)) {
      bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$end = as.numeric(elegans_genome[as.character(bed_out[bed_out$gene %in% elegans_gtf_cur_gene$out_name,]$seqnames),]$end)
    }
  }
}

# Required libraries
library(GenomicRanges)

bed_out$chrom <- gsub("chr", "", bed_out$chrom)
bed_out <- bed_out[bed_out$chrom != "MtDNA",]
# I II III IV V X

# Function to merge intervals for the same gene
merge_intervals <- function(df) {
  gr <- GRanges(seqnames = df$chrom, ranges = IRanges(start = df$start, end = df$end))
  merged_gr <- reduce(gr) # Merge overlapping ranges
  data.frame(chrom = as.character(seqnames(merged_gr)),
             start = start(merged_gr),
             end = end(merged_gr),
             gene = df$gene[1]) # Keep the gene name
}

# Apply merging function to each gene
merged_bed <- bed_out %>%
  arrange(chrom, start, gene) %>%
  group_by(gene) %>%
  do(merge_intervals(.)) 

merged_sorted_bed <- merged_bed %>%
  arrange(chrom, start, gene) %>%
  group_by(gene) %>%
  mutate(gene = paste0(gene, "_", row_number()))

# Convert the result back to data frame
merged_sorted_bed <- as.data.frame(merged_sorted_bed)

merged_sorted_bed$chrom <- paste0("chr", merged_sorted_bed$chrom)

options("scipen"=100, "digits"=4)
colnames(merged_sorted_bed) <- c('#chrom', 'chromStart', 'chromEnd', 'name')
bed = merged_sorted_bed %>%
  mutate(across(everything(), as.character))

write.table(bed, paste0(dir, "Objects/scenic/WS290_track_multipeak.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# also export the info for the chipseq experiments appended onto the end of a file:
# exp TF exp exp exp "1" "0" "None" "None" "1" "None" "None" "gene is directly annotated"
tf_info_out <- data.frame(matrix(ncol = 13, nrow = 0))
for(exp in unique(peaks_plus$experiment)) {
  temp_exp <- peaks_plus[peaks_plus$experiment %in% exp,]
  temp_exp <- temp_exp[! duplicated(temp_exp$experiment),]
  tf_info_out <- rbind(tf_info_out, data.frame(
    temp_exp$exp,
    temp_exp$TF,
    temp_exp$exp,
    temp_exp$exp,
    temp_exp$exp,
    "1",
    "0",
    "None",
    "None",
    "1",
    "None",
    "None",
    "gene is directly annotated"
  ))
}

write.table(tf_info_out, paste0(dir, "Objects/scenic/tf_info_out.tbl"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(tf_info_out[,1], paste0(dir, "Objects/scenic/track_file_names.txt"),
            sep = "\n", quote = FALSE, row.names = FALSE, col.names = FALSE)

