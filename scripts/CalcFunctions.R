library(monocle3)
library(cowplot)
library(ggtext)
library(rtracklayer) # for gene view
library(Gviz) # for 
library(ggplotify)

# get.norm.expr.matrix
# CalcTPM
# CalcTPMRobust
# trinarize
# trinarize_cds
# jaccard

# Normalize expression matrix by Size factor
get.norm.expr.matrix <- function(cds) {
  mat = counts(cds)
  mat@x = mat@x / rep.int(pData(cds)$Size_Factor, diff(mat@p))
  return(mat)
}

# Basic TPM Calculation
CalcTPM <- function(cur_cds, column) {
  norm.expr = get.norm.expr.matrix(cur_cds) # normalize by size factor
  cell.bins <- unique(pData(cur_cds)[[column]])[! is.na(unique(pData(cur_cds)[[column]]))]
  
  bin.norm.means = sapply(cell.bins, function(this.bin) {
    if(sum(pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]) < 2) {
      print("Yes")
      norm.expr[,pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]]
    } else {
      Matrix::rowMeans(norm.expr[,pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]])
    }
  }) #for each cell.bin do the rowMean

  bin.tpm = base::sweep(
    bin.norm.means, 2,
    Matrix::colSums(bin.norm.means), "/") * 1000000
  
  return(bin.tpm)
}

# Basic TPM Calculation without size factor correction
CalcTPM_without_size <- function(cur_cds, column) {
  expr = exprs(cur_cds) # normalize by size factor
  cell.bins <- unique(pData(cur_cds)[[column]])[! is.na(unique(pData(cur_cds)[[column]]))]
  
  bin.norm.means = sapply(cell.bins, function(this.bin) {
    if(sum(pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]) < 2) {
      print("Yes")
      expr[,pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]]
    } else {
      Matrix::rowMeans(expr[,pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]])
    }
  }) #for each cell.bin do the rowMean
  
  bin.tpm = base::sweep(
    bin.norm.means, 2,
    Matrix::colSums(bin.norm.means), "/") * 1000000
  
  return(bin.tpm)
}

# Robust TPM calculation that excludes the top and bottom sorted values
CalcTPMRobust <- function(cur_cds, column, cur_species) {
  norm.expr = get.norm.expr.matrix(cur_cds) # normalize by size factor
  
  cell.bins <- unique(pData(cur_cds)[[column]])[! is.na(unique(pData(cur_cds)[[column]]))]
  
  bin.norm.means.robust = sapply(cell.bins, function(this.bin) {
    message(this.bin)
    if(sum(pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]) < 2) {
      norm.expr[,pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]]
    } else {
      tmp.mat = as.matrix(norm.expr[,pData(cur_cds)[[column]] %in% cell.bins[cell.bins == this.bin]])
      apply(tmp.mat, 1, function(x) {
        x = sort(x);
        mean(x[2:(length(x)-1)]) # try to make estimates more robust against outliers
      })
    }
  }) #for each cell.bin do the rowMean
  colnames(bin.norm.means.robust) = names(cell.bins)
  
  bin.tpm.robust = sweep(
    bin.norm.means.robust, 2,
    Matrix::colSums(bin.norm.means.robust), "/") * 1000000
  
  rownames(bin.tpm.robust) = fData(cur_cds)$gene_short_name
  
  pb_TPM_log <- data.frame(apply(bin.tpm.robust, 2, function(x) {log2(x + 1)}))
  
  return(pb_TPM_log)
}

trinarize <- function(a = 1.5, b = 2, f  = 0.2, n, k) {
  # Calculate confidence in measurements
  # We calculate probability p that at least f of the cells detect (in each group),
  # and compare with pep, setting the binary pattern to 1 if p > pep,
  # -1 if p < (1 - pep) and 0 otherwise.
  # Args:
  #   k (int):	Number of observed positive cells
  #   n (int):	Total number of cells
  #   a (int):  Hyperparameter for beta(a,b)
  #   b (int):  Hyperparameter for beta(a,b)
  
  incb = pbeta(f, a + k, b - k + n, lower.tail = TRUE)
  if(incb == 0) {
    p = 1
  } else {
    p = 1 - exp(log(incb) + lbeta(a + k, b - k + n) + lgamma(a + b + n) - lgamma(a + k) - lgamma(b - k + n))
  }
  return(p)
}

# Out will be similar to TPMList[[species]] wherein
# CellType by gene matrix with values set at
# Expressed = 1, NotExpressed = -1, Intermediate = 0
# PEP FDR threshold
# f Expressed in x% of cells

trinarize_cds <- function(cur_cds, column, PEP = 0.05, f = 0.10) {
  # Wraps trinarize on a cds per cell.type based on a cell.type column
  # Depending on posterior probability, counts as detected, not detected, 
  # or intermediate as 1, -1, or 0 respectively
  #
  # Args:
  #   cur_cds (cds):	Cds to be analyzed
  #   column (char):	Column of cell.type info in cds pData
  #
  cell.bins <- unique(pData(cur_cds)[[column]])[! is.na(unique(pData(cur_cds)[[column]]))]
  
  gene_list <- rownames(cur_cds)
  bin.binarize <- data.frame(matrix(ncol = length(cell.bins), nrow = length(gene_list)))
  
  colnames(bin.binarize) <- cell.bins
  rownames(bin.binarize) <- gene_list
  
  for(cell_type in cell.bins) {
    print(paste0("Starting: ", cell_type))

    temp_counts = counts(cur_cds[,which(colData(cur_cds)[[column]] == cell_type)])
    total_cells = length(temp_counts[1,])
    
    bin.binarize[[cell_type]] <- sapply(rowSums(temp_counts > 0), function(non_zero_counts) {
      p = trinarize(n = total_cells, k = non_zero_counts, f = f)
      
      if(p > 1 - PEP) { # Detected
        return(1)
      } else if(p < PEP) { # Not detected
        return(-1)
      } else { # Else is indeterminate
        return(0)
      }
    })
  }
  return(bin.binarize)
}

# Jaccard Index
# Expression in cell type
# Go by gene, compare two groups
jaccard <- function(gene_a, gene_b) {
  a <- names(gene_a)[gene_a > 0]
  b <- names(gene_b)[gene_b > 0]
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

js_divg <- function(p, q) {
  r <- 0.5 * (p + q)
  # Define each samples component
  P1 <- p * log2(p / r)
  P2 <- q * log2(q / r)
  # In the case of zeroes entries log is undefined, JSD is defined as zero
  P1[! is.finite(P1)] <- 0
  P2[! is.finite(P2)] <- 0
  
  d <- (P1 + P2) / 2
  return(sum(d, na.rm = TRUE))
}

cos_dist <- function(p, q) {
  return(1 - (sum(p*q) / (sqrt(sum(p^2)) * sqrt(sum(q^2)))))
}

## For plotting number in a category on a plot
n_fun <- function(x) {
  return(data.frame(y = 0, label = length(x)))
}

plot_read_track <- function(gene, window = 6000) {
  options(ucscChromosomeNames=FALSE)
  if(! exists("elegans_gff")) {
    print("Loading elegans gff")
    elegans_gff <- import.gff(paste0(ExternalPath, "LaDeana/C_elegans_WS260_3p_UTR_extended_canonical_common_names_jonathan.gff"))
  }
  
  if(! exists("briggsae_gtf")) {
    print("Loading briggsae gtf")
    briggsae_gtf <- rtracklayer::import(paste0(NewBobPath,"BigWigs/Cbriggsae-filtered-ensemble.opt.flattened.plus_MtDNA.format.plusCDS.fixed.gt"))
    
    briggsae_gtf$gene_name_short <- sapply(strsplit(briggsae_gtf$gene_name, "-"), function(x) {
      if(length(x) > 1) {
        paste0(x[2:length(x)], collapse ="-")
      } else {
        x
      }
    }
    )
  }
  
  if(! exists("elegans_embryo_mutant.bottom.bw")) {
    print("Loading big wig files")
    elegans_bulk <- import(paste0(NewBobPath, "BigWigs/elegans_modENCODE_RNAseq.rdpm.WS245.bw"))
    elegans_embryo_mutant.bottom.bw <- import(paste0(NewBobPath, "BigWigs/embryo_mutant.bottom.rdpm.WS245.bw"))
    elegans_embryo_mutant.top.bw <- import(paste0(NewBobPath, "BigWigs/embryo_mutant.top.rdpm.WS245.bw"))
    briggsae_embryo_mutant.bottom.bw <- import(paste0(NewBobPath, "BigWigs/Cb_-_120to500.rdpm.bw"))
    briggsae_embryo_mutant.top.bw <- import(paste0(NewBobPath, "BigWigs/Cb_+_120to500.rdpm.bw"))
    briggsae_embryo_bulk <- import(paste0(NewBobPath, "BigWigs/AF16_AG1109.raw_normalized_covg.cb4.wig.bw"))
  }
  
  # Do elegans first
  ele_gene_query <- data.frame(elegans_gff[elegans_gff$gene_name == gene,])
  
  ele_query_range = GRanges(seqnames = unique(ele_gene_query$seqnames),
                  ranges=IRanges(start = min(c(ele_gene_query$start, ele_gene_query$end)) - window,
                                 end = max(c(ele_gene_query$start, ele_gene_query$end)) + window))
  
  ele_temp <- data.frame(subsetByOverlaps(elegans_gff, ele_query_range))
  colnames(ele_temp) <- c("chromosome", "start", "end", "width", "strand", "source", "feature", "score", "phase", "gene", "gene_source", "symbol", "transcript")
  
  ele_temp <- ele_temp[ele_temp$feature %in% c("exon"),]
  ele_temp$feature <- as.character(factor(ele_temp$feature))
  
  ele_gtrack <- GenomeAxisTrack()
  ele_dTrack_top <- DataTrack(range = subsetByOverlaps(elegans_embryo_mutant.top.bw, ele_query_range), genome = "elegans", name = "Single-cell top")
  ele_dTrack_bottom <- DataTrack(range = subsetByOverlaps(elegans_embryo_mutant.bottom.bw, ele_query_range), genome = "elegans", name = "Single-cell bottom")
  ele_grtrack <- GeneRegionTrack(ele_temp, name = "Gene Model", genome = "elegans", showId = TRUE, fill = "#009E73", col = "black", transcriptAnnotation = "transcript")
  
  ele_hlight <- HighlightTrack(trackList = list(ele_dTrack_top, ele_dTrack_bottom),
                           start = min(c(ele_temp[ele_temp$symbol == gene,]$start, ele_temp[ele_temp$symbol == gene,]$end)),
                           end = max(c(ele_temp[ele_temp$symbol == gene,]$start, ele_temp[ele_temp$symbol == gene,]$end)),
                           chromosome = unique(ele_gene_query$seqnames), fill = "transparent", col = "black", lty = "dashed")
  
  ele_gene_track <- as.grob(~plotTracks(list(ele_gtrack, ele_hlight, ele_grtrack), type = "hist", background.title = "#009E73"))
  
  # Then briggsae
  bri_gene_query <- data.frame(briggsae_gtf[briggsae_gtf$gene_name_short == gene,])
  
  bri_query_range = GRanges(seqnames = as.character(unique(bri_gene_query$seqnames)),
                            ranges=IRanges(start = min(c(bri_gene_query$start,bri_gene_query$end)) - window,
                                           end = max(c(bri_gene_query$start,bri_gene_query$end)) + window))
  
  bri_temp <- data.frame(subsetByOverlaps(briggsae_gtf, bri_query_range))
  colnames(bri_temp) <- c("chromosome", "start", "end", "width", "strand", "source", "type", "score", "phase", "gene_id", "transcript_id", "gene_name", "gene_source", "gene_biotype", "transcript", "transcript_source", "transcript_biotype", "exon_number", "exon_id", "exon_version", "protein_id", "protein_version", "tag", "gene")
  
  bri_temp <- bri_temp[bri_temp$type %in% c("exon"),]
  bri_temp$type <- as.character(factor(bri_temp$type))
  
  bri_gtrack <- GenomeAxisTrack()
  bri_dTrack_top <- DataTrack(range = subsetByOverlaps(briggsae_embryo_mutant.top.bw, bri_query_range), genome = "briggsae", name = "Single-cell top")
  bri_dTrack_bottom <- DataTrack(range = subsetByOverlaps(briggsae_embryo_mutant.bottom.bw, bri_query_range), genome = "briggsae", name = "Single-cell bottom")
  bri_dTrack_bulk <- DataTrack(range = subsetByOverlaps(briggsae_embryo_bulk, bri_query_range), genome = "briggsae", name = "Bulk")
  bri_grtrack <- GeneRegionTrack(bri_temp, name = "Gene Model", genome = "briggse", showId = TRUE, fill = "#56B4E9", col = "black", transcriptAnnotation = "transcript")

  bri_hlight <- HighlightTrack(trackList = list(bri_dTrack_top, bri_dTrack_bottom, bri_dTrack_bulk),
                           start = min(c(bri_temp[bri_temp$gene == gene,]$start, bri_temp[bri_temp$gene == gene,]$end)),
                           end = max(c(bri_temp[bri_temp$gene == gene,]$start, bri_temp[bri_temp$gene == gene,]$end)),
                           chromosome = unique(bri_gene_query$seqnames), fill = "transparent", col = "black", lty = "dashed")
  
  bri_gene_track <- as.grob(~plotTracks(list(bri_gtrack, bri_hlight, bri_grtrack), type = "hist", background.title = "#56B4E9"))

  ## could try to return as a plot_grid
  
  png(paste0(NewBobPath, "test_briggsae_2.png"), height = 500, width = 2000)
  # plotTracks(list(gtrack, hlight, grtrack), type = "hist",
  #            background.title = "#56B4E9", cex = 1.2,
  #            cex.axis = 1.4, cex.main = 2, cex.title = 1.4, cex.id = 1.4, cex.group = 1)
  plot_grid(ele_gene_track, bri_gene_track,nrow=2)
  dev.off()
  
  return(list(ele_gene_track, bri_gene_track))
}

# For plotting by quantile
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

plot_gene_heatmap <- function(TPMList, diag_corr_genes, Category, Category_Name) {
  genes <- rownames(diag_corr_genes[which(diag_corr_genes[, Category_Name] == Category),])
  fold_change <- log2(TPMList[["elegans"]][genes, ] + 1) - log2(TPMList[["C.briggsae"]][genes, ] + 1)
  
  fold_change <- t(scale(t(fold_change)))
  
  return(pheatmap(fold_change,
           #scale = "row",
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           color = inferno(length(quantile_breaks(fold_change, n = 11)) - 1),
           breaks = quantile_breaks(fold_change, n = 11),
           drop_levels = TRUE,
           legend_breaks = quantile_breaks(fold_change, n = 11)
           )
         )
}
