
dir <- "/kimdata/livinlrg/scAnalysis/WS290/"
BobPath <- "/kimdata/livinlrg/scAnalysis/BobData/"
NewBobPath <- "/kimdata/livinlrg/scAnalysis/BobDataComb/"
ExternalPath <- "/kimdata/livinlrg/scAnalysis/ExternalDataset/"
source(paste0(BobPath, "../Scripts/CalcFunctions.R"))
BobPathTreeOut <- paste0(BobPath, "TreeView/")

# Load Orthogroups.txt
OGList_Raw <- scan("/kimdata/livinlrg/GenomeAnalysis/orthofinder_sans_uteleia/fastas/OrthoFinder/Results_Apr19/Orthogroups/Orthogroups.txt", what="", sep="\n")

OGNames <- unlist(lapply(OGList_Raw, function(OG) {
  return(gsub(":", "", strsplit(OG, " ")[[1]][1]))
}))

OGFullGeneList <- lapply(OGList_Raw, function(OG) {
  temp <- strsplit(OG, " ")[[1]]
  return(temp[2:length(temp)])
})

names(OGFullGeneList) <- OGNames

species_list <- data.frame(species = c("briggsae", "brenneri", "tropicalis", "zanzibari", "sulstoni", "elegans", "panamensis", "becei", "tribulationis", "nigoni", "inopinata", "latens", "remanei", "nigoni"),
                           species_short = c("CBG", "CBN", "Csp11", "CSP26", "CSP32", "CELEG", "CSP28", "CSP29", "CSP40", "nigoni", "Sp34", "FL83", "GCK72", "Cnig"))

OGFullSpeciesList <- lapply(OGFullGeneList, function(genes) {
  gene_concat <- paste0(genes)
  
  return(unlist(lapply(genes, function(gene) {
    species_list[unlist(lapply(species_list$species_short, function(species) {
      grepl(species, gene)
    })),]$species
  })))
  # (species_list[,]$species)
})

raw_gff_dir = "/kimdata/livinlrg/GenomeAnalysis/orthofinder_sans_uteleia/raw_genomes/"
gff_list <- list()
gff_list_mrna <- list()
# c_becei.PRJEB28243.WS290.annotations.gff3_longest_isoforms
file_list <- list.files(raw_gff_dir)
file_list <- file_list[file_list != "softmasked_genomic"]
for(species in file_list) {
  files <- list.files(paste0(raw_gff_dir, species, "/"))
  gff_list[[species]] <- rtracklayer::import(paste0(raw_gff_dir, species, "/", files[endsWith(files, ".gff3_longest_isoforms")]), format="gff3")
  gff_list_mrna[[species]] <- data.frame(gff_list[[species]][gff_list[[species]]$type == "mRNA",])
  gff_list_mrna[[species]]$gene_name <- apply(gff_list_mrna[[species]], 1, function(row) {strsplit(row["Parent"][[1]], ":")[[1]][2]})
}

## Now have elegans and briggsae in WBGene names
library(parallel)

gff_list_mrna[["briggsae"]]$short_name <- gsub("[.][1-9]$", "", gff_list_mrna[["briggsae"]]$Name)
gff_list_mrna[["briggsae"]]$short_name <- gsub("[a-z]$", "", gff_list_mrna[["briggsae"]]$short_name)

gff_list_mrna[["elegans"]]$short_name <- unlist(strsplit(gff_list_mrna[["elegans"]]$Name, "[.][1-9]$"))

gff_list_mrna[["elegans"]]$out_name <- ifelse(is.na(gff_list_mrna[["elegans"]]$locus), gff_list_mrna[["elegans"]]$short_name, gff_list_mrna[["elegans"]]$locus)
gff_list_mrna[["elegans"]]$out_name <- gsub("[a-z]$", "", gff_list_mrna[["elegans"]]$out_name)
gff_list_mrna[["briggsae"]]$out_name <- ifelse(is.na(gff_list_mrna[["briggsae"]]$locus), gff_list_mrna[["briggsae"]]$short_name, gff_list_mrna[["briggsae"]]$locus)

gff_list_mrna[["brenneri"]]$Name <- unlist(lapply(strsplit(gff_list_mrna[["brenneri"]]$Name, ".", fixed = TRUE), function(x) x[1]))
gff_list_mrna[["brenneri"]]$out_name <- ifelse(is.na(gff_list_mrna[["brenneri"]]$locus), gff_list_mrna[["brenneri"]]$Name, gff_list_mrna[["brenneri"]]$locus)

OGList <- mclapply(OGNames, function(OG) {
  temp_gene <- OGFullGeneList[[OG]]
  temp_species <- OGFullSpeciesList[[OG]]

  temp_ele <- temp_gene[grep("CELEG", temp_gene)]
  temp_ele <- gsub("CELEG.", "", temp_ele)
  
  #temp_ele <- gsub("[a-z]$", "", temp_ele)
  temp_ele_short <- gff_list_mrna[["elegans"]][which(gff_list_mrna[["elegans"]]$short_name %in% temp_ele),]$gene_name
  
  if(sum(temp_ele %in% gff_list_mrna[["elegans"]]$short_name) != length(temp_ele)) {
    print(temp_ele)
    print(OG)
  }
  
  temp_bri <- temp_gene[grep("CBG", temp_gene)]
  temp_bri <- gsub("[a-z]$", "", temp_bri)
  temp_bri_short <- gff_list_mrna[["briggsae"]][which(gff_list_mrna[["briggsae"]]$short_name %in% temp_bri),]$gene_name
  
  if(sum(temp_bri %in% gff_list_mrna[["briggsae"]]$short_name) != length(temp_bri)) {
    print(temp_bri)
  }
  
  return((c(temp_ele_short, temp_bri_short)))
}, mc.cores = 32)

names(OGList) <- OGNames
names(OGFullSpeciesList) <- OGNames

OGList_transcript_name <- mclapply(OGNames, function(OG) {
  temp_gene <- OGFullGeneList[[OG]]
  temp_species <- OGFullSpeciesList[[OG]]
  
  temp_ele <- temp_gene[grep("CELEG", temp_gene)]
  temp_ele <- gsub("CELEG.", "", temp_ele)
  temp_ele <- gsub("[a-z]$", "", temp_ele)
  
  temp_ele_short <- gff_list_mrna[["elegans"]][which(gsub("[a-z]$", "", gff_list_mrna[["elegans"]]$short_name) %in% temp_ele),]$short_name
  
  if(sum(temp_ele %in% gsub("[a-z]$", "", gff_list_mrna[["elegans"]]$short_name)) != length(temp_ele)) {
    print(temp_ele[! temp_ele %in% gsub("[a-z]$", "", gff_list_mrna[["elegans"]]$short_name)])
  }
  
  temp_bri <- temp_gene[grep("CBG", temp_gene)]
  temp_bri <- gsub("[a-z]$", "", temp_bri)
  temp_bri_short <- gff_list_mrna[["briggsae"]][which(gff_list_mrna[["briggsae"]]$short_name %in% temp_bri),]$short_name
  
  if(sum(temp_bri %in% gff_list_mrna[["briggsae"]]$short_name) != length(temp_bri)) {
    print(temp_bri[! temp_bri %in% gff_list_mrna[["briggsae"]]$short_name])
  }
  
  return((c(temp_ele_short, temp_bri_short)))
}, mc.cores = 32)

names(OGList_transcript_name) <- OGNames

OGList_short_name <- mclapply(OGNames, function(OG) {
  temp_gene <- OGFullGeneList[[OG]]
  temp_species <- OGFullSpeciesList[[OG]]
  
  temp_ele <- temp_gene[grep("CELEG", temp_gene)]
  temp_ele <- gsub("CELEG.", "", temp_ele)
  temp_ele <- gsub("[a-z]$", "", temp_ele)
  
  temp_ele_short <- gff_list_mrna[["elegans"]][which(gsub("[a-z]$", "", gff_list_mrna[["elegans"]]$short_name) %in% temp_ele),]$out_name
  
  if(sum(temp_ele %in% gsub("[a-z]$", "", gff_list_mrna[["elegans"]]$short_name)) != length(temp_ele)) {
    print(temp_ele[! temp_ele %in% gsub("[a-z]$", "", gff_list_mrna[["elegans"]]$short_name)])
  }
  
  temp_bri <- temp_gene[grep("CBG", temp_gene)]
  temp_bri <- gsub("[a-z]$", "", temp_bri)
  temp_bri_short <- gff_list_mrna[["briggsae"]][which(gff_list_mrna[["briggsae"]]$short_name %in% temp_bri),]$out_name
  
  if(sum(temp_bri %in% gff_list_mrna[["briggsae"]]$short_name) != length(temp_bri)) {
    print(temp_bri[! temp_bri %in% gff_list_mrna[["briggsae"]]$short_name])
  }
  
  return((c(temp_ele_short, temp_bri_short)))
}, mc.cores = 32)

names(OGList_short_name) <- OGNames

## fix some weird naming issues
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$out_name == "eif-2Bgamm",]$out_name <- "eif-2Bgamma"
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$out_name == "eif-2Bepsilo",]$out_name <- "eif-2Bepsilon"
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$out_name == "eif-2Bdelt",]$out_name <- "eif-2Bdelta"
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$out_name == "eif-2bet",]$out_name <- "eif-2beta"
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$out_name == "eif-2alph",]$out_name <- "eif-2alpha"
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$out_name == "eif-2gamm",]$out_name <- "eif-2gamma"
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$out_name == "eif-2Bbet",]$out_name <- "eif-2Bbeta"
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$out_name == "eif-2Balph",]$out_name <- "eif-2Balpha"

saveRDS(gff_list_mrna, paste0(NewBobPath, "Objects/WS290/gff_list_mrna.rds"))
saveRDS(OGList, paste0(NewBobPath, "Objects/WS290/OGList.rds"))
saveRDS(OGList_transcript_name, paste0(NewBobPath, "Objects/WS290/OGList_transcript_name.rds"))
saveRDS(OGList_short_name, paste0(NewBobPath, "Objects/WS290/OGList_short_name.rds"))
saveRDS(OGFullSpeciesList, paste0(NewBobPath, "Objects/WS290/OGFullSpeciesList.rds"))

gff_list_mrna <- readRDS(paste0(NewBobPath, "Objects/WS290/gff_list_mrna.rds"))
OGList <- readRDS(paste0(NewBobPath, "Objects/WS290/OGList.rds"))
OGFullSpeciesList <- readRDS(paste0(NewBobPath, "Objects/WS290/OGFullSpeciesList.rds"))
paml <- readRDS(paste0(NewBobPath, "Objects/WS290/paml.rds"))
briggsae_gtf <- data.frame(rtracklayer::import("/kimdata/livinlrg/scAnalysis/Reference/briggsae/C_briggsae.ws290.extendedUTRv3.format.gtf"))

# add mt genes to gff_list_mrna for briggsae
mt_gtf <- briggsae_gtf[briggsae_gtf$seqnames == "MtDNA",]

mt_gtf[,c("exon_number",
          "exon_numbe",
          "gene_source",
          "gene_biotype",
          "transcript_source",
          "transcript_biotype",
          "protein_id",
          "exon_id",
          "transcr")] <- NULL

colnames(mt_gtf) <- c("seqnames", "start", "end", "width", "strand", "source", "type", "score", "phase",
                            "gene_name", "short_name", "out_name")
mt_gtf[,c("ID", "Alias",
  "Name", "biotype",
  "curie", "sequence_name",
  "so_term_name", "Parent",
  "wormpep", "prediction_status",
  "protein_id", "Note",
  "locus", "start_not_found",
  "end_not_found", "matureType")] <- NA

gff_list_mrna[["briggsae"]] <- rbind(gff_list_mrna[["briggsae"]], mt_gtf)

gff_list_mrna[["elegans"]]$omega <- NA
gff_list_mrna[["elegans"]]$OG <- NA
gff_list_mrna[["elegans"]]$num_saturated_branches <- NA
gff_list_mrna[["elegans"]]$total_branches <- NA

gff_list_mrna[["briggsae"]]$omega <- NA
gff_list_mrna[["briggsae"]]$OG <- NA
gff_list_mrna[["briggsae"]]$num_saturated_branches <- NA
gff_list_mrna[["briggsae"]]$total_branches <- NA

rownames(gff_list_mrna[["elegans"]]) <- gff_list_mrna[["elegans"]]$gene_name
rownames(gff_list_mrna[["briggsae"]]) <- gff_list_mrna[["briggsae"]]$gene_name

for(OG in OGNames) {
  if(length(OGList[[OG]]) < 0) {
    print("nope")
  } else {
    cel_genes <- OGList[[OG]][OGList[[OG]] %in% gff_list_mrna[["elegans"]]$gene_name]
    cbr_genes <- OGList[[OG]][OGList[[OG]] %in% gff_list_mrna[["briggsae"]]$gene_name]
    if(length(cel_genes) > 0 & OG %in% paml$OG) {
      gff_list_mrna[["elegans"]][cel_genes,]$omega <- paml[OG,]$omega
      gff_list_mrna[["elegans"]][cel_genes,]$OG <- paml[OG,]$OG
      gff_list_mrna[["elegans"]][cel_genes,]$num_saturated_branches <- paml[OG,]$num_saturated_branches
      gff_list_mrna[["elegans"]][cel_genes,]$total_branches <- paml[OG,]$total_branches
    } else if (length(cel_genes) > 0) {
      gff_list_mrna[["elegans"]][cel_genes,]$OG <- OG
    }
    if(length(cbr_genes) > 0 & OG %in% paml$OG) {
      gff_list_mrna[["briggsae"]][cbr_genes,]$omega <- paml[OG,]$omega
      gff_list_mrna[["briggsae"]][cbr_genes,]$OG <- paml[OG,]$OG
      gff_list_mrna[["briggsae"]][cbr_genes,]$num_saturated_branches <- paml[OG,]$num_saturated_branches
      gff_list_mrna[["briggsae"]][cbr_genes,]$total_branches <- paml[OG,]$total_branches
    } else if (length(cbr_genes) > 0) {
      gff_list_mrna[["briggsae"]][cbr_genes,]$OG <- OG
    }
  }
}

dNdS <- read.csv(paste0(ExternalPath, "Ma_2023/media-5.csv"))
gff_list_mrna[["elegans"]]$Ma_omega <- dNdS[match(gff_list_mrna[["elegans"]]$gene_name, dNdS$WBGene),]$dN.dS.w

CutterKa <- read.csv(paste0(ExternalPath, "Tu-etal-2015-NAR_kaks_forChrisLarge_WS230_ele-bri_aligned_merge-YN_edit.csv"))
CutterKa <- CutterKa[! duplicated(CutterKa$ElegansSequence),]

CutterKa_elegans <- CutterKa[CutterKa$ElegansSequence %in% gff_list_mrna[["elegans"]]$gene_name,]
gff_list_mrna[["elegans"]]$Cutter_Ka <- CutterKa_elegans[match(gff_list_mrna[["elegans"]]$gene_name, CutterKa_elegans$ElegansSequence),]$Ka
gff_list_mrna[["elegans"]]$Cutter_Ks <- CutterKa_elegans[match(gff_list_mrna[["elegans"]]$gene_name, CutterKa_elegans$ElegansSequence),]$Ks
gff_list_mrna[["elegans"]]$Cutter_Ka.Ks <- CutterKa_elegans[match(gff_list_mrna[["elegans"]]$gene_name, CutterKa_elegans$ElegansSequence),]$Ka.Ks

CutterKa_briggsae <- CutterKa[CutterKa$BriggsaeSequence %in% gff_list_mrna[["briggsae"]]$gene_name,]
gff_list_mrna[["briggsae"]]$Cutter_Ka <- CutterKa_briggsae[match(gff_list_mrna[["briggsae"]]$gene_name, CutterKa_briggsae$BriggsaeSequence),]$Ka
gff_list_mrna[["briggsae"]]$Cutter_Ks <- CutterKa_briggsae[match(gff_list_mrna[["briggsae"]]$gene_name, CutterKa_briggsae$BriggsaeSequence),]$Ks
gff_list_mrna[["briggsae"]]$Cutter_Ka.Ks <- CutterKa_briggsae[match(gff_list_mrna[["briggsae"]]$gene_name, CutterKa_briggsae$BriggsaeSequence),]$Ka.Ks

## Need to sort through the OG's and pull out the number of elegans and briggsae genes in each
## Number of species from each OG
EleLength <- lapply(OGFullSpeciesList, function(OG) {
  return(sum(OG == "elegans", na.rm = TRUE))
})
names(EleLength) <- OGNames
EleLength <- unlist(EleLength)

BriLength <- lapply(OGFullSpeciesList, function(OG) {
  return(sum(OG == "briggsae", na.rm = TRUE))
})
names(BriLength) <- OGNames
BriLength <- unlist(BriLength)

gff_list_mrna[["elegans"]]$cel_OG_count <- NA
gff_list_mrna[["elegans"]][! is.na(gff_list_mrna[["elegans"]]$OG),]$cel_OG_count <- EleLength[gff_list_mrna[["elegans"]][! is.na(gff_list_mrna[["elegans"]]$OG),]$OG]
gff_list_mrna[["elegans"]]$cbr_OG_count <- NA
gff_list_mrna[["elegans"]][! is.na(gff_list_mrna[["elegans"]]$OG),]$cbr_OG_count <- BriLength[gff_list_mrna[["elegans"]][! is.na(gff_list_mrna[["elegans"]]$OG),]$OG]

gff_list_mrna[["briggsae"]]$cel_OG_count <- NA
gff_list_mrna[["briggsae"]][! is.na(gff_list_mrna[["briggsae"]]$OG),]$cel_OG_count <- EleLength[gff_list_mrna[["briggsae"]][! is.na(gff_list_mrna[["briggsae"]]$OG),]$OG]
gff_list_mrna[["briggsae"]]$cbr_OG_count <- NA
gff_list_mrna[["briggsae"]][! is.na(gff_list_mrna[["briggsae"]]$OG),]$cbr_OG_count <- BriLength[gff_list_mrna[["briggsae"]][! is.na(gff_list_mrna[["briggsae"]]$OG),]$OG]

# Worm Cat
WormCat <- read.delim(paste0(ExternalPath, "WormCat/", "whole_genome_v2_nov-11-2021.csv"), sep = ",")
WormCat <- WormCat[WormCat$Wormbase.ID %in% gff_list_mrna[["elegans"]]$gene_name,]

gff_list_mrna[["elegans"]]$WormCat.1 <- NA
gff_list_mrna[["elegans"]][match(WormCat$Wormbase.ID, gff_list_mrna[["elegans"]]$gene_name),]$WormCat.1 <- WormCat$Category.1

gff_list_mrna[["elegans"]]$WormCat.2 <- NA
gff_list_mrna[["elegans"]][match(WormCat$Wormbase.ID, gff_list_mrna[["elegans"]]$gene_name),]$WormCat.2 <- WormCat$Category.2

gff_list_mrna[["elegans"]]$WormCat.3 <- NA
gff_list_mrna[["elegans"]][match(WormCat$Wormbase.ID, gff_list_mrna[["elegans"]]$gene_name),]$WormCat.3 <- WormCat$Category.3

# SimpleMine data
SimpleMine_elegans <- read.delim(paste0(ExternalPath, "SimpleMine/", "simplemine_results_elegans.txt"), sep = "\t")
SimpleMine_elegans <- SimpleMine_elegans[SimpleMine_elegans$WormBase.Gene.ID %in% gff_list_mrna[["elegans"]]$gene_name,]

gff_list_mrna[["elegans"]]$operon <- NA
gff_list_mrna[["elegans"]][match(SimpleMine_elegans$WormBase.Gene.ID, gff_list_mrna[["elegans"]]$gene_name),]$operon <- SimpleMine_elegans$Operon

gff_list_mrna[["elegans"]]$RNAi_Phenotype <- NA
gff_list_mrna[["elegans"]][match(SimpleMine_elegans$WormBase.Gene.ID, gff_list_mrna[["elegans"]]$gene_name),]$RNAi_Phenotype <- SimpleMine_elegans$RNAi.Phenotype.Observed

gff_list_mrna[["elegans"]]$Allele_Phenotype <- NA
gff_list_mrna[["elegans"]][match(SimpleMine_elegans$WormBase.Gene.ID, gff_list_mrna[["elegans"]]$gene_name),]$Allele_Phenotype <- SimpleMine_elegans$Allele.Phenotype.Observed

SimpleMine_briggsae <- read.delim(paste0(ExternalPath, "SimpleMine/", "simplemine_results_briggsae.txt"), sep = "\t")
SimpleMine_briggsae <- SimpleMine_briggsae[SimpleMine_briggsae$WormBase.Gene.ID %in% gff_list_mrna[["briggsae"]]$gene_name,]

gff_list_mrna[["briggsae"]]$operon <- NA
gff_list_mrna[["briggsae"]][match(SimpleMine_briggsae$WormBase.Gene.ID, gff_list_mrna[["briggsae"]]$gene_name),]$operon <- SimpleMine_briggsae$Operon

embryonic_lethal_direct_inferred <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "embryonic_lethal_direct_inferred_rnai_variation.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]
embryonic_lethal_direct_inferred_rnai <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "embryonic_lethal_direct_inferred_rnai.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]
embryonic_lethal_direct_inferred_variation <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "embryonic_lethal_direct_inferred_variation.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]
embryonic_lethal_direct <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "embryonic_lethal_direct_rnai_variation.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]

larval_lethal_direct_inferred <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "larval_lethal_direct_inferred_rnai_variation.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]
larval_lethal_direct_inferred_rnai <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "larval_lethal_direct_inferred_rnai.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]
larval_lethal_direct_inferred_variation <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "larval_lethal_direct_inferred_variation.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]
larval_lethal_direct <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "larval_lethal_direct_rnai_variation.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]

lethal_direct_inferred <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "lethal_direct_inferred_rnai_variation.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]
lethal_direct_inferred_rnai <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "lethal_direct_inferred_rnai.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]
lethal_direct_inferred_variation <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "lethal_direct_inferred_variation.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]
lethal_direct <- read.delim(paste0(ExternalPath, "Lethal_Wormbase/", "lethal_direct_rnai_variation.txt"), sep = "\t", header = FALSE, col.names = c("Name", "gene_name_short", "X"))[,1:2]

gff_list_mrna[["elegans"]]$embryonic_lethal_direct_inferred <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% embryonic_lethal_direct_inferred$Name,]$embryonic_lethal_direct_inferred <- TRUE

gff_list_mrna[["elegans"]]$embryonic_lethal_direct_inferred_rnai <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% embryonic_lethal_direct_inferred_rnai$Name,]$embryonic_lethal_direct_inferred_rnai <- TRUE

gff_list_mrna[["elegans"]]$embryonic_lethal_direct_inferred_variation <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% embryonic_lethal_direct_inferred_variation$Name,]$embryonic_lethal_direct_inferred_variation <- TRUE

gff_list_mrna[["elegans"]]$embryonic_lethal_direct <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% embryonic_lethal_direct$Name,]$embryonic_lethal_direct <- TRUE

gff_list_mrna[["elegans"]]$larval_lethal_direct_inferred <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% larval_lethal_direct_inferred$Name,]$larval_lethal_direct_inferred <- TRUE

gff_list_mrna[["elegans"]]$larval_lethal_direct_inferred_rnai <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% larval_lethal_direct_inferred_rnai$Name,]$larval_lethal_direct_inferred_rnai <- TRUE

gff_list_mrna[["elegans"]]$larval_lethal_direct_inferred_variation <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% larval_lethal_direct_inferred_variation$Name,]$larval_lethal_direct_inferred_variation <- TRUE

gff_list_mrna[["elegans"]]$larval_lethal_direct <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% larval_lethal_direct$Name,]$larval_lethal_direct <- TRUE

gff_list_mrna[["elegans"]]$lethal_direct_inferred <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% lethal_direct_inferred$Name,]$lethal_direct_inferred <- TRUE

gff_list_mrna[["elegans"]]$lethal_direct_inferred_rnai <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% lethal_direct_inferred_rnai$Name,]$lethal_direct_inferred_rnai <- TRUE

gff_list_mrna[["elegans"]]$lethal_direct_inferred_variation <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% lethal_direct_inferred_variation$Name,]$lethal_direct_inferred_variation <- TRUE

gff_list_mrna[["elegans"]]$lethal_direct <- FALSE
gff_list_mrna[["elegans"]][gff_list_mrna[["elegans"]]$gene_name %in% lethal_direct$Name,]$lethal_direct <- TRUE

# PhyloStrata
MA_PS <- read.csv(paste0(ExternalPath, "Ma_2023/", "pnas.2216351120.sd01.csv"))

colnames(MA_PS)[4] <- "PS.Value"
MA_PS <- MA_PS[MA_PS$WBGene.ID %in% gff_list_mrna[["elegans"]]$gene_name,]

gff_list_mrna[["elegans"]]$PS.Value <- MA_PS[match(gff_list_mrna[["elegans"]]$gene_name, MA_PS$WBGene.ID),]$PS.Value

PS.species <- c("Mus musculus", "Drosophila melanogaster", "Trichinella spiralis", "Brugia malayi", "Pristionchus pacificus", "Oscheius tipulae", "Caenorhabditis bovis", "Caenorhabditis uteleia", "Caenorhabditis becei", "Caenorhabditis briggsae", "Caenorhabditis inopinata", "Caenorhabditis elegans")
names(PS.species) <- seq(1, 12)

gff_list_mrna[["elegans"]]$PS.Name <- PS.species[gff_list_mrna[["elegans"]]$PS.Value]

# Maternal
Goldstein <- read.table(paste0(ExternalPath, "FromFelicia/", "TableS2_RPKMs_Goldstein_singleCell2016.txt"),
                        header = TRUE)

Goldstein <- Goldstein[,c("Label", "P0_1.cell_st411", "P0_1.cell_st413", "P0_1.cell_st441", "P0_1.cell_st449", "P0_1.cell_st451")]
Goldstein$Median <- apply(Goldstein[,2:6], 1, median)
Goldstein$Maternal <- Goldstein$Median > 0

rownames(Goldstein) <- Goldstein$Label

gff_list_mrna[["elegans"]]$maternal <- ifelse(! is.na(Goldstein[match(gff_list_mrna[["elegans"]]$out_name, Goldstein$Label),]$Maternal),
                                              Goldstein[match(gff_list_mrna[["elegans"]]$out_name, Goldstein$Label),]$Maternal,
                                              ifelse(! is.na(Goldstein[match(gff_list_mrna[["elegans"]]$gene_name, Goldstein$Label),]$Maternal),
                                              Goldstein[match(gff_list_mrna[["elegans"]]$gene_name, Goldstein$Label),]$Maternal,
                                              Goldstein[match(gsub("[a-z]$", "", gff_list_mrna[["elegans"]]$short_name), Goldstein$Label),]$Maternal))

# Recombination arms
# Use the recombination bins from Rockman and Kruglyak 2009
CelRecombArms <- data.frame(chrom = rep(c("I", "II", "III", "IV", "V", "X"), each = 5),
                         local = rep(c("LeftTip", "LeftArm", "Center", "RightArm", "RightTip"), times = 6),
                         start = c(0, 527, 3858, 11040, 14875,
                                   0, 306, 4879, 12020, 14609,
                                   0, 494, 3722, 10340, 13217,
                                   0, 720, 3896, 12970, 16712,
                                   0, 643, 5897, 16550, 20337,
                                   0, 572, 6137, 12480, 16417),
                         end = c(527, 3858, 11040, 14875, 15072,
                                 306, 4879, 12020, 14609, 15279,
                                 494, 3722, 10340, 13217, 13784,
                                 720, 3896, 12970, 16712, 17949,
                                 643, 5897, 16550, 20337, 20920,
                                 572, 6137, 12480, 16417, 17719))

# Use the recombinant arms from C. briggsae Ross et al., 2011
CbrRecombArms <- data.frame(chrom = rep(c("I", "II", "III", "IV", "V", "X"), each = 5),
                            local = rep(c("LeftTip", "LeftArm", "Center", "RightArm", "RightTip"), times = 6),
                            start = c(0, 370, 4210, 10890, 15200,
                                      0, 90, 4530, 11640, 16540,
                                      0, 360, 4700, 10810, 14100,
                                      0, 460, 5730, 13840, 17150,
                                      0, 120, 6170, 14100, 19210,
                                      0, 790, 8490, 15440, 21370),
                            end = c(370, 4210, 10890, 15200, 15440,
                                    90, 4530, 11640, 16540, 16630,
                                    360, 4700, 10810, 14100, 14590,
                                    460, 5730, 13840, 17150, 17490,
                                    120, 6170, 14100, 19210, 19490,
                                    790, 8490, 15440, 21370, 21550))

gff_list_mrna[["elegans"]]$RecombRegion <- NA

#Define regions of chromosome in the gff_list based on CelRecombArms
for(i in seq(1, nrow(gff_list_mrna[["elegans"]]))) {
  print(paste0(i, " ", gff_list_mrna[["elegans"]][i,]$seqnames))
  if(is.na(gff_list_mrna[["elegans"]][i,]$seqnames) | gff_list_mrna[["elegans"]][i,]$seqnames == "MtDNA") next
  
  temp <- CelRecombArms[CelRecombArms$chrom %in% gff_list_mrna[["elegans"]][i,]$seqnames,]
  for(j in seq(1, nrow(temp))) {
    if(gff_list_mrna[["elegans"]][i,]$start/1000 > temp[j,]$start &
       gff_list_mrna[["elegans"]][i,]$end/1000 < temp[j,]$end) {
      gff_list_mrna[["elegans"]][i,]$RecombRegion <- paste0(temp[j,]$chrom, "_", temp[j,]$local)
    }
  }
}

gff_list_mrna[["briggsae"]]$RecombRegion <- NA

#Define regions of chromosome in gff_list based on CbrRecombArms
for(i in seq(1, nrow(gff_list_mrna[["briggsae"]]))) {
  print(paste0(i, " ", gff_list_mrna[["elegans"]][i,]$seqnames))
  if(is.na(gff_list_mrna[["briggsae"]][i,]$seqnames) | (! gff_list_mrna[["briggsae"]][i,]$seqnames %in% CbrRecombArms$chrom)) next
  
  temp <- CbrRecombArms[CbrRecombArms$chrom %in% gff_list_mrna[["briggsae"]][i,]$seqnames,]
  for(j in seq(1, nrow(temp))) {
    if(gff_list_mrna[["briggsae"]][i,]$start/1000 > temp[j,]$start &
       gff_list_mrna[["briggsae"]][i,]$end/1000 < temp[j,]$end) {
      gff_list_mrna[["briggsae"]][i,]$RecombRegion <- paste0(temp[j,]$chrom, "_", temp[j,]$local)
    }
  }
}

gff_list_mrna[["elegans"]][synteny_filt[synteny_filt$canonical == 1,]$cel_WBG_name,]$RecombRegion ==
gff_list_mrna[["briggsae"]][synteny_filt[synteny_filt$canonical == 1,]$cbr_WBG_name,]$RecombRegion

saveRDS(gff_list_mrna, paste0(NewBobPath, "Objects/WS290/gff_list_mrna.rds"))

saveRDS(OGList, paste0(NewBobPath, "Objects/WS290/OGList.rds"))
saveRDS(OGNames, paste0(NewBobPath, "Objects/WS290/OGNames.rds"))
saveRDS(OGFullGeneList, paste0(NewBobPath, "Objects/WS290/OGFullGeneList.rds"))
saveRDS(OGFullSpeciesList, paste0(NewBobPath, "Objects/WS290/OGFullSpeciesList.rds"))
saveRDS(paml, paste0(NewBobPath, "Objects/paml.rds"))
paml <- readRDS(paste0(NewBobPath, "Objects/paml.rds"))

gff_list_mrna <- readRDS(paste0(NewBobPath, "Objects/WS290/gff_list_mrna.rds"))
synteny_filt <- readRDS(paste0(dir, "Objects/synteny_filt.rds"))

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_bg_20240829.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

rm(cds)

gff_list_mrna[["elegans"]] <- left_join(gff_list_mrna[["elegans"]], data.frame(rowData(cds_filt)[,c("elegans_id", "elegans_gene_short_name", "elegans_gene_long_name",
                                                        "briggsae_id", "briggsae_gene_short_name", "briggsae_gene_long_name",
                                                        "gene.type", "orthology_conf")]), by = c("out_name" = "elegans_gene_short_name"))

gff_list_mrna[["briggsae"]] <- left_join(gff_list_mrna[["briggsae"]], data.frame(rowData(cds_filt)[,c("elegans_id", "elegans_gene_short_name", "elegans_gene_long_name",
                                                                                               "briggsae_id", "briggsae_gene_short_name", "briggsae_gene_long_name",
                                                                                               "gene.type", "orthology_conf")]), by = c("out_name" = "briggsae_gene_short_name"))

gff_list_mrna[["briggsae"]] <- left_join(gff_list_mrna[["briggsae"]],
          gff_list_mrna[["elegans"]][,c("elegans_id", "WormCat.1", "WormCat.2", "WormCat.3")],
          by = "elegans_id",
          na_matches = "never")


## relate to the gff files to add orthology relationships
dir <- "/kimdata/livinlrg/scAnalysis/WS290/"
synteny <- readRDS(paste0(dir, "Objects/synteny_filt.rds"))

cds <- readRDS(paste0(dir, "Objects/cds_objects/cds_no_bg_20240828.rds"))
cds_filt <- cds[, which(! colData(cds)$filter %in% c("damaged", "doublet", "mutant", "stressed"))]

rm(cds)

rownames(cds_filt) <- rowData(cds_filt)$id
gff_list_mrna[["elegans"]]$cds_gene_name <- rowData(cds_filt)[gff_list_mrna[["elegans"]]$gene_name,]$gene_short_name
gff_list_mrna[["briggsae"]]$cds_gene_name <- rowData(cds_filt)[ifelse(is.na(gff_list_mrna[["briggsae"]]$elegans_id), gff_list_mrna[["briggsae"]]$briggsae_id, gff_list_mrna[["briggsae"]]$elegans_id),]$gene_short_name

# elegans
gff_list_mrna[["elegans"]][,c("cel_protein_length", "cel_percent_align", "cbr_percent_align",
                              "cbr_protein_length", "percent_identity", "percent_similarity",
                              "syntenic", "mutual_best_blat", "mutual_best_sw",
                              "cel_cbr_sw_best", "cel_cbr_sw_best_score", "cel_cbr_sw_next_best", "cel_cbr_sw_next_best_score",
                              "cbr_cel_sw_best", "cbr_cel_sw_best_score", "cbr_cel_sw_next_best", "cbr_cel_sw_next_best_score")] <- NA
for(gene in gff_list_mrna[["elegans"]][! is.na(gff_list_mrna[["elegans"]]$cds_gene_name),]$gene_name) {
  temp_row <- gff_list_mrna[["elegans"]][which(gff_list_mrna[["elegans"]]$gene_name == gene),]
  
  if(gene %in% synteny$cel_WBG_name & temp_row$briggsae_id %in% synteny$cbr_WBG_name) {
    temp_synteny <- synteny[which(synteny$cel_WBG_name == gene & synteny$cbr_WBG_name == temp_row$briggsae_id),
                            c("cel_protein_length", "cel_percent_align", "cbr_percent_align",
                              "cbr_protein_length", "percent_identity", "percent_similarity",
                              "syntenic", "mutual_best_blat", "mutual_best_sw",
                              "cel_cbr_sw_best", "cel_cbr_sw_best_score", "cel_cbr_sw_next_best", "cel_cbr_sw_next_best_score",
                              "cbr_cel_sw_best", "cbr_cel_sw_best_score", "cbr_cel_sw_next_best", "cbr_cel_sw_next_best_score")]
    
    gff_list_mrna[["elegans"]][which(gff_list_mrna[["elegans"]]$gene_name == gene),
                               c("cel_protein_length", "cel_percent_align", "cbr_percent_align",
                                 "cbr_protein_length", "percent_identity", "percent_similarity",
                                 "syntenic", "mutual_best_blat", "mutual_best_sw",
                                 "cel_cbr_sw_best", "cel_cbr_sw_best_score", "cel_cbr_sw_next_best", "cel_cbr_sw_next_best_score",
                                 "cbr_cel_sw_best", "cbr_cel_sw_best_score", "cbr_cel_sw_next_best", "cbr_cel_sw_next_best_score")] <- temp_synteny
  }
}

# briggsae
gff_list_mrna[["briggsae"]][,c("cel_protein_length", "cel_percent_align", "cbr_percent_align",
                              "cbr_protein_length", "percent_identity", "percent_similarity",
                              "syntenic", "mutual_best_blat", "mutual_best_sw",
                              "cel_cbr_sw_best", "cel_cbr_sw_best_score", "cel_cbr_sw_next_best", "cel_cbr_sw_next_best_score",
                              "cbr_cel_sw_best", "cbr_cel_sw_best_score", "cbr_cel_sw_next_best", "cbr_cel_sw_next_best_score")] <- NA
for(gene in gff_list_mrna[["briggsae"]][! is.na(gff_list_mrna[["briggsae"]]$cds_gene_name),]$gene_name) {
  temp_row <- gff_list_mrna[["briggsae"]][which(gff_list_mrna[["briggsae"]]$gene_name == gene),]
  
  if(gene %in% synteny$cbr_WBG_name & temp_row$elegans_id %in% synteny$cel_WBG_name) {
    temp_synteny <- synteny[which(synteny$cbr_WBG_name == gene & synteny$cel_WBG_name == temp_row$elegans_id),
                            c("cel_protein_length", "cel_percent_align", "cbr_percent_align",
                              "cbr_protein_length", "percent_identity", "percent_similarity",
                              "syntenic", "mutual_best_blat", "mutual_best_sw",
                              "cel_cbr_sw_best", "cel_cbr_sw_best_score", "cel_cbr_sw_next_best", "cel_cbr_sw_next_best_score",
                              "cbr_cel_sw_best", "cbr_cel_sw_best_score", "cbr_cel_sw_next_best", "cbr_cel_sw_next_best_score")]
    
    gff_list_mrna[["briggsae"]][which(gff_list_mrna[["briggsae"]]$gene_name == gene),
                               c("cel_protein_length", "cel_percent_align", "cbr_percent_align",
                                 "cbr_protein_length", "percent_identity", "percent_similarity",
                                 "syntenic", "mutual_best_blat", "mutual_best_sw",
                                 "cel_cbr_sw_best", "cel_cbr_sw_best_score", "cel_cbr_sw_next_best", "cel_cbr_sw_next_best_score",
                                 "cbr_cel_sw_best", "cbr_cel_sw_best_score", "cbr_cel_sw_next_best", "cbr_cel_sw_next_best_score")] <- temp_synteny
  }
}

saveRDS(gff_list_mrna, paste0(dir, "Objects/gff_list_mrna.rds"))