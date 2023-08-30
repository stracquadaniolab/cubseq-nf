#!/usr/bin/env Rscript

"cubseq-get-heg-leg-fasta.R

Usage:
    cubseq-get-heg-leg-fasta.R <input> [--gtf=<gtf>] [--mut-transcriptome-dir=<mut_transcriptome_dir>] [--heg-geneID=<heg_geneID>] [--leg-geneID=<leg_geneID>]
    
Options:
    --gtf=<gtf>                                         GTF file.
    --mut-transcriptome-dir=<mut_transcriptome_dir>     Path to mut-transcriptome directory.
    --heg_geneID=<heg_geneID>                           CSV file containing gene IDs for high expression genes.
    --leg_geneID=<leg_geneID>                           CSV file containing gene IDs for low expression genes.
    -h --help                                           Show this screen.
    --version                                           Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-get-heg-leg-fasta.R")

# load libraries
# suppressMessages(library(seqinr))
suppressMessages(library(Biostrings))
suppressMessages(library(tidyverse))

# load files
metadata <- read.delim(arguments$input)
gtf <- read.table(arguments$gtf, header = FALSE, sep = '\t')
# heg.geneID <- read.table(arguments$heg_geneID, header = TRUE)
# leg.geneID <- read.table(arguments$leg_geneID, header = TRUE)
heg.geneID <- read_csv(arguments$heg_geneID)
leg.geneID <- read_csv(arguments$leg_geneID)

# prepare vector of file names and associated paths (from salmon-quant) for tximport
metadata$files <- file.path(arguments$mut_transcriptome_dir, metadata$run_accession, paste0(metadata$run_accession, ".mut-transcriptome.fa"))

files <- metadata$files
names(files) <- metadata$run_accession

# for debugging purposes only
print(all(file.exists(files)))
print(files)

# # filter each fasta file by HEG and LEG transcript IDs - seqinr method
# for (file in names(files)) {
  
#   fasta <- read.fasta(files[[file]])
#   fasta.heg <- fasta[names(fasta) %in% heg.geneID[[1]]]
#   fasta.leg <- fasta[names(fasta) %in% leg.geneID[[1]]]
  
#   heg.name <- paste0(file, ".heg-mut-transcriptome.fa")
#   leg.name <- paste0(file, ".leg-mut-transcriptome.fa")
  
#   write.fasta(sequences = fasta.heg, names = names(fasta.heg), file.out = heg.name, open = "w")
#   write.fasta(sequences = fasta.leg, names = names(fasta.leg), file.out = leg.name, open = "w")
  
# }

getTrID <- function(gtf) {
  
  # Function creates map linking gene IDs with transcript IDs from gtf file.
  
  trID.map <- gtf %>%
    # filter to keep only transcripts
    filter(., grepl("exon", V3)) %>%
    # extract gene ID and transcript ID from defline
    extract(
      ., V9, into = c("gene_id", "transcript_id"),
      regex = "gene_id\\s*(.*?); transcript_id\\s*(.*?);"
    ) %>%
    dplyr::select(c(gene_id, transcript_id))
  
  return(trID.map)
  
}

# create map of gene IDs to transcript IDs
trID.map <- getTrID(gtf)

# map HEG and LEG to their transcript IDs
heg.geneID <- subset(trID.map, gene_id %in% heg.geneID$genes)
leg.geneID <- subset(trID.map, gene_id %in% leg.geneID$genes)

# filter each fasta file by HEG and LEG transcript IDs - Biostrings method
for (file in names(files)) {
  
  # read file and convert to dataframe for dplyr filtering
  fasta <- readDNAStringSet(files[[file]])
  fasta_df <- data.frame(Defline = names(fasta), Read = paste(fasta))

  # filter for HEG and LEG transcript IDs
  fasta.heg <- fasta_df %>%
    filter(str_detect(Defline, paste(heg.geneID[[2]], collapse = "|")))
  fasta.leg <- fasta_df %>%
    filter(str_detect(Defline, paste(leg.geneID[[2]], collapse = "|")))

  # convert dataframe back to DNAStringSet object
  DNAstr.heg <- DNAStringSet(fasta.heg[["Read"]])
  DNAstr.leg <- DNAStringSet(fasta.leg[["Read"]])
  names(DNAstr.heg) <- fasta.heg[["Defline"]]
  names(DNAstr.leg) <- fasta.leg[["Defline"]]

  # prepare file names (i.e. we name using run accession)
  heg.name <- paste0(file, ".heg-mut-transcriptome.fa")
  leg.name <- paste0(file, ".leg-mut-transcriptome.fa")

  # write to fasta file
  writeXStringSet(DNAstr.heg, heg.name)
  writeXStringSet(DNAstr.leg, leg.name)
  
}
