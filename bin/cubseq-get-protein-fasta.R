#!/usr/bin/env Rscript

"cubseq-get-protein-fasta.R

Usage:
    cubseq-get-protein-fasta.R <input> [--gtf=<gtf>] [--mut-transcriptome-dir=<mut_transcriptome_dir>]
    
Options:
    --gtf=<gtf>                                         GTF file.
    --mut-transcriptome-dir=<mut_transcriptome_dir>     Path to mut-transcriptome directory [default: .].
    -h --help                                           Show this screen.
    --version                                           Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-get-protein-fasta.R")

# load libraries
# suppressMessages(library(seqinr))
suppressMessages(library(Biostrings))
suppressMessages(library(tidyverse))
suppressMessages(library(rtracklayer))

# load files
metadata <- read.delim(arguments$input)
gtf <- import(arguments$gtf)

# filter to keep protein-coding genes from GTF
global_gtf <- as.data.frame(gtf) %>%
  filter(type == "CDS")
print(dim(global_gtf))

# prepare vector of file names and associated paths (from salmon-quant) for tximport
metadata$files <- file.path(arguments$mut_transcriptome_dir, metadata$run_accession, paste0(metadata$run_accession, ".mut-transcriptome.fa"))

files <- metadata$files
names(files) <- metadata$run_accession

# for debugging purposes only
print(all(file.exists(files)))
print(files)

# filter each fasta file by protein transcript IDs - Biostrings method
for (file in names(files)) {
  
  # read file and convert to dataframe for dplyr filtering
  fasta <- readDNAStringSet(files[[file]])
  fasta_df <- data.frame(Defline = names(fasta), Read = paste(fasta))

  # filter fasta file for protein transcript IDs
  fasta.protein <- fasta_df %>%
    filter(str_detect(Defline, paste(global_gtf$transcript_id, collapse = "|")))

  # convert dataframe back to DNAStringSet object
  DNAstr.protein <- DNAStringSet(fasta.protein[["Read"]])
  names(DNAstr.protein) <- fasta.protein[["Defline"]]

  # prepare file names (i.e. we name using run accession)
  protein.name <- paste0(file, ".protein-mut-transcriptome.fa")

  # write to fasta file
  writeXStringSet(DNAstr.protein, protein.name)
  
}
