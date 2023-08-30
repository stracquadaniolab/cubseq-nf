#!/usr/bin/env Rscript

"cubseq-get-featurecounts-tpm.R

Usage:
    cubseq-get-featurecounts-tpm.R <input> [--gtf=<gtf>] [--output=<output>]
    
Options:
    --gtf=<gtf>                             GTF file used to map to gene IDs.
    --output=<output>                       CSV featureCounts output file with additional RPK and TPM counts.
    -h --help                               Show this screen.
    --version                               Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-get-featurecounts-tpm.R")

# load required packages
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(rtracklayer)) # to import GTF

# load data
counts <- read_delim(arguments$input, delim = "\t", skip = 1)
gtf <- rtracklayer::import(arguments$gtf)

# calculate TPMs
counts <- counts %>%
  dplyr::rename(NumReads = 7) %>%
  dplyr::select(Geneid, Length, NumReads) %>%
  # calculate RPK
  mutate(RPK = NumReads/(Length/1000))
# calculate PM scaling factor
PM <- sum(counts$RPK)/1e6
# calculate TPM
counts <- counts %>%
  mutate(TPM = RPK / PM)

# create map of gene IDs to transcript IDs
trID.map <- as.data.frame(gtf) %>%
  dplyr::filter(., type == "exon") %>%
  dplyr::select(., gene_id, transcript_id)

# map transcript IDs to gene IDs
# this is required by tximport so it can correctly map to the 
# first column of its tx2gene object containing transcript IDs
# counts <- merge(counts, trID.map, by.x = "Geneid", by.y = "gene_id", all = TRUE)
counts <- inner_join(counts, trID.map, by = c("Geneid" = "gene_id")) %>%
  distinct(Geneid, .keep_all = TRUE)

# write to file
write_tsv(counts, arguments$output) # need to save as TSV file for tximport
