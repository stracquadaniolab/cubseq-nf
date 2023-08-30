#!/usr/bin/env Rscript

"cubseq-get-strain-tpm.R

Usage:
    cubseq-get-strain-tpm.R <input> [--metadata=<metadata>] [--num_strains=<num_strains>]
    
Options:
    --metadata=<metadata>             CSV file containing sample metadata.
    --num_strains=<num_strains>       Number of strains.
    -h --help                         Show this screen.
    --version                         Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-get-strain-tpm.R")

# load data
txi <- readRDS(arguments$input)
gtf <- read.table(arguments$gtf, header = FALSE, sep = '\t')
metadata <- read_delim(arguments$metadata, delim = "\t", col_names = TRUE)
num_strains <- as.numeric(arguments$num_strains)

# Filter metadata to keep only top strains (does not include NAs)
filtered_metadata <- metadata %>%
  filter(strain %in% names(sort(table(strain), decreasing = TRUE)[1:num_strains]))

# filter TPMs by top 5 strains
tpm <- txi$abundance
filtered_txi <- tpm[, (colnames(tpm) %in% filtered_metadata$run_accession)]

# replace colnames with strain names
colnames(filtered_txi) <- filtered_metadata$strain[match(colnames(filtered_txi), filtered_metadata$run_accession)]

# split TPMs by strain
filtered_txi <- as.data.frame(filtered_txi) # needs to be in df format
strain_txi <- split.default(filtered_txi, colnames(filtered_txi))

