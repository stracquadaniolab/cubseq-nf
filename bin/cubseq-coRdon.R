#!/usr/bin/env Rscript

"cubseq-coRdon.R

Usage:
    cubseq-coRdon.R <metadata> <inputdir> [--cub-statistic=<cub_statistic>] [--cub_expr_statistic=<cub_expr_statistic>] [--len-threshold=<len_threshold>] <codon_counts_outputfile> <cu_measure_outputfile> <cu_expressivity_outputfile>
    
Options:
    --cub-statistic=<cub_statistic>             Statistic to calculate codon usage bias (options: MILC, B, MCB, 
                                                ENCprime, ENC, SCUO) [default: MILC].
    --cub_expr_statistic=<cub_expr_statistic>   Statistic to calculate CU expressivity (options: MELP, E, CAI, 
                                                Fop and GCB) [default: MELP].
    --len-threshold=<len_threshold>             Remove sequences shorter than specified length [default: 80].
    -h --help                                   Show this screen.
    --version                                   Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-coRdon.R")

print(arguments) # for debugging purposes only

# load required packages
suppressMessages(library(coRdon))
suppressMessages(library(ggplot2))
suppressMessages(library(Biostrings))
suppressMessages(library(Biobase))

# reading sample sheet to get run_acc to construct file paths
metadata_file <- read.delim(arguments$metadata)

# prepare vector of file names and associated paths
metadata_file$files <- file.path(arguments$inputdir, metadata_file$run_accession, "*.mut-transcriptome.fa")
files <- metadata_file$files
names(files) <- metadata_file$run_accession

# read directory containing .fasta files
seqs <- readSet(file = files)

# convert sequence to codon table
seq_codon_table <- codonTable(seqs)

# get counts, ID and length and combine into dataframe for later output
seq_codon_counts <- codonCounts(seq_codon_table)
ID <- getID(seq_codon_table)
Length <- getlen(seq_codon_table)
ct_df <- data.frame(ID, Length, seq_codon_counts)

# calculate CU bias for each sequence
# TODO: initially compare to test subset (first 50 sequences) as prototype
# TODO: set ribosomal as default once functional annotation added; give user multiple docopt options
halfcT <- codonTable(codonCounts(seq_codon_table)[1:50,])

# calculate CU measure
cu_measure <- function(cub_statistic) {
  if (cub_statistic == "MILC") {
    CU <- MILC(seq_codon_table, self = TRUE, subsets = list(half = halfcT), filtering = "hard", len.threshold = 80)
  } else if (cub_statistic == "B") {
    CU <- B(seq_codon_table, self = TRUE, subsets = list(half = halfcT), filtering = "hard", len.threshold = 80)
  } else if (cub_statistic == "MCB") {
    CU <- MCB(seq_codon_table, self = TRUE, subsets = list(half = halfcT), filtering = "hard", len.threshold = 80)
  } else if (cub_statistic == "ENCprime") {
    CU <- ENCprime(seq_codon_table, self = TRUE, subsets = list(half = halfcT), filtering = "hard", len.threshold = 80)
  } else if (cub_statistic == "ENC") {
    CU <- ENC(seq_codon_table, filtering = "hard", len.threshold = 80)
  } else if (cub_statistic == "SCUO") {
    SCUO(seq_codon_table, filtering = "hard", len.threshold = 80)
  }
  else {
    stop("Input value should be one of the following: MILC, B, MCB, ENCprime, ENC or SCUO.")
  }
}

cu_measure_table <- cu_measure(cub_statistic = arguments$cub_statistic)

# calculate CU expressivity measure
# TODO: add "seed" and "subsets" as docopt params.
cu_expressivity <- function(cub_expr_statistic) {
  if (cub_expr_statistic == "MELP") {
    CU <- MELP(seq_codon_table, subsets = list(half = halfcT), filtering = "hard", len.threshold = 80) 
  } else if (cub_expr_statistic == "E") {
    CU <- E(seq_codon_table, self = TRUE, subsets = list(half = halfcT), filtering = "hard", len.threshold = 80)
  } else if (cub_expr_statistic == "CAI") {
    CU <- CAI(seq_codon_table, self = TRUE, subsets = list(half = halfcT), filtering = "hard", len.threshold = 80)
  } else if (cub_expr_statistic == "Fop") {
    CU <- Fop(seq_codon_table, self = TRUE, subsets = list(half = halfcT), filtering = "hard", len.threshold = 80)
  } else if (cub_expr_statistic == "GCB") {
    CU <- GCB(seq_codon_table, seed = halfcT, perc = 0.05, filtering = "hard", len.threshold = 80)
  }
  else {
    stop("Input value should be one of the following: MELP, E, CAI, Fop and GCB.")
  }
}

cu_expressivity_table <- cu_expressivity(cub_expr_statistic = arguments$cub_expr_statistic)

# TODO: combine all together into one dataframe.
# TODO: add option to calculate multiple/all CU measures and expressivity measures.

# write to file
write.table(ct_df, file = arguments$codon_counts_outputfile, quote = FALSE, sep = "\t", na = "", row.names = FALSE)
write.table(cu_measure_table, file = arguments$cu_measure_outputfile, quote = FALSE, sep = "\t", na = "", row.names = FALSE)
write.table(cu_expressivity_table, file = arguments$cu_expressivity_outputfile, quote = FALSE, sep = "\t", na = "", row.names = FALSE)
