#!/usr/bin/env Rscript

"cubseq-summarize-to-gene.R

Usage:
    cubseq-summarize-to-gene.R [<inputfile>] [--gtf=<gtf>] [--quant-dir=<quant_dir>] [--counts-from-abundance=<counts_model>] [--output=<output>] 
    
Options:
    --gtf=<gtf>                             GTF file used to count transcripts.
    --quant-dir=<quant_dir>                 Directory containing salmon results [default: .].
    --counts-from-abundance=<counts_model>  Generate estimated counts using abundance estimates either:
                                            no, scaledTPM, lengthScaledTPM, dtuScaledTPM [default: no].
    --output=<output>                       RDS output file [default: txi-summarized-experiment.rds].
    -h --help                               Show this screen.
    --version                               Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-summarize-to-gene.R")

# load required packages
suppressMessages(library(readr)) # tximport uses readr package if available
suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures)) # for converting gtf file into tx2gene object

# read metadata file
metadata_file <- read.delim(arguments$inputfile)

# prepare vector of file names and associated paths (from salmon-quant) for tximport
metadata_file$files <- file.path(arguments$quant_dir, metadata_file$run_accession, "quant.sf")
files <- metadata_file$files
names(files) <- metadata_file$run_accession

# files <- arguments$quant_dir

# for debugging purposes only
print(arguments)
print(all(file.exists(files)))
print(files)

# obtain annotation
txdb <- makeTxDbFromGFF(arguments$gtf)

# Associate transcripts with gene IDs for gene-level summarization
# NB: this must have specific column order: 1) transcript ID and 2) gene ID
k <- keys(txdb, keytype = "TXNAME") # get transcript names
tx2gene <- select(txdb, keys = k, columns="GENEID", keytype = "TXNAME")

# quantify transcripts and summarize counts to gene level
txi <- tximport(files, type = "salmon", countsFromAbundance = arguments$counts_model, tx2gene = tx2gene)
print(dim(txi$abundance)) # for debugging/sanity-check purposes

# save matrix to rds file
saveRDS(txi, file = arguments$output)
