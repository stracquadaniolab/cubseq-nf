#!/usr/bin/env Rscript

"cubseq-summarize-featurecounts-to-gene.R

Usage:
    cubseq-summarize-featurecounts-to-gene.R <input> [--gtf=<gtf>] [--featurecounts-tpm-dir=<featurecounts_tpm_dir>] [--counts-from-abundance=<counts_model>] [--output=<output>] 
    
Options:
    --gtf=<gtf>                                         GTF file used to count transcripts.
    --featurecounts-tpm-dir=<featurecounts_tpm_dir>     Directory containing featureCounts files [default: .].
    --counts-from-abundance=<counts_model>              Generate estimated counts using abundance estimates either:
                                                        no, scaledTPM, lengthScaledTPM, dtuScaledTPM [default: no].
    --output=<output>                                   RDS output file [default: 
                                                        txi-featurecounts-summarized-experiment.rds].
    -h --help                                           Show this screen.
    --version                                           Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-summarize-featurecounts-to-gene.R")

# load required packages
suppressMessages(library(readr)) # tximport uses readr package if available
suppressMessages(library(dplyr))
suppressMessages(library(tximport))
suppressMessages(library(GenomicFeatures)) # for converting gtf file into tx2gene object

# read metadata file
metadata_file <- read.delim(arguments$input)

# prepare vector of file names and associated paths from featureCounts for tximport
metadata_file$files <- file.path(arguments$featurecounts_tpm_dir, 
                                 paste(metadata_file$run_accession, 
                                       "featureCounts-tpm.txt", 
                                       sep = "."))
files <- metadata_file$files
names(files) <- metadata_file$run_accession

# for debugging purposes only
print(arguments)
print(all(file.exists(files)))
print(files)

# obtain annotation
txdb <- makeTxDbFromGFF(arguments$gtf, format = "gtf")

# Associate transcripts with gene IDs for gene-level summarization
# NB: this must have specific column order: 1) transcript ID and 2) gene ID
k <- keys(txdb, keytype = "TXNAME") # get transcript names
tx2gene <- select(txdb, keys = k, columns="GENEID", keytype = "TXNAME")

# quantify transcripts and summarize counts to gene level
txi <- tximport(files,
                type = "none",
                countsFromAbundance = arguments$counts_model,
                tx2gene = tx2gene,
                geneIdCol = "Geneid", # N.B. refers to column in featureCounts file
                txIdCol = "transcript_id", # N.B. refers to column in featureCounts file
                abundanceCol = "TPM",
                countsCol = "NumReads",
                lengthCol = "Length",
                importer=function(x) read_tsv(x)) # specify the importer function on how to read lines

print(dim(txi$abundance)) # for debugging/sanity-check purposes

# save matrix to rds file
saveRDS(txi, file = arguments$output)
