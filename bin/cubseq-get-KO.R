#!/usr/bin/env Rscript

"cubseq-get-KO.R

Usage:
    cubseq-get-KO.R <input> [--org-id=<org_id>] <output> 
    
Options:
    --org-id=<org_id>               KEGG organism code [default: eco].
    -h --help                       Show this screen.
    --version                       Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-get-KO.R")

# load packages
suppressMessages(library(tidyverse))
suppressMessages(library(KEGGREST))
suppressMessages(library(Biostrings))

# get list of all genes for organism
genes <- as.data.frame(keggList(arguments$org_id))
colnames(genes)[1] <- "Gene"

# extract gene names
genes <- genes %>%
  separate("Gene", into = c("Gene", "Definition"), sep = ";", extra = "merge")

# get KO IDs for each gene
ko <- as.data.frame(keggLink(target = "ko", source = "eco"))

# map KO IDs (using row names)
gene2ko <- merge(genes, ko, by=0, all.x = TRUE)
gene2ko$ko.id <- gsub("ko:", "", gene2ko[[4]])
gene2ko <- gene2ko[,-c(4)]

# load fasta file
fasta = readDNAStringSet(arguments$input, format = "fasta")

# extract fasta deflines and sequences
fasta.defline <- names(fasta)
fasta.read <- paste(fasta)
fasta.gene <- str_match(fasta.defline, "gene_name=\\s*(.*?)\\s*;")[,2]
fasta.id <- gsub(" .*", "", fasta.defline)

fastadf <- data.frame(ID=fasta.id, Gene=fasta.gene, Sequence=fasta.read)

# link each gene name to KO ID
fasta2ko <- merge(fastadf, gene2ko, by.x = "Gene", by.y = "Gene.name", all.x = TRUE)

# remove rows where Gene is NA
fasta2ko <- fasta2ko[!is.na(fasta2ko$Gene), ]

# create fasta file with KO IDs (format: KO ID, ID, gene name, definition)
fasta2ko$Defline <- paste0(">", fasta2ko$ko.id, "|", fasta2ko$ID, "|", fasta2ko$Gene, "|", fasta2ko$Definition)
fasta.out <- data.frame(Defline=fasta2ko$Defline, Sequence=fasta2ko$Sequence)
fasta.out <- do.call(rbind, lapply(seq(nrow(fasta.out)), function(i) t(fasta.out[i, ])))

write.table(fasta.out, file = arguments$outputfile, row.names = FALSE, col.names = FALSE, quote = FALSE)
