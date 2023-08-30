#!/usr/bin/env Rscript

"cubseq-summarise-codon-counts.R

Usage:
    cubseq-summarise-codon-counts.R <input>
    
Options:
    -h --help                                   Show this screen.
    --version                                   Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-summarise-codon-counts.R")

# load packages
suppressMessages(library(tidyverse))
suppressMessages(library(Biostrings))

# load raw codon counts dataset - HEG
# counts <- read.table("heg-codon-counts.csv", row.names = NULL)
counts <- read.table(arguments$input, header = TRUE, sep = " ", row.names = NULL)

# combine rows and sum codon counts
count.all <- counts %>%
  group_by(row.names) %>%
  summarise_all(funs(sum)) %>%
  pivot_longer(cols = c(-row.names))

# load codon to AA map - NB: GENETIC_CODE is an object from Biostrings
codon2aa <- data.frame(codon = names(GENETIC_CODE), aa = paste(GENETIC_CODE))

# map codons to AAs
count.all <- merge(count.all, codon2aa, by.x = "name", by.y = "codon")

# rename columns and reorder
colnames(count.all) <- c("codon", "transcript.ID", "codon.count", "aa")
count.all <- count.all[,c("transcript.ID", "codon", "aa", "codon.count")]

# group codons by amino acid and sum respective counts
count.aa <- count.all %>% 
  group_by(transcript.ID, aa) %>% 
  summarise(aa.count = sum(codon.count))

# combine
count.all <- merge(count.all, count.aa, by.x = c("transcript.ID", "aa"), 
                   by.y = c("transcript.ID", "aa"))

# write to csv
write_csv(count.all, "summary-counts.csv")
write_csv(count.aa, "aa-counts.csv")