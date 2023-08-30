#!/usr/bin/env Rscript

"cubseq-strain-cu-heatmap.R

Usage:
    cubseq-strain-cu-heatmap.R <input> [--fasta-dir=<fasta_dir>] [--heatmap-width=<heatmap_width>] [--heatmap-height=<heatmap_height>]
    
Options:
    --fasta-dir=<fasta_dir>                 Directory containing fasta files.
    --heatmap-width=<heatmap_width>         Specified width for heatmap figure [default: 15].
    --heatmap-height=<heatmap_height>       Specified height for heatmap figure [default: 15].
    -h --help                               Show this screen.
    --version                               Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-strain-cu-heatmap.R")

# load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(paletteer))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(Biostrings))

get_relative_freq <- function(.count_table) {
  
  # This function calculates codon relative frequencies for each 
  # strain and stores output in list
  
  # compute total codon counts
  ct <- as.data.frame(colSums(.count_table))
  ct["codon"] <- rownames(ct)
  ct$codon <- gsub("T", "U", ct$codon) # convert to RNA
  colnames(ct) <- c("codon.count", "codon")
  
  # load codon to AA map - NB: RNA_GENETIC_CODE is an object from Biostrings
  codon2aa <- data.frame(codon = names(RNA_GENETIC_CODE), aa = paste(RNA_GENETIC_CODE))
  
  # map codons to AAs
  ct <- merge(ct, codon2aa, by.x = "codon", by.y = "codon")
  
  # calculate relative frequency
  # formula used: Cij(AAi) / sum(Cij(AAi)), so that all codons for ith AA adds up to 1
  sum.aa <- ct %>% 
    group_by(aa) %>% 
    summarise(sum.aa = sum(codon.count))
  ct <- merge(ct, sum.aa, by.x = "aa", by.y = "aa")
  ct <- transform(ct, rf = codon.count / sum.aa)
  
  return(ct)
  
}

get_rf_by_strain <- function(.data) {
  
  # initialise variables
  tax_dna_list <- list()
  tax_ct_list <- list()
  tax_rf_list <- list()
  
  # split df by tax_id
  metadata_taxid <- split(.data, .data$tax_id)
  
  for (taxid in names(metadata_taxid)) {
    
    # load sequences
    seqs <- readDNAStringSet(metadata_taxid[[taxid]]$file_path)
    tax_dna_list[[taxid]] <- seqs
    
    # calculate codon counts
    ct <- trinucleotideFrequency(seqs, step=3)
    tax_ct_list[[taxid]] <- ct
    
    # calculate codon relative frequencies
    rf <- get_relative_freq(ct)
    tax_rf_list[[taxid]] <- rf
    
  }
  
  return(list(ct = tax_ct_list,
              rf = tax_rf_list))
  
}

main <- function() {
    
    # load data
    metadata <- read_delim(arguments$input, col_names = TRUE, delim = "\t")
    inputdir <- arguments$fasta_dir

    # map to associated file path to heg-fasta-dir
    metadata$file_path <- file.path(inputdir, paste0(metadata$run_accession, ".heg-mut-transcriptome.fa"))
    files <- metadata$file_path
    names(files) <- metadata$run_accession
    all(file.exists(files))

    # calculate RF by strain
    tax <- get_rf_by_strain(metadata)

    # reformat data to plot as heatmap
    df <- as.data.frame(tax$rf)
    colnames(df) <- sub("X", "", colnames(df))
    df <- df %>% 
    dplyr::rename("aa" = "562.aa", "codon" = "562.codon") %>%
    select(aa, codon, matches(".rf"))

    mat.codon <- as.matrix(df[ ,!(colnames(df) %in% c("aa", "codon"))])
    rownames(mat.codon) <- df$codon

    # save data
    saveRDS(tax, file = "strain-tax.rds")
    write_csv(df, file = "strain-rf.csv")

    # # plot as heatmap
    # png("heatmap.png", width=as.numeric(arguments$heatmap_width), height=as.numeric(arguments$heatmap_height), units="in", res=300)
    # par(mar=c(4,4,1,1))

    # Heatmap(mat.codon,
    #         name = "Codon\nrelative\nfrequency",
    #         column_title = "Heatmap of strain CU-RFs", 
    #         column_title_gp = gpar(fontsize = 20, fontface = "bold"),
    #         clustering_distance_columns = "euclidean",
    #         heatmap_height = unit(10, "npc"))
            
    # dev.off()

}

main()
