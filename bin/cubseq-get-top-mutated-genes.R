#!/usr/bin/env Rscript

"cubseq-get-top-mutated-genes.R

Usage:
    cubseq-get-top-mutated-genes.R <input> [--annotated-mut-dir=<annotated_mut_dir>]
    
Options:
    --annotated-mut-dir=<annotated_mut_dir>   Directory containing annotated VCF files in GTF format (output from bedtools annotate).
    -h --help                                 Show this screen.
    --version                                 Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-get-top-mutated-genes.R")

# load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

# load data
metadata <- read_delim(arguments$input, col_names = TRUE, delim = "\t")
inputdir <- arguments$annotated_mut_dir

# map to associated annotated vcf-gtf file path
metadata$file_path <- file.path(inputdir, paste0(metadata$run_accession, ".annotated-vcf.gtf"))
files <- metadata$file_path
names(files) <- metadata$run_accession
print("All files exists:")
all(file.exists(files))

# import all gtf files into single df
print("Binding all annotated GTF files into vcf dataframe...")
vcf <- do.call(rbind, lapply(files, read.table, header = FALSE, sep = "\t"))
print("Done.")
print("Dimensions of vcf dataframe:")
dim(vcf)

# tidy up df
gtf_names <- c("seqname", "source", "feature", "start", "end", 
               "score", "strand", "frame", "attribute", "count",
               "coverage_fraction")
colnames(vcf) <- gtf_names

# filter and extract attributes
print("Filtering and extracting attributes...")
vcf <- vcf %>%
  filter(., grepl("CDS", feature)) %>%
  extract(., attribute, into = c("gene_id", "transcript_id", "exon_number",
                                 "gene_name", "gene_source", "gene_biotype", 
                                 "transcript_name", "transcript_source",
                                 "transcript_biotype", "protein_id"),
          regex =  "gene_id\\s*(.*?); transcript_id\\s*(.*?); exon_number\\s*(.*?); gene_name\\s*(.*?); gene_source\\s*(.*?); gene_biotype\\s*(.*?); transcript_name\\s*(.*?); transcript_source\\s*(.*?); transcript_biotype\\s*(.*?); protein_id\\s*(.*?);") %>%
  filter(., grepl("protein_coding", gene_biotype)) %>%
  dplyr::select(c(start, end, gene_id, transcript_id, gene_name, count, coverage_fraction))
print("Done")

# new method: count number of times each gene is mutated per sample
sample_mut_count <- vcf %>%
  group_by(gene_name) %>%
  summarise(count = sum(count > 0))

# plot barplot of top mutated genes
top_mutated_gene_count_barplot_20 <- sample_mut_count %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  ggplot(., aes(y = reorder(gene_name, count), x = count)) +
  geom_bar(stat = "identity") +
  geom_text(aes(x = count, y = gene_name, label = count), hjust = -0.25) +
  labs(x = "Number of samples with mutated gene", y = "Gene") +
  theme_cowplot(12)

# group and sum count columns by gene
print("Grouping and summing count columns by gene...")
mut <- vcf %>%
  mutate(length = end - start) %>%
  mutate(fraction = count / length) %>%
  group_by(gene_name) %>%
  summarise(mut.fraction = sum(fraction))
print("Done.")

# plot barplot of top 100 mutated genes
print("Plotting top 100 and top 50 mutated genes...")
top_mutated_genes_barplot_100 <- mut %>%
  arrange(desc(mut.fraction)) %>%
  slice(1:100) %>%
  ggplot(., aes(y = reorder(gene_name, mut.fraction), x = mut.fraction)) +
  geom_bar(stat = "identity") +
  labs(x = "Fraction of mutations", y = "Gene") +
  theme_cowplot(12)

top_mutated_genes_barplot_50 <- mut %>%
  arrange(desc(mut.fraction)) %>%
  slice(1:50) %>%
  ggplot(., aes(y = reorder(gene_name, mut.fraction), x = mut.fraction)) +
  geom_bar(stat = "identity") +
  labs(x = "Fraction of mutations", y = "Gene") +
  theme_cowplot(12)
print("Done.")

# save data
print("Saving data...")
saveRDS(vcf, "vcf-bedtools-annotate-counts.rds")
write_csv(mut, "gene-mutation-fraction.csv")
write_csv(sample_mut_count, "sample-mutation-counts.csv")
print("Done.")

# save figures
print("Saving plots...")
save_plot(filename = "top_mutated_gene_count_barplot_20.pdf", plot = top_mutated_gene_count_barplot_20, base_height = 15, base_width = 4)
save_plot(filename = "top-mutated-genes-barplot_100.pdf", plot = top_mutated_genes_barplot_100, base_height = 15, base_width = 4)
save_plot(filename = "top-mutated-genes-barplot_50.pdf", plot = top_mutated_genes_barplot_50, base_height = 8, base_width = 4)
print("Done.")
print("Finished!")