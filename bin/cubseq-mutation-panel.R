#!/usr/bin/env Rscript

"cubseq-mutation-panel.R

Usage:
    cubseq-mutation-panel.R <input> [--gtf=<gtf>] [--gene_tpm_scores=<gene_tpm_scores>] [--fdr-threshold=<fdr>] 
    
Options:
    --gtf=<gtf>                               GTF file.
    --gene_tpm_scores=<gene_tpm_scores>       Gene TPM scores in csv format.
    --fdr-threshold=<fdr>                     FDR threshold.
    -h --help                                 Show this screen.
    --version                                 Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-mutation-panel.R")

# load libraries
# suppressMessages(library(vcfR))
suppressMessages(library(dplyr))
suppressMessages(library(readr))
suppressMessages(library(tidyr))
suppressMessages(library(rtracklayer))
suppressMessages(library(circlize))

# # load data
# print("Loading data.")
# vcf <- read.vcfR(arguments$input, verbose = FALSE)
# vcf <- cbind(as.data.frame(getFIX(vcf)), INFO2df(vcf))

# # Select columns POS and AC
# print("Selecting POS and AC columns.")
# var <- vcf[, c("POS", "AC")]

# # Separate rows for multiple values in AC column
# split_AC <- unlist(strsplit(var$AC, ","))
# var <- var[rep(seq_len(nrow(var)), times = lengths(strsplit(var$AC, ","))), ]

# # Update AC column
# var$AC <- as.numeric(split_AC)

# # Calculate pos_kb
# var$pos_kb <- as.numeric(var$POS) / 1000

# # Calculate AC_norm
# print("Calculating AC norm.")
# var$AC_norm <- var$AC / 6763

# # save data
# write_csv(var, "mutation-frequency-data.csv")
# print("Done.")

# load data
print("Loading data.")
var <- read_csv(arguments$input)
gtf <- rtracklayer::import(arguments$gtf)
gene_tpm_scores <- read_csv(arguments$gene_tpm_scores)

# reformat GTF
print("Reformatting GTF file.")
gtf <- as.data.frame(gtf) %>%
  filter(., type == "exon") %>%
  mutate(start = as.numeric(start),
         end = as.numeric(end))

# filter for significant highly expressed genes
message("Filtering for highly expressed genes.")
fdr <- arguments$fdr
meta_sig <- gene_tpm_scores[gene_tpm_scores$fdr <= fdr, ]
meta_high <- gene_tpm_scores[gene_tpm_scores$mean_score >= 1.644, ]
dim(meta_sig)
dim(meta_high)

# map gene_id to positional coordinates
message("Mapping gene_id to POS.")
result <- var %>%
  rowwise() %>%
  mutate(gene_id = if_else(POS >= min(gtf$start) & POS <= max(gtf$end),
                           if (any(POS >= gtf$start & POS <= gtf$end)) {
                             paste(gtf$gene_id[POS >= gtf$start & POS <= gtf$end], collapse = ",")
                           } else {
                             NA
                           },
                           NA_character_)) %>%
  # identify highly expressed genes
  separate_rows(gene_id, sep = ",") %>%
  mutate(heg = ifelse(gene_id %in% meta_high$genes, "TRUE", "FALSE")) %>%
  # map to abundance scores and FDR
  left_join(gene_tpm_scores, by = c("gene_id" = "genes")) %>%
  mutate(
    mean_score = ifelse(!is.na(mean_score), mean_score, score),
    fdr = ifelse(!is.na(fdr), fdr, pvalue)
  ) %>%
  select(-score, -pvalue)

# plot as circos
message("Starting plotting as circos...")
# add sector
result$sectors <- 1

# Set up the PDF device
pdf("circlize_mutations.pdf", width = 10, height = 10)

# set track height
circos.par("track.height" = 0.3, start.degree = 90, gap.degree = 10)

# add track for pos_kb
circos.initialize(result$sectors, x = result$pos_kb)

# add track for AC (frequency of mutation)
message("Adding track for mutation frequency.")
circos.track(result$sectors, y = result$AC_norm,
             panel.fun = function(x, y) {
               # set x-axis font size
               circos.axis(labels.cex = 0.8)
               # set y-axis font size
               par(cex = 0.8)
               circos.yaxis(side = "left", 
                            at = c(0, 0.2, 0.4, 0.6, 0.8, 1), # set y-axis tick values
                            labels = TRUE, 
                            labels.cex = par("cex"))
             })

circos.trackLines(result$sectors, 
                  result$pos_kb, 
                  result$AC_norm, 
                  col = "#6DA2CF", 
                  pch = 16, 
                  cex = 0.5, 
                  type = "h")

# setting any mean_score with NA values to 0 (circlize can't handle NA)
result$mean_score[is.na(result$mean_score)] <- 0

message("Adding track for gene TPM mean scores.")
circos.track(result$sectors, 
             x = result$pos_kb, 
             y = result$mean_score,
             panel.fun = function(x, y) {
               # gap.after = 20
               # N.B. we do not set x-axis labels here as this is in first track
               # set y-axis font size for second track
               par(cex = 0.8)
               circos.yaxis(side = "left", 
                            # at = c(-2, -1, 0, 1, 2, 3.1), # set y-axis tick values
                            labels = TRUE, 
                            labels.cex = par("cex"))
             })

# Define the colours
score_cols <- c("#A8DCB8", "#de2d26")

# Map the result$heg category to the corresponding color
score_cols_mapping <- score_cols[as.numeric(factor(result$heg))]

circos.trackLines(result$sectors, 
                  result$pos_kb, 
                  result$mean_score, 
                  col = score_cols_mapping, 
                  baseline = 0,
                  pch = -16, 
                  cex = 0.5, 
                  type = "h")

# Create the Circos image
message("Creating circos image.")
circos.trackPlotRegion(ylim = NULL)

# Close the PDF device
message("Closing PDF device.")
dev.off()