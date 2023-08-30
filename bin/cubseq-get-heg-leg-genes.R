#!/usr/bin/env Rscript

"cubseq-get-heg-leg-genes.R

Usage:
    cubseq-get-heg-leg-genes.R <input> [--gtf=<gtf>] [--proportion-genes=<proportion_genes>] [--permutations=<permutations>] [--output-heg=<output_heg>] [--output-leg=<output_leg>]
    
Options:
    --gtf=<gtf>                               GTF file.
    --proportion-genes=<proportion_genes>     Proportion of genes to subset [default: 0.1].
    --permutations=<permutations>             Number of permutations to run [default: 100000].
    --output-heg=<output_heg>                 Output file name for high expression gene transcript IDs [default: heg-TrID.csv].
    --output-leg=<output_leg>                 Output file name for low expression gene transcript IDs [default: leg-TrID.csv].
    -h --help                                 Show this screen.
    --version                                 Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-get-heg-leg-genes.R")

# load packages
suppressMessages(library(tidyverse))

# # OLD P-VALUE CODE

# get_p_value <- function(data, n) {

#   # This function calculates the p-value threshold for the top and bottom 5% genes in a given dataset. This can act as a significance threshold to determine genes that are significantly highly and lowly expressed.
#   # Arguments:
#   # data = dataframe containing genes (rows) and abundance/TPM values per sample run (cols).
#   # n = number times to repeat.
#   # Returns: a numeric vector containing p-value for HEG and LEG, respectively.
  
#   heg.vals <- c()
#   leg.vals <- c()
#   set.seed(42) # for reproducibility
  
#   for (i in 1:n) {
    
#     # shuffle column-wise
#     tpm.shuffle <- apply(data, 2, sample)
#     rownames(tpm.shuffle) <- rownames(data)
    
#     # calculate geometric mean of each ranked gene (dim = 1 col)
#     tpm.gm <- apply(tpm.shuffle, 1, function(x) geometric_mean(x))
#     tpm.gm.df <- as.data.frame(tpm.gm)
    
#     # order genes by geometric rank and take top 5% of genes
#     tpm.heg <- tpm.gm.df %>%
#       slice_max(order_by = tpm.gm, prop = 0.05)
#     # take bottom 10% of genes
#     tpm.leg <- tpm.gm.df %>%
#       slice_min(order_by = tpm.gm, prop = 0.05)
    
#     heg.min <- max(tpm.heg) # get min value in top 5%
#     leg.max <- min(tpm.leg) # get max value in bottom 5%
    
#     # append to vals
#     heg.vals <- append(heg.vals, heg.min)
#     leg.vals <- append(leg.vals, leg.max)
    
#   }
  
#   # get p-value (in 0.05th position)
#   heg.vals <-sort(heg.vals)
#   leg.vals <- sort(leg.vals)
#   pos <- 0.05*n
#   p.heg <- heg.vals[pos]
#   p.leg <- leg.vals[pos]
  
#   return(c(p.heg, p.leg))
  
# }

# main <- function() {

#   # load files
#   txi <- readRDS(arguments$input)
#   gtf <- read.table(arguments$gtf, header = FALSE, sep = '\t')
  
#   # get HEG and LEG
#   # rank gene TPMs per sample run (this normalises for differing sample library sizes)
#   tpm.rank <- apply(txi$abundance, 2, rank, ties.method = "average") # TPM information is stored in abundance col
#   # get p-value for HEG and LEG
#   p <- get_p_value(tpm.rank, 100000)
#   # calculate geometric mean of ranked TPMs for each gene
#   tpm.gm <- apply(tpm.rank, 1, function(x) geometric_mean(x)) # '1' to apply over rows
#   tpm.gm.df <- as.data.frame(tpm.gm)
#   # filter HEG and LEG by p-value cut-off
#   heg <- filter(tpm.gm.df, tpm.gm >= p[1])
#   leg <- filter(tpm.gm.df, tpm.gm <= p[2])
  
#   # get HEG and LEG transcript IDs
#   heg.trID <- geneID2trID(gtf, heg)
#   leg.trID <- geneID2trID(gtf, leg)
  
#   # write transcript IDs to .txt file
#   write.table(heg.trID$transcript_id, file = arguments$output_heg, quote = FALSE, row.names = FALSE, col.names = FALSE)
#   write.table(leg.trID$transcript_id, file = arguments$output_leg, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
# }

# main()

# NEW P-VALUE CODE

geometric_mean <- function(x, na.rm = TRUE) {
  # Function calculates and returns geometric mean of each gene
  gm <- exp(mean(log(x[x>0]), na.rm=na.rm))
  return(gm)
}

get_shuffled_gm <- function(data, proportion_genes) {

  # Function shuffles original dataset of ranked TPMs, and calculates geometric mean of the
  # shuffled ranks. Geometric means are then ordered (ascending value) and a set proportion
  # of max and min genes are taken (proportion set by user).
  # 
  # ARGUMENT "data" -- ranked TPMs of original tximport dataset.
  # ARGUMENT "proportion_genes" -- proportion of genes to subset.
  # RETURNS : list containing [1] genes with highest geometric means, and [2] genes with
  # lowest geometric means.
  
  # shuffle original dataset of ranked TPMs (column-wise)
  tpm.shuffle <- apply(data, 2, sample)
  rownames(tpm.shuffle) <- rownames(data)

  # calculate geometric mean of shuffled ranks (row-wise)
  tpm.gm.shuffle <- apply(tpm.shuffle, 1, function(x) geometric_mean(x))
  tpm.gm.shuffle <- as.data.frame(tpm.gm.shuffle)

  # take genes with highest geometric means
  tpm.heg.shuffle <- tpm.gm.shuffle %>% slice_max(order_by = tpm.gm.shuffle, prop = proportion_genes)
  heg.max.shuffle <- max(tpm.heg.shuffle) # take max value if we want to be conservative

  # TODO: we can just directly take min
  # take genes with lowest geometric means
  tpm.leg.shuffle <- tpm.gm.shuffle %>% slice_min(order_by = tpm.gm.shuffle, prop = proportion_genes)
  leg.min.shuffle <- min(tpm.leg.shuffle) # take min value if we want to be conservative

  # store values in vector
  return(c(heg.max.shuffle, leg.min.shuffle))

}

get_p_value <- function(perms, data, proportion_genes) {

  # Function performs permutation testing to generate a null distribution, and returns
  # values in the confidence thresholds for high and low expression genes.
  #
  # ARGUMENT "perms" -- number of permutations to perform.
  # ARGUMENT "data" -- this is passed onto the get_shuffled_gm() function.
  # ARGUMENT "proportion_genes" -- this is passed onto the get_shuffled_gm() function.

  set.seed(42)
  # perform permutations to generate null distribution for HEG and LEG
  heg.vals <- replicate(perms, get_shuffled_gm(data, proportion_genes)[1], simplify = "vector")
  leg.vals <- replicate(perms, get_shuffled_gm(data, proportion_genes)[2], simplify = "vector")
  
  # sort values (ascending order)
  heg.vals <- sort(heg.vals)
  leg.vals <- sort(leg.vals)

  # set confidence intervals for HEG and LEG
  heg.pos <- 0.95*perms
  leg.pos <- 0.05*perms

  # get values in confidence interval
  p.heg <- heg.vals[heg.pos]
  p.leg <- leg.vals[leg.pos]

  return(c(p.heg, p.leg))

}

# TODO: to determine adequate number of replications could do a meta-permutation of get_p_values(), test range of permutations.

getTrID <- function(gtf) {
  
  # Function creates map linking gene IDs with transcript IDs from gtf file.
  
  trID.map <- gtf %>%
    # filter to keep only transcripts
    filter(., grepl("transcript", V3)) %>%
    # extract gene ID and transcript ID from defline
    extract(
      ., V9, into = c("gene_id", "transcript_id"),
      regex = "gene_id\\s*(.*?); transcript_id\\s*(.*?);"
    ) %>%
    dplyr::select(c(gene_id, transcript_id))
  
  return(trID.map)
  
}

main <- function() {

  # load data
  data <- readRDS(arguments$input)
  gtf <- read.table(arguments$gtf, header = FALSE, sep = '\t')
  
  # initialise variables
  proportion_genes <- as.numeric(arguments$proportion_genes)
  perms <- as.numeric(arguments$permutations)

  # get ranked TPMs (stored in "abundance" column), N.B. genes with higher TPMs = larger rank value
  tpm.rank <- apply(data$abundance, 2, rank, ties.method = "average")

  # calculate confidence intervals for HEG and LEG
  p.heg <- get_p_value(perms, tpm.rank, proportion_genes)[1]
  print('Computed p.heg: ')
  print(p.heg)
  p.leg <- get_p_value(perms, tpm.rank, proportion_genes)[2]
  print('Computed p.leg: ')
  print(p.leg)
  
  # calculate observed geometric means of original ranked TPMs
  obs.gm <- apply(tpm.rank, 1, function(x) geometric_mean(x))
  print('Computed obs.gm: ')
  print(obs.gm)
  obs.gm <- as.data.frame(obs.gm)
  obs.heg <- obs.gm %>% slice_max(order_by = obs.gm, prop = proportion_genes)
  print('Computed obs.heg: ')
  print(obs.heg)
  obs.leg <- obs.gm %>% slice_min(order_by = obs.gm, prop = proportion_genes)
  print('Computed obs.leg: ')
  print(obs.leg)
  
  # filter HEG and LEG by p-value cut-off
  sig.heg <- filter(obs.heg, obs.gm >= p.heg)
  print('Computed sig.heg: ')
  print(sig.heg)
  sig.leg <- filter(obs.leg, obs.gm <= p.leg)
  print('Computed sig.leg: ')
  print(sig.leg)
  
  # create map of gene IDs to transcript IDs
  trID.map <- getTrID(gtf)
  print('Computed trID.map: ')
  print(head(trID.map))
  
  # map HEG and LEG to their transcript IDs
  heg.trID <- subset(trID.map, gene_id %in% rownames(sig.heg))
  print('heg.trID mapping done: ')
  print(head(heg.trID))
  leg.trID <- subset(trID.map, gene_id %in% rownames(sig.leg))
  print('leg.trID mapping done: ')
  print(head(leg.trID))
  
  # write transcript IDs to .txt file
  write.table(heg.trID$transcript_id, file = arguments$output_heg, quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(leg.trID$transcript_id, file = arguments$output_leg, quote = FALSE, row.names = FALSE, col.names = FALSE)

}

main()
