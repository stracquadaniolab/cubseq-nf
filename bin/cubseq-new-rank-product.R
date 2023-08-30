#!/usr/bin/env Rscript

"cubseq-new-rank-product.R

Usage:
    cubseq-new-rank-product.R <input> [--gtf=<gtf>] [--permutations=<permutations>] [--fdr-threshold=<fdr_threshold>] [--output-real-rp=<output_real_rp>] [--output-sim-rp=<output_sim_rp>] [--output-summary-rp=<output_summary_rp>]
    
Options:
    --gtf=<gtf>                               GTF file.
    --permutations=<permutations>             Number of permutations to run [default: 5000].
    --fdr-threshold=<fdr_threshold>           Set FDR threshold [default: 0.05].
    --output-real-rp=<output_real_rp>         Output file name for observed RP.
    --output-sim-rp=<output_sim_rp>           Output file name for simulated RP.
    --output-summary-rp=<output_summary_rp>   Final output file name containing genes, real RP, average expected RP and % FP.
    -h --help                                 Show this screen.
    --version                                 Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-new-rank-product.R")

# load required packages
suppressMessages(library(tidyverse))

geometric_mean <- function(x, na.rm = TRUE) {

  # Function calculates and returns geometric mean
  gm <- exp(mean(log(x[x>0]), na.rm=na.rm))

  return(gm)
  
}

compute_rp <- function(ranked_tpm) {
  
  # compute rank product via geometric mean (row-wise)
  rp <- apply(ranked_tpm, 1, function(x) geometric_mean(x))
  rp <- as.matrix(rp)
  
  return(rp)
  
}

compute_simulated_rp <- function(ranked_tpm) {
  
  # shuffle original ranked TPMs (column-wise) to generate simulated dataset
  sim.tpm.rank <- apply(ranked_tpm, 2, sample)
  rownames(sim.tpm.rank) <- rownames(ranked_tpm)
  
  # compute rank product via geometric mean
  sim.rp <- compute_rp(sim.tpm.rank)
  
  return(sim.rp)
  
}

compute_expected_rp <- function(real_rp, sim_rp, permutations) {
  
  # merge df by rowname
  all.rp <- merge(real_rp, sim_rp, by="row.names", all.y=TRUE)
  
  # calculate how many simulated genes have RP >= real genes
  all.rp$exp.rp <- rowSums(all.rp[,-1:-2] >= all.rp$V1.x)

  # calculate average expected RP per gene
  all.rp$avg.exp.rp <- all.rp$exp.rp / permutations
  
  # select rows for final output
  all.rp <- all.rp %>%
    select(Row.names, V1.x, exp.rp, avg.exp.rp) %>%
    rename(real.rp = V1.x, genes = Row.names)
  
  return(all.rp)
  
}

main <- function() {

    # load data
    data <- readRDS(arguments$input)
    gtf <- read.table(arguments$gtf, header = FALSE, sep = '\t')

    # initialise variables
    permutations <- as.numeric(arguments$permutations)

    print(arguments) # logs
    cat("\nDimensions of data:") # logs
    print(dim(data$abundance)) # logs

    # rank each gene by its TPM (stored in "abundance" column)
    # N.B. gene with lowest TPM = rank 1
    tpm.rank <- apply(data$abundance, 2, rank, ties.method = "average")
    cat("\nDimensions of ranked TPMs:")
    print(dim(tpm.rank)) # logs

    # compute rank product via geometric mean on real and simulated datasets
    real.rp <- compute_rp(tpm.rank)
    cat("\nDimensions and head of real RP:\n") # logs
    print(dim(real.rp)) # logs
    print(head(real.rp)) # logs
    set.seed(42)
    sim.rp <- replicate(permutations, compute_simulated_rp(tpm.rank), simplify = "matrix")
    rownames(sim.rp) <- rownames(tpm.rank)
    cat("\nDimensions and head of simulated RP:\n") # logs
    print(dim(sim.rp)) # logs
    print(head(sim.rp)) # logs

    # compute average expected rank product
    summary.rp <- compute_expected_rp(real.rp, sim.rp, permutations)
    cat("\nHead of summary RP with average expected RP:\n")
    print(head(summary.rp))

    # compute percentage false positives
    summary.rp$pfp <- summary.rp$avg.exp.rp / summary.rp$real.rp
    cat("\nHead of summary RP with % FP:\n")
    print(head(summary.rp))

    # outputs
    write.table(real.rp, file = arguments$output_real_rp, quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(sim.rp, file = arguments$output_sim_rp, quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(summary.rp, file = arguments$output_summary_rp, quote = FALSE, row.names = FALSE, col.names = TRUE)
    
    # TODO: add histogram plots for p-value, simulated RP and % FP.

}

main()
