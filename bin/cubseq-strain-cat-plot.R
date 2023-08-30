#!/usr/bin/env Rscript

"cubseq-strain-cat-plot.R

Usage:
    cubseq-strain-cat-plot.R <input> [--metadata=<metadata>] [--tax-names=<tax_names>] [--heg-scores=<heg_scores>] [--leg-scores=<leg_scores>] [--num-strains=<num_strains>] [--fdr=<fdr>]
    
Options:
    --metadata=<metadata>             CSV file containing sample metadata.
    --tax-names=<tax_names>           CSV file containing tax_id and tax_name columns.
    --heg-scores=<heg_scores>         CSV file containing gene rankings for highly expressed genes.
    --leg-scores=<leg_scores>         CSV file containing gene rankings for lowly expressed genes.
    --num-strains=<num_strains>       Set number of strains [default: 5].
    --fdr=<fdr>                       Set FDR threshold [default: 0.01].
    -h --help                         Show this screen.
    --version                         Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-strain-cat-plot.R")

# load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(paletteer))

chi_meta_analysis <- function(.dataset) {
  
  # This function identifies significantly highly expressed genes
  # from TPM values.
  # Input: TPM matrix (output from tximport).
  # Output: list containing 1) original ranked data and 2) results
  # of genes with corresponding p-value and FDR.
  
  # extract number of rows and columns from dataset
  .genes <- nrow(.dataset)
  .samples <- ncol(.dataset)
  
  # 1) rank genes per sample; 2) normalise ranks; 3) inverse normal transform
  .rank_data <- apply(.dataset, 2, rank, ties.method = "random")
  .rank_data <- .rank_data / (.genes + 1)
  .norm_data <- qnorm(.rank_data)
  
  # 1) sum the INTs; 2) get mean score; 3) get chi scores
  .norm_sum <- rowSums(.norm_data)
  .mean_score <- rowMeans(.norm_data)
  .chi_scores <- rowSums(.norm_data**2)
  
  # 1) calculate chisq p-value; 2) calculate fdr
  .pvalues <- pchisq(.chi_scores, df = .samples, lower.tail = F)
  .fdr <- p.adjust(.pvalues, method = "fdr")
  
  # summarise results into dataframe
  .results <- data.frame(genes = rownames(.dataset), 
                         score = .norm_sum, 
                         mean_score = .mean_score, 
                         chisq = .chi_scores, 
                         pvalue = .pvalues, 
                         fdr = .fdr)
  
  # return list containing: 1) original ranked data; and 2) results with scores
  return(list(data = .rank_data, results = .results))
  
}


calculate_heg_leg_scores <- function(.meta_results, .fdr) {
  
  # get significant genes
  meta_sig <- .meta_results[.meta_results$fdr <= .fdr, ]
  
  # get highly expressed genes
  meta_high <- .meta_results[.meta_results$mean_score >= 1.644, ]
  
  # reorder by increasing mean score
  meta_high <- meta_high[order(meta_high$mean_score), ]
  
  # repeat to get lowly expressed genes and reorder by increasing mean score
  meta_low <- .meta_results[.meta_results$mean_score <= -1.644, ]
  meta_low <- meta_low[order(meta_low$mean_score), ]
  
  return(list(heg = meta_high, leg = meta_low))
  
}


concordance_analysis <- function(.metadata, .gene_ranks, .tpm, .num_taxid, .fdr) {
  
  # extract top n frequently occurring tax_ids from metadata
  tax_metadata <- .metadata %>%
    filter(tax_id %in% names(sort(table(tax_id), decreasing = TRUE)[1:.num_taxid]))
  
  # filter TPM df samples by top n tax_ids
  tax_tpm <- .tpm[, (colnames(.tpm) %in% tax_metadata$run_accession)]
  
  # replace TPM df colnames with tax_id names
  colnames(tax_tpm) <- tax_metadata$tax_id[match(colnames(tax_tpm), tax_metadata$run_accession)]
  
  # split TPMs by strain
  tax_tpm <- as.data.frame(tax_tpm) # needs to be in df format
  tax_tpm_list <- split.default(tax_tpm, colnames(tax_tpm))
  
  # perform chi meta-analysis per strain
  # initialise empty list
  tax_chisq_list <- list()
  tax_scores_list <- list()
  meta_tax_rank_list <- list()
  cat_results_list <- c()
  rank_max <- nrow(.gene_ranks)
  rank_list <- 1:rank_max
  
  for (tax_name in names(tax_tpm_list)) {
    
    # looping over names of the list, each name has associated dataframe
    tax_df <- tax_tpm_list[[tax_name]]
    
    # calculate chisq meta-analysis for meta-HEG
    tax_chisq <- chi_meta_analysis(tax_df)
    # assign chisq results dataframe per strain to list
    tax_chisq_list[[tax_name]] <- tax_chisq
    
    # get scores for HEG and LEG
    tax_scores <- calculate_heg_leg_scores(tax_chisq$results, .fdr)
    # append scores dataframe per strain to list
    tax_scores_list[[tax_name]] <- tax_scores
    
    # rank each strain HEG by chisq
    tax_scores_list[[tax_name]]$heg$chisq_rank <- rank(-tax_scores_list[[tax_name]]$heg$chisq)
    #tax_scores_list[[tax_name]][["heg"]][[paste0(tax_name, ".chisq_rank")]] <- rank(-tax_scores_list[[tax_name]]$heg$chisq)
    
    # join strain rankings with meta
    meta_tax_rank <- full_join(.gene_ranks, tax_scores_list[[tax_name]][["heg"]][, c("genes", "chisq_rank")], by = c("genes"))
    #meta_tax_rank <- full_join(meta_rank, tax_scores_list[[tax_name]][["heg"]][[, c("genes", paste0(tax_name, ".chisq_rank"))]], by = c("genes"))
    
    meta_tax_rank_list[[tax_name]] <- meta_tax_rank
    
    # calculate matches and concordance per strain with meta-HEG
    cat_results <- tibble(rank = rank_list,
                          match = rep(0, rank_max),
                          concordance = rep(0, rank_max))
    
    for (i in 1:rank_max) {
      
      meta_rank <- meta_tax_rank$meta.chisq_rank[meta_tax_rank$meta.chisq_rank <= i]
      tax_rank <- meta_tax_rank$meta.chisq_rank[meta_tax_rank$chisq_rank <= i]
      
      cat_results$match[i] <- length(intersect(meta_rank, tax_rank))-1
      cat_results$concordance[i] <- cat_results$match[i] / i
      
    }
    
    cat_results_list[[tax_name]] <- cat_results
    
    # join ranks and concordance per strain into one dataframe
    cat_results_df <- do.call(cbind, cat_results_list)
    
    # tidy up CAT plot df
    cat_results_df <- cat_results_df %>%
      dplyr::rename("rank" = colnames(cat_results_df[1])) %>%
      dplyr::select(., starts_with("rank"), ends_with(".concordance")) %>%
      dplyr::rename_with(~str_remove(., '.concordance'))
    
  }
  
  return(list(chisq = tax_chisq_list, 
              scores = tax_scores_list, 
              ranks = meta_tax_rank_list,
              cat_results = cat_results_list,
              concordance = cat_results_df))
  
}

get_cat_plot <- function(.data, .ec_names) {

    # TODO: fix later, tax_names match to tax_id but cols get mixed up
    # # map and paste NCBI taxdump tax names along with tax IDs for figure legend
    # colnames(.data$concordance)[na.omit(match(.ec_names$tax_id, names(.data$concordance)))] <- paste0(names(.data$concordance)[-1], " (", (.ec_names$tax_name[na.omit(match(names(.data$concordance)[-1], .ec_names$tax_id))]), ")")

    cat_plot <- .data$concordance %>%
    as_tibble() %>%
    pivot_longer(-1) %>%
    ggplot(aes(rank, value, color = name)) +
    geom_line() +
    labs(color = "Taxonomy ID") +
    xlab("Rank") +
    ylab("Concordance") +
    theme_cowplot() +
    scale_color_paletteer_d("ggsci::nrc_npg")
    
}

main <- function() {
  
  # load data
  txi <- readRDS(arguments$txi)
  metadata <- read_delim(arguments$metadata, delim = "\t", col_names = TRUE)
  ec_taxdump <- read_csv(arguments$tax_names) # NCBI taxdump file for figure legend
  heg_scores <- read_delim(arguments$heg_scores)
  leg_scores <- read_delim(arguments$leg_scores)
  
  # initialise variables
  fdr <- as.numeric(arguments$fdr)
  num_taxid <- as.numeric(arguments$num_taxid)
  
  # get TPMs
  tpm <- txi$abundance
  
  # get CAT plot for highly and lowly expressed genes
  heg_tax_concordance <- concordance_analysis(metadata, heg_scores, tpm, num_taxid, fdr)
  leg_tax_concordance <- concordance_analysis(metadata, leg_scores, tpm, num_taxid, fdr)
  
  # # plot CAT plot
  # cat_plot <- heg_tax_concordance$concordance %>%
  #   as_tibble() %>%
  #   pivot_longer(-1) %>%
  #   ggplot(aes(rank, value, color = name)) +
  #   geom_line() +
  #   labs(color = "Taxonomy ID") +
  #   xlab("Rank") +
  #   ylab("Concordance") +
  #   theme_cowplot() +
  #   scale_color_paletteer_d("ggsci::nrc_npg")

  # generate figures
  cat_plot_heg <- get_cat_plot(heg_tax_concordance, ec_taxdump)
  cat_plot_leg <- get_cat_plot(leg_tax_concordance, ec_taxdump)

  # save data
  saveRDS("heg_tax_concordance.rds", heg_tax_concordance)
  saveRDS("leg_tax_concordance.rds", leg_tax_concordance)
  
  # save plots
  save_plot(filename = "cat-plot-HEG.pdf", plot = cat_plot_heg, base_width = 12, base_height = 9)
  save_plot(filename = "cat-plot-LEG.pdf", plot = cat_plot_leg, base_width = 12, base_height = 9)
  
}

main()
