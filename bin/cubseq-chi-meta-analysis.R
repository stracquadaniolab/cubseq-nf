#!/usr/bin/env Rscript

"cubseq-chi-meta-analysis.R

Usage:
    cubseq-chi-meta-analysis.R <input> [--metadata=<metadata>] [--gtf=<gtf>] [--fdr-threshold=<fdr_threshold>] [--num-strains=<num_strains>]
    
Options:
    --metadata=<metadata>                 Metadata file.
    --gtf=<gtf>                           GTF file.
    --fdr-threshold=<fdr_threshold>       Set FDR threshold [default: 0.01].
    --num-strains=<num_strains>           Set number of strains [default: 5].
    -h --help                             Show this screen.
    --version                             Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-chi-meta-analysis.R")

# load required packages
suppressMessages(library(tidyverse))
suppressMessages(library(rtracklayer)) # to import GTF
suppressMessages(library(ggplot2))
suppressMessages(library(paletteer))
suppressMessages(library(cowplot))

extract_protein_coding_genes <- function(gtf, tpm) {

    # This function filters to keep only protein-coding genes.
    # Input: 1) GTF file, and 2) TPM matrix.
    # Output: TPM matrix containing only protein-coding genes.

    # filter GTF for protein coding genes
    proteins_gtf <- as.data.frame(gtf) %>%
        filter(gene_biotype == "protein_coding", type == "exon")

    # keep only protein coding genes in TPMs
    proteins_tpm <- subset(tpm, rownames(tpm) %in% proteins_gtf$gene_id)

    # TODO: remove ribosomal proteins

    return(proteins_tpm)

}

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

calculate_heg_leg_scores <- function(meta_results, fdr) {

    # get significant genes
    meta_sig <- meta_results[meta_results$fdr <= fdr, ]

    # get highly expressed genes
    meta_high <- meta_results[meta_results$mean_score >= 1.644, ]

    # reorder by increasing mean score
    meta_high <- meta_high[order(meta_high$mean_score), ]

    # repeat to get lowly expressed genes and reorder by increasing mean score
    meta_low <- meta_results[meta_results$mean_score <= -1.644, ]
    meta_low <- meta_low[order(meta_low$mean_score), ]

    return(list(heg = meta_high, leg = meta_low))

}

concordance_analysis <- function(metadata, meta_chi_rank, tpm, num_strains, fdr) {
  
  # extract top 5 frequently occurring strains from metadata
  strain_metadata <- metadata %>%
    filter(strain %in% names(sort(table(strain), decreasing = TRUE)[1:num_strains]))
  
  # filter TPM df samples by top 5 strains
  strain_tpm <- tpm[, (colnames(tpm) %in% strain_metadata$run_accession)]
  
  # replace TPM df colnames with strain names
  colnames(strain_tpm) <- strain_metadata$strain[match(colnames(strain_tpm), strain_metadata$run_accession)]
  
  # split TPMs by strain
  strain_tpm <- as.data.frame(strain_tpm) # needs to be in df format
  strain_tpm_list <- split.default(strain_tpm, colnames(strain_tpm))
  
  # perform chi meta-analysis per strain
  # initialise empty list
  strain_chisq_list <- list()
  strain_scores_list <- list()
  meta_strain_rank_list <- list()
  cat_results_list <- c()
  rank_max <- 117
  rank_list <- 1:rank_max
  
  for(strain_name in names(strain_tpm_list)) {
    
    # looping over names of the list, each name has associated dataframe
    strain_df <- strain_tpm_list[[strain_name]]
    
    # calculate chisq meta-analysis for meta-HEG
    strain_chisq <- chi_meta_analysis(strain_df)
    # assign chisq results dataframe per strain to list
    strain_chisq_list[[strain_name]] <- strain_chisq
    
    # get scores for HEG and LEG
    strain_scores <- calculate_heg_leg_scores(strain_chisq$results, fdr)
    # append scores dataframe per strain to list
    strain_scores_list[[strain_name]] <- strain_scores
    
    # rank each strain HEG by chisq
    strain_scores_list[[strain_name]]$heg$chisq_rank <- rank(-strain_scores_list[[strain_name]]$heg$chisq)
    #strain_scores_list[[strain_name]][["heg"]][[paste0(strain_name, ".chisq_rank")]] <- rank(-strain_scores_list[[strain_name]]$heg$chisq)
    
    # join strain rankings with meta
    meta_strain_rank <- full_join(meta_chi_rank, strain_scores_list[[strain_name]][["heg"]][, c("genes", "chisq_rank")], by = c("genes"))
    #meta_strain_rank <- full_join(meta_rank, strain_scores_list[[strain_name]][["heg"]][[, c("genes", paste0(strain_name, ".chisq_rank"))]], by = c("genes"))
    
    meta_strain_rank_list[[strain_name]] <- meta_strain_rank
    
    # calculate matches and concordance per strain with meta-HEG
    cat_results <- tibble(rank = rank_list,
                          match = rep(0, rank_max),
                          concordance = rep(0, rank_max))

    for (i in 1:rank_max) {

      meta_rank <- meta_strain_rank$meta.chisq_rank[meta_strain_rank$meta.chisq_rank <= i]
      strain_rank <- meta_strain_rank$meta.chisq_rank[meta_strain_rank$chisq_rank <= i]

      cat_results$match[i] <- length(intersect(meta_rank, strain_rank))-1
      cat_results$concordance[i] <- cat_results$match[i] / i

    }

    cat_results_list[[strain_name]] <- cat_results
    
    # join ranks and concordance per strain into one dataframe
    cat_results_df <- do.call(cbind, cat_results_list)
    
    # tidy up CAT plot df
    cat_results_df <- cat_results_df %>%
      dplyr::rename("rank" = colnames(cat_results_df[1])) %>%
      dplyr::select(., starts_with("rank"), ends_with(".concordance")) %>%
      dplyr::rename_with(~str_remove(., '.concordance'))
    
  }

  return(list(chisq = strain_chisq_list, 
              scores = strain_scores_list, 
              ranks = meta_strain_rank_list,
              cat_results = cat_results_list,
              concordance = cat_results_df))

}

main <- function() {
  
  print(arguments)

  # load data
  txi <- readRDS(arguments$input)
  metadata <- read_delim(arguments$metadata, delim = "\t", col_names = TRUE)
  gtf <- import(arguments$gtf)
  fdr <- as.numeric(arguments$fdr_threshold)
  num_strains <- as.numeric(arguments$num_strains)

  # filter TPMs to keep only protein-coding genes
  tpm <- txi$abundance
  proteins_tpm <- extract_protein_coding_genes(gtf, tpm)

  # Meta-HEG analysis
  # calculate chisq meta-analysis for meta-HEG (all strains)
  meta_chisq <- chi_meta_analysis(proteins_tpm)
  # get meta-HEG scores
  meta_scores <- calculate_heg_leg_scores(meta_chisq$results, fdr)
  # rank meta-HEG by chisq
  meta_scores$heg$meta.chisq_rank <- rank(-meta_scores$heg$chisq)
  meta_chisq_rank <- meta_scores$heg[, c("genes", "meta.chisq_rank")]
  
  # strain HEG analysis - get CAT plot
  strain_concordance_results <- concordance_analysis(metadata, meta_chisq_rank, proteins_tpm, num_strains, fdr)

  # plot CAT plot
  cat_plot <- strain_concordance_results$concordance %>%
    as_tibble() %>%
    pivot_longer(-1) %>%
    ggplot(aes(rank, value, color = name)) +
    geom_line() +
    labs(color = "Species") +
    xlab("Rank") +
    ylab("Concordance") +
    theme_cowplot() +
    scale_color_paletteer_d("ggsci::nrc_npg")

  save_plot(filename = "cat-plot.pdf", plot = cat_plot, base_width = 12, base_height = 9)

  # # statistics for log file
  # print("--- GTF STATISTICS ---")
  # print("Tally of gene biotype: ")
  # gtf_gene_stats <- as.data.frame(gtf) %>% count(gene_biotype)
  # print(gtf_gene_stats)
  
  # print("--- TPM STATISTICS ---")
  # print("TPM dimensions (original): ")
  # print(dim(tpm))
  # print("TPM dimensions (after filtering for protein-coding genes": )
  # # TODO: TPM dimensions after filtering

  # print("--- STATISTICS FOR STRAIN ANALYSIS ---")
  # print("Top strains tally: ")
  # strain_metadata_stats <- as.data.frame(filtered_metadata) %>% count(strain)
  # print(strain_metadata_stats)

}

main()
