#!/usr/bin/env Rscript

"cubseq-fisher-transform.R

Usage:
    cubseq-fisher-transform.R <input> [--metadata=<metadata>] [--gtf=<gtf>] [--gff=<gff>] [--fdr-threshold=<fdr_threshold>]
    
Options:
    --metadata=<metadata>                           Metadata file.
    --gtf=<gtf>                                     GTF file.
    --gff=<gff>                                     GFF file.
    --fdr-threshold=<fdr_threshold>                 Set FDR threshold [default: 0.01].
    -h --help                                       Show this screen.
    --version                                       Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-fisher-transform.R")

# load packages
suppressMessages(library(readr))
suppressMessages(library(dplyr))
suppressMessages(library(rtracklayer))

extract_protein_coding_genes <- function(.gff, .gtf, .tpm) {
  
  # This function filters to keep only protein-coding genes,
  # filters out genes encoding ribosomal proteins and removes
  # any genes with a CDS length not divisible by 3.
  # Input: 1) GTF file, and 2) TPM matrix.
  # Output: TPM matrix containing only protein-coding genes.

  # Get list of ribosomal proteins (RPs) from GFF file.
  # We use GFF because it contains a "description field"
  # that allows us to identify RPs.
  # Here, we consider RPs to be any proteins that make up the
  # core ribosomal structure (i.e. not co-factors). This should
  # give a final list of 57 RPs (55 RPs plus 2 putative RPs).
  rp_gff <- as.data.frame(.gff) %>%
  filter(., type == "gene") %>%
  filter(., ifelse("biotype" %in% names(.), biotype == "protein_coding", gene_biotype == "protein_coding")) %>%
  filter(grepl("ribosomal subunit protein", description) | 
           grepl("putative ribosomal protein", description)) %>%
  # filter out any enzymes (e.g. methyltransferases)
  filter(!grepl("ase", description))
  
  # filter GTF for protein coding genes
  proteins_gtf <- as.data.frame(.gtf) %>%
    filter(gene_biotype == "protein_coding", type == "exon") %>%
    # filter to keep only genes with CDS width divisible by 3
    # N.B. these should be filtered out anyway (but we perform this
    # operation in case some proteins have isoforms with odd CDS)
    mutate(is_divisible_by_3 = ifelse(width %% 3 == 0, "yes", "no")) %>%
    filter(is_divisible_by_3 == "yes")
  
  # keep only protein coding genes in TPMs
  proteins_tpm <- subset(.tpm, rownames(.tpm) %in% proteins_gtf$gene_id)

  # remove RP genes in TPMs
  proteins_tpm_non_rp <- proteins_tpm[!rownames(proteins_tpm) %in% rp_gff$gene_id, ]
  
  return(proteins_tpm_non_rp)
  
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

main <- function() {
  
  # load data
  txi <- readRDS(arguments$input)
  metadata <- read_delim(arguments$metadata, delim = "\t", col_names = TRUE)
  gtf <- rtracklayer::import(arguments$gtf)
  gff <- rtracklayer::import(arguments$gff)

  print(arguments)
  
  # initialise variables
  fdr <- as.numeric(arguments$fdr_threshold)
  
  # filter TPMs to keep only protein-coding genes
  tpm <- txi$abundance
  print(dim(tpm)) # debugging purposes only
  proteins_tpm <- extract_protein_coding_genes(gff, gtf, tpm)
  print(dim(proteins_tpm)) # debugging purposes only (should be 4141)
  
  # calculate highly and lowly expressed genes (all strains)
  meta_chisq <- chi_meta_analysis(proteins_tpm)
  print(dim(meta_chisq$data)) # debugging purposes only
  print(dim(meta_chisq$results)) # debugging purposes only

  meta_scores <- calculate_heg_leg_scores(meta_chisq$results, fdr)
  print(dim(meta_scores$heg)) # debugging purposes only
  print(dim(meta_scores$leg)) # debugging purposes only
  
  # rank highly and lowly expressed genes by chisq
  meta_scores$heg$meta.chisq_rank <- rank(-meta_scores$heg$chisq)
  meta_scores$leg$meta.chisq_rank <- rank(meta_scores$leg$chisq)
  
  # save results and scores to csv
  write_csv(as.data.frame(meta_chisq$data), file = "rank-tpm.csv")
  write_csv(as.data.frame(meta_chisq$results), file = "chisq-results.csv")
  write_csv(as.data.frame(meta_scores$heg), file = "heg-scores.csv")
  write_csv(as.data.frame(meta_scores$leg), file = "leg-scores.csv")
  saveRDS(proteins_tpm, "proteins_tpm.rds")
  
}

main()