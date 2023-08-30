#!/usr/bin/env Rscript

"cubseq-pca.R

Usage:
    cubseq-pca.R <input> [--metadata=<metadata>]
    
Options:
    --metadata=<metadata>             CSV file containing sample metadata.
    -h --help                         Show this screen.
    --version                         Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-pca.R")

# load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(paletteer))

figure_tpm_pca <- function(pca_data, scores_data, colour_by_col, legend, title) {
  
  # add sample number ("n=42") to each legend variable
  pca_fig <- ggplot(data = pca_data) + 
    geom_point(aes(x = PC1, y = PC2, col = colour_by_col)) + 
    labs(col = legend) +
    xlab(paste0("PC1: ", sprintf("%0.1f%%", scores_data[1,2]))) +
    ylab(paste0("PC2: ", sprintf("%0.1f%%", scores_data[2,2]))) +
    theme_cowplot(12) +
    ggtitle(paste0(title, " (n=", nrow(pca_data), ")")) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.line = element_blank(),
          panel.border = element_rect(colour = "#737373", fill=NA, size=0.3))
  
  return(pca_fig)
  
}

main <- function() {
  
  # load datasets
  txi <- readRDS(arguments$input)
  metadata <- read_delim(arguments$metadata, delim = "\t", col_names = TRUE)
  print(arguments)
  
  # preprocess data
  tpm <- t(txi$abundance) # transpose to see sample variance (i.e. set samples as rows)
  tpm <- tpm[ , which(apply(tpm, 2, var) != 0)] # remove zero-variance genes

  # TODO: run on real dataset and see whether to use X% high variance genes
  
  # perform pca
  pca <- prcomp(tpm, center = TRUE, scale. = TRUE)
  pca.df <- data.frame('Samples' = rownames(tpm), pca$x[,1:2])
  
  # get relative importance of each PC
  pca.summ = summary(pca)
  pca.exp_var = pca.summ$importance[2,] * 100
  pca.exp_var <- data.frame(PC = c(names(pca.exp_var)),
                            perc.variance = c(as.numeric(paste(pca.exp_var))))
  
  # map variables of interest to pca.df by run_accession
  vars <- metadata %>%
    select('tax_id', 'strain', 'study_accession', 'sample_accession', 
           'run_accession', 'library_selection', 'instrument_model', 
           'read_count', 'base_count')
  pca.df <- merge(pca.df, vars, by.x = "Samples", by.y = "run_accession")
  
  # get figures
  plot1 <- ggplot(data=pca.exp_var, aes(x=reorder(PC, -perc.variance), y=perc.variance)) +
    geom_bar(stat="identity", width=0.9) +
    xlab("PC") +
    ylab("Relative Importance (%)") +
    ggtitle("Relative importance of each PC") +
    geom_text(aes(label = sprintf("%0.1f%%", perc.variance)), 
              position=position_dodge(width=2), vjust=-0.5) +
    theme_cowplot(12) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.line = element_blank(),
          panel.border = element_rect(colour = "#737373",
          fill=NA,
          size=0.3))         

  plot2 <- figure_tpm_pca(pca.df, pca.exp_var, as.factor(pca.df$tax_id), "Taxonomy ID", "Taxonomy")
  plot3 <- figure_tpm_pca(pca.df, pca.exp_var, pca.df$strain, "Strain", "Strain")
  plot4 <- figure_tpm_pca(pca.df, pca.exp_var, pca.df$study_accession, "Study Accession", "Study")
  plot5 <- figure_tpm_pca(pca.df, pca.exp_var, pca.df$library_selection, "Library Selection", "Library Selection")
  plot6 <- figure_tpm_pca(pca.df, pca.exp_var, pca.df$instrument_model, "Instrument Model", "Instrument Model")
  plot7 <- figure_tpm_pca(pca.df, pca.exp_var, log2(pca.df$read_count), 
                           (expression(log[2]*"(Library Size)")), "Library Size")
  plot8 <- figure_tpm_pca(pca.df, pca.exp_var, log2(pca.df$base_count), (expression(log[2]*"(Base Count)")), "Base Count")

  # plot as grid
  plot_pca <- plot_grid(plot2, plot3, plot5, plot6, plot7, plot8,
                        align = 'vh',
                        labels = c('B', 'C', 'D', 'E', 'F', 'G'),
                        hjust = -1,
                        nrow = 3,
                        ncol = 2,
                        label_size = 12)
  
  save_plot(filename = "pca.pdf", plot = plot_pca, base_width = 12, base_height = 9)
  message("Figure PCA done.")
  save_plot(filename = "pca-study.pdf", plot = plot4, base_width = 12, base_height = 9)
  message("Figure PCA-study done.")
  save_plot(filename = "pca-pc-barplot.pdf", plot = plot1, base_width = 12, base_height = 9)
  message("Figure PC barplot done.")

  # TODO: PCA faceted by strain/tax_id (other species greyed out)
  
}

main()
