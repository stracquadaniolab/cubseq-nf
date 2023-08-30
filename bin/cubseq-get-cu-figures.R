#!/usr/bin/env Rscript

"cubseq-get-figures.R

Usage:
    cubseq-get-figures.R [--cu-table=<cu_table>] [--aa-property=<aa_property>]
    
Options:
    --cu-table=<cu_table>             CSV file containing CU relative frequency, paired difference (etc) data.
    --aa-property=<aa_property>       CSV file containing amino acids and corresponding property value.
    -h --help                         Show this screen.
    --version                         Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-get-figures.R")

# load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(paletteer))
suppressMessages(library(ggrepel))

figure_cu_scatter <- function(data, x, y, aa_label, col, legend, xlab, ylab) {
  
  cu_scatter <- ggplot(data, aes(x=x, y=y, color=col)) + 
    geom_point() + geom_abline(intercept = 0, slope = 1) +
    theme_cowplot(12) +
    scale_color_paletteer_d("ggsci::nrc_npg") +
    labs(color=legend) +
    xlab(xlab) +
    ylab(ylab) +
    geom_text_repel(label = aa_label, size=2.5, show.legend = FALSE)
  
}

figure_aa_barplot <- function(data, x, y, fill, legend, xlab, ylab, title) {
  
  barplot_aa <- ggplot(data, aes(x=x, y=y, group=factor(x), fill=fill)) + # NB: here we group by "aa" column
    geom_bar(stat = "identity", position = "dodge") +
    labs(fill=legend) +
    xlab(xlab) +
    ylab(ylab) +
    theme_cowplot(12) +
    scale_fill_paletteer_d("ggsci::nrc_npg") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  
}

figure_codon_barplot <- function(data, x, y, fill, legend, xlab, ylab, title) {
  
  barplot_codon <- ggplot(data, aes(x=x, y=y, group=factor(aa), fill=fill)) + # NB: "aa" column is hard-coded here
    geom_bar(stat = "identity") +
    labs(fill=legend) +
    xlab(xlab) +
    ylab(ylab) +
    scale_y_discrete(drop = TRUE, expand = c(0, 0)) +
    facet_grid(aa~., scales = "free", space = "free_y", switch = "y") + # NB: "aa" column is hard-coded here
    theme(strip.placement = "outside",
          panel.spacing = unit(0, "in"),
          strip.background.y = element_rect(fill = "white", color = "gray75"),
          panel.background = element_rect(fill = "white"),
          panel.grid = element_line(color = "gray95")) +
    scale_fill_paletteer_d("ggsci::nrc_npg") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) # cannot use cowplot as it messes up x-label, so have to set title theme manually
  
}

figure_aa_boxplot <- function(data, x, y, fill, legend, xlab, ylab, title) {
  
  # NB: x-axis label colours are set manually following ggsci::nrc_npg paletteer colours for aa.property column
  # This was done to increase readability and interpretation of the graph - it was the only way (I think)
  aa.axis.cols <- c("#8391B4", "#4DBBD5", "#F39B7E", "#E64B35", "#E64B35", "#00A087", "#4DBBD5", "#3C5387", "#4DBBD5", "#3C5387", "#4DBBD5", "#4DBBD5", "#F39B7E", "#4DBBD5", "#F39B7E", "#3C5387", "#F39B7E", "#F39B7E", "#4DBBD5", "#00A087", "#00A087")
  
  boxplot_aa <- ggplot(data, aes(x=x, y=y, fill=fill)) +
    geom_boxplot(lwd=0.3) +
    labs(fill=legend) +
    xlab(xlab) +
    ylab(ylab) +
    theme_cowplot(12) +
    scale_fill_paletteer_d("ggsci::nrc_npg") +
    theme(axis.text.x = element_text(color = aa.axis.cols)) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))
  
}

# TODO: add log2 aa and codon bar plots

# TODO: add figure_cu_heatmap()

# TODO: add figure_aa_heatmap()

# TODO: add figure_strain_heatmap()

main <- function() {
  
  # load datasets
  cu <- read.table(arguments$cu_table, header = TRUE)
  aa.property <- read.table(arguments$aa_property, header = TRUE) # to colour amino acids by property (user-specified)
  cu <- merge(cu, aa.property, by.x = "aa", by.y = "aa")
  
  # figure 1
  plot1A <- figure_cu_scatter(cu, cu$heg.rf, cu$kazusa.rf, cu$aa, cu$aa.property, "AA property", "cubseq-HEG", "Kazusa E. coli K-12")
  plot1B <- figure_cu_scatter(cu, cu$leg.rf, cu$kazusa.rf, cu$aa, cu$aa.property, "AA property", "cubseq-LEG", "Kazusa E. coli K-12")
  plot1C <- figure_cu_scatter(cu, cu$leg.rf, cu$heg.rf, cu$aa, cu$aa.property, "AA property", "cubseq-LEG", "cubseq-HEG")
  plot1 <- plot_grid(plot1A + theme(legend.position="none"), 
                     plot1B + theme(legend.position="none"),
                     plot1C + theme(legend.position="none"),
                     align = 'v',
                     labels = c('A', 'B', 'C'), 
                     hjust = -1,
                     ncol = 1,
                     label_size = 12)
                     
  legend <- get_legend(plot1A + theme(legend.box.margin = margin(0, 0, 0, 12)))
  plot1 <- plot_grid(plot1, legend, rel_widths = c(3, .4))

  save_plot(filename = "figure1.pdf", plot = plot1, base_height = 9, base_width = 9)
  message("Figure 1 done.")

  # figure 2
  plot2A <- figure_aa_barplot(cu, cu$aa, cu$sq.diff.HEG_kazusa, cu$aa.property, "AA property", "", "Squared difference in relative frequency", "Kazusa vs. cubseq-HEG")
  plot2B <- figure_aa_barplot(cu, cu$aa, cu$sq.diff.LEG_kazusa, cu$aa.property, "AA property", "Amino Acids", "", "Kazusa vs. cubseq-LEG")
  plot2C <- figure_aa_barplot(cu, cu$aa, cu$sq.diff.HEG_LEG, cu$aa.property, "AA property", "", "", "cubseq-HEG vs. cubseq-LEG")
  plot2 <- plot_grid(plot2A + theme(legend.position="none"), 
                     plot2B + theme(legend.position="none"),
                     plot2C + theme(legend.position="none"),
                     align = 'vh',
                     labels = c('A', 'B', 'C'), 
                     hjust = -1,
                     nrow = 1,
                     label_size = 12)
  
  legend <- get_legend(plot2A + theme(legend.box.margin = margin(0, 0, 0, 10)))
  plot2 <- plot_grid(plot2, legend, rel_widths = c(3, .4))

  save_plot(filename = "figure2.pdf", plot = plot2, base_width = 12)
  message("Figure 2 done.")

  # figure 3
  plot3A <- figure_codon_barplot(cu, cu$sq.diff.HEG_kazusa, cu$codon, cu$aa.property, "AA property", "", "Amino Acid (with corresponding codon)", "Kazusa vs. cubseq-HEG")
  plot3B <- figure_codon_barplot(cu, cu$sq.diff.LEG_kazusa, cu$codon, cu$aa.property, "AA property", "Squared difference in relative frequency", "", "Kazusa vs. cubseq-LEG")
  plot3C <- figure_codon_barplot(cu, cu$sq.diff.HEG_LEG, cu$codon, cu$aa.property, "AA property", "", "", "cubseq-HEG vs. cubseq-LEG")
  plot3 <- plot_grid(plot3A + theme(legend.position="none"), 
                     plot3B + theme(legend.position="none"),
                     plot3C + theme(legend.position="none"),
                     align = 'vh',
                     labels = c('A', 'B', 'C'), 
                     hjust = -1,
                     nrow = 1,
                     label_size = 12)
  
  legend <- get_legend(plot3A + theme(legend.box.margin = margin(0, 0, 0, 10)))
  plot3 <- plot_grid(plot3, legend, rel_widths = c(3, .4))

  save_plot(filename = "figure3.pdf", plot = plot3, base_height = 9, base_width = 12)
  message("Figure 3 done.")

  # figure 4
  plot4A <- figure_aa_boxplot(cu, cu$aa, cu$sq.diff.HEG_kazusa, cu$aa.property, "AA Property", "", "Squared difference in relative frequency", "Kazusa vs.cubseq-HEG")
  plot4B <- figure_aa_boxplot(cu, cu$aa, cu$sq.diff.LEG_kazusa, cu$aa.property, "AA Property", "Amino Acids", "", "Kazusa vs.cubseq-LEG")
  plot4C <- figure_aa_boxplot(cu, cu$aa, cu$sq.diff.HEG_LEG, cu$aa.property, "AA Property", "", "", "cubseq-HEG vs. cubseq-LEG")
  plot4 <- plot_grid(plot4A + theme(legend.position="none"), 
                     plot4B + theme(legend.position="none"),
                     plot4C + theme(legend.position="none"),
                     align = 'vh',
                     labels = c('A', 'B', 'C'), 
                     hjust = -1,
                     nrow = 1,
                     label_size = 12)
  
  legend <- get_legend(plot4A + theme(legend.box.margin = margin(0, 0, 0, 10)))
  plot4 <- plot_grid(plot4, legend, rel_widths = c(3, .4))
  save_plot(filename = "figure4.pdf", plot = plot4, base_width = 12)
  message("Figure 4 done.")

}

main()
