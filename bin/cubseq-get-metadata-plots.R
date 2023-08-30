#!/usr/bin/env Rscript

"cubseq-get-metadata-plots.R

Usage:
    cubseq-get-metadata-plots.R <input>
    
Options:
    -h --help                         Show this screen.
    --version                         Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-get-metadata-plots.R")

# load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(paletteer))
suppressMessages(library(ggrepel))

# load metadata
metadata <- read_delim(arguments$input, delim = "\t", col_names = TRUE)

# plot <- function(.data, .variable, .row_start, .row_end, .x_lab) {
  
#   # count samples
#   data_count <- .data %>%
#     count(.variable, sort = TRUE)
  
#   # plot data
#   var_plot <- data_count[as.numeric(.row_start),as.numeric(.row_end)] %>%
#     ggplot(., aes(x = reorder(data_count[,1], data_count[,2]))) +
#     geom_bar(aes(y = data_count[,2]), stat = "identity", position = "dodge") +
#     geom_text(aes(x = data_count[,1], 
#                   y = data_count[,2], 
#                   label = data_count[,2]), 
#               hjust = -0.25) +
#     labs(x = .x_lab, y = "Total sequencing runs") +
#     coord_flip() +
#     theme_cowplot(12)
  
# }

# bar plot of total samples per strain
top_strains <- metadata %>% 
  count(strain, sort = TRUE)

top_strains_plot <- top_strains[1:25,] %>%
  ggplot(., aes(y = reorder(strain, n))) +
  geom_bar(aes(x = n), stat = "identity", position = "dodge") +
  geom_text(aes(x = n, y = strain, label = n), hjust = -0.25) +
  labs(x = "Total sequencing runs", y = "Strain") +
  theme_cowplot(12)

# bar plot of total samples per taxonomy ID
taxid_count <- metadata %>% 
  count(as.character(tax_id), sort = TRUE)
colnames(taxid_count) <- c("tax_id", "n")

top_taxid_plot <- taxid_count[1:25,] %>%
  ggplot(., aes(y = reorder(tax_id, n))) +
  geom_bar(aes(x = n), stat = "identity", position = "dodge") +
  geom_text(aes(x = n, y = tax_id, label = n), hjust = -0.25) +
  labs(x = "Total sequencing runs", y = "Taxonomy ID") +
  theme_cowplot(12)

# bar plot of total samples per year
year_count <- metadata %>% 
  count(format(as.Date(first_public, format="%Y/%m/%d"),"%Y"))
colnames(year_count) <- c("Year", "Count")

year_plot <- year_count %>%
  ggplot(., aes(x = Year)) +
  geom_bar(aes(y=Count), stat = "identity", position = "dodge") +
  geom_text(aes(x=Year, y=Count, label=Count), hjust=-0.25) +
  labs(x = "Year", y = "Total sequencing runs") +
  coord_flip() +
  theme_cowplot(12)

# bar plot of total samples per instrument model
instrument_model_plot <- metadata %>%
  count(instrument_model, sort = TRUE) %>%
  ggplot(., aes(y = reorder(instrument_model, n))) +
  geom_bar(aes(x=n), stat = "identity", position = "dodge") +
  geom_text(aes(x = n, y=instrument_model, label=n), hjust=-0.25) +
  labs(x = "Total sequencing runs", y = "Instrument Model") +
  theme_cowplot(12)

# bar plot of total samples per library selection method
library_selection_plot <- metadata %>%
  count(library_selection, sort = TRUE) %>%
  ggplot(., aes(y = reorder(library_selection, n))) +
  geom_bar(aes(x = n), stat = "identity", position = "dodge") +
  geom_text(aes(x = n, y = library_selection, label = n), hjust = -0.25) +
  labs(x = "Total sequencing runs", y = "Library Selection") +
  theme_cowplot(12)

# histogram and box plot of read count by sequencing run and taxonomy ID, respectively
read_count_plot <- metadata %>%
  mutate_at(vars(read_count), funs(./ 1000000)) %>%
  ggplot(., aes(x = read_count)) +
  geom_bar(binwidth = 2, stat = "bin") +
  labs(x = bquote("Read Count ("~10^6~")"), y = "Total sequencing runs") +
  theme_cowplot(12)

taxid_read_count_plot <- metadata %>%
  mutate_at(vars(read_count), funs(./ 1000000)) %>%
  ggplot(., aes(x = read_count, 
                y = reorder(as.character(tax_id), read_count),
                group=factor(tax_id))) +
  geom_boxplot(lwd=0.3) +
  labs(x = bquote("Read Count ("~10^6~")"), y = "Taxonomy ID") +
  theme_cowplot(12)

# histogram and box plot of base count by sequencing run and taxonomy ID, respectively
base_count_plot <- metadata %>%
  mutate_at(vars(base_count), funs(./ 1000000000)) %>%
  ggplot(., aes(x = base_count)) +
  geom_bar(binwidth = 0.5, stat = "bin") +
  labs(x = bquote("Base Count ("~10^9~")"), y = "Total sequencing runs") +
  theme_cowplot(12)

taxid_base_count_plot <- metadata %>%
  mutate_at(vars(base_count), funs(./ 1000000000)) %>%
  ggplot(., aes(x = base_count, 
                y = reorder(as.character(tax_id), base_count),
                group=factor(tax_id))) +
  geom_boxplot(lwd=0.3) +
  labs(x = bquote("Base Count ("~10^9~")"), y = "Taxonomy ID") +
  theme_cowplot(12)

# histogram and box plot of sra bytes by sequencing run and taxonomy ID, respectively
sra_bytes_plot <- metadata %>%
  mutate_at(vars(sra_bytes), funs(./ 1000000000)) %>%
  ggplot(., aes(x = sra_bytes)) +
  geom_bar(binwidth = 0.2, stat = "bin") +
  labs(x = "Read file size (GB)", y = "Total sequencing runs") +
  theme_cowplot(12)

taxid_sra_bytes_plot <- metadata %>%
  mutate_at(vars(sra_bytes), funs(./ 1000000000)) %>%
  ggplot(., aes(x = sra_bytes, 
                y = reorder(as.character(tax_id), sra_bytes),
                group=factor(tax_id))) +
  geom_boxplot(lwd=0.3) +
  labs(x = "Read file size (GB)", y = "Taxonomy ID") +
  theme_cowplot(12)

# plot as grid
plot1 <- plot_grid(top_strains_plot,
                   top_taxid_plot, 
                   year_plot, 
                   machines_count_plot, 
                   lib_count_plot,
                   align = 'vh',
                   labels = c('A', 'B', 'C', 'D', 'E'), 
                   hjust = -1, 
                   nrow = 3, 
                   ncol = 2, 
                   label_size = 12)

plot2 <- plot_grid(read_count_plot,
                   base_count_plot, 
                   sra_bytes_plot,
                   align = 'h',
                   labels = c('A', 'B', 'C'),
                   hjust = -1,
                   nrow = 1,
                   ncol = 3,
                   label_size = 12)

plot3 <- plot_grid(taxid_read_count_plot,
                   taxid_base_count_plot, 
                   taxid_sra_bytes_plot,
                   align = 'h',
                   labels = c('A', 'B', 'C'),
                   hjust = -1,
                   nrow = 1,
                   ncol = 3,
                   label_size = 12)

save_plot(filename = "descriptive-panel.pdf", plot = plot1, base_width = 23, base_height = 15)
message("Plot1 done.")

save_plot(filename = "read-base-file-panel.pdf", plot = plot2, base_width = 15, base_height = 7)
message("Plot2 done.")

save_plot(filename = "taxid-read-base-file-panel.pdf", plot = plot3, base_width = 15, base_height = 15)
message("Plot3 done.")