#!/usr/bin/env Rscript

# load required packages
suppressMessages(library(vcfR))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(paletteer))
library(circlize)

# load data
vcf <- read.vcfR("merged.vcf", verbose = FALSE)
vcf_df <- cbind(as.data.frame(getFIX(vcf)), INFO2df(vcf))

plot_variants <- vcf_df %>%
  select(., POS, AC) %>%
  separate_rows(AC, sep = ",") %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  mutate(., AC = as.numeric(AC)) %>%      # AC is total number of alt alleles in called genotypes.        
  mutate(., AC_norm = AC / 6763) %>%    # divide by total number of sequencing runs.
  ggplot(., aes(pos_kb, AC_norm)) +
  geom_bar(stat = "identity", colour = "#8491B4FF") +
  labs(x = "Position (kb)", 
       y = "Mutation frequency") +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "#737373",
                                    fill=NA,
                                    linewidth=0.3))

plot_variants

# plot as circos plot
circos.genomicInitialize(plot_variants)



plot_variants <- vcf_df %>%
  select(., POS, AC) %>%
  separate_rows(AC, sep = ",") %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  mutate(., AC = as.numeric(AC)) %>%    # AC is total number of alt alleles in called genotypes.        
  mutate(., AC_norm = AC / 6763) %>%    # divide by total number of sequencing runs.
  arrange(pos_kb) %>%                   # arrange by position for circular plot
  select(., pos_kb, AC_norm)

# TODO: need to fix
circos.par("start.degree" = 90)    # start at the top of the circle
circlize(plot_variants$pos_kb, plot_variants$AC_norm, track.height = 0.4)

%>%    # create circular plot with track height of 0.4
  circos.trackPlotRegion(ylim = c(0, max(.$AC_norm)), track.height = 0.4) %>%
  circos.axis(h = "bottom") %>%          # add axis labels at bottom of circle
  circos.text(sect = 1:length(unique(.$pos_kb)), labels = unique(.$pos_kb), 
              facing = "clockwise", niceFacing = TRUE) %>%  # add tick marks at position labels
  circos.barplot(.$AC_norm, 
                 col = "#8491B4FF", 
                 border = NA, 
                 track.height = 0.4, 
                 add = TRUE) %>%          # add bars to circular plot
  circos.clear()
