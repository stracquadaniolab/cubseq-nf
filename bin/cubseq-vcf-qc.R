#!/usr/bin/env Rscript

"cubseq-vcf-qc.R

Usage:
    cubseq-vcf-qc.R <input>
    
Options:
    -h --help                               Show this screen.
    --version                               Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-vcf-qc.R")
print(arguments)

# load required packages
suppressMessages(library(vcfR))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(paletteer))

# # load FASTA file
# message("Loading fasta file.")
# dna <- ape::read.dna(arguments$fasta, format = "fasta")

# # load GFF file
# message("Loading gff file.")
# gff <- read.table(arguments$gff, sep = "\t", quote = "")

# load VCF file
message("Loading vcf file.")
vcf <- read.vcfR(arguments$input, verbose = FALSE)
message("Total samples:")
dim(vcf@gt)
message("Total variants recorded:")
dim(vcf@fix)

# # create chromR objects
# message("Creating chromR objects.")
# chrom <- create.chromR(name='', vcf=vcf, seq=dna, ann=gff)
# chrom <- proc.chromR(chrom)

# # plot
# message("Plotting chrom object.")
# pdf(file = "vcfr-chrom-plot.pdf")
# plot(chrom)
# dev.off()

# # visualise chromo data
# message("Plotting chromo QC.")
# pdf(file = "vcfr-chromoqc-plot.pdf")
# chromoqc(chrom, dp.alpha = 22)
# dev.off()

# convert to dataframe
vcf_df <- cbind(as.data.frame(getFIX(vcf)), INFO2df(vcf))                                 

# plot alternative allele observation (AO)
message("Plotting alternative allele observation (AO) - raw values.")
ao <- vcf_df %>%
  select(., POS, AO) %>%
  separate_rows(AO)

# th_ao_95 <- quantile(as.numeric(ao$AO), prob=0.95)

plot_ao <- vcf_df %>%
  select(., POS, AO) %>%
  separate_rows(AO) %>%
  # mutate(NAO = pmin(as.numeric(AO), th_ao_95)) %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  ggplot(., aes(pos_kb, as.numeric(AO))) +
  geom_point(alpha = 0.2, colour = "#982A86", size = 0.4) +
  # labs(x = "Position (kb)", y = "Alternate Allele Observation Count") +
  labs(x = "Position (kb)", y = "AO") +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "#737373",
                                    fill=NA,
                                    size=0.2))
message("Plot_ao done.")

th_ao_99 <- quantile(as.numeric(ao$AO), prob=0.99)
plot_ao_99 <- vcf_df %>%
  select(., POS, AO) %>%
  separate_rows(AO) %>%
  mutate(NAO = pmin(as.numeric(AO), th_ao_99)) %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  ggplot(., aes(pos_kb, as.numeric(NAO))) +
  geom_point(alpha = 0.2, colour = "#982A86", size = 0.4) +
  # labs(x = "Position (kb)", y = "Alternate Allele Observation Count") +
  labs(x = "Position (kb)", y = "AO") +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "#737373",
                                    fill=NA,
                                    size=0.2))
message("Plot_ao_99 done.")

th_ao_999 <- quantile(as.numeric(ao$AO), prob=0.999)
plot_ao_999 <- vcf_df %>%
  select(., POS, AO) %>%
  separate_rows(AO) %>%
  mutate(NAO = pmin(as.numeric(AO), th_ao_999)) %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  ggplot(., aes(pos_kb, as.numeric(NAO))) +
  geom_point(alpha = 0.2, colour = "#982A86", size = 0.4) +
  # labs(x = "Position (kb)", y = "Alternate Allele Observation Count") +
  labs(x = "Position (kb)", y = "AO") +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "#737373",
                                    fill=NA,
                                    size=0.2))
message("Plot_ao_999 done.")

th_ao_95 <- quantile(as.numeric(ao$AO), prob=0.999)
plot_ao_95_log10 <- vcf_df %>%
  select(., POS, AO) %>%
  separate_rows(AO) %>%
  mutate(NAO = pmin(as.numeric(AO), th_ao_95)) %>%
  mutate(log_NAO = log10(NAO)) %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  ggplot(., aes(pos_kb, as.numeric(log_NAO))) +
  geom_point(alpha = 0.2, colour = "#982A86", size = 0.4) +
  # labs(x = "Position (kb)", y = "Alternate Allele Observation Count") +
  labs(x = "Position (kb)", y = bquote('AO ('~log[10]*')')) +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "#737373",
                                    fill=NA,
                                    size=0.2))
message("Plot_ao_95_log10_ao_999 done.")


# plot read depth along position
message("Plotting read depth.")
plot_dp <- vcf_df %>%
  mutate(., DP.1e6 = DP / 1e6) %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  ggplot(., aes(pos_kb, DP.1e6)) +
  geom_point(alpha = 0.2, colour = "#001889FF", size = 0.4) +
  # labs(x = "Position (kb)", y = "Read Depth (1e6)") +
  labs(x = "Position (kb)", y = "DP (1e6)") +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "#737373",
                                    fill=NA,
                                    size=0.3))
message("Plot DP done.")

message("Plotting read depth (log10-scaled.)")
plot_dp_log10 <- vcf_df %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  mutate(., log_DP = log10(DP)) %>%
  ggplot(., aes(pos_kb, log_DP)) +
  geom_point(alpha = 0.2, colour = "#001889FF", size = 0.4) +
  # labs(x = "Position (kb)", y = "Read Depth (1e6)") +
  labs(x = "Position (kb)", y = bquote('DP ('~log[10]*')')) +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "#737373",
                                    fill=NA,
                                    size=0.3))
message("Plot dp_log10 done.")

# plot mapping quality (MQ)
message("Plotting mapping quality (MQM).")
plot_mqm <- vcf_df %>%
  select(., POS, MQM) %>%
  separate_rows(MQM, sep = ",") %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  ggplot(., aes(pos_kb, as.numeric(MQM))) +
  geom_point(alpha = 0.2, colour = "#DC0000FF", size = 0.4) +
  # labs(x = "Position (kb)", 
       # y = "Mean mapping quality of observed alternate alleles") +
  labs(x = "Position (kb)", y = "MQM") +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "#737373",
                                    fill=NA,
                                    size=0.3))
message("Done.")

# plot Phred scale quality (QUAL)
message("Plotting Phred scale quality (QUAL).")
plot_qual <- vcf_df %>%
  select(., POS, QUAL) %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  mutate(., QUAL.1e6 = as.numeric(QUAL) / 1e6) %>%
  ggplot(., aes(pos_kb, QUAL.1e6)) +
  geom_point(alpha = 0.2, colour = "#00A087FF", size = 0.4) +
  # labs(x = "Position (kb)", 
  #      y = "Phred-scaled quality (1e6)") +
  labs(x = "Position (kb)", y = "QUAL (1e6)") +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "#737373",
                                    fill=NA,
                                    size=0.3))
message("Done.")

message("Plotting Phred scale quality (QUAL), log10-scaled.")
plot_qual_log10 <- vcf_df %>%
  select(., POS, QUAL) %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  mutate(., log_QUAL = log10(as.numeric(QUAL))) %>%
  ggplot(., aes(pos_kb, log_QUAL)) +
  geom_point(alpha = 0.2, colour = "#00A087FF", size = 0.4) +
  # labs(x = "Position (kb)", 
  #      y = "Phred-scaled quality (1e6)") +
  labs(x = "Position (kb)", y = bquote('QUAL ('~log[10]*')')) +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "#737373",
                                    fill=NA,
                                    size=0.3))
message("Done.")

# plot variants per site
message("Plotting variants per site.")
plot_variants <- vcf_df %>%
  select(., POS, AC) %>%
  separate_rows(AC, sep = ",") %>%
  mutate(., pos_kb = as.numeric(POS) / 1000) %>%
  mutate(., AC = as.numeric(AC)) %>%                            # AC is total number of alt alleles in called genotypes.
  mutate(., AC_norm = AC / 6763) %>%                            # divide by total number of sequencing runs.
  ggplot(., aes(pos_kb, AC_norm)) +
  geom_bar(stat = "identity", colour = "#8491B4FF") +
  labs(x = "Position (kb)", 
       y = "Mutation frequency") +
  theme_cowplot(12) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "#737373",
                                    fill=NA,
                                    size=0.3))
message("Done.")

# # plot nucleotide content
# message("Plotting nucleotide content.")

# plot grid
message("Plotting grid with AO raw values.")
plot1 <- plot_grid(plot_ao,
                   plot_dp,
                   plot_mqm,
                   plot_qual,
                   align = 'vh',
                   labels = c('A', 'B', 'C', 'D'),
                   ncol = 2,
                   nrow = 2,
                   label_size = 12)

bottom_plot <- plot(plot_variants, labels = c('E'), label_size = 12)

vcfqc_plot <- plot_grid(plot1, 
          bottom_plot, 
          labels = c('', 'E'), 
          label_size = 12, 
          ncol = 1,
          rel_heights = c(3, 2))

save_plot(filename = "vcf-qc.pdf", plot = vcfqc_plot, base_height = 6, base_width = 9)
message("Finished!")


message("Plotting grid with AO 99% threshold + log10.")
plot2 <- plot_grid(plot_ao_99,
                   plot_dp_log10,
                   plot_mqm,
                   plot_qual_log10,
                   align = 'vh',
                   labels = c('A', 'B', 'C', 'D'),
                   ncol = 2,
                   nrow = 2,
                   label_size = 12)

vcfqc_plot2 <- plot_grid(plot2, 
          bottom_plot, 
          labels = c('', 'E'), 
          label_size = 12, 
          ncol = 1,
          rel_heights = c(3, 2))

save_plot(filename = "vcf-qc.pdf", plot = vcfqc_plot, base_height = 6, base_width = 9)
message("Finished!")


message("Plotting grid with AO 99.9% threshold + log10.")
plot3 <- plot_grid(plot_ao_999,
                   plot_dp_log10,
                   plot_mqm,
                   plot_qual_log10,
                   align = 'vh',
                   labels = c('A', 'B', 'C', 'D'),
                   ncol = 2,
                   nrow = 2,
                   label_size = 12)

vcfqc_plot3 <- plot_grid(plot3, 
          bottom_plot, 
          labels = c('', 'E'), 
          label_size = 12, 
          ncol = 1,
          rel_heights = c(3, 2))


message("Plotting grid with AO log10-scaled 95% threshold + log10.")
plot4 <- plot_grid(plot_ao_95_log10,
                   plot_dp_log10,
                   plot_mqm,
                   plot_qual_log10,
                   align = 'vh',
                   labels = c('A', 'B', 'C', 'D'),
                   ncol = 2,
                   nrow = 2,
                   label_size = 12)

vcfqc_plot4 <- plot_grid(plot4, 
          bottom_plot, 
          labels = c('', 'E'), 
          label_size = 12, 
          ncol = 1,
          rel_heights = c(3, 2))

save_plot(filename = "vcf-qc-raw.pdf", plot = vcfqc_plot, base_height = 6, base_width = 9)
save_plot(filename = "vcf-qc-99_log10.pdf", plot = vcfqc_plot2, base_height = 6, base_width = 9)
save_plot(filename = "vcf-qc-999_log10.pdf", plot = vcfqc_plot3, base_height = 6, base_width = 9)
save_plot(filename = "vcf-qc-log10-95_log10.pdf", plot = vcfqc_plot4, base_height = 6, base_width = 9)
message("Finished!")
