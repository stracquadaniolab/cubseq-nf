# load packages
library(tidyverse)
library(ggplot2)

# load datasets
cu <- read.table("scripts/datasets/kazusa-cubseq-CU.csv", header = TRUE)

# GET SCATTER PLOTS
# plot of cubseq-HE vs. Kazusa
ggplot(cu, aes(x=cubseq.HE, y=Kazusa, color=AA)) + 
  geom_point() + geom_abline(yintercept = 0)

# plot of cubseq-LE vs. Kazusa
ggplot(cu, aes(x=cubseq.LE, y=Kazusa, color=AA)) + 
  geom_point() + geom_abline(yintercept = 0)

# plot of cubseq-HE vs. cubseq-LE
ggplot(cu, aes(x=cubseq.LE, y=cubseq.HE, color=AA)) + 
  geom_point() + geom_abline(yintercept = 0)


# CALCULATE PAIRWISE DIFFERENCE
# calculate pairwise difference (absolute)
cu$abs.diff.HE <- abs(cu$Kazusa - cu$cubseq.HE)
cu$abs.diff.LE <- abs(cu$Kazusa - cu$cubseq.LE)
cu$abs.diff.HE_LE <- abs(cu$cubseq.HE - cu$cubseq.LE)

# calculate pairwise difference (squared)
cu$sq.diff.HE <- (cu$Kazusa - cu$cubseq.HE)^2
cu$sq.diff.LE <- (cu$Kazusa - cu$cubseq.LE)^2
cu$sq.diff.HE_LE <- (cu$cubseq.HE - cu$cubseq.LE)^2

# group data by AA
df <- cu %>% group_by(AA, abs.diff.HE)
df <- arrange(df, AA, abs.diff.HE)

# plot bar plot (grouped by AA)
ggplot(df, aes(x=codon, y=abs.diff.HE, group=factor(AA), fill=AA)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("Codons") +
  ylab("Kazusa vs. cubseq-HE (absolute difference)")

ggplot(df, aes(x=codon, y=abs.diff.LE, group=factor(AA), fill=AA)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("Codons") +
  ylab("Kazusa vs. cubseq-LE (absolute difference)")


# plot bar plot (grouped by AA) - absolute values
#TODO: plot all three as one bar plot, colour by x_lab type.
ggplot(df, aes(x=abs.diff.HE, y=AA, group=factor(AA), fill=AA)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("Kazusa vs. cubseq-HE (absolute difference)") +
  ylab("Amino acids")

ggplot(df, aes(x=abs.diff.LE, y=AA, group=factor(AA), fill=AA)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("Kazusa vs. cubseq-LE (absolute difference)") +
  ylab("Amino acids")

ggplot(cu, aes(x=abs.diff.HE_LE, y=AA, group=factor(AA), fill=AA)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("cubseq-HE vs. cubseq-LE (absolute difference)") +
  ylab("Amino acids")

# plot bar plot (grouped by AA) -- squared difference:
#TODO: plot all three as one bar plot, colour by x_lab type.
ggplot(df, aes(x=sq.diff.HE, y=AA, group=factor(AA), fill=AA)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("Kazusa vs. cubseq-HE (squared difference)") +
  ylab("Amino acids")

ggplot(df, aes(x=sq.diff.LE, y=AA, group=factor(AA), fill=AA)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("Kazusa vs. cubseq-LE (squared difference)") +
  ylab("Amino acids")

ggplot(cu, aes(x=sq.diff.HE_LE, y=AA, group=factor(AA), fill=AA)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("cubseq-HE vs. cubseq-LE (squared difference)") +
  ylab("Amino acids")

# plot bar plot (grouped by codon)
ggplot(df, aes(x=codon, y=abs.diff.LE, group=factor(AA), fill=AA)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("Codons") +
  ylab("Kazusa vs. cubseq-LE (absolute difference)")

ggplot(df, aes(x=codon, y=sq.diff.HE, group=factor(AA), fill=AA)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("Codons") +
  ylab("Kazusa vs. cubseq-HE (squared difference)")

ggplot(cu, aes(x=codon, y=sq.diff.HE_LE, group=factor(AA), fill=AA)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  xlab("Codons") +
  ylab("cubseq-HE vs. cubseq-LE (squared difference)")

# boxplot per AA - absolute values
#TODO: plot all three as one bar plot, colour by x_lab type.
ggplot(df, aes(x=AA, y=abs.diff.HE, fill=AA)) + 
  geom_boxplot() +
  xlab("Amino acids") +
  ylab("Absolute difference (Kazusa vs. cubseq-HE)")

ggplot(df, aes(x=AA, y=abs.diff.LE, fill=AA)) + 
  geom_boxplot() +
  xlab("Amino acids") +
  ylab("Absolute difference (Kazusa vs. cubseq-LE)")

  ggplot(cu, aes(x=AA, y=abs.diff.HE_LE, fill=AA)) + 
  geom_boxplot() +
  xlab("Amino acids") +
  ylab("Absolute difference (cubseq-HE vs.cubseq-LE)")

# boxplot per AA - squared difference values
#TODO: plot all three as one bar plot, colour by x_lab type.
ggplot(df, aes(x=AA, y=sq.diff.HE, fill=AA)) + 
  geom_boxplot() +
  xlab("Amino acids") +
  ylab("Squared difference (Kazusa vs.cubseq-HE)")

ggplot(df, aes(x=AA, y=sq.diff.LE, fill=AA)) + 
  geom_boxplot() +
  xlab("Amino acids") +
  ylab("Squared difference (Kazusa vs.cubseq-LE)")

ggplot(cu, aes(x=AA, y=sq.diff.HE_LE, fill=AA)) + 
  geom_boxplot() +
  xlab("Amino acids") +
  ylab("Squared difference (cubseq-HE vs.cubseq-LE)")

# plot heatmap
cu.mat <- as.matrix(cu[3:5])
rownames(cu.mat) <- cu$codon
heatmap(cu.mat, cexCol = 1.2)
