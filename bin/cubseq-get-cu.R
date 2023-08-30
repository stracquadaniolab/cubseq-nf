#!/usr/bin/env Rscript

"cubseq-get-cu.R

Usage:
    cubseq-get-cu.R [--heg-fasta-dir=<heg_fasta_dir>] [--leg-fasta-dir=<leg_fasta_dir>] [--protein-fasta-dir=<protein_fasta_dir>] [--tax-id=<tax_id>]
    
Options:
    --heg-fasta-dir=<heg_fasta_dir>             Directory containing fasta files for high expression genes [default: .].
    --leg-fasta-dir=<leg_fasta_dir>             Directory containing fasta files for low expression genes [default: .].
    --protein-fasta-dir=<protein_fasta_dir>     Directory containing fasta files for protein-coding genes [default: .].
    --tax-id=<tax_id>                           Taxonomy ID to filter Kazusa and CoCoPUTs databases [default: 83333].
    -h --help                                   Show this screen.
    --version                                   Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-get-cu.R")

# load packages
suppressMessages(library(Biostrings))
suppressMessages(library(tidyverse))

get_kazusa_codons <- function(tax_id) {
    
    # initialise variables
    query.spsum <- URLencode("ftp://ftp.ebi.ac.uk/pub/databases//cutg/gbbct.spsum.gz") # contains codon counts
    query.spsum_label <- URLencode("ftp://ftp.ebi.ac.uk/pub/databases//cutg/SPSUM_LABEL.gz") # contains column headers

    # get codon counts and label
    spsum <- read_tsv(query.spsum, col_names = FALSE)
    spsum_label <- read_tsv(query.spsum_label, col_names = FALSE)
    
    # reformat into codon count table
    spsum <- spsum %>%
        mutate(ind = rep(c(1, 2),length.out = n())) %>%
        group_by(ind) %>%
        mutate(id = row_number()) %>%
        spread(ind, X1) %>%
        select(-id) %>%
        separate(`1`, into = c("organism", "cds"), sep = ": ", extra = "warn", fill = "warn") %>%
        separate("organism", into = c("tax.id", "organism"), sep = ":", extra = "merge", fill = "warn") %>%
        separate(`2`, into = str_split(spsum_label[2,], pattern = " ")[[1]], sep = " ", extra = "warn")
    
    # filter by tax_id
    filtered.spsum <- spsum %>%
        filter(tax.id == tax_id)
    
    # prepare codon counts dataset (take only codons)
    kazusa.ct <- as.data.frame(t(filtered.spsum[, 4:67]))
    kazusa.ct[] <- sapply(kazusa.ct, as.numeric) # counts are originally in character so convert to numeric
    kazusa.ct <- t(kazusa.ct) # transpose because sapply converts it to single col of values (we want in row-format)
    
    return(kazusa.ct)
  
}

get_cocoputs_codons <- function(tax_id) {

    # initialise variables
    query <- URLencode(paste0("https://dnahive.fda.gov/dna.cgi?cmd=ionTaxidCollapseExt&svcType=svc-codon-usage&objId=537&fileSource=Refseq_species.tsv&plen=3&taxid=", tax_id, "&filterInColName=%5B%22Organelle%22%5D&filterIn=%5B%22genomic%22%5D&searchDeep=true&raw=1&raw=1&bust=1674829688569"))

    cocoputs.cc <- read_csv(query)

    cocoputs.ct <- cocoputs.cc %>%
        filter(!row_number() %in% c(1:8)) %>%
        pivot_wider(names_from = id, values_from = value)

    return(cocoputs.ct)

}

get_cubseq_codon_counts <- function(fasta_dir) {

    # get fasta file paths
    files <- list.files(fasta_dir, pattern = "*-mut-transcriptome.fa", recursive = TRUE, full.names = TRUE)

    # count codons per sequence
    reads <- readDNAStringSet(files) # NB: fasta files are in DNA format, so we use readDNAStringSet()
    cubseq.ct <- trinucleotideFrequency(reads, step=3)

    # rename row names to transcript IDs (to use for future analysis)
    # rownames(cubseq.ct <- names(reads)) # fasta defline
    rownames(cubseq.ct) <- sub(" .*", "", names(reads)) # transcript ID
    # rownames(cubseq.ct) <- str_match(names(reads), "gene_name=\\s*(.*?)\\s*;")[,2] # gene name
    
    return(cubseq.ct)
}

# # compute frequencies of codons per thousand (N.B. this is the Kazusa method)
# getFreqPerThousand <- function(x) {
#   freqPerThousand = (x / totalCodon)*1000
#   return(freqPerThousand)
# }

get_relative_freq <- function(count_table, prefix) {
  
  # compute total codon counts
  ct <- as.data.frame(colSums(count_table))
  ct["codon"] <- rownames(ct)
  ct$codon <- gsub("T", "U", ct$codon) # convert to RNA
  colnames(ct) <- c("codon.count", "codon")

  # load codon to AA map - NB: RNA_GENETIC_CODE is an object from Biostrings
  codon2aa <- data.frame(codon = names(RNA_GENETIC_CODE), aa = paste(RNA_GENETIC_CODE))
  
  # map codons to AAs
  ct <- merge(ct, codon2aa, by.x = "codon", by.y = "codon")
  
  # calculate relative frequency
  # formula used: Cij(AAi) / sum(Cij(AAi)), so that all codons for ith AA adds up to 1
  sum.aa <- ct %>% 
    group_by(aa) %>% 
    summarise(sum.aa = sum(codon.count))
  ct <- merge(ct, sum.aa, by.x = "aa", by.y = "aa")
  ct <- transform(ct, relative.freq = codon.count / sum.aa)
  ct <- ct %>% rename(., !! paste0(prefix, ".count") := codon.count, !! paste0(prefix, ".sum.aa") := sum.aa, !! paste0(prefix, ".rf") := relative.freq)
  
  return(ct)
  
}

# calculate_difference <- function(df) {
  
#   # calculate pairwise difference (absolute)
#   df$abs.diff.HEG_kazusa <- abs(df$kazusa.rf - df$heg.rf)
#   df$abs.diff.LEG_kazusa <- abs(df$kazusa.rf - df$leg.rf)
#   df$abs.diff.HEG_LEG <- abs(df$heg.rf - df$leg.rf)
  
#   # calculate pairwise difference (squared)
#   df$sq.diff.HEG_kazusa <- (df$kazusa.rf - df$heg.rf)^2
#   df$sq.diff.LEG_kazusa <- (df$kazusa.rf - df$leg.rf)^2
#   df$sq.diff.HEG_LEG <- (df$heg.rf - df$leg.rf)^2
  
#   return(df)
  
# }

main <- function() {

    # # load CoCoPUTs K-12 counts data
    # test.cocoputs.cc <- read_csv(arguments$cocoputs_counts)

    # get codon counts for Kazusa
    kazusa.cc <- get_kazusa_codons(arguments$tax_id)

    # get codon counts for CoCoPUTs
    cocoputs.cc <- get_cocoputs_codons(arguments$tax_id)

    # count codons for HEG and LEG and protein genes
    heg.cc <- get_cubseq_codon_counts(arguments$heg_fasta_dir)
    leg.cc <- get_cubseq_codon_counts(arguments$leg_fasta_dir)
    protein.cc <- get_cubseq_codon_counts(arguments$protein_fasta_dir)

    # calculate relative frequency of codons
    heg.rf <- get_relative_freq(heg.cc, prefix = "heg")
    leg.rf <- get_relative_freq(leg.cc, prefix = "leg")
    protein.rf <- get_relative_freq(protein.cc, prefix = "protein")
    kazusa.rf <- get_relative_freq(kazusa.cc, prefix = "kazusa")
    cocoputs.rf <- get_relative_freq(cocoputs.cc, prefix = "cocoputs")

    # merge dataframes into one
    all.rf <- merge(heg.rf, leg.rf) %>%
        merge(protein.rf) %>%
        merge(kazusa.rf) %>%
        merge(cocoputs.rf) %>%
        select(., aa, codon, heg.rf, leg.rf, protein.rf, kazusa.rf, cocoputs.rf)

    # # calculate absolute and squared difference
    # all.rf <- calculate_difference(all.rf)

    # calculate log2 comparison
    # cubseq vs. kazusa
    all.rf$log2.heg.kazusa <- log2(all.rf$heg.rf / all.rf$kazusa.rf)
    all.rf$log2.leg.kazusa <- log2(all.rf$leg.rf / all.rf$kazusa.rf)
    all.rf$log2.protein.kazusa <- log2(all.rf$protein.rf / all.rf$kazusa.rf)
    # cubseq vs. cubseq
    all.rf$log2.heg.leg <- log2(all.rf$heg.rf / all.rf$leg.rf)
    all.rf$log2.heg.protein <- log2(all.rf$heg.rf / all.rf$protein.rf)
    # cubseq vs. cocoputs
    all.rf$log2.heg.cocoputs <- log2(all.rf$heg.rf / all.rf$cocoputs.rf)
    all.rf$log2.leg.cocoputs <- log2(all.rf$leg.rf / all.rf$cocoputs.rf)
    all.rf$log2.protein.cocoputs <- log2(all.rf$protein.rf / all.rf$cocoputs.rf)
    # cocoputs vs. kazusa
    all.rf$log2.cocoputs.kazusa <- log2(all.rf$cocoputs.rf / all.rf$kazusa.rf)

    # write count dataframes to csv
    write.table(heg.cc, file = "heg-codon-counts.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(leg.cc, file = "leg-codon-counts.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)
    write.table(protein.cc, file = "protein-codon-counts.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)
    write_csv(as.data.frame(kazusa.cc), file = paste0("kazusa-", arguments$tax_id, "-codon-counts.csv"))
    write_csv(cocoputs.cc, file = paste0("cocoputs-", arguments$tax_id, "-codon-counts.csv"))

    # write summary CU frequency dataframe to csv
    write_csv(all.rf, file = paste0("all-kazusa-cocoputs-", arguments$tax_id, "-rf.csv"))

}

main()
