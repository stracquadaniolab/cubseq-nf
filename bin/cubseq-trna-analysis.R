#!/usr/bin/env Rscript

"cubseq-trna-analysis.R

Usage:
    cubseq-trna-analysis.R <input> [--metadata=<metadata>] [--gtf=<gtf>] [--fdr-threshold=<fdr_threshold>] 
    [--aa-property=<aa_property>]
    
Options:
    --txi=<txi>                                     Tximport file.
    --gtf=<gtf>                                     GTF file.
    --fdr-threshold=<fdr_threshold>                 Set FDR threshold [default: 0.01].
    --aa-property=<aa_property>                     AA properties csv.
    -h --help                                       Show this screen.
    --version                                       Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-trna-analysis.R")

# load libraries
suppressMessages(library(tidyverse))
suppressMessages(library(rtracklayer))
suppressMessages(library(data.table))
suppressMessages(library(tximport))
suppressMessages(library(Biostrings))
suppressMessages(library(stringi))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(paletteer))
suppressMessages(library(ggrepel))
suppressMessages(library(viridis))


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


anticodon2codon <- function(.data) {

    # function to convert anticodon to codon
  
    # first reverse anticodon (5'->3')
    .data$rev.anticodon <- stri_reverse(.data$anticodon)
    
    # replace reversed anticodon with complementary base pair
    # use chartr to map characters to characters
    .data$codon = chartr("ATGC", "TACG", .data$rev.anticodon)
    
    # convert to RNA format (replace T with U)
    .data$codon <- gsub("T", "U", .data$codon)
    
    return(.data)
  
}

# figure_trna_scatter <- function(data, x, y, aa_label, col, legend, xlab, ylab) {
  
#   cu_scatter <- ggplot(data, aes(x=x, y=y, color=col)) + 
#     geom_point() +
#     theme_cowplot(12) +
#     labs(color=legend) +
#     xlab(xlab) +
#     ylab(ylab) +
#     geom_text_repel(label = aa_label, size=2.5, show.legend = FALSE)
  
# }

figure_trna_scatter <- function(.data, .y, .title) {
    
    trna_plot_scatter <- .data %>%
        ggplot(., aes(x=trna_mean_score, y=.data[[.y]], color=aa.property)) + 
        geom_point() +
        geom_smooth(method="lm", fill = "transparent") +
        theme_cowplot(12) +
        labs(color="Amino acid\nproperty") +
        xlab("tRNA mean abundance score") +
        ylab("Corresponding tRNA codon relative frequency") +
        geom_text_repel(label = trna_cu_unique$codon, size=2.5, show.legend = FALSE) +
        facet_wrap(~tRNA_type, scales = "free") +
        theme_cowplot(12) +
        ggtitle("Relationship between tRNA abundance and corresponding codon frequency") +
        scale_color_paletteer_d("ggsci::nrc_npg") +
        theme(plot.title = element_text(hjust = 0.5),
                strip.background =element_rect(fill="white"),
                strip.text = element_text(face="bold"),
                axis.line = element_blank(),
                panel.border = element_rect(colour = "#737373",
                                            fill=NA,
                                            linewidth=0.3))

}



plot_trna_facet <- function(.data, .y, .title) {
  
  trna_plot <- .data %>%
    group_by(tRNA_type) %>%
    filter(n_distinct(codon) > 1) %>%
    ggplot(., aes(x = reorder(gene_id, trna_rank), y=.data[[.y]], fill=aa.property)) +
    geom_bar(aes(x = reorder(AC, trna_rank)), stat = "identity") +
    labs(fill="Amino acid\nproperty") +
    xlab("tRNA anticodon abundance (ordered by mean normalised rank)") +
    ylab("Corresponding tRNA codon relative frequency") +
    facet_wrap(vars(tRNA_type), scales = "free") +
    theme_cowplot(12) +
    # theme(axis.text.x = element_text(angle = 45, vjust = 0.8)) +
    ggtitle(.title) +
    scale_fill_paletteer_d("ggsci::nrc_npg") +
    theme(plot.title = element_text(hjust = 0.5),
          strip.background =element_rect(fill="white"),
          strip.text = element_text(face="bold"),
          axis.line = element_blank(),
          panel.border = element_rect(colour = "#737373",
                                      fill=NA,
                                      linewidth=0.3))
  
}


main <- function() {
    
    # load datasets
    gtf <- import(arguments$gtf)
    txi <- readRDS(arguments$txi)
    cu <- read_csv(arguments$cu)
    aa.property <- read_csv(arguments$aa_property)


    # get tRNA codons from GtRNAdb
    gtrnadb <- read_delim(
        "eschColi_K_12_MG1655-tRNAs/eschColi_K_12_MG1655-tRNAs.out", 
        delim = "\t", escape_double = FALSE,
        col_names = FALSE, 
        trim_ws = TRUE, 
        skip = 3)

    col_names <- c("Seq", "tRNA_num", "start", "end", 
                 "tRNA_type", "anticodon", "intron_bounds_begin",
                 "intron_bounds_end", "Inf_HMM_score", "2Str_score",
                 "Isotype_score", "Isotype_CM", "Score", "Note")
  
    setnames(gtrnadb, names(gtrnadb), col_names)


    # filter GTF for genes encoding tRNA
    trna_gtf <- as.data.frame(gtf) %>%
        filter(gene_biotype == "tRNA", type == "exon")


    # tidy dataframe
    trna_gtf <- trna_gtf %>% 
        mutate(start1 = ifelse(strand == "-", end, start), 
            end = ifelse(strand =="-", start, end)) %>% 
        select(seqnames, start = start1, end, width, strand, gene_id, transcript_id) %>%
        mutate(start = ifelse(gene_id == "hisR", 3982510, start))


    # # get tRNA codons from GtRNAdb
    # eco_trna <- read_tsv("http://gtrnadb.ucsc.edu/download/GtRNAdb/search/gtrnadb-search153044.out")
    # # tidy dataframe
    # eco_trna <- eco_trna %>%
    #     separate(Locus, into = c("Locus", "Strand1"), sep = " ") %>%
    #     separate(Locus, into = c("start", "end"), sep = "-") %>%
    #     mutate_all(~gsub("chr:", "", .))
    # # fix column names
    # names(eco_trna)[8:(ncol(eco_trna)-1)] <- names(eco_trna)[9:ncol(eco_trna)]
    # # delete last column
    # eco_trna <- subset(eco_trna, select=-c(15))
    
    # map tRNA anticodons to GTF by locus start and end
    trna_gtf <- merge(trna_gtf, gtrnadb, by=c("start", "end"))

    # filter TPMs for tRNA genes only
    tpm <- txi$abundance
    trna_tpm <- subset(tpm, rownames(tpm) %in% trna_gtf$gene_id)

    # get mean normalised rank score - we'll use "mean_score"
    trna_chisq <- chi_meta_analysis(trna_tpm)

    # # get df --> gene_id, anticodon, codon, TPM score
    # trna_df <- as.data.frame(trna_chisq$results$mean_score)
    # trna_df$gene_id <- rownames(trna_df)
    # colnames(trna_df)[1] <- "mean_norm_rank"

    # trna_df <- merge(trna_df, trna_gtf, by=c("gene_id"))

    # trna_df <- trna_gtf %>%
    #     select(gene_id, Isotype, Anticodon, tRNAscanID) %>%
    #     separate(tRNAscanID, into = c("tRNAscanID", "tRNA_name"), sep = "-") %>%
    #     select(-tRNAscanID) %>%
    #     merge(., mean_score_df, by = c("gene_id")) %>%
    #     distinct()

    # get df --> gene_id, AA, AC, codon, avg.rank
    trna_df <- as.data.frame(trna_chisq$results) %>%
        select(., gene_id = genes, mean_score)

    trna_df <- merge(trna_df, trna_gtf, by=c("gene_id"))

    trna_df <- anticodon2codon(trna_df)

    # map df codons to CU R.F. codons
    trna_cu <- cu %>%
        select(aa, codon, heg.rf, leg.rf, kazusa.rf, cocoputs.rf) %>%
        merge(., trna_df, by = c("codon")) %>%
        # remove Met, STOP and Trp
        filter(!str_detect(aa, 'M|W')) %>%
        filter(!str_detect(tRNA_type, 'SeC')) %>%
        # rank tRNAs by mean_score (highest score = rank 1)
        mutate(trna_rank = rank(-mean_score, ties.method = "average")) %>%
        select(., gene_id, tRNA_type, AA = aa, AC = anticodon, codon, 
            trna_mean_score = mean_score, trna_rank, heg.rf, leg.rf, 
            kazusa.rf, cocoputs.rf)

    # map amino acid properties to df
    trna_cu <- merge(trna_cu, aa.property, by.x = "AA", by.y = "aa")

    # plot scatter
    trna.cubseq.heg <- plot_trna_facet(trna_cu, "heg.rf", "CUBseq highly expressed genes")
    trna.cubseq.leg <- plot_trna_facet(trna_cu, "leg.rf", "CUBseq lowly expressed genes")
    trna.kazusa <- plot_trna_facet(trna_cu, "kazusa.rf", "Kazusa genes")
    trna.cocoputs <- plot_trna_facet(trna_cu, "cocoputs.rf", "CoCoPUTs genes")

    # save figures
    save_plot(filename = "trna-fig1-cubseq-heg.png", plot = trna.cubseq.heg, base_height = 6, base_width = 9)
    save_plot(filename = "trna-fig2-cubseq-leg.png", plot = trna.cubseq.leg, base_height = 6, base_width = 9)
    save_plot(filename = "trna-fig3-kazusa.png", plot = trna.kazusa, base_height = 6, base_width = 9)
    save_plot(filename = "trna-fig4-cocoputs.png", plot = trna.cocoputs, base_height = 6, base_width = 9)

    # save files
    write_csv(trna_gtf, "trna-gtf.csv")
    write_csv(eco_trna, "gtrnadb-eco_k12_mg1655.csv")
    write_csv(trna_chisq$results, "trna-scores.csv")
    write_csv(trna_df, "trna-all.csv")
    saveRDS(trna_tpm, "trna-tpm.rds")

}

main()
