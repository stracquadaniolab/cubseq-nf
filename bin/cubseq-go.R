#!/usr/bin/env Rscript

"cubseq-go.R

Usage:
  cubseq-go.R <inputfile> <outputfile> [--fdr=<alpha>] [--db-organism=<db_org>] [--gene-id=<gene_id>] [--remove-gencode-version=<gv>]
  
Options:
  -f --fdr=<alpha>                      False Discovery Rate threshold to filter differentially
                                        expressed genes [default: 0.01].
  -d --db-organism=<db_org>             Organism database for gene annotation [default: org.Hs.eg.db].
  -g --gene-id=<gene_id>                Gene ID type [default: ensembl].
  --remove-gencode-version=<gv>         Remove Genecode gene version [default: no].
  -h --help                             Show this screen.
  --version                             Show version.
" -> doc

# parsing command line arguments
library(docopt)
arguments <- docopt(doc, version = "cubseq-go.R")

# loading data processing libraries
suppressMessages(library(tidyverse))
suppressMessages(library(GO.db))
suppressMessages(library(topGO))

# reading deseq object
res <- read_csv(arguments$inputfile)
res <- res[!is.na(res$padj), ]

if (arguments$remove_gencode_version == "yes") {
  print("removing gencode version")
  res$gene_id <- str_replace(res$gene_id, "\\.[0-9]*", "")
}

# gene list
gene_list <- res$padj
names(gene_list) <- res$gene_id

# GO data
go_data <- new("topGOdata",
  ontology = "BP",
  allGenes = gene_list,
  geneSelectionFun = function(.x) {
    .x <= as.numeric(arguments$fdr)
  },
  annot = annFUN.org, mapping = arguments$db_organism, ID = arguments$gene_id
)

# perform classical fisher test
fisher_test <- runTest(go_data, algorithm = "classic", statistic = "fisher")
fisher_results <- GenTable(go_data, classicFisher = fisher_test, topNodes = 1000)

# write results to file
write_csv(fisher_results, arguments$outputfile)