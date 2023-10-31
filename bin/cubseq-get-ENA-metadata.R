#!/usr/bin/env Rscript

"cubseq-get-ENA-metadata.R

Usage:
    cubseq-get-ENA-metadata.R [--taxon-id=<taxon_id>] [--limit-search=<limit_search>] [--remove-run=<remove_run>] [--max-sra-bytes=<max_sra_bytes>] [--date-min=<date_min>] [--date-max=<date_max>] <outputfile> 
    
Options:
    --taxon-id=<taxon_id>           Taxon ID to filter by [default: 562].
    --limit-search=<limit_search>   Limit number of records output from search query [default: 0].
    --remove-run=<remove_run>       Remove run by specifying its run accession [default: NULL].
    --max-sra-bytes=<max_sra_bytes> Specify runs to remove if they exceed size of sra_bytes [default: 55000000000].
    --date-min=<date_min>           Set minimum date (YYYY/MM/DD) to filter runs by (inclusive) [default: 1950-01-01].
    --date-max=<date_max>           Set maximum date (YYYY/MM/DD) to filter runs by (inclusive), uses current date by default [default: 2100-01-01].
    -h --help                       Show this screen.
    --version                       Show version.
" -> doc

# parsing command line arguments
suppressMessages(library(docopt))
arguments <- docopt(doc, version = "cubseq-get-ENA-metadata.R")

# load packages
suppressMessages(library(tidyverse))
print(arguments)
# query ENA repository
api_url <- "https://www.ebi.ac.uk/ena/portal/api/search"
tax_id <- arguments$taxon_id
lib_src <- "TRANSCRIPTOMIC"
limit_search <- arguments$limit_search
remove_run <- c(arguments$remove_run)

criteria <- sprintf("tax_tree(%s) AND library_layout=PAIRED AND library_source=%s AND instrument_platform=ILLUMINA &limit=%s", tax_id, lib_src, limit_search)
fields <- "tax_id,scientific_name,strain,first_public,last_updated,study_accession,study_title,experiment_accession,experiment_title,sample_accession,sample_alias,sample_description,run_accession,run_alias,description,collection_date,library_source,library_selection,instrument_platform,instrument_model,library_strategy,library_layout,read_count,base_count,sra_bytes,fastq_ftp"

query <- URLencode(sprintf("%s?result=read_run&query=%s&fields=%s", api_url, criteria, fields))

metadata <- read_tsv(query)

# initialise variables for metadata cleaning
max_sra_bytes <- as.numeric(arguments$max_sra_bytes)

# clean metadata
metadata <- metadata %>%
  filter(library_strategy == "RNA-Seq") %>%
  drop_na(fastq_ftp) %>%
  # remove single fastq files (could be interleaved PE data, or mislabelled SE data)
  filter(str_detect(fastq_ftp, ";")) %>%
  # filter out problem run_acc (e.g. for dev testing purposes)
  filter(!run_accession %in% remove_run) %>%

  # extract only fastq1 and fastq2 ftp links
  mutate(fastq_ftp_2 = strsplit(fastq_ftp, ";")) %>%
  unnest(fastq_ftp_2) %>%
  group_by(run_accession) %>%
  mutate(row = row_number()) %>%
  filter(grepl("_1.fastq.gz|_2.fastq.gz", fastq_ftp_2)) %>%
  spread(row, fastq_ftp_2) %>%
  ungroup() %>%
  unite(., col = "fastq_ftp_2", `1`:last_col(), na.rm = TRUE, sep = ";") %>%
  separate(fastq_ftp_2, into = c("fastq1", "fastq2"), sep = ";") %>%

  # filter out abnormal files (e.g. large file size but keep NA)
  # N.B. these are files that could potentially halt pipeline
  filter(., sra_bytes <= max_sra_bytes | is.na(sra_bytes)) %>%

  # filter by specified date range
  filter(first_public >= arguments$date_min & first_public <= if_else(is.null(arguments$date_max), paste0(Sys.Date()), arguments$date_max))

# write to csv
write.table(metadata, file = arguments$outputfile, quote = FALSE, sep = "\t", na = "", row.names = FALSE)
