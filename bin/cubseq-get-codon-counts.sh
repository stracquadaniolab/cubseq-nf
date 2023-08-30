#!/bin/bash

# WHY DO I NEED THIS SCRIPT? My protein-codon-counts.csv file was
# too big to load in R. Then had issues loading tidyverse, so HPC
# was no longer an option. Hence I am using Bash. This script takes
# the large codon counts file and summarises counts at codon level.
# Resulting output file can be used for downstream analysis, such as
# calculating total codon counts, and codon frequency per thousand
# for generating the codon table.

# first we need to remove transcript IDs (first column, consider these as row names)
awk 'NR==1 {print $0} NR>1 { for(i=2;i<=NF;i++) printf "%s%s", $i, (i==NF?ORS:OFS)}' protein-codon-counts.csv > protein-codon-counts-new.csv

# now we can sum the codon counts
# initialise file names
filename="protein-codon-counts-new.csv"
output_file="protein-codon-summary-counts.csv"

# get total number of columns in csv file
num_columns=$(head -n 1 "$filename" | tr ' ' '\n' | wc -l)

# sum each column
awk -F ' ' -v num_columns="$num_columns" 'BEGIN {
    for (i=1; i<=num_columns; i++) {
        sum[i] = 0
    }
}
NR == 1 {
    for (i=1; i<=num_columns; i++) {
        headers[i] = $i
    }
}
{
    for (i=1; i<=num_columns; i++) {
        sum[i] += $i
    }
}
END {
    printf("codon,sum\n")
    for (i=1; i<=num_columns; i++) {
        printf("%s,%s\n", headers[i], sum[i])
    }
}' "$filename" > "$output_file"
