# load libraries
library(tidyverse)

# load metadata
metadata <- read.delim("final-data/metadata/metadata.csv")

filtered_metadata <- metadata %>%
  mutate(scientific_name = str_squish(metadata$scientific_name)) %>%
  mutate(strain = str_squish(metadata$strain)) %>%
  mutate(scientific_name = toupper(metadata$scientific_name)) %>%
  mutate(strain = toupper(metadata$strain))

# for debugging purposes
print(dim(metadata))
print(dim(filtered_metadata))
print(all_equal(filtered_metadata))