setwd("~/Dropbox/QIB/BEAM R Analysis/Shotgun data/")

library(dplyr)
library(tidyverse)
library(ggplot2)

# Data formatting- GENUS
genus <- read_tsv("combined_genus.tsv", col_names = TRUE)

# Seperate fraction and numbers, we will use frac (%)
genus_frac <- genus %>% 
  select(name, ends_with("genus_frac"))
genus_num <- genus %>% 
  select(name, ends_with("genus_num"))

# Format the data table, so sample ID is the first column 
# known as wide format, we actually prefer long format
genus_transposed <- column_to_rownames(genus_frac, "name")
genus_transposed <- as.data.frame(t(genus_transposed))
genus_transposed <- rownames_to_column(genus_transposed, var = "Sample")

genus_transposed$Sample <- gsub("-nonhost-genus_frac", "", 
                                genus_transposed$Sample)

SampleID <- genus_transposed %>%
  select(Sample)

# Check normalisation of the genus_frac
abundances <- genus_transposed %>%
  select("Cloacibacillus":ncol(genus_transposed)) 
rowSums(abundances) # To check if it's 100%
abundances_norm <- abundances/rowSums(abundances)*100
rowSums(abundances_norm)

genus_norm_sampleID <- cbind(SampleID, abundances_norm)

write.csv(genus_norm_sampleID, "BEAM_genus_normalised.csv")

# Dataformatting -SPECIES
species <- read_tsv("combined_species.tsv", col_names = TRUE)

# Seperate fraction and numbers, we will use frac (%)
species_frac <- species %>% 
  select(name, ends_with("species_frac"))
species_num <- species %>% 
  select(name, ends_with("species_num"))

# Format the data table, so sample ID is the first column 
# known as wide format, we actually prefer long format
species_transposed <- column_to_rownames(species_frac, "name")
species_transposed <- as.data.frame(t(species_transposed))
species_transposed <- rownames_to_column(species_transposed, var = "Sample")

species_transposed$Sample <- gsub("-nonhost-species_frac", "", 
                                species_transposed$Sample)

SampleID <- species_transposed %>%
  select(Sample)

# Check normalisation of the species_frac
abundances <- species_transposed %>%
  select("Synechococcus sp. Minos11":ncol(species_transposed)) 
rowSums(abundances) # To check if it's 100%
abundances_norm <- abundances/rowSums(abundances)*100
rowSums(abundances_norm)

species_norm_sampleID <- cbind(SampleID, abundances_norm)

write.csv(species_norm_sampleID, "BEAM_species_normalised.csv")
