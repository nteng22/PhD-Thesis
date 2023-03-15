setwd("~/Dropbox/QIB/BEAM R Analysis/Shotgun data/")

library(dplyr)
library(tidyverse)
library(ggplot2)

# Data formatting- GENUS
genus_JPUH <- read_tsv("combined_genus.tsv", col_names = TRUE)
genus_NNUH <- read_tsv("combined_genus-Nancy.tsv", col_names = TRUE)

# Separate fraction and numbers, we will use frac (%)
genus_frac_JPUH <- genus_JPUH %>% 
  select(name, ends_with("genus_frac"))

# Format the data table, so sample ID is the first column 
# known as wide format, we actually prefer long format
genus_transposed_JPUH <- column_to_rownames(genus_frac_JPUH, "name")
genus_transposed_JPUH <- as.data.frame(t(genus_transposed_JPUH))
genus_transposed_JPUH <- rownames_to_column(genus_transposed_JPUH, var = "Sample")

genus_transposed_JPUH$Sample <- gsub("-nonhost-genus_frac", "", 
                                genus_transposed_JPUH$Sample)

SampleID_JPUH <- genus_transposed_JPUH %>%
  select(Sample)

# Check normalisation of the genus_frac
abundances_JPUH <- genus_transposed_JPUH %>%
  select("Cloacibacillus":ncol(genus_transposed_JPUH)) 
rowSums(abundances_JPUH) # To check if it's 100%
abundances_norm_JPUH <- abundances_JPUH/rowSums(abundances_JPUH)*100
rowSums(abundances_norm_JPUH)

genus_sample_JPUH <- cbind(SampleID_JPUH, abundances_norm_JPUH)
  
# Separate fraction and numbers, we will use frac (%)
genus_frac_NNUH <- genus_NNUH %>% 
  select(name, ends_with("genus_frac"))

# Format the data table, so sample ID is the first column 
# known as wide format, we actually prefer long format
genus_transposed_NNUH <- column_to_rownames(genus_frac_NNUH, "name")
genus_transposed_NNUH <- as.data.frame(t(genus_transposed_NNUH))
genus_transposed_NNUH <- rownames_to_column(genus_transposed_NNUH, var = "Sample")

genus_transposed_NNUH$Sample <- gsub("-nonhost-genus_frac", "", 
                                     genus_transposed_NNUH$Sample)

SampleID_NNUH <- genus_transposed_NNUH %>%
  select(Sample)

# Check normalisation of the genus_frac
abundances_NNUH <- genus_transposed_NNUH %>%
  select("Cloacibacillus":ncol(genus_transposed_NNUH)) 
rowSums(abundances_NNUH) # To check if it's 100%
abundances_norm_NNUH <- abundances_NNUH/rowSums(abundances_NNUH)*100
rowSums(abundances_norm_NNUH)

genus_sample_NNUH <- cbind(SampleID_NNUH, abundances_norm_NNUH)

genus_frac_BEAM <- full_join(genus_sample_JPUH, genus_sample_NNUH)
genus_frac_BEAM[is.na(genus_frac_BEAM)] <- 0

write.csv(genus_frac_BEAM, "BEAM_genus_normalised.csv")

# Dataformatting -SPECIES
species_JPUH <- read_tsv("combined_species.tsv", col_names = TRUE)
species_NNUH <- read_tsv("combined_species-Nancy.tsv", col_names = TRUE)

# Separate fraction and numbers, we will use frac (%)
species_frac_JPUH <- species_JPUH %>% 
  select(name, ends_with("species_frac"))

# Format the data table, so sample ID is the first column 
# known as wide format, we actually prefer long format
species_transposed_JPUH <- column_to_rownames(species_frac_JPUH, "name")
species_transposed_JPUH <- as.data.frame(t(species_transposed_JPUH))
species_transposed_JPUH <- rownames_to_column(species_transposed_JPUH, var = "Sample")

species_transposed_JPUH$Sample <- gsub("-nonhost-species_frac", "", 
                                     species_transposed_JPUH$Sample)

SampleID_JPUH <- species_transposed_JPUH %>%
  select(Sample)

# Check normalisation of the species_frac
abundances_JPUH <- species_transposed_JPUH %>%
  select(`Synechococcus sp. Minos11`:ncol(species_transposed_JPUH)) 
rowSums(abundances_JPUH) # To check if it's 100%
abundances_norm_JPUH <- abundances_JPUH/rowSums(abundances_JPUH)*100
rowSums(abundances_norm_JPUH)

species_sample_JPUH <- cbind(SampleID_JPUH, abundances_norm_JPUH)

# Separate fraction and numbers, we will use frac (%)
species_frac_NNUH <- species_NNUH %>% 
  select(name, ends_with("species_frac"))

# Format the data table, so sample ID is the first column 
# known as wide format, we actually prefer long format
species_transposed_NNUH <- column_to_rownames(species_frac_NNUH, "name")
species_transposed_NNUH <- as.data.frame(t(species_transposed_NNUH))
species_transposed_NNUH <- rownames_to_column(species_transposed_NNUH, var = "Sample")

species_transposed_NNUH$Sample <- gsub("-nonhost-species_frac", "", 
                                     species_transposed_NNUH$Sample)

SampleID_NNUH <- species_transposed_NNUH %>%
  select(Sample)

# Check normalisation of the species_frac
abundances_NNUH <- species_transposed_NNUH %>%
  select("Bacteroides sp. CBA7301":ncol(species_transposed_NNUH)) 
rowSums(abundances_NNUH) # To check if it's 100%
abundances_norm_NNUH <- abundances_NNUH/rowSums(abundances_NNUH)*100
rowSums(abundances_norm_NNUH)

species_sample_NNUH <- cbind(SampleID_NNUH, abundances_norm_NNUH)

species_frac_BEAM <- full_join(species_sample_JPUH, species_sample_NNUH)
species_frac_BEAM[is.na(species_frac_BEAM)] <- 0

write.csv(species_frac_BEAM, "BEAM_species_normalised.csv")
