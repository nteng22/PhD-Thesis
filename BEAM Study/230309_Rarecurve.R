setwd("~/Dropbox/QIB/BEAM R Analysis/Shotgun data/")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(vegan)

species <- read_tsv("combined_species.tsv", col_names = TRUE)
species_num <- species %>% 
  select(name, ends_with("species_num"))
species_transposed <- column_to_rownames(species_num, "name")
species_transposed <- as.data.frame(t(species_transposed))
species_transposed <- rownames_to_column(species_transposed, var = "Sample")

species_transposed$Sample <- gsub("-nonhost-species_num", "", 
                                  species_transposed$Sample)

species_transposed <- column_to_rownames(species_transposed, var = "Sample")

# Make the Genus names into row names so when the dataframe is all numbers.
raw_data <- species_transposed

rarecurve((raw_data),
          step = 20,
          sample = 30000,
          col = "blue",
          cex = 0.4,
          label = TRUE,
          #label = FALSE,
          lwd = 0.001,
          ylab = "No. of species",
          xlab = "No. of sequences"
          #xlim = c(0, 250000)
)

# Genus
genus <- read_tsv("combined_genus.tsv", col_names = TRUE)
genus_num <- genus %>% 
  select(name, ends_with("genus_num"))
genus_transposed <- column_to_rownames(genus_num, "name")
genus_transposed <- as.data.frame(t(genus_transposed))
genus_transposed <- rownames_to_column(genus_transposed, var = "Sample")

genus_transposed$Sample <- gsub("-nonhost-genus_num", "", 
                                  genus_transposed$Sample)

genus_transposed <- column_to_rownames(genus_transposed, var = "Sample")

NNUH_genus <- read_tsv("combined_genus-Nancy.tsv", col_names = TRUE)
NNUH_genus_num <- NNUH_genus %>% 
  select(name, ends_with("genus_num"))
NNUH_genus_transposed <- column_to_rownames(NNUH_genus_num, "name")
NNUH_genus_transposed <- as.data.frame(t(NNUH_genus_transposed))
NNUH_genus_transposed <- rownames_to_column(NNUH_genus_transposed, var = "Sample")

NNUH_genus_transposed$Sample <- gsub("-nonhost-genus_num", "", 
                                NNUH_genus_transposed$Sample)
NNUH_genus_transposed$NNUH <- "NNUH00"
NNUH_genus_transposed <- unite(NNUH_genus_transposed,
                               NNUH, Sample, col = "Sample", sep = "")

NNUH_genus_transposed <- column_to_rownames(NNUH_genus_transposed, var = "Sample")

# Make the Genus names into row names so when the dataframe is all numbers.
raw_data <- rbind(genus_transposed, NNUH_genus_transposed)
raw_data <- full_join()
rarecurve((raw_data),
          step = 200,
          sample = 30000,
          col = "blue",
          cex = 0.4,
          label = TRUE,
          #label = FALSE,
          lwd = 0.001,
          ylab = "No. of genus",
          xlab = "No. of sequences"
          #xlim = c(0, 250000)
)

