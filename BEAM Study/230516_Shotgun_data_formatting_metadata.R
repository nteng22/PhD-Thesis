#### HOUSEKEEPING ####
setwd("~/Dropbox/QIB/BEAM R Analysis/Shotgun data/")

library("dplyr")
library("tidyverse")
library("vegan")
library("dplyr")
library("phyloseq")

#### METADATA ####
library(readxl)
metadata <- read_xlsx("../BEAM_metadata.xlsx")

metadata$BrCa <- ifelse(metadata$BrCa== "NEUROENDOCRINE", "Unclassified", metadata$BrCa)
metadata$BrCa[is.na(metadata$BrCa)] <- "Unclassified"
metadata$chemotherapy[is.na(metadata$chemotherapy)] <- "Unclassified"
metadata$age[is.na(metadata$age)] <- "Unclassified"

final_data_species <- read_csv("BEAM_species_normalised.csv")
final_data_genus <- read_csv("BEAM_genus_normalised.csv")

metadata_finaldata <- final_data_species %>%
  select(Sample)
NNUH <- metadata_finaldata[17:34,] 
NNUH$prefix <- c("NNUH")
NNUH <- unite(NNUH, prefix, Sample,
              sep = "00",
              col = Sample)
metadata_finaldata[17:34,] <- NNUH

sampleID_separated <- separate(metadata_finaldata, Sample,
                               sep = "(?<=[0-9])(?=[A-Z])",
                               into = c("patient", "time_point")) # Need to remove the P in NNUH
final_metadata <- cbind(sampleID_separated, metadata_finaldata)

metadata_final <- inner_join(final_metadata, metadata, by = "patient")

final_data_species$Sample <- final_metadata$Sample

abundance_metaadata_final <- inner_join(metadata_final, final_data_species, 
                                        by = "Sample")

abundance_metaadata_final$time_point <- gsub("D", "One-Year", abundance_metaadata_final$time_point)
abundance_metaadata_final$time_point <- gsub("C", "6-Months", abundance_metaadata_final$time_point)
abundance_metaadata_final$time_point <- gsub("B", "Post-Surgery", abundance_metaadata_final$time_point)
abundance_metaadata_final$time_point <- gsub("A", "Baseline", abundance_metaadata_final$time_point)

abundance_metaadata_final <- abundance_metaadata_final[,-(7)]

write_csv(abundance_metaadata_final, "BEAM_final_metadata_normalised_species_shotgun.csv")

###
abundance_metaadata_final <- inner_join(metadata_final, final_data_genus, 
                                        by = "Sample")

abundance_metaadata_final$time_point <- gsub("D", "One-Year", abundance_metaadata_final$time_point)
abundance_metaadata_final$time_point <- gsub("C", "6-Months", abundance_metaadata_final$time_point)
abundance_metaadata_final$time_point <- gsub("B", "Post-Surgery", abundance_metaadata_final$time_point)
abundance_metaadata_final$time_point <- gsub("A", "Baseline", abundance_metaadata_final$time_point)

abundance_metaadata_final <- abundance_metaadata_final[,-(7)]

write_csv(abundance_metaadata_final, "BEAM_final_metadata_normalised_genus_shotgun.csv")
