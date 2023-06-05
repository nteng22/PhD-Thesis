#### Housekeeping ####
setwd("~/Dropbox/QIB/KELLY Study/Clinical data/")
setwd("C:\\Users/tengn/Dropbox/QIB/KELLY Study/Clinical data/")

#### Download libraries ####
#install.packages("dplyr")
library("dplyr")

#install.packages("tidyverse")
library("tidyverse")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("readxl")
library("readxl")

# Import the spreadsheet
clinical <- read_excel("../BBDD_Outcomes_Kelly_NLR_20211015 (1).xlsx")

clinical <- clinical %>%
  select(ID_Patient, PDL1, NLR, CBR, `Resp_RECIST Calculated`, tPFS)

clinical$ID_Patient <- gsub("-", "_", clinical$ID_Patient)
clinical <- dplyr::rename(clinical, RECIST = "Resp_RECIST Calculated")

#Change tPFS to less than 6 months and after 6 months
# 1 tPFS < 6
# 0 tPFS > 6
clinical$tPFS <- ifelse(clinical$tPFS > 6.0000 , "1", clinical$tPFS)
clinical$tPFS <- ifelse(!clinical$tPFS == 1 , "0", clinical$tPFS)

#For genus data
# Try to make a master spread sheet with all clincal data and abundance info together
normalised_data <- read.csv("../genus_sample_merged.csv")
normalised_data <- normalised_data %>%
  select(ID_Patient:ncol(normalised_data))

metadata_clinical_merged <- inner_join(clinical, normalised_data, by = "ID_Patient") 
metadata_clinical_merged <- drop_na(metadata_clinical_merged)

# Only clinical data
metadata_clinical <- metadata_clinical_merged %>%
  select(ID_Patient:Time_point)
write.csv(metadata_clinical, "KELLY clinical metadata.csv")

# clinical with reads
write.csv(metadata_clinical_merged, "KELLY clinical metadata data-genus.csv")

#For species data
normalised_data <- read.csv("species_sample_merged.csv")
normalised_data <- normalised_data %>%
  select(ID_Patient:ncol(normalised_data))
metadata_clinical_merged <- inner_join(clinical, normalised_data, by = "ID_Patient") 
metadata_clinical_merged <- drop_na(metadata_clinical_merged)

write.csv(metadata_clinical_merged, "KELLY clinical metadata data-species.csv")
# data is not in percentages!