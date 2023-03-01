#### Housekeeping ####
setwd("~/Dropbox/QIB/KELLY Study/")

#### Download libraries ####
#install.packages("dplyr")
library("dplyr")

#install.packages("tidyverse")
library("tidyverse")

#install.packages("ggplot2")
library("ggplot2")

genus <- read_tsv("combined_genus.tsv") # Data would have been given as a .tsv from the metagenomic pipeline given by Raymond Kiu.
samples <- read_csv("KELLY_samples.csv")

# gsub is subs the 'pattern' with your 'sub'
samples$`ID Patient` <- gsub("-", "_", samples$`ID Patient`)
samples$`Time point` <- gsub(" ", "_", samples$`Time point`)
samples$`Time point` <- gsub("Eot", "EoT", samples$`Time point`)

# rename columns by piping (from tidyverse)
# Remove the spaces in the column headers
samples <- samples %>%
            rename(ID_Patient = "ID Patient", Stool_sample = "Stool Sample",
                   Time_point = "Time point")

# Select only the columns we want i.e. fractions or numbers
genus_frac <- genus %>% 
          select(name, ends_with("genus_frac"))
genus_num <- genus %>% 
          select(name, ends_with("genus_num"))

# make genus name into row names (remove the: "1,2,3,4,...")
genus_transposed <- column_to_rownames(genus_frac, "name")
# t() simply transposes i.e. swaps rows and columns, but it'll be a matrix
genus_transposed <- t(genus_transposed)
# Need to change to data frame to use downstream
genus_transposed <- as.data.frame(genus_transposed)

# Now switch the rownames to columns, by Stool Sample
# So we can then identify by stool samples ID
genus_transposed <- rownames_to_column(genus_transposed, var = "Stool_sample")

# remove "genus_frac" of the column
genus_transposed$Stool_sample <- gsub("-genus_frac", "", genus_transposed$Stool_sample)

# merge samples and genus frac together
# inner_join only keeps samples that are in both files you're trying to merge
  # full_join does everything
genus_sample_merge <- full_join(samples, genus_transposed, by = "Stool_sample")
genus_sample_merge_inner <- inner_join(samples, genus_transposed, by = "Stool_sample")  

#### Example for using filtering ####
# Use select, to specify what you want to analyse
# filtered_samples <- genus_sample_merge_inner %>% 
#                    filter(Time_point == "Baseline")

# Example:  Multiple selections
# filtered_samples <- genus_sample_merge_inner %>% 
#  filter(Time_point %in% c("Baseline", "Eot"))

# by_patient <- sample_data_long %>%
#  group_by(ID_Patient)

# Save as final csv file (just in case)
# write.csv(genus_sample_merge_inner, "genus_sample_merged.csv")

# Function to make long data
library("reshape2")

# Melt, melts a variable into the data set- making it long data vs wide
genus_sample_long <- melt(genus_sample_merge_inner, 
                              id.vars = c("ID_Patient",
                                          "Stool_sample",
                                          "Site",
                                          "Time_point"), 
                              variable.name = "genus")

#### Rarecurve plot ####
genus <- read_tsv("combined_genus.tsv")
normalised_data <- read_csv("genus_sample_merged.csv")

genus_num <- genus %>% 
  select(name, ends_with("genus_num"))

# make genus name into row names (remove the: "1,2,3,4,...")
data <- column_to_rownames(genus_num, "name")
# t() simply transposes i.e. swaps rows and columns, but it'll be a matrix
data <- t(data)
# Need to change to data frame to use downstream
data <- as.data.frame(data)

# Now switch the rownames to columns, by Stool Sample
# So we can then identify by stool samples ID
data <- rownames_to_column(data, var = "Stool_sample")

# remove "genus_frac" of the column
data$Stool_sample <- gsub("-genus_num", "", data$Stool_sample)

data <- column_to_rownames(data, var = "Stool_sample")

rarecurve(data,
          step = 200,
          sample = 30000,
          col = "blue",
          cex = 0.3,
          label = TRUE,
          #label = FALSE,
          lwd = 0.001,
          ylab = "No. of genus",
          xlab = "No. of sequences"
)

#### Make stacked bar plot ####
#### Make heat map ####
