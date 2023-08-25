#### HOUSEKEEPING ####
setwd("~/Dropbox/QIB/BEAM R Analysis/MiSeq/Final/")
setwd("C:\\Users/tengn/Dropbox/QIB/BEAM R Analysis/")

library("dplyr")
library("tidyverse")
library("vegan")
library("dplyr")
library("phyloseq")

# Import the biom file.
biom <- import_biom("BEAM_Nextseq.biom")
# Check the contents of the file.
biom

# Extract the otu table from the phyloseq object.
otu_data <- otu_table(biom)
# Convert the otu table into a dataframe.
otu_data <- as.data.frame(otu_data)
otu_data

# Extract the taxonomic table from the phyloseq object.
tax_data <- tax_table(biom)
tax_data
# Convert the taxonomic table into a dataframe.
tax_data <- as.data.frame(tax_data)
tax_data

# Rename the columns in the taxonomic data to be more meaningful using dplyr.
tax_data <- tax_data %>%
  dplyr::rename(
    Kingdom = Rank1,
    Phylum = Rank2,
    Class = Rank3,
    Order = Rank4,
    Family = Rank5,
    Genus = Rank6,
    Species = Rank7
  )
# Check the renamed dataframe.
tax_data

# Combine the taxonomic data and the otu data together into one dataframe.
# The rows are all in the same order so this works okay.
taxa_otu_data <- cbind(tax_data, otu_data)
# Check the new combined dataframe.
taxa_otu_data

################################################################
# Sum the count data into genus categories.
# Example: All the reads classifed as g__[Eubacterium] go into a single row.

# Select the Genus column and count data into a new dataframe.
genus_data_subset <- dplyr::select(taxa_otu_data, c(Genus, "JPUH001A":"NNUHP006D"))

# Rename the unclassified data ("g__", NA) to "Unclassified".
genus_data_subset$Genus <- ifelse(genus_data_subset$Genus== "g__", "Unclassified", genus_data_subset$Genus)
genus_data_subset$Genus[is.na(genus_data_subset$Genus)] <- "Unclassified"
# Rename two of the genus names that have square brackets added.
genus_data_subset$Genus <- ifelse(genus_data_subset$Genus== "g__[Ruminococcus]", "Ruminococcus", genus_data_subset$Genus)
genus_data_subset$Genus <- ifelse(genus_data_subset$Genus== "g__[Eubacterium]", "Eubacterium", genus_data_subset$Genus)
genus_data_subset$Genus <- ifelse(genus_data_subset$Genus== "g__[Prevotella]", "Prevotella", genus_data_subset$Genus)
# Replace the "g__" in the Genus names with nothing.
genus_data_subset$Genus <- gsub("g__", "" , genus_data_subset$Genus)

# Make the Genus column into a factor.
genus_data_subset$Genus <- as.factor(genus_data_subset$Genus)

# Sum all the count data by the genus groups.
genus_aggregated <- aggregate(.~Genus, genus_data_subset, sum)
genus_aggregated

# You can save this into a excel .csv file.
write.csv(genus_aggregated, "BEAM_Final Genus otu reads.csv")

################################################################
# Check how many sequences you have in your genus catagory after removing the unclassified OTUs.

# Select every row except the "Unclassified" row. The "!" means those rows are removed.
genus_aggregated <- dplyr::filter(genus_aggregated,
                                  !Genus == "Unclassified")

# Make the Genus names into row names so when the dataframe is all numbers.
raw_data <- genus_aggregated
raw_data <- column_to_rownames(raw_data, var = "Genus")
raw_data <- raw_data %>%
  select(!`...1`)

# Draw a rarecure plot to compare the depth of sequencing of the samples.
# The "t(raw_data)" transposes the dataframe so the genus go from rows to columns within this function.
rarecurve(t(raw_data),
          step = 200,
          sample = 30000,
          col = "blue",
          cex = 0.4,
          label = TRUE,
          #label = FALSE,
          lwd = 0.1,
          ylab = "No. of genus",
          xlab = "No. of sequences"
          #xlim = c(0, 250000)
)

# Find the sample with the smallest number of reads.
col_sums <- colSums(raw_data)
col_sums <- as.data.frame(col_sums)
col_sums
# JPUH002C only 59
# Alternatively for more samples you can see which samples have counts less than a certain number i.e 20000 in this case.
sample_less_than <- dplyr::filter(col_sums, col_sums < 20000)
sample_less_than

raw_data <- raw_data%>%
  select("JPUH001A":"JPUH002B",
         "JPUH003A":ncol(raw_data))

################################################################
# Variance Stabilizing Transformation to normalize the sequence data.

#install.packages("DESeq2")
library("DESeq2")

# Transpose the dataframe.
raw_data_t <- t(raw_data)
# Turn it back into a dataframe.
raw_data_t <- as.data.frame(raw_data_t)
# Make the rownames into a sample_ID column.
raw_data_t <- rownames_to_column(raw_data_t, var = "Sample_ID")

# Select the Sample_ID column into its own dataframe.
coldata <- dplyr::select(raw_data_t,
                         Sample_ID)

# Select all the sequence data columns into their own dataframe.
count_data <- dplyr::select(raw_data_t,
                            Sample_ID,
                            Actinomyces :ncol(raw_data_t))
# Make the Sample_ID column into row names.
count_data <- column_to_rownames(count_data, var = "Sample_ID")
# Add 1 to every cell in the dataframe (Variance Stabilization can't cope with zeros).

count_data <- count_data + 1  
count_data_t <- t(count_data)
countdata <- as.matrix(count_data_t)

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = coldata,
                              design = ~ 1)

# Variance Stabilizing Transformation
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds, fitType = "local")

varianceStabilizingTransformation(dds, blind = TRUE, fitType = "local")

getVarianceStabilizedData(dds)

# Turn it back into counts.
normalised_data <- counts(dds, normalized = TRUE)
# Turn it back into a dataframe.
normalised_data_df <- data.frame(normalised_data)

# Transpose the dataframe.
normalised_data_df_t <- t(normalised_data_df)
# Turn it back into a dataframe.
normalised_data_df_t <- as.data.frame(normalised_data_df_t)

# Make the row names into a column.
normalised_data <- rownames_to_column(normalised_data_df_t, var = "sample_ID")

abundance_table <- normalised_data %>%
  select(Actinomyces:ncol(normalised_data))
abundance_normalised <- abundance_table/rowSums(abundance_table)*100
rowSums(abundance_normalised)

sample_ID <- normalised_data %>%
  select(sample_ID)

final_data <- cbind(sample_ID, abundance_normalised)

write_csv(final_data, "BEAM_final_normalised_data_16S.csv")

#final_data <- read_csv("BEAM_final_normalised_data_16S.csv")
#### METADATA ####
library(readxl)
metadata <- read_xlsx("../../BEAM_metadata.xlsx")
final_data <- read_csv("BEAM_final_normalised_data_16S.csv")

metadata$BrCa <- ifelse(metadata$BrCa== "NEUROENDOCRINE", "Unclassified", metadata$BrCa)
metadata$BrCa[is.na(metadata$BrCa)] <- "Unclassified"
metadata$chemotherapy[is.na(metadata$chemotherapy)] <- "Unclassified"
metadata$age[is.na(metadata$age)] <- "Unclassified"

metadata$radiotherapy
metadata$radiotherapy[is.na(metadata$radiotherapy) == FALSE] <- "yes"
metadata$radiotherapy[is.na(metadata$radiotherapy)] <- "no"

# P002A is actually P001D and P001D is actually P002A. Remove the Ar and Dr
metadata_finaldata <- final_data %>%
  select(sample_ID)
sampleID_separated <- separate(metadata_finaldata, sample_ID,
                     sep = "(?<=[0-9])(?=[A-Z])",
                     into = c("patient", "time_point")) # Need to remove the P in NNUH
sampleID_separated$patient <- gsub("NNUHP00", "NNUH00", sampleID_separated$patient)
final_metadata <- cbind(sampleID_separated, metadata_finaldata)

final_metadata <- final_metadata[-c(37, 39),]
final_metadata$time_point <- gsub("Ar", "A", final_metadata$time_point)
final_metadata$time_point <- gsub("Dr", "D", final_metadata$time_point)
final_metadata$sample_ID <- gsub("NNUHP001Dr", "NNUHP001D", final_metadata$sample_ID)
final_metadata$sample_ID <- gsub("NNUHP002Ar", "NNUHP002A", final_metadata$sample_ID)

metadata_final <- inner_join(final_metadata, metadata, by = "patient")

abundance_metaadata_final <- inner_join(metadata_final, final_data, 
                                        by = "sample_ID")
abundance_metaadata_final$sample_ID <- gsub("NNUHP00", "NNUH00", abundance_metaadata_final$sample_ID)

abundance_metaadata_final$time_point <- gsub("D", "One-Year", abundance_metaadata_final$time_point)
abundance_metaadata_final$time_point <- gsub("C", "6-Months", abundance_metaadata_final$time_point)
abundance_metaadata_final$time_point <- gsub("B", "Post-Surgery", abundance_metaadata_final$time_point)
abundance_metaadata_final$time_point <- gsub("A", "Baseline", abundance_metaadata_final$time_point)

write_csv(abundance_metaadata_final, "BEAM_final_metadata_normalised_16S.csv")
