#install.packages("phyloseq")
library("phyloseq")
#install.packages("dplyr")
library("dplyr")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("vegan")
library("vegan")
#install.packages("reshape2")
library("reshape2")

# Set the file location
setwd("~/Dropbox/QIB/BEAM R Analysis/")
setwd("C:\\Users/tengn/Dropbox/QIB/BEAM R Analysis/")

#### Make the biom file into a readable format ####
# Import the biom file.
biom <- import_biom("230111_merged_otu_table.biom")
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


# You can save this into a excel .csv file.
#write.csv(taxa_otu_data, "BEAM Combined taxanomic names and subsampled otu reads.csv")

#### Clean up file ####
# Sum the count data into genus categories.
# Example: All the reads classifed as g__[Eubacterium] go into a single row.

# Select the Genus column and count data into a new dataframe.
genus_data_subset <- dplyr::select(taxa_otu_data, c(Genus, "P001A":"P006D"))

# Rename the unclassified data ("g__", NA) to "Unclassified".
genus_data_subset$Genus <- ifelse(genus_data_subset$Genus== "g__", "Unclassified", genus_data_subset$Genus)
genus_data_subset$Genus[is.na(genus_data_subset$Genus)] <- "Unclassified"
# Rename two of the genus names that have square brackets added.
genus_data_subset$Genus <- ifelse(genus_data_subset$Genus== "g__[Ruminococcus]", "Ruminococcus", genus_data_subset$Genus)
genus_data_subset$Genus <- ifelse(genus_data_subset$Genus== "g__[Eubacterium]", "Eubacterium", genus_data_subset$Genus)
# Replace the "g__" in the Genus names with nothing.
genus_data_subset$Genus <- gsub("g__", "" , genus_data_subset$Genus)

# Make the Genus column into a factor.
genus_data_subset$Genus <- as.factor(genus_data_subset$Genus)

# Sum all the count data by the genus groups.
genus_aggregated <- aggregate(.~Genus, genus_data_subset, sum)
genus_aggregated

# You can save this into a excel .csv file.
#write.csv(genus_aggregated, "BEAM Genus otu reads.csv")

# Check how many sequences you have in your genus catagory after removing the unclassified OTUs.

# Select every row except the "Unclassified" row. The "!" means those rows are removed.
genus_aggregated <- dplyr::filter(genus_aggregated,
                                  !Genus == "Unclassified")

# Make the Genus names into row names so when the dataframe is all numbers.
raw_data <- genus_aggregated
raw_data <- column_to_rownames(raw_data, var = "Genus")

# Draw a rarecure plot to compare the depth of sequencing of the samples.
# The "t(raw_data)" transposes the dataframe so the genus go from rows to columns within this function.
rarecurve(t(raw_data),
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

# Find the sample with the smallest number of reads.
col_sums <- colSums(raw_data)
col_sums <- as.data.frame(col_sums)
col_sums
# samples P005B has 17183, which is a good high number.

# Alternatively for more samples you can see which samples have counts less than a certain number i.e 20000 in this case.
sample_less_than <- dplyr::filter(col_sums, col_sums < 20000)
sample_less_than
# 5 samples less than, P1A P1B P2D P3A P5B

#### Normalise the sequence data ####
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

# Save this into a excel .csv file.
#write.csv(normalised_data, "BEAM_normalised_data_Genus_data.csv")
