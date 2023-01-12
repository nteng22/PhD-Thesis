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
setwd("C:\\Users/tengn/Dropbox/QIB/BEAM R Analysis"/)

# Import the biom file.
biom <- import_biom("211008_merged_otu_table.biom")
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

################################################################
# Sum the count data into genus categories.
# Example: All the reads classifed as g__[Eubacterium] go into a single row.

# Select the Genus column and count data into a new dataframe.
genus_data_subset <- dplyr::select(taxa_otu_data, c(Genus, "P001D":"P3A"))

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

################################################################
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

# Save this into a excel .csv file.
#write.csv(normalised_data, "BEAM_normalised_data_Genus_data.csv")

################################################################
# Make a stacked box plot

#Import the recoded metadata.
normalised_data <- read.csv("BEAM_normalised_data_Genus_data.csv")
BrCa <- read_csv("BEAM 16S normalised%+metadata.csv") %>%
  select(BrCa)

normalised_data$sample_ID <- factor(normalised_data$sample_ID, 
                                    levels = c("P1A", "P1B", "P001D", 
                                               "P002A", "P002B", "P002D",
                                               "P3A", "P003B", "P003D",
                                               "P004A", "P004B", "P004D",
                                               "P005B", "P005D",
                                               "P006A", "P006B"))
##Select data
data <- dplyr::select(normalised_data,
                      c(Actinomyces:ncol(normalised_data)))

metadata <- normalised_data %>%
  select(sample_ID)

metadata <- cbind(metadata, BrCa)

# Make the data into percentages
data <- data/rowSums(data)*100

# Just check all the rows sum to 100.
rowsums <- rowSums(data)
rowsums <- as.data.frame(rowsums)
rowsums

# Make a dataframe "top" that contains the top 10 most abundant genera. You can alter these as you want.
top <- data[,order(colSums(data),decreasing=TRUE)]
#Extract list of top N Taxa
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- data[,order(colSums(data),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(normalised_data)[2] # This needs to be the total number of colums in the dataframe.
# dim(normalised_data)[2], calls the second element of the vector dim(normalised_data)
taxa_list2 <- colnames(other)[11:N]
N <- length(taxa_list2)
Others <- data.frame(other[,colnames(other) %in% taxa_list2])
Others
# Sum all the other columns into one.
Others <- rowSums(Others)
Others <- as.data.frame(Others)
Others

# Combine the top 10 genus and the other column and check the rows all still add up to 100.
top_other <- cbind(top, Others)
top_other
top_other_rowsums <- rowSums(top_other)
top_other_rowsums <- as.data.frame(top_other_rowsums)
top_other_rowsums

# Select sample_ID
sample_ID <- dplyr::select(normalised_data)

# Combine the data and the sample names back together again.
sample_data <- cbind(metadata, top_other)

# Order the samples by increasing proportion of one genus.
sample_data <- sample_data %>%
  arrange(Faecalibacterium
          #desc(Bacteroides) # Using desc() If you want to arrange in descending order.
  )

# Fix the order of sample IDs in the order of genus proportion.
sample_data$sample_ID <- factor(sample_data$sample_ID,
                                levels=unique(sample_data$sample_ID))


# Use melt to turn the data from wide format into long format. This puts all the genus data into a single column.
sample_data_long <- melt(sample_data, id.vars = c("sample_ID", "BrCa"), variable.name = "Genus")
sample_data_long

# Melt "summarizes" each data point i.e. for sample ID P001D
# Genus is x, value is y. 
# Easier to plot? 

# Make a palette of colours for your top genus. There are lots of colours to use in R.
taxa_list # This shows you what those top ones are.
genusPalette <- c( Bacteroides = "cadetblue1",
                   Blautia = "red",
                   Faecalibacterium = "purple",
                   Coprococcus = "blue",
                   Bifidobacterium = "green3",
                   Ruminococcus = "yellow",
                   Dialister = "red4",
                   Eubacterium = "green",
                   Prevotella = "orange",
                   Sutterella = "pink",
                   
                   Others = "grey")


# Use ggplot2 to make the stacked box plot.
ggplot(sample_data_long,
       aes(x = sample_ID,
           y = value,
           fill = Genus)) +
  #Set the width of the bars in the plot
  geom_bar(stat = "identity",
           width = 0.7) +
  #facet_grid(. ~ BrCa)+
  #Set the colors to use the genusPalette already created
  scale_fill_manual(values = genusPalette) +
  #Remove extra space at the top and bottom of the plot
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100.1)) +
  #Set the axis labels and title
  labs(title="Ten most abundant genera",
       x = "Sample ID",
       y = "Relative abundance (%)") +
  theme(#Set the title font size
    plot.title = element_text(size=8),
    #Set the legend title position
    legend.position = "right",
    #Set the legend title font size
    legend.title = element_text(size=8),
    #Define the size of the legend text and make it italic
    legend.text = element_text(size=8,
                               face = "italic"),
    #Remove the grey background from the legend
    legend.background = element_blank(),
    #Remove the box from the legend
    legend.key = element_blank(),
    #Remove the grey background
    panel.background = element_blank(),
    #Remove the plot border
    panel.border = element_blank(),
    #Remove the major plot grid lines
    panel.grid.major = element_blank(),
    #Remove the minor plot grid lines
    panel.grid.minor = element_blank(),
    #Change orientation of x axis labels
    axis.text.x = element_text(angle=90,
                               hjust=1,
                               vjust=0.5,
                               size=8),
    axis.title.x = element_text(size=8),
    #Define the axis title text size
    axis.title = element_text(size=8),
    #Define the axis label text size
    axis.text.y = element_text(size=8),
    #Add back the x and y axis lines and define thickness (size), line type, and color
    axis.line = element_line(size = 0.5,
                             linetype = "solid",
                             colour = "black"),
    #Set the aspect ratio of the plot
    aspect.ratio = 1
  )


#Save as a pdf for size to go into Inkspace figure
ggsave("BEAM Ten most abundant genus stacked box plots ~ BrCa.pdf", width = 150, height = 150, units = c("mm"), dpi = 300)

##############################################################################
# Make a complex heatmap clustered by overall similarity of the genus microbiota using a Bray-Curtis matrix.

#install.packages("ComplexHeatmap")
# if (!requireNamespace("BiocManager", quietly = TRUE))

#install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")
#library("Cairo")
#install.packages("Cairo")

library("ComplexHeatmap")

#Import the recoded metadata.
normalised_data <- read_csv("BEAM_normalised_data_Genus_data.csv")

# Make the sample_ID into row names.
normalised_data <- column_to_rownames(normalised_data, var = "sample_ID")
BrCa <- read_csv("BEAM 16S normalised%+metadata.csv") %>%
  select(Patient, BrCa, Time_point)

normalised_data <- cbind(BrCa, normalised_data)

# Select abundance data and put it into a new variable
data <- dplyr::select(normalised_data,
                      c(Actinomyces:ncol(normalised_data)))
metadata <- dplyr::select(normalised_data,
                          c(Patient:Time_point))

# Make the data into percentages
data <- data/rowSums(data)*100

# Just check all the rows sum to 100.
rowsums <- rowSums(data)
rowsums <- as.data.frame(rowsums)
rowsums

# Use the data to calculate the clustering for the heatmap
# Calculate the Bray-Curtis dissimilarity matrix on the sample rows
# You can change the method e.g. Jacard, see manual for more information
data_dist_rows <- vegdist(data, method = "bray")
# Cluster the rows by hierarchical clustering
row_clustering <- hclust(data_dist_rows, "average")
# Calculate the Bray-Curtis dissimilarity matrix on the genus columns
data_dist_columns <- vegdist(t(data), method = "bray")
#data_dist_columns <- vegdist(t(top), method = "bray") # this calculates the 'distance' based on the "top"
# N abundances, but you need to have done the stacked bar plot previously so it's saved.
# Cluster the columns by hierarchical clustering
col_clustering <- hclust(data_dist_columns, "average")

# Define the color palette for the heatmap colours
colour_palette <- colorRampPalette(colors=c("white",
                                            "orange",
                                            "darkorange",
                                            "red",
                                            "darkred",
                                            "black"))(100)

# This is to add annotation to the heatmap
# Make a dataframe of the variable you want to have 'annotated' with the heatmap
annot_df1 <- data.frame(Time_point = metadata$Time_point)
annot_df2 <- data.frame(BrCa = metadata$BrCa)

# Make a list, where the levels of the dataframe are allocated a colour
# colour list can be found here: http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
col1 = list(Time_point = c("baseline" = "lightgreen", 
                           "post_surgery" = "lightslateblue",
                           "one_year" = "plum3"))
col2 = list(BrCa = c("lobular" = "hotpink4",
                     "ductal" = "lightpink3",
                     "control" = "indianred1"))

# Add everything together, change the df, col and title accordingly. 
sidebar_annotation1 <- rowAnnotation(df = annot_df1, # Dataframe containing treatment groups
                                     col = col1, # The list of treatment groups and their assigned colours
                                     show_annotation_name = FALSE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Time point", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size
sidebar_annotation2 <- rowAnnotation(df = annot_df2, # Dataframe containing treatment groups
                                     col = col2, # The list of treatment groups and their assigned colours
                                     show_annotation_name = FALSE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "BrCa type", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size

# Create the heatmap
heatmap <- Heatmap(as.matrix(data), # The dataframe containing the heatmap data
                   name = "Proportion", # Name is used as the title of the heatmap legend if shown
                   col = colour_palette, # The predefined colour palette
                   cluster_rows = row_clustering, # Cluster the rows using the predefined clustering
                   # cluster_columns = col_clustering, # Cluster the columns using the predefined clustering
                   show_row_names = TRUE, # Show or hide the row names, TRUE to show rownames
                   row_names_gp = gpar(fontsize = 6), # Row name font size
                   column_names_gp = gpar(fontsize = 8, # Column name font size
                                          fontface = "italic"), # Column names in italics
                   column_title_gp = gpar(fontsize = 8), # Column title font size
                   row_title = "Samples", # Set row title
                   row_title_gp = gpar(fontsize = 8), # Set row title font size
                   column_title = "Genus", # Set column title
                   column_title_side = "bottom", # Set column title font size
                   heatmap_legend_param = list(title = "Relative\nabundance\n(%)", # Set legend title
                                               at = c(0,20,40,60,80,100), # Set legend scale breaks
                                               labels = c("0","20","40","60","80","100"), # Set legend scale labels
                                               title_gp = gpar(fontsize = 8), # Set legend title font size
                                               labels_gp = gpar(fontsize = 8))) # Set legend label font size

# Combine the heatmap and the annotation together in the order in which they are to appear
p <- heatmap # + sidebar_annotation1 + sidebar_annotation2
p

# Set the saved pdf file name and define the size of the plot
# This isn't ggplot so you need to save using pdf. The width and height are in inches because of Americans.
pdf(file = "~/Dropbox/QIB/Heatmap of BEAM samples.pdf", width = 8, height = 8)
# print(p) saves the figure into a file
print(heatmap)
dev.off()
