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

setwd("~/Dropbox/QIB/BEAM R Analysis/Shotgun data/")
#install.packages("ComplexHeatmap")
# if (!requireNamespace("BiocManager", quietly = TRUE))

library("ComplexHeatmap")

normalised_data <- read_csv("BEAM_species_normalised.csv", col_names =  TRUE)

# Make the sample_ID into row names.
#normalised_data <- column_to_rownames(normalised_data, var = "sample_ID")
#BrCa <- normalised_data %>%
#  select(Patient, BrCa, Time_point)

# Select abundance data and put it into a new variable
data <- normalised_data %>%
  select(`Synechococcus sp. Minos11`:ncol(normalised_data))
metadata <- dplyr::select(normalised_data,
                          Sample)

# Just check all the rows sum to 100.
rowsums <- rowSums(data)
rowsums <- as.data.frame(rowsums)
rowsums

# Data formatting for top 20
# Make a dataframe "top" that contains the top 10 most abundant genera. You can alter these as you want.
top <- data[,order(colSums(data),decreasing=TRUE)]
#Extract list of top N Taxa
N <- 20
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- data[,order(colSums(data),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(normalised_data)[2] # This needs to be the total number of colums in the dataframe.
# dim(normalised_data)[2], calls the second element of the vector dim(normalised_data)
taxa_list2 <- colnames(other)[21:N]
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

# Use the data to calculate the clustering for the heatmap
# Calculate the Bray-Curtis dissimilarity matrix on the sample rows
# You can change the method e.g. Jacard, see manual for more information
data_dist_rows <- vegdist(top_other, method = "bray")
# Cluster the rows by hierarchical clustering
row_clustering <- hclust(data_dist_rows, "average")
# Calculate the Bray-Curtis dissimilarity matrix on the genus columns
data_dist_columns <- vegdist(t(top_other), method = "bray")
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

# Create the heatmap
heatmap <- Heatmap(as.matrix(top_other), # The dataframe containing the heatmap data
                   name = "Proportion", # Name is used as the title of the heatmap legend if shown
                   col = colour_palette, # The predefined colour palette
                   cluster_rows = row_clustering, # Cluster the rows using the predefined clustering
                   #cluster_columns = col_clustering, # Cluster the columns using the predefined clustering
                   show_row_names = TRUE, # Show or hide the row names, TRUE to show rownames
                   row_names_gp = gpar(fontsize = 6), # Row name font size
                   column_names_gp = gpar(fontsize = 8, # Column name font size
                                          fontface = "italic"), # Column names in italics
                   column_title_gp = gpar(fontsize = 8), # Column title font size
                   row_title = "Samples", # Set row title
                   row_title_gp = gpar(fontsize = 8), # Set row title font size
                   column_title = "Species", # Set column title
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
#pdf(file = "~/Dropbox/QIB/Heatmap of BEAM samples.pdf", width = 8, height = 8)
#pdf(file = "230228_Heatmap of BEAM samples ALL.pdf", width = 8, height = 8)
# print(p) saves the figure into a file
#print(p)
#dev.off()
