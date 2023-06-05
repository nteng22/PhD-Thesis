#### Housekeeping ####
setwd("~/Dropbox/QIB/KELLY Study/initial analysis/")

#### Download libraries ####
#install.packages("dplyr")
library("dplyr")

#install.packages("tidyverse")
library("tidyverse")

#install.packages("ggplot2")
library("ggplot2")

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
# Make a complex heatmap clustered by overall similarity of the genus microbiota using a Bray-Curtis matrix.

#install.packages("ComplexHeatmap")
# if (!requireNamespace("BiocManager", quietly = TRUE))

#install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap")

library("ComplexHeatmap")

# https://www.datacamp.com/community/tutorials/hierarchical-clustering-R
#Import the recoded metadata.
normalised_data <- read_csv("../KELLY clinical metadata data-genus.csv")

Baseline <- (normalised_data %>%
               filter(Time_point == "Baseline"))$ID_Patient

EoT <- (normalised_data %>%
          filter(Time_point == "EoT"))$ID_Patient

test <- intersect(EoT, Baseline)

metadata <- normalised_data %>%
  select(ID_Patient:Time_point)
metadata$CBR <- gsub("1", "yes", metadata$CBR)
metadata$CBR <- gsub("0", "no", metadata$CBR)


data <- normalised_data %>%
  # filter(ID_Patient %in% test) %>%
  select(Cloacibacillus:ncol(normalised_data))

# Make the data into percentages
data <- data/rowSums(data)*100

rowsums <- rowSums(data)
rowsums <- as.data.frame(rowsums)
rowsums

# remove the reads that are less than 0 
coldata <- normalised_data[,colSums(normalised_data != 0) > 0]

#### Only the top x number
top <- data[,order(colSums(data),decreasing=TRUE)]

# Extract list of top N Taxa
N <- 30
# Select the columns, 1:N and put them in a new object
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
# Create a data frame
# with the 'ordered' columns
# and only selecting those present in 'taxa_list', which is top 20.
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

# Just check all the rows sum to 100.
rowsums <- rowSums(data)
rowsums <- as.data.frame(rowsums)
rowsums

# Use the data to calculate the clustering for the heatmap
# Calculate the Bray-Curtis dissimilarity matrix on the sample rows
# Alternative vegdist methods:
# "manhattan", "euclidean", "canberra", "bray", "kulczynski", "jaccard", 
# "gower", "altGower", "morisita", "horn", "mountford", "raup" , "binomial", 
# "chao", "cao" or "mahalanobis".
data_dist_rows <- vegdist(data, method = "bray")
# Cluster the rows by hierarchical clustering
# hclust how to compare the clusters
# hclust alternatives:
# "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA),
# "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
row_clustering <- hclust(data_dist_rows, "average")
# Calculate the Bray-Curtis dissimilarity matrix on the genus columns
# Bray-Curtis dissimalirty matrix on genus columns
data_dist_columns <- vegdist(t(top), method = "bray")
# Cluster the columns by hierarchical clustering
col_clustering <- hclust(data_dist_columns, "average")

# Define the color palette for the heatmap colours
colour_palette <- colorRampPalette(colors=c("white",
                                            "orange",
                                            "darkorange",
                                            "red",
                                            "darkred",
                                            "black"))(100)
# top <- cbind(metadata, top)
# top <- column_to_rownames(top, var = "ID_Patient")

# ~ could change name of patient ID to suggest time point
# or add column

# Define the groups for the colour bar
annot_df1 <- data.frame(Time_point = metadata$Time_point)
annot_df2 <- data.frame(Patient_ID = metadata$ID_Patient)
annot_df3 <- data.frame(CBR = metadata$CBR)

# Define colors for each groups
col1 = list(Time_point = c("Baseline" = "lightgreen", 
                           "Week_9" = "lightslateblue", 
                           "EoT" = "plum3"))
col2 = list(Patient_ID = c("101_06" = "aquamarine",
                           "101_07" = "aquamarine1",
                           "101_08" = "aquamarine2",
                           "101_09" = "aquamarine3",
                           "101_11" = "aquamarine4",
                           "102_04" = "chartreuse",
                           "102_05" = "chartreuse1",
                           "102_06" = "chartreuse2",
                           "102_07" = "chartreuse3",
                           "102_08" = "chartreuse4",
                           "103_01" = "aliceblue",
                           "104_03" = "steelblue",
                           "104_05" = "steelblue2",
                           "104_06" = "steelblue4",
                           "106_06" = "darkorchid",
                           "106_07" = "darkorchid1",
                           "107_03" = "yellow1",
                           "107_04" = "yellow2",
                           "107_05" = "yellow3",
                           "107_06" = "yellow4",
                           "107_07" = "yellowgreen",
                           "108_01" = "mediumturquoise",
                           "109_02" = "olivedrab",
                           "109_03" = "olivedrab1",
                           "109_04" = "olivedrab2",
                           "110_02" = "lightcyan1",
                           "110_03" = "lightcyan2"
))
col3 = list(CBR = c("no" = "red",
                    "yes" = "green"))

# Create the coloured side bar group annotation for the heatmap
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
                                     annotation_legend_param = list(title = "Patient ID", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size
sidebar_annotation3 <- rowAnnotation(df = annot_df3, # Dataframe containing treatment groups
                                     col = col3, # The list of treatment groups and their assigned colours
                                     show_annotation_name = FALSE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "CBR status", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size

# Create the heatmap
heatmap <- Heatmap(as.matrix(top), # The dataframe containing the heatmap data
                   name = "Proportion", # Name is used as the title of the heatmap legend if shown
                   col = colour_palette, # The predefined colour palette
                   cluster_rows = row_clustering, # Cluster the rows using the predefined clustering
                   # cluster_columns = col_clustering, # Cluster the columns using the predefined clustering
                   show_row_names = TRUE, # Show or hide the row names, TRUE to show rownames
                   row_names_gp = gpar(fontsize = 8), # Row name font size
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
# It will add it by the order you give it
p <- heatmap # + sidebar_annotation1 + sidebar_annotation2
p

# Set the saved pdf file name and define the size of the plot
# This isn't ggplot so you need to save using pdf. The width and height are in inches because of Americans.
pdf(file = "~/Dropbox/QIB/KELLY Study/Heatmap_filtered_genus.pdf", width = 11, height = 8)
# print(p) saves the figure into a file
print(p)
dev.off()
