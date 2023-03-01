setwd("~/Dropbox/QIB/KELLY Study/AEs/")

library("tidyverse")
library("ggplot2")
library("ComplexHeatmap")
library("readxl")
library("dplyr")
library("vegan")

####
data <- read_xlsx("20210211_Listado_AEs_KELLY.xlsx")

CTCAE <- data %>%
  filter(`CTCAE Grade` == "Grade 3" |`CTCAE Grade` == "Grade 4" ) #| `CTCAE Grade` == "Grade 2")

CTCAE <- dplyr::rename(CTCAE, "ID_Patient" = "Patient")
CTCAE$ID_Patient <- gsub("-", "_", CTCAE$ID_Patient)

CTCAE_related <- CTCAE %>%
  filter(`Relation Eribulin` == "Possibly Related" | `Relation Pembrolizumab` == "Possibly Related")

# Clean up data :')
CTCAE_related$AE <- gsub("neutrophyl count decreased",
                         "Neutrophil count decreased",
                         CTCAE_related$AE)
CTCAE_related$AE <- gsub("NEUTROPHIL COUNT DECREASED",
                         "Neutrophil count decreased",
                         CTCAE_related$AE)
CTCAE_related$AE <- gsub("NEUTROPENIA",
                         "Neutropenia",
                         CTCAE_related$AE)
CTCAE_related$AE <- gsub("neutropenia",
                         "Neutropenia",
                         CTCAE_related$AE)
CTCAE_related$AE <- gsub("HYPERTRIGLYCERIDEMIA",
                         "hypertriglyceridemia",
                         CTCAE_related$AE)
CTCAE_related$AE <- gsub("NEUROTOXICITY",
                         "Neurotoxicity",
                         CTCAE_related$AE)
CTCAE_related$AE <- gsub("alopecia",
                         "Alopecia",
                         CTCAE_related$AE)
CTCAE_related$AE <- gsub("ALOPECIA",
                         "Alopecia",
                         CTCAE_related$AE)
CTCAE_related$AE <- gsub("ANEMIA",
                         "Anemia",
                         CTCAE_related$AE)
CTCAE_related$AE <- gsub("anemia",
                         "Anemia",
                         CTCAE_related$AE)
CTCAE_related <- CTCAE_related %>%
  select(ID_Patient, AE, `CTCAE Grade`)

### Merge with genus data ###
MCB_data <- read_csv("../Saliva_data/KELLY_saliva_normalised%_with_metadata.csv")
MCB_data <- dplyr::rename(MCB_data, "ID_Patient" = "patient_ID")
MCB_data$ID_Patient <- gsub("-", "_", MCB_data$ID_Patient)

merged_data <- inner_join(CTCAE_related, MCB_data, by = "ID_Patient")

#### HEATMAP
reads <- merged_data %>%
  select("5-7N15":ncol(merged_data))

metadata <- merged_data %>%
  select(ID_Patient:NLR)

rowsums <- rowSums(reads)
rowsums <- as.data.frame(rowsums)
rowsums # Data is normalised, no need for the next step

normalised_data <- reads
#Normalize reads
# normalised_data <- reads/rowSums(reads)*100

# remove the reads that are less than 0 
coldata <- normalised_data[,colSums(normalised_data != 0) > 0] # none

#### Only the top x number
top <- normalised_data[,order(colSums(normalised_data),decreasing=TRUE)]

# Extract list of top N Taxa
N <- 10
# Select the columns, 1:N and put them in a new object
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)

# Create a data frame
# with the 'ordered' columns
# and only selecting those present in 'taxa_list', which is top 20.
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

data_dist_rows <- vegdist(top, method = "bray")
row_clustering <- hclust(data_dist_rows, "average")
data_dist_columns <- vegdist(t(top), method = "bray")
col_clustering <- hclust(data_dist_columns, "average")
colour_palette <- colorRampPalette(colors=c("white",
                                            "orange",
                                            "darkorange",
                                            "red",
                                            "darkred",
                                            "black"))(100)

annot_df1 <- data.frame("CTCAE Grade" = metadata$`CTCAE Grade`)
col1 = list(CTCAE.Grade= c("Grade 2" = "darkgreen",
                           "Grade 3" = "darkorange",
                           "Grade 4" = "darkred"))
sidebar_annotation1 <- rowAnnotation(df = annot_df1, # Dataframe containing treatment groups
                                     col = col1, # The list of  groups and their assigned colours
                                     show_annotation_name = FALSE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "CTCAE Grade", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size

heatmap <- Heatmap(as.matrix(top), # The dataframe containing the heatmap data
                   name = "Proportion", # Name is used as the title of the heatmap legend if shown
                   col = colour_palette, # The predefined colour palette
                   cluster_rows = row_clustering, # Cluster the rows using the predefined clustering
                   cluster_columns = col_clustering, # Cluster the columns using the predefined clustering
                   show_row_names = TRUE, # Show or hide the row names, TRUE to show rownames
                   row_names_gp = gpar(fontsize = 8), # Row name font size
                   column_names_gp = gpar(fontsize = 8, # Column name font size
                                          fontface = "italic"), # Column names in italics
                   column_title_gp = gpar(fontsize = 8), # Column title font size
                   row_title = "Samples by AE", # Set row title
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
p <- heatmap + sidebar_annotation1
p

# Set the saved pdf file name and define the size of the plot
# This isn't ggplot so you need to save using pdf. The width and height are in inches because of Americans.
pdf(file = "~/Dropbox/QIB/KELLY Study/AEs/Heatmap_AEs_saliva.pdf", width = 11, height = 8)
# print(p) saves the figure into a file
print(p)
dev.off()
