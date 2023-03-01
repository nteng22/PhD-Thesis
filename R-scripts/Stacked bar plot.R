#### Stacked bar plot using 16S rRNA gene amplicon sequencing reads ####
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

#Import the recoded metadata.
normalised_data <- read.csv("BEAM_normalised_data_Genus_data.csv")

normalised_data$sample_ID <- factor(normalised_data$sample_ID,    # To order the samples on the plot
                                    levels = c("P001A", "P001B", "P001D", 
                                               "P002A", "P002B", "P002D",
                                               "P003A", "P003B", "P003D",
                                               "P004A", "P004B", "P004D",
                                               "P005B", "P005D",
                                               "P006A", "P006B", "P006C", "P006D"))
##Select data
data <- dplyr::select(normalised_data,
                      c(Actinomyces:ncol(normalised_data)))

metadata <- normalised_data %>%
  select(sample_ID)

#metadata <- cbind(metadata, data)

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

# OPTIONAL: Order the samples by increasing proportion of one genus.
sample_data <- sample_data %>%
  arrange(Faecalibacterium
          #desc(Bacteroides) # Using desc() If you want to arrange in descending order.
  )

# OPTIONAL: Fix the order of sample IDs in the order of genus proportion.
sample_data$sample_ID <- factor(sample_data$sample_ID,
                                levels=unique(sample_data$sample_ID))


# Use melt to turn the data from wide format into long format. This puts all the genus data into a single column.
sample_data_long <- melt(sample_data, id.vars = c("sample_ID"), variable.name = "Genus")
sample_data_long

# Melt "summarizes" each data point i.e. for sample ID P001D

# Make a palette of colours for your top genus. There are lots of colours to use in R.
taxa_list # This shows you what those top ones are.
genusPalette <- c( Faecalibacterium = "cadetblue1",
                   Bacteroides = "red",
                   Blautia = "purple",
                   Coprococcus = "blue",
                   Bifidobacterium = "green3",
                   Ruminococcus = "yellow",
                   Dialister = "red4",
                   Eubacterium = "green",
                   Prevotella = "orange",
                   Collinsella = "pink",
                   
                   Others = "grey")


# Use ggplot2 to make the stacked box plot.
ggplot(sample_data_long,
       aes(x = sample_ID,
           y = value,
           fill = Genus)) +
  #Set the width of the bars in the plot
  geom_bar(stat = "identity",
           width = 0.7) +
  #facet_grid(. ~ BrCa)+ # Choose to show per category of metadata
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
#ggsave("BEAM Ten most abundant genus stacked box plots.pdf", width = 150, height = 150, units = c("mm"), dpi = 300)
