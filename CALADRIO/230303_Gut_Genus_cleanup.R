#### Housekeeping ####
setwd("~/Dropbox/QIB/KELLY Study/")

#### Download libraries ####
#install.packages("dplyr")
library("dplyr")

#install.packages("tidyverse")
library("tidyverse")

#install.packages("ggplot2")
library("ggplot2")

genus <- read_tsv("Raw data/combined_genus.tsv")
samples <- read_csv("Raw data/KELLY_samples.csv")

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
filtered_samples <- genus_sample_merge_inner %>% 
                    filter(Time_point == "Baseline")

# Example:  Multiple selections
filtered_samples <- genus_sample_merge_inner %>% 
  filter(Time_point %in% c("Baseline", "Eot"))

by_patient <- sample_data_long %>%
  group_by(ID_Patient)

# Save as final csv file (just in case)
write.csv(genus_sample_merge_inner, "genus_sample_merged.csv")

# Function to make long data
library("reshape2")

# Melt, melts a variable into the data set- making it long data vs wide
genus_sample_long <- melt(genus_sample_merge_inner, 
                              id.vars = c("ID_Patient",
                                          "Stool_sample",
                                          "Site",
                                          "Time_point"), 
                              variable.name = "genus")

#### Make stacked bar plot ####
# Taken from 16S metax file
normalised_data <- read.csv("genus_sample_merged.csv")

# Select what I need
data <- dplyr::select(normalised_data,
                      c(Cloacibacillus:ncol(normalised_data)))

# Make the data into percentages- already in %
# data <- data/rowSums(data)*100

# Just check all the rows sum to 100 or 1
rowsums <- rowSums(data)
rowsums <- as.data.frame(rowsums)
rowsums

# Make a dataframe "top" that contains the top to least most abundant genera.
# Taking the total 'sum' of columns
# i.e. per genus which ones have the most reads. 
# Then order by decreasing so that the top are first
top <- data[,order(colSums(data),decreasing=TRUE)]

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

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- data[,order(colSums(data),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(normalised_data)[2] # This needs to be the total number of colums in the dataframe.
# dim(normalised_data)[2], calls the second element of the vector dim(normalised_data)
# this is the number of columns. 
# Select the columns after top 20
taxa_list2 <- colnames(other)[11:N]
N <- length(taxa_list2)
Others <- data.frame(other[,colnames(other) %in% taxa_list2])
Others

# Sum all the other columns into one.
Others <- rowSums(Others)
Others <- as.data.frame(Others)
Others

# Combine the top genus and the other column and check the rows all still add up to 100.
top_other <- cbind(top, Others)
top_other
top_other_rowsums <- rowSums(top_other)
top_other_rowsums <- as.data.frame(top_other_rowsums)
top_other_rowsums

# Select sample_ID
Sample_info <- dplyr::select(normalised_data,
                           c(ID_Patient, Stool_sample, Site, Time_point))

# Combine the data and the sample names back together again.
sample_data <- cbind(Sample_info, top_other)

# Use melt to turn the data from wide format into long format. This puts all the genus data into a single column.
sample_data_long <- melt(sample_data, 
                         id.vars = c("ID_Patient",
                                     "Stool_sample",
                                     "Site",
                                     "Time_point"), 
                         variable.name = "genus")
sample_data_long

# Order the 'time point' so that it's chronological
sample_data_long$Time_point <- factor(sample_data_long$Time_point,
                                      levels = c("Baseline", "Week_9", "EoT"))

# Make a palette of colours for your top genus. There are lots of colours to use in R.
taxa_list # This shows you what those top ones are.
genusPalette <- c( Bacteroides = "cadetblue1",
                   Phocaeicola = "red",
                   Faecalibacterium = "purple",
                   Alistipes = "blue",
                   Prevotella = "green3",
                   Parabacteroides = "yellow",
                   Roseburia = "red4",
                   Lachnospira = "green",
                   Streptococcus = "orange",
                   Mediterraneibacter = "pink",
                   # Bifidobacterium = "darkorchid",
                   # Ruminococcus = "darkgreen",
                   # Klebsiella = "blueviolet",
                   # Blautia = "darkmagenta",
                   # Coprococcus = "chartreuse2",
                   # Escherichia = "darkolivegreen3",
                   # Anaerostipes = "cornflowerblue",
                   # Enterocloster = "brown",
                   # Akkermansia = "darkorange",
                   # Odoribacter = "aquamarine",
                   Others = "grey")


# Use ggplot2 to make the stacked box plot.
ggplot(sample_data_long,
       aes(x = Stool_sample,
           y = value,
           fill = genus)) +
  facet_grid(. ~ Time_point,
             drop = TRUE,
             # Free_x removes empty bars
             scale = "free_x") +
  #Set the width of the bars in the plot
  geom_bar(stat = "identity",
           width = 0.7) +
  #Set the colors to use the genusPalette already created
  scale_fill_manual(values = genusPalette) +
  #Remove extra space at the top and bottom of the plot
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100.1)) +
  #Set the axis labels and title
  labs(title="Ten most abundant genera",
       x = "Stool Sample ID",
       y = "Relative abundance (%)") +
  theme(#Set the title font size
    plot.title = element_text(size=8),
    #Set the legend title position
    legend.position = "bottom",
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
                               hjust=0,
                               vjust=0.5, # "justifying" text
                               size=5),
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

# Save plot
ggsave("KELLY Ten most abundant genus stacked box plots.png", width = 200, height = 200, units = c("mm"), dpi = 400)

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
