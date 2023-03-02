#### Housekeeping ####
setwd("~/Dropbox/QIB/KELLY Study/initial analysis/Pasta party/")

#### Download libraries ####
#install.packages("dplyr")
library("dplyr")

#install.packages("tidyverse")
library("tidyverse")

#install.packages("ggplot2")
library("ggplot2")

#install.packages("reshape2")
library("reshape2")

#install.packages("vegan")
library("vegan")

#install.packages("rstatix")
library("rstatix")

# Read in the saved csv file from previously
normalised_data <- read.csv("../../Saliva_data/KELLY_saliva_normalised%_genus+species_with_metadata.csv")

abundance_table <- normalised_data %>%
  select("Desulfovibrio_D168":ncol(normalised_data))

metadata <- normalised_data %>%
  select(patient_ID:NLR)

# is the data normalised
rowSums(abundance_table)

Baseline <- (normalised_data %>%
               filter(Timepoint == "Baseline"))$patient_ID

EoT <- (normalised_data %>%
          filter(Timepoint == "EoT"))$patient_ID

W9 <- (normalised_data %>%
         filter(Timepoint == "Week9"))$patient_ID

metadata$PDL1 <- factor(metadata$PDL1)
metadata$CBR <- factor(metadata$CBR)
metadata$NLR <- factor(metadata$NLR)
metadata$tPFS <- factor(metadata$tPFS)

#Check assumptions for stats
library("dplyr")
#install.package("ggpubr")
library("ggpubr")

# Full course patients
#full_course <- intersect(intersect(Baseline, W9), EoT)

abundance_table <- normalised_data %>%
  #filter(patient_ID %in% full_course) %>%
  select("Desulfovibrio_D168":ncol(normalised_data))

metadata <- normalised_data %>%
  #filter(patient_ID %in% full_course) %>%
  select(patient_ID:NLR)

#### Species Richness ####
# Calculate species number 
# Margin = 1, by row
species_number <- specnumber(abundance_table, MARGIN = 1)
species_number <- as.data.frame(species_number)
species_metadata <- cbind(metadata, species_number)

# Shannon Diversity
shannon_diversity <- diversity(abundance_table,
                               index = "shannon")
shannon_diversity <- as.data.frame(shannon_diversity)
shannon_metadata <- cbind(metadata, shannon_diversity)

#write_csv(shannon_metadata, file = "Shannon_metadata_genus.csv")

# Simpsons index
simpson_diversity <- diversity(abundance_table,
                               index = "invsimpson")
simpson_diversity <- as.data.frame(simpson_diversity)
simpson_metadata <- cbind(metadata, simpson_diversity)

#write_csv(simpson_metadata, file = "Simpson_metadata_genus.csv")

# Plot diversity indexes
# Order the choronological order of time points
shannon_metadata$Timepoint <- factor(shannon_metadata$Timepoint, 
                                      levels = c("Baseline", "Week9", "EoT"))
# rename time point
shannon_metadata$Timepoint <- gsub("Week9", "C4D1", shannon_metadata$Timepoint)
simpson_metadata$Timepoint <- gsub("Week9", "C4D1", simpson_metadata$Timepoint)

simpson_metadata$Timepoint <- factor(simpson_metadata$Timepoint,
                                      levels = c("Baseline", "C4D1", "EoT"))
shannon_metadata$CBR <- factor(shannon_metadata$CBR)
shannon_metadata$CBR <- gsub("1", "CB", shannon_metadata$CBR)
shannon_metadata$CBR <- gsub("0", "No CB", shannon_metadata$CBR)

simpson_metadata$CBR <- factor(simpson_metadata$CBR)
simpson_metadata$CBR <- gsub("1", "yes", simpson_metadata$CBR)
simpson_metadata$CBR <- gsub("0", "no", simpson_metadata$CBR)

shannon_metadata$NLR <- factor(shannon_metadata$NLR)
shannon_metadata$NLR <- gsub("1", "bad", shannon_metadata$NLR)
shannon_metadata$NLR <- gsub("0", "good", shannon_metadata$NLR)

simpson_metadata$NLR <- factor(simpson_metadata$NLR)
simpson_metadata$NLR <- gsub("1", "bad", simpson_metadata$NLR)
simpson_metadata$NLR <- gsub("0", "good", simpson_metadata$NLR)

#Pasta plot 
plot <- ggplot(shannon_metadata, aes(x = Timepoint, y = shannon_diversity)) +
  #facet_grid(. ~ CBR) +
  geom_boxplot(aes(fill = Timepoint, shape = Timepoint, colour = Timepoint)) + #, outlier.shape = NA) +
  geom_line(aes(group = patient_ID),
            linetype = "solid",
            linewidth = 0.05) +
  geom_point(position = "identity", aes(fill = Timepoint, 
                                        shape = Timepoint, 
                                        colour = Timepoint,
                                        group = patient_ID), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("red",
                              "blue",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  #Set the axis labels and title
  labs(title="SALIVA genus Shannon Diversity all samples across Timepoints by CBR",
       y = "Shannon Diversity Index",
       x = "Time point") +
  theme(#Set the title font size
    plot.title = element_text(size=8,
                              hjust = 0.5),
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
    axis.text.x = element_text(angle=0,
                               hjust=0.5,
                               vjust=0, # "justifying" text
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
    aspect.ratio = 1)
# to add annotations to the plot
#annotate("text", x = 1.5, y = 16, label = "p = 1", size = 2.5) +
#annotate("segment", x = 1, xend = 2, y = 15, yend=15, lwd = 0.5)

#ggsave("CALADRIO Saliva genus shannon diversity all samples across Tp~CBR.png", width = 200, height = 200, units = c("mm"), dpi = 400)

# Save by pdf
plot

pdf(file = "~/Dropbox/QIB/KELLY Study/initial analysis/Pasta party/CALADRIO Saliva genus SHANNON.pdf", width = 11, height = 8)

print(plot)
dev.off()
