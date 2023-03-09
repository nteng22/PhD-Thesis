setwd("~/Dropbox/QIB/BEAM R Analysis/Shotgun data/")

library(dplyr)
library(tidyverse)
library(ggplot2)
library(vegan)

data <- read_csv("BEAM_species_normalised.csv", col_names = TRUE)

# Select only abundance data
abundance_table <- data %>%
  select("Synechococcus sp. Minos11":ncol(data))
rowSums(abundance_table)  # Check is 100%

# Select and format metadata
metadata <- data %>%
  select(Sample)

metadata <- separate(metadata, Sample,
                     sep = "(?<=[0-9])(?=[A-Z])",
                     into = c("PatientID", "time_point"))
metadata$time_point <- gsub("A", "baseline", metadata$time_point)
metadata$time_point <- gsub("B", "post-surgery", metadata$time_point)

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

wilcox.test(shannon_diversity ~ time_point, data = shannon_metadata)
#write_csv(shannon_metadata, file = "Shannon_metadata_genus.csv")

# Simpsons index
simpson_diversity <- diversity(abundance_table,
                               index = "invsimpson")
simpson_diversity <- as.data.frame(simpson_diversity)
simpson_metadata <- cbind(metadata, simpson_diversity)

wilcox.test(simpson_diversity ~ time_point, data = simpson_metadata)
#write_csv(simpson_metadata, file = "Simpson_metadata_genus.csv")

# Plot diversity indexes
# Order the choronological order of time points
shannon_metadata$Time_point <- factor(shannon_metadata$time_point, 
                                      levels = c("Baseline", "post-surgery"))

#Alpha diversity plot 
ggplot(shannon_metadata, aes(x = time_point, y = shannon_diversity),
       fill = time_point, shape = time_point) +
  geom_boxplot(aes(fill = time_point, shape = time_point, colour = time_point)) + #, outlier.shape = NA) +
  geom_point(position = "identity", aes(fill = time_point, shape = time_point, colour = time_point), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("blue",
                              "red",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  #Set the axis labels and title
  labs(title="Simpson Diversity filtered species"  ) +
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

#ggsave("KELLY Simpson Diversity filtered baseline_EoT genus.png", width = 200, height = 200, units = c("mm"), dpi = 400)

# Beta diversity
adonis2(abundance_table ~ time_point,
       data = metadata, permutations = 999, 
       method = "bray")

nmds_calc = metaMDS(abundance_table, distance =  "bray", k = 2, 
                    plot = TRUE,
                    try = 100,
                    engine = "monoMDS")

# adding columns from data set to put it into context
data.scores = as.data.frame(scores(nmds_calc)$sites)

plot.data <- cbind(metadata, data.scores)

# Plot using ggplot2
ggplot(data = plot.data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = time_point, shape = time_point, fill = time_point),
             size = 3) + 
  stat_ellipse(aes(colour = time_point)) +
    # scale_"" is used to design the plot
  scale_fill_manual(values = c("white", "white", "white")) + 
  scale_colour_manual(values = c("blue", "red", "darkgreen")) + 
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  #Add label to points
  geom_text(label=(plot.data$PatientID), 
            nudge_x = 0.01, nudge_y = 0.02, 
            check_overlap = T, size = 2) +
  #scale_shape_manual(values = c(21, 22, 23)) + 
  labs(title = "Species NMDS BEAM") +
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

#ggsave("KELLY NMDS Baseline_EoT genus.png", width = 200, height = 200, units = c("mm"), dpi = 400)
