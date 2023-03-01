#### Housekeeping ####
setwd("~/Dropbox/QIB/KELLY Study/initial analysis/")

#### Download libraries ####
#install.packages("dplyr")
library("dplyr")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("ggplot2")
library("ggplot2")
#install.packages("reshape2")
library("reshape2")
#install.paclages("vegan")
library("vegan")
#install.packages("rstatix")
library("rstatix")

# Read in the saved csv file from previously
normalised_data <- read.csv("../KELLY clinical metadata data-genus.csv")

abundance_table <- normalised_data %>%
  select(Cloacibacillus:ncol(normalised_data))
rowSums(abundance_table) #check if in percentages
abundance_table <- abundance_table/rowSums(abundance_table)*100

metadata <- normalised_data %>%
    select(ID_Patient:Time_point)
# Filter out the samples that made it to EoT 
# Filter by time point (stats test)
# assumption that they're independent of each other

# Variables for analysis by time points
# Baseline <- (normalised_data %>%
#               filter(Time_point == "Baseline"))$ID_Patient

# EoT <- (normalised_data %>%
#          filter(Time_point == "EoT"))$ID_Patient

# W9 <- (normalised_data %>%
#         filter(Time_point == "Week_9"))$ID_Patient

# full_course <- intersect(intersect(EoT, Baseline), W9)

# Only select the OTU abundance columns
#abundance_table <- data %>%
#  filter(ID_Patient %in% full_course) %>%
  select("Cloacibacillus":ncol(data))

#metadata <- normalised_data %>%
#  filter(ID_Patient %in% full_course) %>%
  select(ID_Patient:Time_point)

metadata$PDL1 <- factor(metadata$PDL1)
metadata$CBR <- factor(metadata$CBR)
metadata$NLR <- factor(metadata$NLR)
metadata$tPFS <- factor(metadata$tPFS)

#Check assumptions for stats
library("dplyr")
#install.package("ggpubr")
library("ggpubr")

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

wilcox.test(shannon_diversity ~ CBR, data = shannon_metadata)
#CBR For all samples W = 418, p = 0.4
#CBR For "full course" samples, W = 196, p = 0.08

#write_csv(shannon_metadata, file = "Shannon_metadata_genus.csv")

# Simpsons index
simpson_diversity <- diversity(abundance_table,
                               index = "invsimpson")
simpson_diversity <- as.data.frame(simpson_diversity)
simpson_metadata <- cbind(metadata, simpson_diversity)

wilcox.test(simpson_diversity ~ CBR, data = simpson_metadata)
# CBR "full course", W = 197, p = 0.07
# CBR all samples, W = 417, p = 0.45

#write_csv(simpson_metadata, file = "Simpson_metadata_genus.csv")

# Plot diversity indexes
# Order the choronological order of time points
shannon_metadata$Time_point <- factor(shannon_metadata$Time_point, 
                                      levels = c("Baseline", "Week_9", "EoT"))
simpson_metadata$Time_point <- factor(simpson_metadata$Time_point,
                                      levels = c("Baseline", "Week_9", "EoT"))
                                      
#Alpha diversity plot 
ggplot(simpson_metadata, aes(x = Time_point, y = simpson_diversity),
       fill = Time_point, shape = Time_point) +
  geom_boxplot(aes(fill = Time_point, shape = Time_point, colour = Time_point)) + #, outlier.shape = NA) +
  geom_point(position = "identity", aes(fill = Time_point, shape = Time_point, colour = Time_point), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("blue",
                              "red",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  #Set the axis labels and title
  labs(title="Simpson Diversity filtered genus"  ) +
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

#### Beta Diversity ####
# Test groupings of data
# Similar to ANOVA. It can tell you potential significant variables
# * is there an interaction between the listed variables
# R2 % effect explained by the variable

# Look at Time point vs. Site
adonis(abundance_table ~ Time_point + Site + Time_point*Site,
       data = metadata, permutations = 999, 
       method = "bray")

# Look at significance of Site
adonis(abundance_table ~ Site, 
       data = metadata, permutations = 999,
       method = "bray")

# Anova of all factors
# not tested independently of each other, so do each variable on its own
adonis(abundance_table ~ tPFS, 
       data = metadata, permutations = 999,
       method = "bray")
# PDL1, p = 0.232
# NLR, p = 0.658 
# CBR, p = 0.269
# RECIST, p = 0.096
# tPFS, p = 0.052

adonis(abundance_table ~ tPFS*RECIST, 
       data = metadata, permutations = 999,
       method = "bray")
# No significant differences for two-variables

# k = 2, specify dimensions
# Look at final stress number, ideally should be below 0.2
# Stress number: roughly the goodness of fit
nmds_calc = metaMDS(abundance_table, distance =  "bray", k = 2, 
                    plot = TRUE,
                    try = 100,
                    engine = "monoMDS")

# adding columns from data set to put it into context
data.scores <- as.data.frame(scores(nmds_calc))
plot.data <- cbind(metadata, data.scores)

plot.data$Time_point <- factor(plot.data$Time_point, levels = c("Baseline", "Week_9", "EoT"))
plot.data$Site <- factor(plot.data$Site) # Site was seen as a string of numbers

# Change class of the metadata
plot.data$PDL1 <- factor(plot.data$PDL1)
plot.data$CBR <- factor(plot.data$CBR)
plot.data$RECIST <- factor(plot.data$RECIST)
plot.data$NLR <- factor(plot.data$NLR)
plot.data$tPFS <- factor(plot.data$tPFS)

# Plot using ggplot2
ggplot(data = plot.data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = Time_point, shape = PDL1, fill = Time_point),
             size = 3) + 
  # scale_"" is used to design the plot
  scale_fill_manual(values = c("white", "white", "white")) + 
  scale_colour_manual(values = c("blue", "red", "darkgreen")) + 
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  #Add label to points
  #geom_text(label=(plot.data$Stool_sample), 
  #          nudge_x = 0.1, nudge_y = 0.1, 
  #          check_overlap = T, size = 3) +
  #scale_shape_manual(values = c(21, 22, 23)) + 
  labs(title = "Baseline_EoT genus PDL1") +
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
