#install.packages("dplyr")
library("dplyr")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("vegan")
library("vegan")
#install.packages("reshape2")
library("reshape2")
library("ggplot2")

#### Data cleanup genus####
# Note 10/01/222: Add metadata order by time point of samples.
setwd("~/Dropbox/QIB/KELLY Study/Saliva_data/Metadata plots/")
normalised_data <- read_csv("../KELLY_saliva_normalised_data_Genus_data.csv")
normalised_data <- normalised_data %>%
  select(sample_ID : ncol(normalised_data))

metadata <- read_csv("../KELLY_Saliva_metadata.csv")
metadata <- metadata %>%
  dplyr::rename(patient_ID = PatientID)
metadata <- metadata %>%
  dplyr::rename(sample_ID = SampleID)

#Change tPFS to less than 6 months and after 6 months
# 1 tPFS < 6
# 0 tPFS > 6
metadata$tPFS <- ifelse(metadata$tPFS > 6.0000 , "1", metadata$tPFS)
metadata$tPFS <- ifelse(!metadata$tPFS == 1 , "0", metadata$tPFS)

# need [] as . alone is a special character
metadata$PDL1 <- gsub("[.]", "2", metadata$PDL1)
metadata$PDL1[metadata$PDL1 == 2] <- NA
metadata <- drop_na(metadata)

count_data <- normalised_data %>%
  select("5-7N15":ncol(normalised_data)) - 1
patient <- normalised_data %>%
  select(sample_ID)

# Make the data into percentages
count_data <- count_data/rowSums(count_data)*100

normalised_data <- cbind(patient, count_data)

data <- inner_join(metadata, normalised_data, var = sample_ID)

#write_csv(data, "../KELLY_saliva_normalised%_with_metadata.csv")

#### Data selection genus ####
data <- read_csv("../KELLY_saliva_normalised%_with_metadata.csv") # %>%
  # filter(Timepoint == "EoT")

metadata <- data %>%
  select(patient_ID:NLR)

metadata$PDL1 <- factor(metadata$PDL1)
metadata$CBR <- factor(metadata$CBR)
metadata$NLR <- factor(metadata$NLR)
metadata$tPFS <- factor(metadata$tPFS)

metadata$CBR <- gsub("1", "CB", metadata$CBR)
metadata$CBR <- gsub("0", "no CB", metadata$CBR)

abundance_table <- data %>%
  select(`5-7N15`: ncol(data))

#### Alpha diversity calcs genus ####
species_number <- specnumber(abundance_table, MARGIN = 1)
species_number <- as.data.frame(species_number)
species_metadata <- cbind(metadata, species_number)

shannon_diversity <- diversity(abundance_table,
                               index = "shannon")
shannon_diversity <- as.data.frame(shannon_diversity)
shannon_metadata <- cbind(metadata, shannon_diversity)

shannon_metadata$Timepoint <- factor(shannon_metadata$Timepoint,
                                     levels = c("Baseline", "Week9", "EoT"))

simpson_diversity <- diversity(abundance_table,
                               index = "invsimpson")
simpson_diversity <- as.data.frame(simpson_diversity)
simpson_metadata <- cbind(metadata, simpson_diversity)
simpson_metadata$Timepoint <- factor(simpson_metadata$Timepoint,
                                     levels = c("Baseline", "Week9", "EoT"))

#### Alpha diversity plots genus ####
ggplot(shannon_metadata, aes(x = CBR, y = shannon_diversity),
       fill = CBR, shape = CBR) +
  geom_boxplot(aes(fill = CBR, shape = CBR, colour = CBR), outlier.shape = NA) +
  #Visualise by time point
  #facet_grid(. ~ Timepoint,
  #           drop = TRUE, # Free_x removes empty bars
  #           scale = "free_x") +
  geom_point(position = "identity", aes(fill = CBR, shape = CBR, colour = CBR), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("blue",
                              "red",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  #Set the axis labels and title
  labs(title="CALADRIO saliva shannon genus CB") +
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

ggsave("KELLY_saliva Shannon CBR genus.png", width = 200, height = 200, units = c("mm"), dpi = 400)

ggplot(simpson_metadata, aes(x = CBR, y = simpson_diversity),
       fill = CBR, shape = CBR) +
  geom_boxplot(aes(fill = CBR, shape = CBR, colour = CBR), outlier.shape = NA) +
  #facet_grid(. ~ CBR,
  #           drop = TRUE, # Free_x removes empty bars
  #           scale = "free_x") +
  geom_point(position = "identity", aes(fill = CBR, shape = CBR, colour = CBR), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("blue",
                              "red",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  #Set the axis labels and title
  labs(title="KELLY saliva Simpson genus CB") +
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

ggsave("KELLY_saliva Simpson genus CB.png", width = 200, height = 200, units = c("mm"), dpi = 400)

ggplot(simpson_metadata, aes(x = NLR, y = simpson_diversity),
       fill = NLR, shape = NLR) +
  geom_boxplot(aes(fill = NLR, shape = NLR, colour = NLR), outlier.shape = NA) +
  geom_point(position = "identity", aes(fill = NLR, shape = NLR, colour = NLR), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("blue",
                              "red",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  #Set the axis labels and title
  labs(title="KELLY saliva Simpson genus EoT") +
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

ggsave("KELLY_saliva Simpson NLR genus EoT.png", width = 200, height = 200, units = c("mm"), dpi = 400)

#### Adonis Genus ####
adonis2(abundance_table ~ CBR*NLR,
       data = metadata, permutations = 999,
       method = "bray")

multi_all <-adonis2(abundance_table ~ CBR + NLR + PDL1 + tPFS,
                    data = metadata, permutations = 999,
                    method = "bray")
multi_all

CB_NLR <- adonis2(abundance_table ~ CBR*NLR,
                  data = metadata, permutations = 999,
                  method = "bray")
CB_NLR
# CBR=0.2, NLR=0.45, PDL1 = 0.09, tPFS=0.52
# CBR*PDL1=PDL1 0.84, CBR:NLR=0.19

#### NMDS Plot Genus ####
# k = 2, specify dimensions
# Look at final stress number, ideally should be below 0.2
# Stress number: roughly the goodness of fit
nmds_calc = metaMDS(abundance_table, distance =  "bray", k = 2, 
                    plot = TRUE,
                    try = 100,
                    engine = "monoMDS")

# adding columns from data set to put it into context
data.scores <- as.data.frame(scores(nmds_calc))
data.scores = as.data.frame(scores(nmds_calc)$sites)

plot.data <- cbind(metadata, data.scores)

plot.data$Timepoint <- factor(plot.data$Timepoint, 
                              levels = c("Baseline", "Week9", "EoT"))

# Change class of the metadata
plot.data$PDL1 <- factor(plot.data$PDL1)
plot.data$CBR <- factor(plot.data$CBR)
plot.data$NLR <- factor(plot.data$NLR)
# plot.data$tPFS <- factor(plot.data$tPFS)

plot.data$CBR <- gsub("1", "CB", plot.data$CBR)
plot.data$CBR <- gsub("0", "no CB", plot.data$CBR)
plot.data <- plot.data %>%
  rename("CB Status" = "CBR")

plot.data$PDL1 <- gsub("1", "present", plot.data$PDL1)
plot.data$PDL1 <- gsub("0", "absent", plot.data$PDL1)

# Plot using ggplot2
ggplot(data = plot.data, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(aes(colour = PDL1, shape = PDL1, fill = PDL1),
             size = 2) + 
  stat_ellipse(aes(colour = PDL1)) +
  # scale_"" is used to design the plot
  scale_fill_manual(values = c("white", "white", "white")) + 
  scale_colour_manual(values = c("blue", "red", "darkgreen")) + 
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)) +
  #Add label to points
  #geom_text(label=(plot.data$Stool_sample), 
  #          nudge_x = 0.1, nudge_y = 0.1, 
  #          check_overlap = T, size = 3) +
  #scale_shape_manual(values = c(21, 22, 23)) + 
  labs(title = "KELLY Saliva microbiota Genus level by PDL1") +
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

ggsave("KELLY saliva microbiota NMDS genus by PDL1.png", width = 200, height = 200, units = c("mm"), dpi = 400)

