#install.packages("dplyr")
library("dplyr")
#install.packages("tidyverse")
library("tidyverse")
#install.packages("vegan")
library("vegan")
#install.packages("reshape2")
library("reshape2")
library("ggplot2")
setwd("~/Dropbox/QIB/KELLY Study/Saliva_data/Metadata plots/")

#### Data cleanup species####
normalised_data <- read_csv("../KELLY_saliva_normalised_data_Genus+Species_data.csv")

data <- normalised_data %>%
  select("Desulfovibrio_D168":ncol(normalised_data))
patient <- normalised_data %>%
  select(sample_ID)

count_data <- data %>%
  select("Desulfovibrio_D168":ncol(data)) - 1
# Data has been rarefied which is why the numbers are so high

# Make the data into percentages
count_data <- count_data/rowSums(count_data)*100

# check if correct
rowsums <- rowSums(count_data)
rowsums <- as.data.frame(rowsums)
rowsums

# Combine metadata with abundances
metadata <- read_csv("..//KELLY_Saliva_metadata.csv")
metadata <- metadata %>% # Rename to merge
  dplyr::rename(patient_ID = PatientID)
metadata <- metadata %>% # Rename to merge
  dplyr::rename(sample_ID = SampleID)

# need [] as . alone is a special character
metadata$PDL1 <- gsub("[.]", "2", metadata$PDL1) #No metadata for one sample
metadata$PDL1[metadata$PDL1 == 2] <- NA
metadata <- drop_na(metadata)

normalised_data <- cbind(patient, count_data) # Combine patient ID with abundances

data <- inner_join(metadata, normalised_data, var = sample_ID)

#write_csv(data, "KELLY_saliva_normalised%_species_with_metadata.csv")

#### Data selection species ####
data <- read_csv("..//KELLY_saliva_normalised%_species_with_metadata.csv") #%>%
# filter(Timepoint == "EoT")

metadata <- data %>%
  select(patient_ID:NLR)

metadata$PDL1 <- factor(metadata$PDL1)
metadata$CBR <- factor(metadata$CBR)
metadata$NLR <- factor(metadata$NLR)
metadata$tPFS <- factor(metadata$tPFS)

metadata$CBR <- gsub("1", "CB", metadata$CBR)
metadata$CBR <- gsub("0", "no CB", metadata$CBR)
metadata <- metadata %>%
  rename("CB Status"="CBR")

metadata$PDL1 <- gsub("1", "present", metadata$PDL1)
metadata$PDL1 <- gsub("0", "absent", metadata$PDL1)

abundance_table <- data %>%
  select("Desulfovibrio_D168":ncol(data))

#### Alpha diversity calcs species ####
species_number <- specnumber(abundance_table, MARGIN = 1)
species_number <- as.data.frame(species_number)
species_metadata <- cbind(metadata, species_number)

shannon_diversity <- diversity(abundance_table,
                               index = "shannon")
shannon_diversity <- as.data.frame(shannon_diversity)
shannon_metadata <- cbind(metadata, shannon_diversity)

simpson_diversity <- diversity(abundance_table,
                               index = "invsimpson")
simpson_diversity <- as.data.frame(simpson_diversity)
simpson_metadata <- cbind(metadata, simpson_diversity)

#### Alpha diversity plots species ####
ggplot(shannon_metadata, aes(x = PDL1, y = shannon_diversity),
       fill = PDL1, shape = PDL1) +
  geom_boxplot(aes(fill = PDL1, shape = PDL1, colour = PDL1), outlier.shape = NA) +
  geom_point(position = "identity", aes(fill = PDL1, shape = PDL1, colour = PDL1), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("blue",
                              "red",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  #Set the axis labels and title
  labs(title="KELLY Saliva Shannon Index on species by CB ") +
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

ggsave("KELLY saliva Shannon PDL1 species.png", width = 200, height = 200, units = c("mm"), dpi = 400)

ggplot(simpson_metadata, aes(x = `CB Status`, y = simpson_diversity),
       fill = `CB Status`, shape = `CB Status`) +
  geom_boxplot(aes(fill = `CB Status`, shape = `CB Status`, colour = `CB Status`), outlier.shape = NA) +
  geom_point(position = "identity", aes(fill = `CB Status`, shape = `CB Status`, colour = `CB Status`), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("blue",
                              "red",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  #Set the axis labels and title
  labs(title="KELLY saliva Simpson species by CB") +
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

ggsave("KELLY_saliva Simpson CB species.png", width = 200, height = 200, units = c("mm"), dpi = 400)

#### Adonis species ####
adonis2(abundance_table ~ PDL1,
       data = metadata, permutations = 999, 
       method = "bray")

# CBR and NLR at ~0.3, not significant but could be due to small number
# PDL1 is significant?? p < 0.05 border line p < 0.01

Shannon_PDL1 <- wilcox.test(shannon_diversity ~ PDL1, data = shannon_metadata)
Shannon_PDL1

Simpson_PDL1 <- wilcox.test(simpson_diversity ~ PDL1, data = simpson_metadata)
Simpson_PDL1

#### NMDS Plot Species ####
# k = 2, specify dimensions
# Look at final stress number, ideally should be below 0.2
# Stress number: roughly the goodness of fit
nmds_calc = metaMDS(abundance_table, distance =  "bray", k = 2, 
                    plot = TRUE,
                    try = 150,
                    engine = "monoMDS")

# The stress number is greater than 0.2. Not ideal. 
# adding columns from data set to put it into context
#data.scores <- as.data.frame(scores(nmds_calc))
data.scores = as.data.frame(scores(nmds_calc)$sites)

plot.data <- cbind(metadata, data.scores)

plot.data$Time_point <- factor(plot.data$Timepoint, levels = c("Baseline", "Week9", "EoT"))

# Change class of the metadata
plot.data$PDL1 <- factor(plot.data$PDL1)
plot.data$`CB Status` <- factor(plot.data$`CB Status`)
plot.data$NLR <- factor(plot.data$NLR)
# plot.data$tPFS <- factor(plot.data$tPFS)

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
  labs(title = "KELLY Saliva species by PDL1") +
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

ggsave("KELLY Saliva microbiota NMDS species by PDL1.png", width = 200, height = 200, units = c("mm"), dpi = 400)
