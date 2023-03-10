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

#install.packages("ggalluvial")
library("ggalluvial")
#### Data ####
data <- read_csv("BEAM_genus_normalised.csv")

abundance_table <- data %>% # Select ID and abundances
  select(Stool_sample,
         Cloacibacillus:ncol(data)
  )

abundance_table <- column_to_rownames(abundance_table, var = "Stool_sample")

# Remove any columns with a sum total of 0.
abundance_table <- abundance_table %>%
  select(where(is.numeric)) %>% 
  select(where(~sum(.) >0))

# check is normalised
rowsums <- rowSums(abundance_table)
rowsums <- as.data.frame(rowsums)
rowsums
# it's not so normalise
abundance_table <- abundance_table/rowSums(abundance_table)*100

metadata <- data %>%
  select(Stool_sample,
         CBR,
         Time_point)

metadata$CBR <- gsub("1", "CB", metadata$CBR)
metadata$CBR <- gsub("0", "no CB", metadata$CBR)
metadata$Time_point <- gsub("Week_9", "C4D1", metadata$Time_point)

top <- abundance_table[,order(colSums(abundance_table),decreasing=TRUE)]

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
other <- abundance_table[,order(colSums(abundance_table),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(abundance_table)[2] # This needs to be the total number of colunms in the dataframe.
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

#Bind dataset together
Data_20 <- cbind(metadata, top_other)

# Long data
Data_long <- melt(Data_20, 
                  id.vars = c("Stool_sample",
                              "Time_point",
                              "CBR"),
                  #"PDL1",
                  #"tPFS",
                  #"NLR"), 
                  variable.name = "Genus")

Data_long$Time_point <- factor(Data_long$Time_point,
                               levels = c("Baseline",
                                          "C4D1",
                                          "EoT"))

# alluvial needs one value per level
Data_long_summarised <- Data_long %>%
  group_by(Time_point,CBR,Genus) %>%
  summarize(value=mean(value))

Data_long_summarised$CBR <- as.factor(Data_long_summarised$CBR)

Data_long_summarised$Time_point <- factor(Data_long_summarised$Time_point,
                                          levels = c("Baseline",
                                                     "C4D1",
                                                     "EoT"))

# Plot
plot <- ggplot(data = Data_long_summarised,
               aes(x = Time_point,
                   stratum = value,
                   alluvium = Genus,
                   y = value,
                   fill = Genus, 
                   label = Genus)) +
  geom_alluvium(aes(fill = Genus), decreasing = FALSE) +
  geom_flow(decreasing = FALSE) +
  geom_stratum(alpha = 0.5, decreasing = FALSE) +
  geom_text(stat = "alluvium",
            aes(label = Genus), size= 1.75,
            decreasing = FALSE) +
  facet_wrap(. ~ CBR ) +
  theme_classic()  

plot <- plot + labs(x = "Time point", y = "Relative abundance (%)",
                    title = "Relative abundance of gut genera") + 
  theme(legend.text = element_text(size=8,
                                   face = "italic"))
plot
ggsave("CALADRIO Genus Alluvial gut final.pdf",
       width = 20,
       height = 10,
       units = "cm",
       dpi = 300)
