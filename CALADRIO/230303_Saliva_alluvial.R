#### Housekeeping ####
setwd("~/Dropbox/QIB/KELLY Study/Saliva_data/")

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

#setwd("C:/Users/dalby/Documents/Hall Lab/Biom R files for Nancy/Alluvial plots")

#### Data ####
Data <- read_csv("KELLY_saliva_normalised%_with_metadata.csv")

Data <- Data %>%
  select(sample_ID,
         #PDL1,
         #NLR,
         #tPFS
         #CBR,
         Timepoint,
         "5-7N15":ncol(Data)
  )

genus <- Data %>%
  select(sample_ID,
         "5-7N15":ncol(Data)
  )

genus <- column_to_rownames(genus, var = "sample_ID")

# Remove any columns with a sum total of 0.
genus <- genus %>%
  select(where(is.numeric)) %>% 
  select(where(~sum(.) >0))

# check is normalised
rowsums <- rowSums(genus)
rowsums <- as.data.frame(rowsums)
rowsums

metadata <- Data %>%
  select(sample_ID,
        #CBR,
         Timepoint)

metadata$CBR <- gsub("1", "CB", metadata$CBR)
metadata$CBR <- gsub("0", "no CB", metadata$CBR)
metadata$Timepoint <- gsub("Week9", "C4D1", metadata$Timepoint)

top <- genus[,order(colSums(genus),decreasing=TRUE)]

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
other <- genus[,order(colSums(genus),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(genus)[2] # This needs to be the total number of colunms in the dataframe.
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
                  id.vars = c("sample_ID",
                              "Timepoint"),
                              #"CBR"),
                  #"PDL1",
                  #"tPFS",
                  #"NLR"), 
                  variable.name = "Genus")

Data_long$Timepoint <- factor(Data_long$Timepoint,
                               levels = c("Baseline",
                                          "C4D1",
                                          "EoT"))

Data_long_summarised <- Data_long %>%
  group_by(Timepoint,#CBR,
           Genus) %>%
  summarize(value=mean(value))

Data_long_summarised$CBR <- as.factor(Data_long_summarised$CBR)

Data_long_summarised$Timepoint <- factor(Data_long_summarised$Timepoint,
                                          levels = c("Baseline",
                                                     "C4D1",
                                                     "EoT"))


# Plot
plot <- ggplot(data = Data_long_summarised,
               aes(x = Timepoint,
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
  #facet_wrap(. ~ CBR ) +
  theme_classic()  

plot <- plot + labs(x = "Time point", y = "Relative abundance (%)",
                    title = "Relative abundances of oral genera") + 
  theme(legend.text = element_text(size=8,
                                   face = "italic"))
plot

ggsave("CALADRIO Genus Alluvial Oral microbiome.pdf",
       width = 20,
       height = 10,
       units = "cm",
       dpi = 300)
