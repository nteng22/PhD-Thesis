setwd("~/Dropbox/QIB/KELLY Study/Venn Diagram/")

library("tidyverse")
library("ggplot2")

gut <- read_csv("../Genus_filtered_reads.csv") %>%
  select(name)
oral <- read_tsv("../Saliva_data/Genus.txt")
oral <- oral[,1] #only select first column

#install.packages("VennDiagram")
library(VennDiagram)

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

colnames(gut) <- "Gut Genera"
colnames(oral) <- "Oral Genera"

# Need to have the data as a list
  # subset of the list needs to be the 'category name'
df <- c(gut, oral)

venn.diagram(df, filename = "221216_VD_gut_oral.png")
             
display_venn(df,
             category.names = c("Gut Genera", "Oral Genera"),
             fill = c("#56B4E9", "#009E73"),
             cex = 1,
             fontface = "italic",)

intersect <- intersect(df$`Gut Genera`, df$`Oral Genera`)
tb <- as.data.frame(intersect)
write_tsv(tb, "Intersect of gut and oral genera.tsv")
