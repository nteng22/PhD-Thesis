plot <- ggplot(shannon_metadata, aes(x = Time_point, y = shannon_diversity)) +
  facet_grid(. ~ CBR) +
  geom_boxplot(aes(fill = Time_point, shape = Time_point, colour = Time_point)) + #, outlier.shape = NA) +
  geom_line(aes(group = ID_Patient),
            linetype = "solid",
            size = 0.05) +
  geom_point(position = "identity", aes(fill = Time_point, 
                                        shape = Time_point, 
                                        colour = Time_point,
                                        group = ID_Patient), # change the shape 
             size = 2,
             stroke = .5) +
  scale_shape_manual(values=c(21, 21, 21)) +
  scale_color_manual(values=c("red",
                              "blue",
                              "goldenrod")) +
  scale_fill_manual(values=c("white", "white", "white")) + # inner colour can change
  #Set the axis labels and title
  labs(title="GUT genus Shannon Diversity all samples across Timepoints by CBR",
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

#ggsave("CALADRIO Gut genus shannon diversity all samples across Tp~CBR.png", width = 200, height = 200, units = c("mm"), dpi = 400)

# Save by pdf
plot

pdf(file = "~/Dropbox/QIB/KELLY Study/initial analysis/Pasta party/CALADRIO Gut genus all samples by CB.pdf", width = 11, height = 8)

print(plot)
dev.off()
