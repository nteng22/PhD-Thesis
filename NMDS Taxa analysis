#install.packages("ggplot2")
library(ggplot2)
#install.packages("vegan")
library(vegan)
#install.packages("dplyr")
library(dplyr)

#################################################################################################
# Script for NMDS plot for time point 2 (samples from 10-29 days of age)

# Set the working directory to import the data.
setwd("~/Dropbox/QIB/KELLY Study/Genus analysis/")
# Import the data.
data <- read.csv("example BAMBI genus data with metadata.csv")

# Subset data using dplyr to select only the columns required.
data_subset <- dplyr::select(filter(data),
                       c(Treatment,
                         birthweight,
                         gestageweeks,
                         antibiotics,
                         Acinetobacter:ncol(data)))


# Remove any rows that contain any missing values (NMDS will not work with any NA's in the data).
data_subset <- na.omit(data_subset, cols=Acinetobacter:ncol(data_subset))

# Remove any columns that only contain zeros (with a sum total of 0).
# Columns with only zeros do not contribute to the analysis but seriously slow it down.
# Check the totals of each column.
colsums <- colSums(data_subset[5:ncol(data_subset)])
colsums <- as.data.frame(colsums)
colsums
# Remove any columns that sum to zero.
data_subset_no_zero_columns <- data_subset[, colSums(data_subset != 0) > 0]
# Check that there are no longer any columns that sum to zero.
colsums2 <- colSums(data_subset_no_zero_columns[5:ncol(data_subset_no_zero_columns)])
colsums2 <- as.data.frame(colsums2)
colsums2

# Select genus sequence data
abund_table <- dplyr::select(filter(data_subset_no_zero_columns),
                             c(Acinetobacter:ncol(data_subset_no_zero_columns)))

# Select metadata
meta_table <- dplyr::select(filter(data_subset_no_zero_columns),
                            c(Treatment,
                              antibiotics,
                              birthweight,
                              gestageweeks
                              ))

# Make any factor metadata variable into a factor.
meta_table$Treatment <- as.factor(meta_table$Treatment)
meta_table$antibiotics <- as.factor(meta_table$antibiotics)

# Make any numeric metadata variable into a numeric variable.
meta_table$birthweight <- as.numeric(meta_table$birthweight)
meta_table$gestageweeks <- as.numeric(meta_table$gestageweeks)

# PERMANOVA test in vegan to test multiple variables against the microbiota composition.
adonis(abund_table ~
         Treatment +
         birthweight +
         gestageweeks +
         antibiotics,
       data = meta_table,
       permutations = 999,
       method = "bray")



#===============================================================================================================
# Create bv.step function to use later
# Code sourced from: userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html

# This R script is an extension of vegan library's bioenv()
# function and uses the bio.env() and bio.step() of
#	http://menugget.blogspot.co.uk/2011/06/clarke-and-ainsworths-bioenv-and-bvstep.html
#	The original author suggested these functions to overcome
#	the inflexibility of the bioenv() function which uses
#	a similarity matrix based on normalized "euclidean" distance.
# The new functions are given below and implement the following algorithms:
# Clarke, K. R & Ainsworth, M. 1993. A method of linking multivariate community structure to environmental variables. Marine Ecology Progress Series, 92, 205-219.
# Clarke, K. R., Gorley, R. N., 2001. PRIMER v5: User Manual/Tutorial. PRIMER-E, Plymouth, UK.
# Clarke, K. R., Warwick, R. M., 2001. Changes in Marine Communities: An Approach to Statistical Analysis and Interpretation, 2nd edition. PRIMER-E Ltd, Plymouth, UK.
# Clarke, K. R., Warwick, R. M., 1998. Quantifying structural redundancy in ecological communities. Oecologia, 113:278-289.

# Create bv.step function.
bv.step <- function(fix.mat, var.mat,
                    fix.dist.method="gower", var.dist.method="euclidean", correlation.method="spearman",
                    scale.fix=FALSE, scale.var=TRUE,
                    max.rho=0.95,
                    min.delta.rho=0.001,
                    random.selection=TRUE,
                    prop.selected.var=0.2,
                    num.restarts=10,
                    var.always.include=NULL,
                    var.exclude=NULL,
                    output.best=10
){

  if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
  if(sum(var.always.include %in% var.exclude) > 0){stop("var.always.include and var.exclude share a variable")}
  require(vegan)

  if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
  if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}

  fix.dist <- vegdist(as.matrix(fix.mat), method=fix.dist.method)

  #an initial removal phase
  var.dist.full <- vegdist(as.matrix(var.mat), method=var.dist.method)
  full.cor <- suppressWarnings(cor.test(fix.dist, var.dist.full, method=correlation.method))$estimate
  var.comb <- combn(1:ncol(var.mat), ncol(var.mat)-1)
  RES <- data.frame(var.excl=rep(NA,ncol(var.comb)), n.var=ncol(var.mat)-1, rho=NA)
  for(i in 1:dim(var.comb)[2]){
    var.dist <- vegdist(as.matrix(var.mat[,var.comb[,i]]), method=var.dist.method)
    temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
    RES$var.excl[i] <- c(1:ncol(var.mat))[-var.comb[,i]]
    RES$rho[i] <- temp$estimate
  }
  delta.rho <- RES$rho - full.cor
  exclude <- sort(unique(c(RES$var.excl[which(abs(delta.rho) < min.delta.rho)], var.exclude)))

  if(random.selection){
    num.restarts=num.restarts
    prop.selected.var=prop.selected.var
    prob<-rep(1,ncol(var.mat))
    if(prop.selected.var< 1){
      prob[exclude]<-0
    }
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  } else {
    num.restarts=1
    prop.selected.var=1
    prob<-rep(1,ncol(var.mat))
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  }

  RES_TOT <- c()
  for(i in 1:num.restarts){
    step=1
    RES <- data.frame(step=step, step.dir="F", var.incl=NA, n.var=0, rho=0)
    attr(RES$step.dir, "levels") <- c("F","B")
    best.comb <- which.max(RES$rho)
    best.rho <- RES$rho[best.comb]
    delta.rho <- Inf
    selected.var <- sort(unique(c(sample(1:dim(var.mat)[2], n.selected.var, prob=prob), var.always.include)))
    while(best.rho < max.rho & delta.rho > min.delta.rho & RES$n.var[best.comb] < length(selected.var)){
      #forward step
      step.dir="F"
      step=step+1
      var.comb <- combn(selected.var, RES$n.var[best.comb]+1, simplify=FALSE)
      if(RES$n.var[best.comb] == 0){
        var.comb.incl<-1:length(var.comb)
      } else {
        var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
        temp <- NA*1:length(var.comb)
        for(j in 1:length(temp)){
          temp[j] <- all(var.keep %in% var.comb[[j]])
        }
        var.comb.incl <- which(temp==1)
      }

      RES.f <- data.frame(step=rep(step, length(var.comb.incl)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]+1, rho=NA)
      for(f in 1:length(var.comb.incl)){
        var.incl <- var.comb[[var.comb.incl[f]]]
        var.incl <- var.incl[order(var.incl)]
        var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
        RES.f$var.incl[f] <- paste(var.incl, collapse=",")
        RES.f$rho[f] <- temp$estimate
      }

      last.F <- max(which(RES$step.dir=="F"))
      RES <- rbind(RES, RES.f[which.max(RES.f$rho),])
      best.comb <- which.max(RES$rho)
      delta.rho <- RES$rho[best.comb] - best.rho
      best.rho <- RES$rho[best.comb]

      if(best.comb == step){
        while(best.comb == step & RES$n.var[best.comb] > 1){
          #backward step
          step.dir="B"
          step <- step+1
          var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
          var.comb <- combn(var.keep, RES$n.var[best.comb]-1, simplify=FALSE)
          RES.b <- data.frame(step=rep(step, length(var.comb)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]-1, rho=NA)
          for(b in 1:length(var.comb)){
            var.incl <- var.comb[[b]]
            var.incl <- var.incl[order(var.incl)]
            var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
            temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
            RES.b$var.incl[b] <- paste(var.incl, collapse=",")
            RES.b$rho[b] <- temp$estimate
          }
          RES <- rbind(RES, RES.b[which.max(RES.b$rho),])
          best.comb <- which.max(RES$rho)
          best.rho<- RES$rho[best.comb]
        }
      } else {
        break()
      }

    }

    RES_TOT <- rbind(RES_TOT, RES[2:dim(RES)[1],])
    print(paste(round((i/num.restarts)*100,3), "% finished"))
  }

  RES_TOT <- unique(RES_TOT[,3:5])


  if(dim(RES_TOT)[1] > output.best){
    order.by.best <- RES_TOT[order(RES_TOT$rho, decreasing=TRUE)[1:output.best],]
  } else {
    order.by.best <-  RES_TOT[order(RES_TOT$rho, decreasing=TRUE), ]
  }
  rownames(order.by.best)<-NULL

  order.by.i.comb <- c()
  for(i in 1:length(selected.var)){
    f1 <- which(RES_TOT$n.var==i)
    f2 <- which.max(RES_TOT$rho[f1])
    order.by.i.comb <- rbind(order.by.i.comb, RES_TOT[f1[f2],])
  }
  rownames(order.by.i.comb)<-NULL

  if(length(exclude)<1){var.exclude=NULL} else {var.exclude=exclude}
  out <- list(
    order.by.best=order.by.best,
    order.by.i.comb=order.by.i.comb,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(order.by.best$var.incl[1], ",")))], collapse=","),
    best.model.rho=order.by.best$rho[1],
    var.always.include=var.always.include,
    var.exclude=var.exclude
  )
  out

}
#===============================================================================================================

# Define parameter commands (other options are available for each one).
# Bray (Bray-Curtis) is the most commonly used one and in usually best to use.
# You can change these without changing the command in the functions.
cmethod<-"pearson" # Correlation method to use: pearson, spearman, kendall
fmethod<-"bray" # Fixed distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
vmethod<-"bray" # Variable distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
nmethod<-"bray" # NMDS distance method:  euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao

# Create the command to find the top 10 genus driving the separation of samples.
res.bv.step.biobio <- bv.step(wisconsin(abund_table),
                              wisconsin(abund_table),
                              fix.dist.method=fmethod,
                              var.dist.method=vmethod,
                              correlation.method=cmethod,
                              scale.fix=FALSE,
                              scale.var=FALSE,
                              max.rho=0.95,
                              min.delta.rho=0.001,
                              random.selection=TRUE,
                              prop.selected.var=0.3,
                              num.restarts=10,
                              output.best=10,
                              var.always.include=NULL)

res.bv.step.biobio
# WARNINGS: These refer to empty rows and can be ignored.
# Sub-samplING the community matrix in bv.step to only include a subset of genus
# means there will be cases when the "abund_table" will have empty rows,
# especially when we are selecting only 10 genus.


# Find the 10 best subset of genus driving the separation of samples
taxaNames<-colnames(abund_table)
bestTaxaFit <-""
for(i in (1:length(res.bv.step.biobio$order.by.best$var.incl)))
{
  bestTaxaFit[i] <- paste(paste(taxaNames[as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[i],
                                                                     split = ",")))],
                                collapse = ' + '),
                          " = ",
                          res.bv.step.biobio$order.by.best$rho[i],
                          sep = "")
}

# Select the column (genus) names of the top ten genus
bestTaxaFit <- data.frame(bestTaxaFit)
colnames(bestTaxaFit) <- "Best combination of taxa with similarity score"
bestTaxaFit

# Generate the NMDS plot
MDS_res = metaMDS(abund_table,
                  distance = nmethod,
                  k = 3, # Set the number of dimensions
                  smin = 0.0001, # 0.0001 is the default
                  sfgrmin = 0.0000001, # 0.0000001 is the default
                  sratmax = 0.99999, # 0.99999 is the default
                  try = 200, # Minimum number of random starts
                  maxit = 300, # Maximum number of random starts
                  trymax = 300, # Number of attempts at calculating the NMDS plot to find the best one
                  plot = TRUE, # TRUE shows the plots developing as the code is running.
                  engine = "monoMDS")
MDS_res


# Test the significance of the taxa in the NMDS plot
bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[1], ",")))
bio.fit <- envfit(MDS_res, abund_table[,bio.keep,drop=F], perm = 999)
bio.fit

# Get the vectors for bio.fit
# Limit the taxa shown to only those with significant p-values.
# This can be useful to avoid overcrowding the NMDS plot with overlapping names.
# Extract the data for the arrows.
vectors <- as.list(bio.fit$vectors)
vectors
# r2 is the degree to which the taxa explain differences in the NMDS plot.
# Scale the length of the arrows by the r2.
arrows <- as.data.frame(vectors$arrows*sqrt(vectors$r))
arrows
# Extract the p-values for the arrows.
pvals <- as.data.frame(vectors$pvals)
pvals
# Create a dataframe with the arrows and p-values
significant_arrows <- cbind(arrows, pvals)
significant_arrows
# Subset the dataframe to only those with p<0.01
significant_taxa <- subset(significant_arrows, pvals < 0.05)
significant_taxa


# Test the significance of the metadata variables in the NMDS plot
en = envfit(MDS_res, meta_table, permutations = 9999, na.rm = TRUE)
en

# Test the nMDS plot and factors using base R.
#plot(MDS_res)
#plot(en)


# Extract the continuous numeric metadata variables that explain NMDS plot variance.
B <- as.list(en$vectors)
# Create a dataframe with the arrows and p-values
B_pvals <- as.data.frame(B$pvals)
B_arrows <- as.data.frame(B$arrows * sqrt(B$r))
B_sig_arrows <- cbind(B_arrows, B_pvals)
# Subset the dataframe to only those with p<0.01
en_coord_cont_sig <- subset(B_sig_arrows, B$pvals < 0.05)
# Use en_coord_cont_sig to annotate the ggplot.
en_coord_cont_sig

# Extract the categorical factor metadata variables that explain NMDS plot variance.
en_factor <- as.list(en$factors)
en_factor
# Create a dataframe with the centroids (the center point of each factor level) and p-values.
en_factor_pvals <- as.data.frame(en_factor$pvals)
en_factor_pvals <- rownames_to_column(en_factor_pvals, "group")
en_factor_var_id <- as.data.frame(en_factor$var.id)
en_factor_centroids <- as.data.frame(en_factor$centroids)
en_factor_centroids <- rownames_to_column(en_factor_centroids, "centroid")
bind <- cbind(en_factor_var_id, en_factor_centroids)
names(bind)[names(bind) == 'en_factor$var.id'] <- 'group'
merged <- merge(bind, en_factor_pvals, by = "group")
names(merged)[names(merged) == 'en_factor$pvals'] <- 'pvals'
merged
# Limit metadata factors shown to only those with significant p-values
# This is to avoid overcrowding the NMDS plot with overlapping names
# Subset the dataframe to only those with p < 0.05
en_coord_cat_sig <- subset(merged, merged$pvals < 0.05)
en_coord_cat_sig <- en_coord_cat_sig[,2:4]
en_coord_cat_sig <- rownames_to_column(en_coord_cat_sig, "none")
en_coord_cat_sig <- column_to_rownames(en_coord_cat_sig, "centroid")
en_coord_cat_sig <- en_coord_cat_sig[,-1]
# Use en_coord_cat_sig to annotate the ggplot.
en_coord_cat_sig

# Prepare the dataframe to make the ggplot.
# Extract the coordinates to draw the NMDS plot.
dataframe <- scores(MDS_res)
dataframe

# Convert into a dataframe.
dataframe <- as.data.frame(dataframe)
dataframe

# Add Treatment group information to the dataframe.
dataframe <- data.frame(dataframe,
                        Treatment = meta_table[,1])
head(dataframe)

# Set the levels of the metadata factor that is used to color the plot.
# This will set the order groups will appear in the legend and when assigning colors and symbols.
dataframe$Treatment <- factor(dataframe$Treatment,
                              levels = c("Control",
                                         "Bif/Lacto"
                                         ))
head(dataframe)


#Create the ggplot
p <- ggplot(data = dataframe,
            aes(x = NMDS1,
                y = NMDS2)) +
  # Set the metadata to be used to change the colour and shape of the plot symbols.
     geom_point(aes(colour = Treatment,
                    shape = Treatment)) +
  # Set the color of the symbols used for each group.
     scale_colour_manual(name = "Treatment",
                         labels = c("Control",
                                    "Bif/Lacto"),
                         values = c("red1",
                                    "steelblue3")) +
  # Set the shape of the symbols used for each group.
     scale_shape_manual(name = "Treatment",
                        labels = c("Control",
                                   "Bif/Lacto"),
                        # Define the shapes of the symbols used in the plot (each number is a different R symbol).
                        values = c(22, 24)) +
  # Defining the x axis limits helps stop taxa names being cut off.
     scale_x_continuous(limits = c(-1.4, 1.1)) +
  # Set the plot title.
  labs(title = "Infant microbiota at 10-29 days of age") +
  theme(# Define the plot title.
        plot.title = element_text(size=8),
        # Define the plot margin size.
        plot.margin = margin(t = 4, r = 4, b = 4, l = 4),
        # Define the legend position.
        legend.position = "right",
        # Define the x axis legend spacing.
        legend.spacing.x = unit(.1, 'cm'),
        # Define the y axis legend spacing.
        legend.spacing.y = unit(.1, 'cm'),
        # Format the legend title.
        legend.title = element_text(size=8),
        # Format the size of the legend text.
        legend.text = element_text(size=8),
        # Remove the grey background.
        legend.background = element_blank(),
        # Remove rectangle around the legend
        legend.box.background = element_blank(),
        # Remove the box from the legend.
        legend.key = element_blank(),
        # Remove the grey plot background.
        panel.background = element_blank(),
        # Remove the plot border.
        panel.border = element_blank(),
        # Remove the major plot grid lines.
        panel.grid.major = element_blank(),
        # Remove the minor plot grid lines.
        panel.grid.minor = element_blank(),
        # Format x axis title.
        axis.title.x = element_text(size=8, colour = "black"),
        # Format x axis labels.
        axis.text.x = element_text(size=8, colour = "black"),
        # Format the y axis title text size.
        axis.title.y = element_text(size=8, colour = "black"),
        # Format the axis label text size.
        axis.text.y = element_text(size=8, colour = "black"),
        # Format x and y axis ticks.
        axis.ticks = element_line(size = 0.35, colour = "black"),
        # Format x and y axis lines.
        axis.line = element_line(size = 0.35, colour = "black"),
        # Define the plot aspect ratio.
        aspect.ratio = 1) +

  # Add the taxa arrows to the ggplot.
  geom_segment(data = significant_taxa,
               aes(x = 0,
                   y = 0,
                   xend = NMDS1,
                   yend = NMDS2),
               # Format the size of the arrow head.
               arrow = arrow(length = unit(0.2, "cm")),
               # Format the color of the arrows.
               color = "black",
               # Format the thickness of the arrows.
               alpha = .8) +

  # Add the taxa names to the ggplot.
  geom_text(data = as.data.frame(significant_taxa * 1.3), # "significant_taxa*1.5" moves the names 1.5 x the length of the arrow.
            aes(NMDS1,
                NMDS2,
                label = rownames(significant_taxa)),
            # Format genus font.
            fontface = "italic",
            # Format genus names.
            color = "black",
            # Format the size of the genus names.
            size = 3) +

  # Add the continuous metadata variable arrows to the ggplot.
  geom_segment(data = as.data.frame(en_coord_cont_sig),
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               # Format the size of the arrow head.
               arrow = arrow(length = unit(0.2, "cm")),
               # Format the color of the arrows.
               colour = "purple",
               # Format the thickness of the arrows.
               alpha = 0.8) +

  # Add the continuous metadata names to the ggplot.
  geom_text(data = en_coord_cont_sig * 1.7, # Make the name labels 1.6 x away from the length of the arrows.
            aes(x = NMDS1, y = NMDS2),
            label = row.names(en_coord_cont_sig),
            # Format the color of the variable names.
            colour = "purple",
            # Format the font of the variable names.
            fontface = "bold",
            # Format the size of the variable names.
            size = 3) +

  # Add the categorical metadata variable points to the ggplot.
  geom_point(data = en_coord_cat_sig,
             aes(x = NMDS1, y = NMDS2),
             # Format the shape of the symbol.
             shape = "diamond",
             # Format the color of the symbol.
             colour = "purple",
             # Format the size of the symbol.
             size = 5) +

  # Add the categorical metadata variable names to the ggplot.
  geom_text(data = en_coord_cat_sig,
            aes(x = NMDS1, y = NMDS2 + 0.08), # Add +0.08 to move the names above the diamonds.
            label = row.names(en_coord_cat_sig),
            # Format the color of the variable names.
            colour = "purple",
            # Format the font of the variable names.
            fontface = "bold",
            # Format the size of the variable names.
            size = 3,)

# Print the plot within R.
p



# This will print into the current working directory file.
# Set the saved pdf file name and define the size of the plot (height and width are in inches).
pdf("NMDS example BAMBI data with Taxa and Metadata overlay 10.29 days.pdf", width = 6, height = 6)

# print(p) saves the figure as a pdf file.
print(p)
dev.off()
