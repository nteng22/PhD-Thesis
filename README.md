# PhD-Thesis
Scripts used during the PhD journey. It's a mix of R and bash. R script given by Matthew Dalby. Bash scripts given by Raymond Kiu. 

## [16S rRNA gene amplicon QC.R](https://github.com/nteng22/PhD-Thesis/blob/main/R-scripts/16S%20rRNA%20gene%20amplicon%20QC.R)
Scripts given by Matthew Dalby. Format the biom data file to be used in R. <br>
QC reads includes rarefying using DESeq2 to create a .csv file with normalised reads

## [Create stacked bar plot.R](https://github.com/nteng22/PhD-Thesis/blob/main/R-scripts/Stacked%20bar%20plot.R)
R script to make a stacked bar plot. Used in the thesis to show the 16S rRNA gene amplicon sequencing reads, but can be used for anything if needed. Notes from Matthew.

## [Create heatmap.R](https://github.com/nteng22/PhD-Thesis/blob/main/R-scripts/Heatmap.R)
R script to make a heatmap. Used in the thesis to show 16S rRNA gene amplicon sequencing reads. Notes from Matthew.

## [Metagenomic shotgun QC.R](https://github.com/nteng22/PhD-Thesis/blob/main/R-scripts/Metagenomic%20shotgun%20QC.R)
R script to format the data .tsv file into the format used to make the 16S rRNA plots. 

## [Alpha and beta diversity analysis and plots.R](https://github.com/nteng22/PhD-Thesis/blob/main/R-scripts/Alpha%20and%20beta%20diversity%20analysis%20and%20plots.R)
R script to do alpha and beta microbiota analysis and visualisation

## [Spaghetti plot.R](https://github.com/nteng22/PhD-Thesis/blob/main/R-scripts/Spaghetti%20plot.R)
ggplot R script to make a spaghetti plot using Shannon diveristy 

## [NMDS Taxa analysis.R](https://github.com/nteng22/PhD-Thesis/blob/main/R-scripts/NMDS%20Taxa%20analysis.R)
Not used in the thesis, but R script given by Matthew to do investigative taxa analysis based on NMDS data

## [Venn diagram.R](https://github.com/nteng22/PhD-Thesis/blob/main/R-scripts/Venn%20Diagram.R)
R script and package to use a Venn Diagram

## [16S Taxaonomic assignment](https://github.com/nteng22/PhD-Thesis/blob/main/Bioinformatics/16S%20Taxaonomic%20assignment.md)
Notes on 16S taxonomic assignment. Using Qiime2 pipeline. It creates a biom file. Consequently used [R.script](https://github.com/nteng22/PhD-Thesis/blob/main/R-scripts/16S%20rRNA%20gene%20amplicon%20QC.R) to do the QC. Used bar plot and heatmap scripts for visualisation. Notes are from Raymond.

## [Creating 16S rRNA phylogenetic trees](https://github.com/nteng22/PhD-Thesis/blob/main/Bioinformatics/16S%20rRNA%20phylogenetic%20trees.md)
Notes on 16S rRNA gene phylogenetic trees. This was used to create 16S trees for the paper published in ISJEM describing LH1062 and LH1063. Notes are from Raymond.

## [Using Checkm for contamination](https://github.com/nteng22/PhD-Thesis/blob/main/Bioinformatics/Checkm.md)
Notes on using checkm to check for contamination. This was used for the novel isolates and characterising the BEAM faecal culturing isolates.

## [Check if a bacterial chromosome is circular](https://github.com/nteng22/PhD-Thesis/blob/main/Bioinformatics/Circular%20chromosome.md)
Notes on using the software circulator on the HPC. This was used to check if LH1062's chromosome was circular. Notes are from Raymond.

## [GDTB-Tk, to determine taxonomies of isolates](https://github.com/nteng22/PhD-Thesis/blob/main/Bioinformatics/GTDB-Tk%20Notes.md)
Notes on using GDTB-Tk which is the gold-standard to classify bacteria. This also includes notes on how to create a phylogenomic tree using Mashtree (kmer based). Notes are from Raymond.

## [Functional analysis using Humann3](https://github.com/nteng22/PhD-Thesis/blob/main/Bioinformatics/Humann3.md)
Notes on using Humann3 on metagenomic data. I haven't run this personally, Raymond did but the notes are available to try. Notes are from Raymond.

## [Long read assembly](https://github.com/nteng22/PhD-Thesis/blob/main/Bioinformatics/Long%20read%20assembly.md)
Notes on using Unicycler or Flye for long-read assembly. This was used for the paper published in ISJEM.

## [PhyloPhlAn 3](https://github.com/nteng22/PhD-Thesis/blob/main/Bioinformatics/PhyloPhlAn3.md)
Notes on using PhyloPhlAn3 locally. Used to create the genomic tree used in the paper published in ISJEM.

## [Screening genomes using Abricate](https://github.com/nteng22/PhD-Thesis/blob/main/Bioinformatics/Screening%20genomes%20with%20Abricate.md)
Notes on using Abricate to screen for genes in genomes. Used in the BEAM and CALADRIO analysis. Notes from Raymond.

## [WGS analysis](https://github.com/nteng22/PhD-Thesis/blob/main/Bioinformatics/Whole%20genome%20sequencing%20(analysis).md)
Notes on the 'pipeline' for whole genome sequencing. Mix of Raymond and my notes. 
