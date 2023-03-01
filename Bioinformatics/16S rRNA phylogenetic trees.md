## Creating 16S rRNA phylogenetic trees

> <p> Need to align the 16S rRNA sequences. <br>
> Muscle to create alignment. <br>
> It's a software that can align up to 500 sequences of 1MB size. So only rRNA or conserved markers. </p>
> 
```
source muscle-3.8.31
```
To make a multi alignment file. This aligns all the fasta sequences relative to each other into one file
- Input file needs to be a concatenated file with all sequences in one 
- i.e. cat *.fasta > sequences.fasta
```
muscle -in [input file] -out [output file] 
```

## IQtree to create the tree
```
source iqtree-2.0.5 
```

> <p> This one is faster </p>
```
sbatch --wrap "iqtree2 -s core_gene_snp.aln -m TEST -alrt 1000 -B 1000 -T AUTO" -c 30 --mem=25GB -p qib-medium -o iqtree2-UFBS1000.out -J iqtree
```

> <p> This one is standard but slow </p>
```
sbatch --wrap "iqtree2 -s Enterococcus-16S-2.aln -m TEST -b 100 -T AUTO" -c 30 --mem=15GB -p qib-long -o iqtree2.out
```
