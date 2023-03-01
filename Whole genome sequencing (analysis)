Based on the NBI HPC software and packages 
Need to download the files from IRIDA
 To make it easier rename the files so PID-XXX is removed 
 See anyname.sh
 Can also rename the individual files
```
do mv *_R1_*.fastq "$f"_R1.fastq
do mv *_R2_*.fastq "$f"_R2.fastq
```

Note
``` 
$(cat list1.txt)
```
Provided the list has the correct sample name will batch submit. 
#### Step 1 ####
# Then you need to use fastp to remove primers and rubbish reads
# This needs to be run on interactive partition. 

```
source package fastp-0.20.0
sbatch --wrap "fastp -i "$f"_R1.fastq -I "$f"_R2.fastq -o "$f"_R1_filtered.fastq -O "$f"_R2_filtered.fastq -q 20" -c 4 --mem=6GB -p qib-short -J fastp
```
#### Step 2 ####
# De novo assembly using Spades
# --careful flag will make the alignment more stringent.

```
source package spades-3.11
sbatch --wrap "spades.py --careful -1 "$f"_R1_filtered.fastq -2 "$f"_R2_filtered.fastq -o Spades_"$f"" -c 6 -p qib-short --mem=6GB -J Spades
```

#### Step 3 #### You can leave interactive node now
# filter out short contigs
# need the path to the script "filter-contig.pl" 
```
filter-contig.pl 500 contigs.fasta > "$f".fasta
```

```
do PATH_TO/filter-contig.pl 500 contigs.fasta > "$f".fasta; done
```

#### Step 4 ####
# Check contamination by using BactSpecies. 
# For some reason loading dependency packages doesn't work automatically
	bactspeciesID -l 
# Therefore need to load it manually: 

```
source package bactspeciesID-1.2;
source package /nbi/software/testing/bin/samtools-1.7;
source package /nbi/software/testing/bin/bedtools-2.24.0; 
source package /nbi/software/testing/bin/abricate-1.0.1; 
source package /tgac/software/testing/bin/barrnap-0.7; 
source package /nbi/software/testing/bin/any2fasta-0.4.2; 
source package /nbi/software/testing/bin/perl-5.22.1; 
source package /nbi/software/production/bin/emboss-6.5.7;
source package /nbi/software/testing/bin/blast+-2.2.30
```

#load and check dependencies by doing -l, -b respectively
```
sbatch --wrap "bactspeciesID -m TRUE "$f".fasta" -c 4 --mem=4GB -p qib-short -J BactspeciesID
```

# File output will be "$f".fasta.species. 
# To summarize, copy all the *.fasta.species files into one folder
# make a list with the file names: list1
# use anyname.sh to remove the .fasta.species part in the name, name this list 2
```
for j in $(cat list2); do cp ../$j/Spades-$j/$j.fasta.species >> $j.fasta.species; done
```
# >> tells linux to append document i.e. add not overwrite 
	# Want the file name first, before appending the output
# change it to results to columns 
```
for j in $(cat list1.txt); do paste -s $j > ${j%.fasta.species}.id; done
```

#### Step 5 ####
# Generate assembly stats

```
assembly-stats.pl "$f".fasta > "$f".out 
```

# Use the BactspeciesID output to look up the GC% and genome size.
# GC% is more important than genomes size for faecal culturing, as there may be plasmids and MGE that can increase size
# Look at number of contigs to see if its a good quality genome. This is usually short read so it should ideally be less than 100. 
# If happy with stats continue to fastANI

#### Step 6 ####
# FastANI will check the average nucleotide identity to your reference genome
#do not run on interactive node as software won't run properly

```
source package fastANI-1.3  
source package gcc-5.2.0

sbatch --wrap "fastANI -q "$f".fasta -r "reference-genome".fasta -o "$f"-ANI" -p qib-short --mem=4GB
```

# Usually a fastANI > 95% confirms a good identity and a good quality genome
# For the reference genomes, you need to find the type strain genome
