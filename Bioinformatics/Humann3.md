## Functional analysis using Humann3

<p> Concatenate paired-end reads as the inputs for Humann3 (note from github: <br>
  As a result, the best way to use paired-end sequencing data with HUMAnN 3.0 is simply to 
  concatenate all reads into a single FASTA or FASTQ file.) </p>

```
source humann-3.0.0
source metaphlan-3.0.13
source bowtie2-2.3.4.1
source samtools-1.7
```
<p> metaphlan3 is required for humann3 as humann3 is using metaphlan3 output for estimation. <br>
  Both latest version of programs were installed as singularity image on HPC on 11/8/2021. </p>
  
Database Chocophlan (full database) required to start analysing using filtered FASTQ reads:
```
/qib/research-groups/Lindsay-Hall/Raymond/Database/chocophlan/
```
<p> This database was downloaded (16GB) on 24/8/2021 using wget on software node and tar -zxvf in the 
  folder as recommended by author in the Biobakery forum. <br>
Another uniref database is also required (19GB) as a translated search database (uniref90 diamond - recommended) downloaded on 24/8/2021. </p>

```
/qib/research-groups/Lindsay-Hall/Raymond/Database/uniref/
```

Utility mapping database:

```
/qib/research-groups/Lindsay-Hall/Raymond/Database/uniref//utility_mapping
```

Metaphlan database (use in --metaphlan-options):

```
/qib/research-groups/Lindsay-Hall/Raymond/Database/metaphlan
```
<p> To build the original bowtie2 metaphlan database (1.4GB; fasta) after downloading.
  (ALL files are required otherwise humann3 may fail to run) from https://zenodo.org/record/3957592#.YSZDVS1Q3_o </p>
 
```
sbatch --wrap " bowtie2-build -f --threads 30 mpa_v30_CHOCOPhlAn_201901.fna mpa_v30_CHOCOPhlAn_201901" -c 30 --mem=50GB -p qib-short -J bt2
```

### Step 1: QC - using kneaddata
<p> Use fastp to remove adapters and remove low-quality reads. <br<
  Use kneaddata to remove host sequences. <br>
  Please refer to Metagenome analysis for details and notes on how to do this. <br>
  Outputs: filtered fastq reads </p>

### Step 2: Concatenate paired-end reads
<p> Use the filtered fastq reads from above, concatenate them into one file. Humann3 only runs on one fastq file. </p>

```
cat read1.fastq read2.fastq > all.fastq
```

### Step 3: Run Humann3

```
sbatch --wrap "humann -i "$j".fastq -o "$j-humannanalysis" --threads 12 --metaphlan=/nbi/software/testing/humann/3.0.0/x86_64/bin/metaphlan --metaphlan-options="-x=mpa_v30_CHOCOPhlAn_201901" --metaphlan-options="--bowtie2db=/qib/research-groups/Lindsay-Hall/Raymond/Database/metaphlan/" --nucleotide-database /qib/research-groups/Lindsay-Hall/Raymond/Database/chocophlan/ --protein-database /qib/research-groups/Lindsay-Hall/Raymond/Database/uniref/" -c 12 --mem=50GB -p qib-medium -J $j-humann -o $j.humann.out -e $j.humann.err; done
```
>Note: at 30 cpus, 100GB RAM, fastq file size 28GB, it took around 8 hours to finish running humann3 program. The syntaxes are problematic with variable $j, but this version eventually worked!

<p> Output files: <br>
  You should get three output files as below as it is run successfully: <br>
  Pure_reads_genefamilies.tsv <br>
  Pure_reads_pathabundance.tsv <br>
  Pure_reads_pathcoverage.tsv </p>

### Step 4: Group data - pathabundance + genefamilies only

```
humann_join_tables -i . -o humann_all_pathabundance.tsv --file_name pathabundance
```
<p> Copy all pathabundance files into one directory (*pathabundance.tsv) and execute the command above. <br>
  Do the same for genefamilies, put all genefamilies.tsv under one sub-directory: </p>

```
sbatch --wrap "humann_join_tables -i . -o humann_all_genefamilies.tsv --file_name genefamilies" -c 10 --mem=5GB -p nbi-short -J jointables
```
<p> Submit the job as it can take longer as input files are bigger in size.</p>

### Step 5: Normalise data

```
humann_renorm_table --input humann_all_pathabundance.tsv --units cpm --output humann_all_pathabundance_cpm.tsv
```
<p> Normalise humann_all_pathabundance_cpm.tsv generated in step 4.<br>
  Submit the job as it takes longer than you think. </p>

For gene families:

```
sbatch --wrap "humann_renorm_table --input humann_all_genefamilies.tsv --output humann_all_genefamilies_cpm.tsv --units cpm" -c 10 --mem=12GB -p qib-medium -o cpm.out
```

### Step 6: Stratify data

```
humann_split_stratified_table --input humann_all_pathabundance_cpm.tsv --output ./
```

Outputs: 
-
```
rwxrwx--- 1 kiur QIB_FR009 736K Sep 3 14:09 humann_all_pathabundance_cpm_stratified.tsv
-rwxrwx--- 1 kiur QIB_FR009 140K Sep 3 14:09 humann_all_pathabundance_cpm_unstratified.tsv
```
<p> Stratified.tsv provides species/genus estimated contribution. </p>
 
For gene families:

```
humann_split_stratified_table --input humann_all_genefamilies_cpm.tsv --output ./
```

Path abundance is more useful to visualise the functionalities of the microbiome.

---
28 June 2022
Get it unstratified and stratified tsv files. 
Ignore UNITERGRATED and UNMAPPED
