## 16S rRNA gene amplicon sequencing processing
<p> No need to use interactive.
  Use 'ls' to look into the file size of your reads. You're expective 70M+ in size. <br>
  Also look at the number of reads, need a minimum of 10k. </p>
  
  ```
  wc -l "filename". 
  ```
  > Divide this number by 4 (reads in fastq files are 4 lines).

General format for these commands (for loops):

```
for j in "file_name"; do sbatch --wrap " -command- "; done
```

#### To loop across files in a folder
>Make a list with the file names (without extension depending on the input command). <br>
>You need to declare the list a variable, which you do by e.g. <br>
> ```
> for j in $(cat list); do ...
> ```
> Now you're telling the terminal, that it needs to treat every items inside 'list' as a variable.
---

<p> FASTQ files via Illumina sequencing are usually demultiplexed if you have separate samples already <br>
  (R1 and R2) as these were performed on the sequencing machine. </p>

Illumina fastq reads are usually 300bp insert size, on Miseq platform for 16S amplicon sequencing as reads are longer.

**Note: Rename fastq now if you can if not !** <br>
>When renaming, you can only use a maximum of one underscore. Keep it simple i.e. x_1.fastq <br>
>Do not use special characters, only letters and numbers. <br>
>When Qiime merges the files the special characters can mess up the final file name. <br>

### ~~STEP 1:~~ <br>
~~QC: The illumina FASTQ reads are not <10000 reads (pair-end) for better quality. <br>
  May use assembly-stats to generate the stats for read count. <br>
  Then 2 steps of quality filtering is followed. 4 major steps as below highlighted in red:~~

### STEP 2:
Use PEAR to merge reads (300bp MiSeq) to get longer reads <br>
for better OTU assignment (also quality filtering)

```
source pear-0.9.6 
pear -f "$j"_1.fastq -r "$j"_2.fastq -q 30 -o $j-merged
```

> PEAR 0.9.6: previously Shab used -q 40 I think too stringent so <br>
> phred 30 is considered OK for most 16S data analysis. <br>
> Use 2 cpus + 4GB RAM should be enough for each job submission.

PEAR will output stats of assembly and trimming - so specify -o when submit jobs to clusters.

**IMPORTANT**: <br>
assembled fastq must be renamed now, <br>
**do not use e.g. 30122_1.fastq as qiime will consider the name before _ as the sample ID***, <br>
  if you have file names like 30122_1, 30122_2, <br>
  they will be named as 30122 for both this will mess up the final BIOM files merging you can see the id: in biom file.

> Use the slurm.out file to look at the quality score. More than 90%+ assembly is good.

### STEP 3:
Use split_libraries (qiime) to do quality filtering + convert into fasta for downstream analysis

```
source qiime-1.9.1
export PATH=$PATH:/tgac/software/testing/usearch61/6.1/x86_64/bin/
```
This tells the environment that usearch is now in your pathway. <br>

May need to use this to instruct the HPC to use UTF-8 as the encoding language

```
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8
vi ~/.bash_profile 

# Next time it will run everything inside, bash_profile everytime you log on to HPC.
```

```
split_libraries_fastq.py -i $j-merged.assembled.fastq --sample_ids $j -o $j-filtered -q 29 --barcode_type 'not-barcoded' --phred_offset 33
```

> --phred_offset must be included else will fail. <br>
> For illumina sequences, this is only used for quality filtering, <br>
> for 454 data this can be used for demultiplexing. <br>
> -q can be 19 (>20) or -q 29 (for phred 30 quality and above), <br>
> personally prefer -q 29 (to be consistent with PEAR) but this also means more reads will be discarded. <br>
> However high quality reads will remain,  decide this based upon your data. If your sequencing quality is not perfect maybe do -q 19.

filtered fasta sequences: <br>
- seqs.fna <br>

stats:
- split_library_log.txt <br>

histograms:
- histograms.txt <br>

now seqs.fna is ready for chimera detection and deletion.

#### Move the seqs.fna into another folder and rename them accordingly using loop.

```
for j in "file_name"; do sbatch --wrap "cp seqs.fna ../$j.fna"
```

### STEP 4:
Identify and delete chimeras (Illumina amplicon sequencing): <br>
To identify chimeras:

```
identify_chimeric_seqs.py -i $j.fna -m usearch61 -o $j-usearch_checked_chimeras/ -r /qib/research-groups/Lindsay-Hall/Raymond/Database/Chimera_db/uchime_chimera.fa
```

Database location for Chimeras:
```
/qib/research-groups/Lindsay-Hall/Raymond/Database/Chimera_db/uchime_chimera.fa
```

To delete chimeras from sequences:
```
filter_fasta.py -f $j.fna -o $j-seqs_chimeras_filtered.fna -s $j-usearch_checked_chimeras/chimeras.txt -n
```
> -n is important (negate). <br>
> This is to discard passed seq ids (chimeras ids) rather than keep passed seq ids. <br> 
> -s is the seq ids, chimeras id to be deleted in chimeras.txt. <br>

The output file you're interested in is "chimeras_filtered.fna"

### STEP 5:
Pick OTUs: <br>
> Use open reference - will also cluster unassigned OTUs this is a more unbiased method, <br>
> recommend this. This is using 16S silva 132 database, <br>
> updated in 2019. This method is based on uclust by default.
 
Database:

```
/qib/research-groups/Lindsay-Hall/Raymond/Database/silva_132_97_16S.fna
```

To pick open reference OTUs:

```
pick_open_reference_otus.py -o $j-pickOpenOtus/ -i $j-seqs_chimeras_filtered.fna -f -r /qib/research-groups/Lindsay-Hall/Raymond/Database/silva_132_97_16S.fna -s 0.1 --prefilter_percent_id 0.95
```

> -s is the percentage of failed samples to be included in de novo clustering <br>
> --prefilter_percent_id 0.95 is to screen the sequences (95% similarity and above <br>
> will remain and the rest discarded), a quality filtering. <br>
> "pick_open_references_otus.py" is the most unbiased assignment. </p>

To assign at species level (more stringent at >98% similarity):

```
for j in FILE_NAMES; do sbatch --wrap "pick_open_reference_otus.py -o $j-pickOpenOtus99/ -i $j-seqs_chimeras_filtered.fna -f -r /hpc-home/kiur/Raymond/Database/silva_132_99_16S.fna -s 0.1 --prefilter_percent_id 0.98" -c 4 --mem=4GB -p nbi-short;done
```

<p> Prefilter at 98.5%, sequences with similarity below that will be discarded. <br>
  Then assign taxa using 99% clustered SILVA database. </p>

The output file you're interested in is: "no_pynast_failures_otus.py" file output. <br>

### STEP 6:
<p> This step to be done once you've applied the previous steps to all files you're interested in <br>
  i.e. generated all the biom files for your samples. </p>

Copy all the "no_pynast_failures_otus.py" into one folder and rename with your sample name. This essentially becomes your SampleID.biom.

```
for j in $(cat list.txt); do cp $j-pickOpenOtus/otu_table_mc2_w_tax_no_pynast_failures.biom ../biom_files/$j.biom; done
```

To merge multiple OTU BIOM tables into one:

```
merge_otu_tables.py -i otu_table1.biom,otu_table2.biom -o merged_otu_table.biom
```
> -i use comma as separator. <br>
> > Use basename1_comma.sh to do that for multiple files.

```
for j in *$1 ; do basename $j $1 ;done |sed ':a;N;$!ba;s/\n/,/g'
```

### STEP 7: (not necessary if you're using MEGAN to visualise)
To summarise the OTU BIOM files on genus and species level for raw counts (http://qiime.org/scripts/summarize_taxa.html):

```
summarize_taxa.py -i merged.biom -a -L 6,7 -o ./tax
```
> -a raw counts, don’t use -a if you want relative abundance <br>
> L level 6 is genus 7 is species <br>
> -o output directory <br>

<p> You get sth like this from both BIOM and TXT files (raw counts): <br>
  Constructed      from     biom  
  
 | OTU                                                                                                              | ID  | 1A10S17 | 1A11S1  |
 | -----------------------------------------------------------------------------------------------------------------|-----|---------|---------|
 | Unassigned;Other;Other;Other;Other;Other                                                                         | 0.0 |   0.0   |   0.0   |
 | k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanobrevibacter| 0.0 |   0.0   |   0.0   |
 | k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__Methanosphaera    | 0.0 |   0.0   |   0.0   |
 |k__Archaea;p__Euryarchaeota;c__Thermoplasmata;o__E2;f__[Methanomassiliicoccaceae];g__vadinCA11                    | 0.0 |   0.0   |   0.0   | 

</p>

Preferred way: You may also use MEGAN to visualise BIOM file as an alternatives (my preference) - <br>
it will give a cleaner view with only genus name/ species name.
