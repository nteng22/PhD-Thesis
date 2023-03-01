Unicycler is a different de novo assembler.
It works well on short and long reads.
Can also do a hybrid using both short and long reads.

### To run you need to load the dependencies (on the NBI HPC)

```
source package python-3.5.1
source package unicycler-0.4.9b
source package pilon-1.22
source package spades-3.11
source package jre-1.8.0_45
source package gcc-5.2.0
source pakcage racon-1.3.1
source package bowtie2-2.1.0 
source package samtools-1.4.1  
source package blast+-2.9.0
```

>Don't run on interactive mode
---
### Illumina-only assembly:
```
unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -o output_dir
```

### Long-read-only assembly:
```
unicycler -l long_reads.fastq.gz -o output_dir
```

### Hybrid assembly:

```
unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -l long_reads.fastq.gz -o output_dir
```

---
Step 1: filter out poor quality reads
source filtlong-0.2.1 
 
```
sbatch --wrap "filtlong --min_length 1000 --keep_percent 90 Bpseudolongum44.fastq > Bpseudolongum44.filtlong.fastq" -c 4 --mem=4GB -p qib-short -J filtlong -o filtlong.out
```
Basically discard the worst 10% of reads which is reasonable, or 5%. Only keep reads > 1000bp.
Then use Flye next to assemble.
 
Step 2: run assembler
```
source flye-2.9
 
sbatch --wrap "flye --nano-raw Bpseudolongum44.fastq -o flye-assembly --scaffold -g 2.0m -i 5 -t 12" -c 12 --mem=30GB -p qib-medium -o flye.out
```
--scaffold to enable scaffolding, in v2.9 by default it will not scaffold. -i polishing iteration -t cpus -g is the estimated genome size.

---
### Flye assembly
```

flye --nano-hq
```
Best ways to assemble Nanopore sequences with high coverage (>100x):
1. Use filtlong to filter worst 10% and <1000bp reads in the fastq (called light QC)
2. Use flye --nano-hq is guppy basecalling version 5+, polishing with 5 iterations - usually generates very good genome assemblies.
