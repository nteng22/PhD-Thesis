## Screening genomes using Abricate
> Definitely there is a better way - use Abricate v1.0.1 which is on HPC. <br>
  Essentially a Blast tool but more advanced and much easier to use, written in perl. <br>
  You can specify parameters for coverage --mincov and identity --minid and get outcomes within seconds! <br>
  The only thing you need to do is to supply a database of your genes of interest which you can look for in either NCBI or 
  uniprot - has to be nucleotides not amino acid sequences as it is a blastn based tool. <br>
  You can now look for your genes of interest, get those nucleotide sequences and compile a multi fasta file <br>

```
source abricate-1.0.1
source perl-5.22.1
source emboss-6.5.7
source any2fasta-0.4.2
source package blast+-2.6.0

abricate --db [specify database] [fasta file]
```

https://github.com/tseemann/abricate

Need to set up a database if you're making your own:
```
cd /path/to/abricate/db     # this is the --datadir default option
mkdir tinyamr
cp /path/to/tinyamr.fa sequences
head -n 1 sequences
>tinyamr~~~GENE_ID~~~GENE_ACC~~RESISTANCES some description here
abricate --setupdb

# or just do this: 
makeblastdb -in sequences -title tinyamr -dbtype nucl -hash_index

% abricate --list
DATABASE  SEQUENCES  DBTYPE  DATE
tinyamr   173        nucl    2019-Aug-28

% abricate --db tinyamr screen_this.fasta
```
> makeblastdb has worked for me while --setupdb hasn't in the past.
