## Check if a bacterial chromosome is circular
> This was needed for LH1062. NBI HPC resources used.

```
source circlator-1.5.5_RK

sbatch --wrap "circlator all --merge_min_id 85 --merge_breaklen 1000 --verbose --threads 32 consensus.fasta Bpseudolongum44.filtlong.fastq circlator-nano-Bpseudolongum44" -c 32 --mem=20GB -p qib-medium -o circlator-nano.out
```
