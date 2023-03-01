### GDTB-Tk, to determine taxonomies of isolates

```
source gtdbtk-1.5.1
```

<p> Because of a dependency called pplacer, you may only use 1 cpu unfortunately for pplacer, o<br>
therwise it will miscalculate memory as time so no matter how much memory you put it won't be enough even 250GB! some 2TB! </p>

#### To run complete workflow:
```
sbatch --wrap "gtdbtk classify_wf --cpus 20 -x fasta --pplacer_cpus 1 --genome_dir ./gtdb-analysis --out_dir gtdb-outcomes" -c 20 --mem=500GB -p qib-long -o gtdb.out
```

Run it with 500GB and see how it goes. If qib-long is not responding, use nbi-largemem (use snoderes to see which nodes are free). <br>
Put all genome assemblies in one folder in the above command it is gtdb-analysis <br>
--out_dir is the name of output directory will mkdir itself.  <br>
--cpus 20 is to assign 20 cpus to the software, --pplacer_cpus 1 is to assign 1 cpu to pplacer only.
 
#### Next steps after GDTB-Tk
Under the classify directory, there is a file called gtdbtk.bac120.summary.tsv <br>
In this tsv file, use 'less -S gtdbtk.bac120.summary.tsv' <br>
You can see the assigned taxon of each genome. <br>

### To make a phylogenomic tree based on Kmers (Using Mashtree)

```
source mashtree-1.2.0
export LC_ALL=C
```
> This version was installed by Raymond on 11/8/2021 previously was 0.30.0 by A.Page <br>

```
export LC_ALL=C 
# to avoid perl warning setting locale - not sure what is happening but it is working!

sbatch --wrap "mashtree --outmatrix mash --numcpus 12 [ALL FASTA/FNA FILES] > mash.tree" -c 12 --mem=15GB -p nbi-short -o mash-test.out
```

With bootstrapping:
```
sbatch --wrap "mashtree_bootstrap.pl --reps 100 --numcpus 12 ../*.fasta -- --min-depth 0 > BT_mashtree.tree" -c 12 --mem=15GB -p nbi-short -o BT_mashtree_230222.out
```
To generate the list of files above in the same directory, do this: <br>

```
/qib/research-groups/Lindsay-Hall/Raymond/script/basename1.sh
```
