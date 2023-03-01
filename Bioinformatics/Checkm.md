## Using Checkm for contamination 

```
source checkm-1.1.3
```
 
To run workflow and generate marker file:

```
sbatch --wrap "checkm lineage_wf checkm/ checkm-output/ -t 20 -x fna" -c 20 --mem=400GB -p qib-medium -o checkm.out
```

checkm/ is the folder with all the bins/ genome assemblies.
-x is the suffix of the bins/genomes in the checkm/ folder.
This workflow will generate lineage.ms as a marker file for other analyses including quality assessment as below.
To run quality assessment:

```
sbatch --wrap "checkm qa -t 20 -a lineage.aln -f qa-table ./checkm-output/lineage.ms ./checkm-output/" -c 20 --mem=50GB -p qib-medium -o qa.out
```

checkm-output is the folder with all the outputs from the lineage_wf 
-a is to generate alignment file (not very useful?)
-f output is in qa-table
