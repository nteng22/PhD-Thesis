> Downloaded PhyloPhlAn locally as the dependencies didn't work on the HPC. <br>
> Used bioconda and conda env to set up correctly <br>
> Notes written are for my sake in understanding what I did O.o

  ![PhyloPhlan workflow:](https://media.springernature.com/lw685/springer-static/image/art%3A10.1038%2Fs41467-020-16366-7/MediaObjects/41467_2020_16366_Fig1_HTML.png?as=webp)
  
### Step 1: Get reference genomes
Need to grab reference genomes, useful when you need to place your isolate amongst others.
```
phylophlan_get_reference \
    -g g__Coprobacter \
    -o input_isolates/ \
    --verbose 2>&1 | tee logs/phylophlan_get_reference.log # This is to get an up-to-date log <br>
    of what the computer is doing, if not run on HPC
```
> -g specifies what you want in the format "-g s__<species_name>" <br>
> if you do -g all, you can specify how many per X by -n

### Step 2: Generate configuration file
Generate configuration file prior to phylogeny tree... <br>

```
sbatch --wrap "phylophlan_write_default_configs.sh -i /qib/research-groups/Lindsay-Hall/Nancy/Novel_isolate/BP5B-C_tree/Genome_tree/Genus_genomes -d a -o ./phylophlan/LH1062.cfg --db_dna makeblastdb --map_dna blastn --msa mafft --tree1 fasttree" -c 8 -p qib-medium --mem=8GB
```
> Use the tutorial: https://github.com/biobakery/biobakery/wiki/PhyloPhlAn-3.0:-Example-05:-Proteobacteria 

```
sbatch --wrap "phylophlan -i ../ -d a --diversity medium --accurate -f supertree_nt.cfg -o ./LH1062_output --nproc 4 -t a" -c 4 -p qib-medium --mem=4GB 
```

> Ran phylophlan locally instead as the HPC didn't have the correct dependencies paths written in the code
> - https://github.com/biobakery/phylophlan/wiki#requirements
> - Uses python (conda) to set up the installation and required dependencies. So either install locally or re-install on the HPC. 

_You can use the default configuration file that PhyloPhlan has by executing:_

```
phylophlan_write_default_configs.sh [output_folder]
```

_Or you can write your specific configuration file:_
```
phylophlan_write_config_file \
    -d a \
    -o LH1063.cfg \
    --db_aa diamond \
    --map_dna diamond \
    --map_aa diamond \
    --msa mafft \
    --trim trimal \
    --tree1 iqtree \
    --verbose 2>&1 | tee logs/phylophlan_write_config_file.log

 -h, --help            Show this help message and exit
  -o OUTPUT, --output OUTPUT
                        Specify the output file where to write the
                        configurations (default: None)
  -d {n,a}, --db_type {n,a}
                        Specify the type of the database, where "n" stands for
                        nucleotides and "a" for amino acids (default: None)
  --db_dna {makeblastdb}
                        Add the "db_dna" section of the selected software that
                        will be used for building the indexed database
                        (default: None)
  --db_aa {usearch,diamond}
                        Add the "db_aa" section of the selected software that
                        will be used for building the indexed database
                        (default: None)
  --map_dna {blastn,tblastn,diamond}
                        Add the "map_dna" section of the selected software
                        that will be used for mapping the database against the
                        input genomes (default: None)
  --map_aa {usearch,diamond}
                        Add the "map_aa" section of the selected software that
                        will be used for mapping the database against the
                        input proteomes (default: None)
  --msa {muscle,mafft,opal,upp}
                        Add the "msa" section of the selected software that
                        will be used for producing the MSAs (default: None)
  --trim {trimal}       Add the "trim" section of the selected software that
                        will be used for the removal of the gappy regions of
                        the MSAs (default: None)
  --gene_tree1 {fasttree,raxml,iqtree}
                        Add the "gene_tree1" section of the selected software
                        that will be used for building the phylogenies for the
                        markers in the database (default: None)
  --gene_tree2 {raxml}  Add the "gene_tree2" section of the selected software
                        that will be used for refining the phylogenies
                        previously built with what specified in the
                        "gene_tree1" section (default: None)
  --tree1 {fasttree,raxml,iqtree,astral,astrid}
                        Add the "tree1" section of the selected software that
                        will be used for building the first phylogeny
                        (default: None)
  --tree2 {raxml}       Add the "tree2" section of the selected software that
                        will be used for refining the phylogeny previously
                        built with what specified in the "tree1" section
                        (default: None)
  -a, --absolute_path   Write the absolute path to the executable instead of
                        the executable name as found in the system path
                        environment (default: False)
  --force_nucleotides   If specified, sets parameters for phylogenetic analysis
                        software so that they use nucleotide sequences, even
                        in the case of a database of amino acids (default:
                        None)
  --overwrite           Overwrite output file if it exists (default: False)
  --verbose             Print more stuff (default: False)
  -v, --version         Print the current phylophlan_write_config_file.py
                        version and exit
```
### Step 3: Make the phylogenomic tree
Once the configuration file has been made you can proceed to making the tree
```
phylophlan \
    -i input_isolates/ \
    -d phylophlan \
    -f LH1063_CP.cfg \
    --diversity medium \
    --accurate \
    -o output_LH1063_CP \
    --nproc 1 \
    --verbose 2>&1 | tee logs/phylophlan_CP.log
```

> Oddly couldn't run the software with the "fasta" extension, <br>
> despite specifying with "--genome_extension", so change to "fna" for it to be recognised. <br>
> Also noted that if you messed up and want to re run you may need to clear out the "tmp" file. <br>
> Otherwise it will not mark everything properly and look in the "tmp" file and run steps from there. <br>
> Depending on --accurate or --fast it will have implications on the --diversity settings too. See manual.

#### Notes on Diversity parameter:
The --diversity parameter allows for three pre-defined options to set several <br>
parameters at once (e.g., trimming, subsampling, fragmentary removal, etc.) in accordance 
<br> with the expected diversity of the phylogeny to be built. </p>

The user can choose among three values: <br>
**Diversity	Description** <br>
- low	for species- and strain-level phylogenies
- medium	for genus- and family-level phylogenies
- high	for tree-of-life and higher-ranked taxonomic levels phylogenies

