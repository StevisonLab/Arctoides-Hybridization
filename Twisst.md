# Topology Weighting Analysis

This file details the topology weighting analysis. This work relied heavily on this github repository, which desribes the method in depth: https://github.com/simonhmartin/twisst. 

1. "Spruce" trees

The topolgy weigths require the data to be in newick tree format. It is recommended to estimate trees in sliding windows of 50 SNPs at a time. But, this requires at least one sample from each lineage to have data for every SNP. Therefore, the first step was to "spruce" the raw VCF files to make sure that for every site, there was at least one representative from each defined population. This was done using a custom script "[spruce_trees.pl](https://github.com/StevisonLab/Arctoides-Hybridization/blob/main/scripts/spruce_trees.pl)" (see scripts folder). The input file was in the geno format from the [hybridization analysis](https://github.com/StevisonLab/Arctoides-Hybridization/blob/main/Hybridization_Analysis.md). The output was in the same format, since that input was needed for the next step. 

Code used was as follows:

```sh
perl spruce_trees.pl pop_list.txt geno_files/Macaque_merged.chr${i}.geno
```

2. Make newick formatted tree files

The python script "phyml_sliding_windows.py" from the [Genomics General repository](https://github.com/simonhmartin/genomics_general) was used to convert the geno formatted data into newick formatted files. See [here](https://github.com/simonhmartin/genomics_general#trees-for-sliding-windows) for details.

Code used was as follows:

```sh
python ${path}/phyml_sliding_windows.py -g geno_files/Macaque_merged.chr${i}.geno.spruced --prefix tree_files/chr${i} -w 50 --windType sites --model GTR --optimise n
```

3. Run Twisst to output topology weights. The python script "twisst.py" from the Twisst github repo linked above was used to generate weights for sample groupings. See [here](https://github.com/simonhmartin/twisst#weighting-method) for details.

Code used was as follows:

```sh
python ${path}/twisst.py -t tree_files/chr${i}.trees -w tree_files/chr${i}.output.weights.csv -g A -g B -g C -g D --method complete --groupsFile ${path}/sample_ids.txt --outputTopos tree_files/chr${i}.toplogies.trees
```

4. Merge tree files from separate chromosomes. The scripts above were run individually on each chromosome. The following code was used to merge the genome into a single file for downstream analysis.

Code used was as follows:

```sh
awk 'NR>3' chr1.output.weights.csv > merged.output.weights.tsv
for i in {2..20} 'X'; do awk 'NR>4' chr${i}.output.weights.csv;done >> merged.output.weights.tsv

cat chr1.data.tsv > merged.data.tsv
for i in {2..20} 'X'; do awk 'NR>1' chr${i}.data.tsv;done >> merged.data.tsv

wc -l merged.data.tsv merged.output.weights.tsv

paste merged.data.tsv merged.output.weights.tsv >merged.twisst_final.tsv
```

5. Assign top topology based on weights. In various sections of the manuscript, rather than using the raw topology weights, we only applied a topology label if the weight was greater than 50%, 2/3 majority or 100%. This was done using a custom script "[top_topos.pl](https://github.com/StevisonLab/Arctoides-Hybridization/blob/main/scripts/top_topos.pl)" (see scripts folder). This script was hard coded to the output from step #4. To adjust weight percentage required for the labels to be applied, the variable $cutoff was used in the script. If it is set at 1, it corresponds to 100% weights must be in a specific topology for the label to be applied.  
