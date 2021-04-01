# Topology Weighting Analysis

This file details the topology weighting analysis. This work relied heavily on this github repository, which desribes the method in depth: https://github.com/simonhmartin/twisst. 

1. "Spruce" trees

The topolgy weigths require the data to be in newick tree format. It is recommended to estimate trees in sliding windows of 50 SNPs at a time. But, this requires at least one sample from each lineage to have data for every SNP. Therefore, the first step was to "spruce" the raw VCF files to make sure that for every site, there was at least one representative from each defined population. This was done using a custom script "spruce_trees.pl" (see scripts folder). The input file was in the geno format from the [hybridization analysis](https://github.com/StevisonLab/Arctoides-Hybridization/blob/main/Hybridization_Analysis.md). The output was in the same format, since that input was needed for the next step. 

Code used was as follows:

```sh
perl spruce_trees.pl pop_list.txt geno_files/Macaque_merged.chr${i}.geno
```

2. Make newick formatted tree files

The python script "phyml_sliding_windows.py" from the Twisst repository linked above was used to convert the geno formatted data into newick formatted files.

3. 
