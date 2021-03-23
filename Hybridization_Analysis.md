# Hybridization Analysis
This file details the steps used to calculate D and fdm on the macaque genome data. For the most part, the following github repo was used heavily: https://github.com/simonhmartin/genomics_general

1. Using the python script from Simon Martin, the VCF file was first converted to a "geno" formatted file. This was done per chromosome.

Code used was as follows:

```sh
for i in 'X' {1..20} 
do
    echo "#!/bin/sh

python ${path}/VCF_processing/parseVCF.py -i ${path}/chr${i}.filtered.wBaboon.vcf.gz --skipIndels --skipMono | gzip > geno_files/Macaque_merged.chr${i}.geno.gz
" >chr${i}.geno.sh
    chmod +x chr${i}.geno.sh
done
```

2. Next, the program `ABBABABAwindows.py` was used to compute D and fdm in sliding windows of various size scales. This was done in a nested loop.

Code used was as follows:

```sh
#minimum
m=10

#Run ABBA-BABA stats and F-stats for five window sizes; Figure 4, Supplementary Figures 2-5
for size in '5000' '25000' '50000' '100000' '500000' '1000000'
do

   kb=`echo "scale=0 ; $size / 1000" | bc -l`
   min=`echo "( $m * ( $kb / 5 ) ) /1" | bc`
   step=`echo "$size / 5" | bc`
   echo "Now running bin sizes $kb kb; minimum number of sites will be $min; step size is $step"

   for i in 'X' {1..20} 
   do
       echo Now running on chr${i}...
echo "#!/bin/sh

echo \"Now running bin sizes $kb kb; minimum number of sites will be $min; step size is $step\"

python ${path}/ABBABABAwindows.py -g geno_files/Macaque_merged.chr${i}.geno.gz -f phased -P1 Sinica -P2 Fascicularis -P3 Arctoides -O Baboon -w $size -m $min -s $step -o csv_outputs/Macaque_merged.chr${i}.output.${kb}kb.csv -T 5 --writeFailedWindows --popsFile pop_list.txt
" >chr${i}.${kb}.fstats.sh
	
	chmod +x chr${i}.${kb}.fstats.sh
	run_script chr${i}.${kb}.fstats.sh
    done
done
```

3. The csv outputs were combined into a single output for each scale across all chromosomes. This combined output was then used to plot the resulting values in a heatmat across the genome using the template R script "plot_heatmap_new.r" (see scripts folder). The resulting image is in the main text as Figure 2 for the 50,000 bp or 50kb scale. Additional scale plots are found in the supplement.

Code used was as follows:

```sh
for size in '5000' '25000' '50000' '100000' '500000' '1000000'
do
size=50000
kb=`echo "scale=0 ; $size / 1000" | bc -l`
echo Now running bin sizes $kb kb
   
    cat Macaque_merged.chr1.output.${kb}kb.csv >combined.${kb}kb.txt
    for i in {2..20} 'X'; do awk 'NR>1' Macaque_merged.chr$i.output.${kb}kb.csv >>combined.${kb}kb.txt; done
    sed -e "s/50kb/${kb}kb/g" plot_heatmap_new.r >plot_heatmap.${kb}kb.r

    echo "#!/bin/sh                                                               
echo Now making ${kb} plot...

R --vanilla <plot_heatmap.${kb}kb.r

echo Done.
" >plot.${kb}.sh
    chmod +x plot.${kb}.sh
    run_script plot.${kb}.sh
done
```

