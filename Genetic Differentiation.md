# Genetic Differentiation Analysis

This file details the steps to generate DXY in sliding windows across the genome to compare the divergence of regions with shared ancestry with the fascicularis species group to regions with shared ancestry with the sinica species group.

1. Following the Step 5 in the [Twisst Analysis](https://github.com/StevisonLab/Arctoides-Hybridization/blob/main/Twisst.md), awk was used to extract the coordinates of the regions with 100% weights in each of the three topologies. 

Code used was as follows:

```sh
awk '$4=="topo1" {OFS="\t"; print $1, $2, $3}' haplo_summary.tsv >topo1.bed
awk '$4=="topo2" {OFS="\t"; print $1, $2, $3}' haplo_summary.tsv >topo2.bed
awk '$4=="topo3" {OFS="\t"; print $1, $2, $3}' haplo_summary.tsv >topo3.bed
```
2. The python script "popgenWindows.py" from the [Genomics General repository](https://github.com/simonhmartin/genomics_general) was used to estimate divergence as DXY in each region defined in the bed files from Step 1. This was done separately for each chromosome and then combined across the genome.

Code used was as follows:

```sh
for i in {1..20}
do
    grep "^chr${i}" topo1.bed >chr${i}.topo1.bed
    grep "^chr${i}" topo2.bed >chr${i}.topo2.bed
    grep "^chr${i}" topo3.bed >chr${i}.topo3.bed

    input="${path}/Macaque_merged.chr${i}.geno.gz"

    echo "Running pop gen windows for chr${i}..."
    python $path/popgenWindows.py --popsFile $pops --windCoords chr${i}.topo1.bed -g $input -o outputs/chr${i}.topo1.windows.csv.gz -f phased -m 50 -T 4 --windType predefined --writeFailedWindows -p Arctoides -p Sinica -p Baboon -p Fascicularis
    python $path/popgenWindows.py --popsFile $pops --windCoords chr${i}.topo2.bed -g $input -o outputs/chr${i}.topo2.windows.csv.gz -f phased -m 50 -T 4 --windType predefined --writeFailedWindows -p Arctoides -p Sinica -p Baboon -p Fascicularis
    python $path/popgenWindows.py --popsFile $pops --windCoords chr${i}.topo3.bed -g $input -o outputs/chr${i}.topo3.windows.csv.gz -f phased -m 50 -T 4 --windType predefined --writeFailedWindows -p Arctoides -p Sinica -p Baboon -p Fascicularis
done
```

3. A combined file for all regions across the genome for each topology was generated for downstream analysis.

Code used was as follows:

```sh
#make combined file for each topo
zcat outputs/chr1.topo1.windows.csv.gz >outputs/combined.topo1.csv
for i in {2..20}; do zcat outputs/chr${i}.topo1.windows.csv.gz | awk 'NR>1' >>outputs/combined.topo1.csv; done 

zcat outputs/chr1.topo2.windows.csv.gz >outputs/combined.topo2.csv
for i in {2..20}; do zcat outputs/chr${i}.topo2.windows.csv.gz | awk 'NR>1' >>outputs/combined.topo2.csv; done 

zcat outputs/chr1.topo3.windows.csv.gz >outputs/combined.topo3.csv
for i in {2..20}; do zcat outputs/chr${i}.topo3.windows.csv.gz | awk 'NR>1' >>outputs/combined.topo3.csv; done
```
