# Baboon Reference Genome
This file details the steps to combine the final VCFs with the baboon reference genome.

1. The reference genome file of each individual chromosome from the rheMac8 genome was aligned against the papAnu4 reference genome which was downloaded from UCSC in the same way as the rheMac8 reference genome using the software lastz.

Code used was as follows:

```sh
for i in {1..20}  
do
    echo "#! /usr/bin/sh

module load lastz

lastz ${path}/rheMac8.masked.fa[multiple,subsample=${i}/23] ${path}/papAnu4_subset.fa.masked --identity=75 --notransition --step=10 --gapped --chain --gfextend --format=maf+ >chr${i}.rheMac8.papAnu4.maf
" > chr${i}.lastz.sh

  chmod +x chr${i}.lastz.sh
  run_script chr${i}.lastz.sh                                                                                                                  

done
```

2. A custom perl script (see scripts folder: "[parse_maf_new.pl](https://github.com/StevisonLab/Arctoides-Hybridization/blob/main/scripts/parse_maf_new.pl)"), was used to convert the maf 2-way alignment for each chromosome into a mock-vcf formatted file. The output contained eight columns corresponding similarly to that in a VCF formatted file. The first column was the chromosome and then position along the chromosome using the rheMac8 values. There was then a column with a simple period corresponding to the ID column. The next two columns contained the reference bases for the corresponding position for rheMac8 and papAnu4 respectively. The lastz alignment score was placed into the 6th column and used to pick a top alignment below. The format field was as follows `CHROM.POS:%ID:ORI:LEN` and the next field contained a strong with those values as coded in the format field.

`~/parse_maf_new.pl chr${i}.rheMac8.papAnu4.maf chr${i}.rheMac8.papAnu4.vcf`

3. The resulting vcf-like file was sorted into a tmp file first by position, then reverse sorted by the laszt alignment score. Then awk was used to reduce the file to only one line per position (getting rid of alignments with lower scores), and finally filtered to remove any indels. In the last step, a barebones vcf header was added to the vcf-like file. The text of that file is below.

Code for these clean-up steps are as follows:
```sh
sort -k2n -k6rn chr${i}.rheMac8.papAnu4.vcf >chr${i}.sorted.tmp
awk '{ if($2!=x) {print $0; x=$2} }' chr${i}.sorted.tmp > chr${i}.reduced.tmp
awk '$5!~"-" && $6!~"-"' chr${i}.reduced.tmp | cat vcf_header.txt - >chr${i}.rheMac8.papAnu4.sorted.vcf
```

Text of vcf-like header file is:

>##fileformat=VCFv4.0

>#CHROM  POS     ID      rheMac8 papAnu4 QUAL    FILTER  INFO

4. Resulting baboon VCF file aligned to the rheMac8 reference was intersected with the final filtered vcf from the genome analysis pipeline (see Step #12). This was done using a custom perl script (see scripts folder: "[intersect_parsed_vcf_baboon_data.pl](https://github.com/StevisonLab/Arctoides-Hybridization/blob/main/scripts/intersect_parsed_vcf_baboon_data.pl)"), which read in the baboon VCF-like file and the merged/filtered data from the macaque samples. It then matched up each line of the two VCF files to output an intersected file with an additional column with baboon reference site. This site was coded as a genotype (e.g. alleles matching rheMac8 were coded as 0/0 and alleles not matching were coded as 1/1). Before coding alternate alleles, the value in the alternate allele column of the macaque VCF was compared against the baboon reference site. If they were the same, the same coding was used. Otherwise, it was coded as a new alternative allele matching VCF allele coding conventions.

`~/intersect_parsed_vcf_baboon_data.pl ~/${path}/Macaque_merged.chr${i}.filtered.vcf chr${i}.rheMac8.papAnu4.sorted.vcf ~/${path}/chr${i}.filtered.wBaboon.vcf`





