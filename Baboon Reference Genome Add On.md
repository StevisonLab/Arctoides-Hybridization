#This file details the steps to combine the final VCFs with the baboon reference genome.

1. The reference genome file of each individual chromosome from the rheMac8 genome were compared agains the papAnu4 reference genome which was downloaded from UCSC.

Code used was as follows:

```sh
for i in {1..20}  # 'X' 'M' #must manually edit X and M for lastz subsample and female file for chrX
do
    echo "#! /usr/bin/sh

module load lastz

lastz ${path}/rheMac8.masked.fa[multiple,subsample=${i}/23] ${path}/papAnu4_subset.fa.masked --identity=75 --notransition --step=10 --gapped --chain --gfextend --format=maf+ >chr${i}.rheMac8.papAnu4.maf

#gunzip chr${i}.rheMac8.papAnu4.maf

~/bin/parse_maf_new.pl chr${i}.rheMac8.papAnu4.maf chr${i}.rheMac8.papAnu4.vcf

gzip chr${i}.rheMac8.papAnu4.maf

#gunzip chr${i}.rheMac8.papAnu4.vcf

sort -k2n -k6rn chr${i}.rheMac8.papAnu4.vcf >chr${i}.sorted.tmp

gzip chr${i}.rheMac8.papAnu4.vcf

awk '{ if(\$2!=x) {print \$0; x=\$2} }' chr${i}.sorted.tmp > chr${i}.reduced.tmp

rm chr${i}.sorted.tmp

awk '\$5!~\"-\" && \$6!~\"-\"' chr${i}.reduced.tmp | cat vcf_header.txt - >chr${i}.rheMac8.papAnu4.sorted.vcf

rm chr${i}.reduced.tmp

gunzip ~/Macaque_New/split_by_chr/Macaque_merged.chr${i}.filtered.vcf.gz

~/bin/intersect_parsed_vcf_baboon_data.pl ~/Macaque_New/split_by_chr/Macaque_merged.chr${i}.filtered.vcf chr${i}.rheMac8.papAnu4.sorted.vcf ~/Macaque_New/split_by_chr/intersected_vcf/chr${i}.filtered.wBaboon.vcf

wc -l ~/Macaque_New/split_by_chr/intersected_vcf/chr${i}.filtered.wBaboon.vcf ~/Macaque_New/split_by_chr/Macaque_merged.chr${i}.filtered.vcf chr${i}.rheMac8.papAnu4.sorted.vcf

gzip -f ~/Macaque_New/split_by_chr/Macaque_merged.chr${i}.filtered.vcf ~/Macaque_New/split_by_chr/intersected_vcf/chr${i}.filtered.wBaboon.vcf chr${i}.rheMac8.papAnu4.sorted.vcf

" > chr${i}.lastz.sh

  chmod +x chr${i}.lastz.sh
  run_script chr${i}.lastz.sh                                                                                                                  

done```

