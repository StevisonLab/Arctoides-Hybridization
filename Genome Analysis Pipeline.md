# Genome Analysis Pipeline 

This file details the steps in genome analysis of raw reads. 

A visual overview of the pipeline is shown below:

![](https://github.com/StevisonLab/Arctoides-Hybridization/blob/main/GenomeAnalysisWorkflow.png)

1. Public samples were downloaded from NCBI via the tool `sratoolkit` version 2.8.1 of the Alabama Super Computer. For each sample, the raw sra run table was downloaded from ncbi sra website. This was used to generate a list of each sra that would be downloaded via the command line.

Code used was as follows:
```sh
#make array of sra numbers
sra_list=($(awk 'NR>1 {print $1}' sra_list.txt))

#extract in a loop (each run as a separate job)
fastq-dump --origfmt --split-files --gzip -I ${sra_list[$i]}
```

2. Several public sras were combined data across multiple lanes of sequencing. To conduct downstream analysis that was read group aware, these fastq file were split based on the fastq headers into separate files. 

Code used was as follows:
```sh
#uses same array as above
awk \"BEGIN {FS = \\\":\\\"} {lane_id=\\\$1.\\\$2 ; print > \\\"lane.${sra_list[${i}]}.\\\"lane_id\\\".fastq\\\" ; for (i = 1; i <= 3; i++) {getline ; print > \\\"lane.${sra_list[${i}]}.\\\"lane_id\\\".fastq\\\"}}\" ${sra_list[${i}]}.fastq
```

3. Prior to alignment, the reference was downloaded from UCSC and indexed via bwa, samtools and picard tools. 

Code used was as follows:
```sh
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/rheMac8/bigZips/rheMac8.fa.masked.gz' -O rheMac8.fa.masked.gz

#uncompressed and renamed to rheMac8.masked.fa for picard tools to recognize fasta file extension

#indexing
ref="rheMac8.masked.fa"
bwa index -p rheMac8 -a bwtsw ${path}/${ref}
samtools faidx ${path}/${ref}
java -Xms2g -Xmx10g -jar $dict_jar R=${path}/${ref} OUTPUT=${path}/rheMac8.masked.dict
```

4. Once fastq files were split based on lane information and genome was properly indexed, alignment of each file to the reference genome was conducted using BWA version 0.7.12. 

Code used was as follows:
```sh
#make separate alignment jobs per SRA file
for i in ${sra_list[@]}
do
    echo "$i..."
    echo "#!/bin/sh

    #run alignment to get sorted bam
    bwa mem -M -v 2 -t 8 -R \"@RG\tID:$flowcell.$lane.$i\tSM:$name.$species\tPU:$flowcell.$lane.$sample\tPL:Illumina\tLB:$library\" rheMac8 ${i}.fastq.gz  | samtools view -Sb | samtools sort - -@ 8 -m 2GB  >$i.sorted.bam

    #index sorted bam file
    samtools index $i.sorted.bam
    " >${path}/$i.alignment.sh
    chmod +x ${path}/$i.alignment.sh
done
#run each alignment script in parallel until completed
```

5. Merge alignments at the sample level into single bam file and sort/index.

Code used was as follows:
```sh
#make separate alignment jobs per individual sample
for f in sample1 sample2 .. sample10
do

    echo $f
    echo "#!/bin/sh

files=(\$(awk -v species=$f '\$6==species {print \"I=\" \$1 \".sorted.bam\"}' fastq_list.txt | sed 's/_R1//'))
echo \"Files to merge: \${files[@]}\"

#merge runs separately on 1 core, small queue, 2gb ram
java -Xms2g -Xmx4g -jar $mergeSam_path \${files[@]} OUTPUT=$f.merged.bam VALIDATION_STRINGENCY=SILENT    
echo Done merging.

#run sort/index separatly in large queue, 32gb ram, 8 cores
samtools sort -m 2GB -@ 8 $f.merged.bam >$f.merged.sorted.bam
echo Done sorting, now indexing

samtools index $f.merged.sorted.bam
echo Done.
" >${path}/$f.merge.sh
done
#run each merge script in parallel until completed
```

6. Run GATK to mark duplicates and complete indel realignment

Code used was as follows:
```sh
#make separate GATK jobs per individual sample
sample_list=($(awk '{print $1}' sample_list.txt)) #10 total

for i in ${sample_list[@]}
do
    echo "#!/bin/sh
    echo "echo $i...

#mark duplicates
#java -Xms2g -Xmx4g -jar $md_jar INPUT=${init_bam} OUTPUT=${md_bam} REMOVE_DUPLICATES=false METRICS_FILE=${md_bam}.dup_metrics VALIDATION_STRINGENCY=SILENT
    
#build bam index
#java -Xms2g -Xmx4g -jar $index_jar INPUT=${md_bam} OUTPUT=${md_bam}.bai VALIDATION_STRINGENCY=SILENT
    
#run GATK to do indel realignment, need to target regions that need realignment first
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T RealignerTargetCreator -I ${md_bam} -o ${md_bam}.intervals
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T IndelRealigner -I ${md_bam} -o ${realign_bam} -targetIntervals ${md_bam}.intervals

echo Done.
" >>GATK.$i.sh
    chmod +x GATK.$i.sh
done
#run each GATK script in parallel until completed
```

7. BQSR step required initial round of variant calling to get set of known variants. These were then used to run BQSR and output metrics (see PDFs in QC_stats folder). If more than 1% new variants found, then a second iteration was done. See supplement for full explanation.

```sh
#make separate GATK jobs per individual sample
species_list=($(awk '{print $4}' sample_list.txt | sort | uniq))

#initialize filenames
#bams
merge_bam="$i.species.merged.markdup.realigned.bam"
merge_sorted_bam="$i.species.merged.sorted.markdup.realigned.bam"
recal_1_bam="$i.merged.sorted.markdup.realigned.recalibrated1.bam"
recal_2_bam="$i.merged.sorted.markdup.realigned.recalibrated2.bam"
    
#vcfs
init_vcf="$i.initial.vcf"
filt_1_vcf="$i.filtered1"
sec_vcf="$i.second.vcf"
filt_2_vcf="$i.filtered2"
third_vcf="$i.third.vcf"
filt_3_vcf="$i.filtered3"
final_vcf="$i.erc_mode.g.vcf"

for i in ${species_list[@]}
do
    echo "#!/bin/sh

#merge realigned bams into one for the species
files=(\$(awk -v species=$i '\$4==species {print \"I=\" \$1 \".merged.sorted.markdup.realigned.bam\"}' sample_list.txt))

echo \"Files to merge: \${files[@]}\"
   if [ \${#files[@]} -eq 1 ];then
       name=(\$(awk -v species=$i '\$4==species {print \$1 \".merged.sorted.markdup.realigned.bam\"}' sample_list.txt))
       mv \$name ${merge_sorted_bam}
   else
       java -Xms2g -Xmx4g -jar $mergeSam_path \${files[@]} OUTPUT=${merge_bam} VALIDATION_STRINGENCY=SILENT    
       samtools sort -m 2GB -@ 8 ${merge_bam} >${merge_sorted_bam}
   fi

#index bam
samtools index ${merge_sorted_bam}

#initial variant calling with HC
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T HaplotypeCaller -I ${merge_sorted_bam} -nct 8 -o ${init_vcf}
gzip ${init_vcf}

#hard filtering of variants
vcftools --gzvcf ${init_vcf}.gz --minGQ 40 --minDP 15 --maxDP 160 --minQ 50 --recode --recode-INFO-all --out ${filt_1_vcf}

#BQSR Round 1
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T BaseRecalibrator -I ${merge_sorted_bam} -knownSites ${filt_1_vcf}.recode.vcf -o $i.recal_data1.table
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T BaseRecalibrator -I ${merge_sorted_bam} -knownSites ${filt_1_vcf}.recode.vcf -BQSR $i.recal_data1.table -o $i.post_recal_data1.table
java -Xms2g -Xmx4g -jar $gatk_jar -T AnalyzeCovariates -R $ref -before $i.recal_data1.table -after $i.post_recal_data1.table -plots $i.recalibration_plots1.pdf 
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T PrintReads -I ${merge_sorted_bam} -BQSR $i.recal_data1.table -o ${recal_1_bam}

#2nd round of HC
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T HaplotypeCaller -I ${recal_1_bam} -nct 8 -o ${sec_vcf}
gzip ${sec_vcf}

#hard filter variants again
vcftools --gzvcf ${sec_vcf}.gz --minGQ 40 --minDP 15 --maxDP 160 --minQ 50 --recode --recode-INFO-all --out ${filt_2_vcf}

#compare vcfs
bgzip -f ${filt_1_vcf}.recode.vcf
bgzip -f ${filt_2_vcf}.recode.vcf
tabix -f -p vcf ${filt_1_vcf}.recode.vcf.gz
tabix -f -p vcf ${filt_2_vcf}.recode.vcf.gz
vcf-compare ${filt_1_vcf}.recode.vcf.gz ${filt_2_vcf}.recode.vcf.gz >${i}_BQSR_round1.txt

#BQSR Round 2                                                                                                                        
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T BaseRecalibrator -I ${merge_sorted_bam} -knownSites ${filt_2_vcf}.recode.vcf.gz -o $i.recal_data2.table
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T BaseRecalibrator -I ${merge_sorted_bam} -knownSites ${filt_2_vcf}.recode.vcf.gz -BQSR $i.recal_data2.table -o $i.post_recal_data2.table       
java -Xms2g -Xmx4g -jar $gatk_jar -T AnalyzeCovariates -R $ref -before $i.recal_data2.table -after $i.post_recal_data2.table -plots $i.recalibration_plots2.pdf
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T PrintReads -I ${merge_sorted_bam} -BQSR $i.recal_data2.table -o ${recal_2_bam}

#3rd round of HC
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T HaplotypeCaller -I ${recal_2_bam} -nct 8 -o ${third_vcf}
gzip ${third_vcf}

#hard filter variants again
vcftools --gzvcf ${third_vcf}.gz --minGQ 40 --minDP 15 --maxDP 160 --minQ 50 --recode --recode-INFO-all --out ${filt_3_vcf}

#compare vcfs again
bgzip -f ${filt_3_vcf}.recode.vcf
tabix -f -p vcf ${filt_3_vcf}.recode.vcf.gz
vcf-compare ${filt_2_vcf}.recode.vcf.gz ${filt_3_vcf}.recode.vcf.gz >${i}_BQSR_round2.txt

echo Done.
" >>BQSR.$i.sh
#run each BQSR script in parallel until completed
```

8. Collect QC metrics

Samtools depth, samtools flagstat, and vcftools were used to collect various metrics throughout the pipeline and output into centralized table called qc.summary.txt (see folder [QC_stats](https://github.com/StevisonLab/Arctoides-Hybridization/tree/main/QC_stats)). For each individual sample, the raw read count from the combined SRR files and read length were computed. These were used to calculate an initial coverage estimate using the genome size from the fai index file (see below). Using the flagstat output, the read count and mapped read count were used to calculate a percentage of reads that mapped. This too was used to get initial coverage following mapping. The flagstat output from the bam file with duplicate reads marked was used to calculate a percentage duplicated for each sample. A flagstat report after indel realignment was generated to compare before and after. Any differences were reported as an error. Vcf coverges was calculated using vcftools option `--site-mean-depth`. 

Genome size was calculated as follows:
```sh
genome_size=`awk '{sum+=$2} END {print sum}' ${path}/rheMac8.masked.fa.fai`
```

9. VQSR: The output of the third round of HC following BQSR (see #7 above) was used to conduct VQSR. See supplementary methods for specifics. VQSR like BQSR was conducted at the species level rather than the sample level.


10. Haplotype Caller for each sample for downstream analysis.

Code used was as follows:
```sh
#HC for sample $f
java -Xms2g -Xmx4g -jar $gatk_jar -R $ref -T HaplotypeCaller -I ${recal_2_bam} -ERC GVCF -nct 8 -o $f.${final_vcf} --sample_name $sample_ID
bgzip $f.${final_vcf}
tabix -p vcf $f.${final_vcf}.gz
```

11. Genotype GVCFs. 

In conducting this step, we found the GATK program to be extremely slow due to the large number of scaffolds in the whole genome file. Therefore, we used this [script](https://github.com/lstevison/random-scripts/blob/master/split_by_chr.pl) to split the reference genome fasta file into chromosomes. We then combined only the chromosome level scaffolds (chr1-chr20 and chrX and chrY) into a new fasta file that was used as the reference moving forward. The same index steps in #3 above were conducted on the new fasta file. Similarly, the additional scaffold labels were removed from the final vcf in Step 9 above. 

Code used was as follows:

```sh
species=($(awk '{print $4}' sample_list.txt))
sample=($(awk '{print $2}' sample_list.txt))

zcat ${sample[$s]}.${species[$s]}.erc_mode.g.vcf.gz | grep -v "^##contig=<ID=chrUn_" > ${sample[$s]}.${species[$s]}.fixed.erc_mode.g.vcf
bgzip ${sample[$s]}.${species[$s]}.fixed.erc_mode.g.vcf
tabix -p vcf ${sample[$s]}.${species[$s]}.fixed.erc_mode.g.vcf.gz

#to further speed up, genotypeGVCF was run per chromosome for each sample
for i in {1..20} 'X' 'Y' 'M'
do

    echo "#!/bin/sh
    #genotype GVCF

java -Xms2g -Xmx20g -jar $gatk_jar -R $ref -T GenotypeGVCFs --variant ${sample[0]}.${species[0]}.fixed.erc_mode.g.vcf.gz --variant ${sample[1]}.${species[1]}.fixed.erc_mode.g.vcf.gz --variant ${sample[2]}.${species[2]}.fixed.erc_mode.g.vcf.gz --variant ${sample[3]}.${species[3]}.fixed.erc_mode.g.vcf.gz --variant ${sample[4]}.${species[4]}.fixed.erc_mode.g.vcf.gz --variant ${sample[5]}.${species[5]}.fixed.erc_mode.g.vcf.gz --variant ${sample[6]}.${species[6]}.fixed.erc_mode.g.vcf.gz --variant ${sample[7]}.${species[7]}.fixed.erc_mode.g.vcf.gz --variant ${sample[8]}.${species[8]}.fixed.erc_mode.g.vcf.gz --variant ${sample[9]}.${species[9]}.fixed.erc_mode.g.vcf.gz -allSites -nt 8 -L chr${i} -o split_by_chr/Macaque_merged.chr${i}.final.g.vcf

echo \"Compressing chr${i} file...\"

bgzip -f split_by_chr/Macaque_merged.chr${i}.final.g.vcf
tabix -p vcf split_by_chr/Macaque_merged.chr${i}.final.g.vcf.gz
    
echo \"Finished with chr${i}.\" " >GVCF.$i.sh

    chmod +x GVCF.$i.sh

done
```

12. Apply mask. Combining the VCF files from Steps 9 and 11, a near final callset was generated. The file $v below is the output from VQSR from Step 9. These positions were extracted and filtered out of the combined VCF file from Step 11.

Code used was as follows:
```sh
sample_list=($(awk '{print $1}' sample_list.txt)) #10 total
species_list=($(awk '{print $4}' sample_list.txt | sort | uniq))

#extract variant positions per individual per chromosome
zcat $v | grep \"^chr${i}\" | awk 'BEGIN {OFS=\"\t\"} {print \$1,\$2}' >${path}/${species[$s]}.chr${i}.variant_positions.txt

#extract those postions from vcf
vcftools --gzvcf ${path}/Macaque_merged.chr${i}.final.g.vcf.gz --min-alleles 2 --max-alleles 2 --positions ${path}/${species[$s]}.chr${i}.variant_positions.txt --keep ${species[$s]}.sample_list.txt --recode --recode-INFO-all --out ${path}/${species[$s]}.chr${i}.variant

bgzip -f ${path}/${species[$s]}.chr${i}.variant.recode.vcf
tabix -p vcf ${path}/${species[$s]}.chr${i}.variant.recode.vcf.gz

#positions were used to generate a multiIntersect bed file
zcat $v | grep \"^chr${i}\" | awk 'BEGIN {OFS=\"\t\"} {print \$1, \$2-1, \$2+length(\$4)-1}' >${path}/${species[$s]}.chr${i}.variant_positions.bed
multiIntersectBed -header -i ${species_list[0]}.chr${i}.variant_positions.bed ${species_list[1]}.chr${i}.variant_positions.bed ${species_list[2]}.chr${i}.variant_positions.bed ${species_list[3]}.chr${i}.variant_positions.bed ${species_list[4]}.chr${i}.variant_positions.bed -names ${species_list[0]} ${species_list[1]} ${species_list[2]} ${species_list[3]} ${species_list[4]} > Macaque.chr${i}.variant_positions.intersected.bed

#use resulting bed file to extract invariant sites from VCF
grep -v \"$${species_list[$s]}\" Macaque.chr${i}.variant_positions.intersected.bed | awk 'BEGIN {OFS=\"\t\"} {print \$1, \$2, \$3}' >${species_list[$s]}.chr${i}.invariant_positions.bed

bcftools view -T ${species_list[$s]}.chr${i}.invariant_positions.bed -S ${species_list[$s]}.sample_list.txt -i 'FMT/DP>5 & FMT/GQ>20' Macaque_merged.chr${i}.final.g.vcf.gz | bcftools view -e 'GT=\"0/1\" | GT=\"1/0\" | GT=\"1/1\"' -o ${species_list[$s]}.chr${i}.invariant.vcf

bgzip -f ${species_list[$s]}.chr${i}.invariant.vcf 
tabix -p vcf ${species_list[$s]}.chr${i}.invariant.vcf.gz

#merge variant and invariant sites for ${species_list[$s]}
vcf-concat ${species_list[$s]}.chr${i}.invariant.vcf.gz ${species_list[$s]}.chr${i}.variant.recode.vcf.gz | vcf-sort | bgzip -c > ${species_list[$s]}.chr${i}.combined.vcf.gz
tabix -p vcf ${species_list[$s]}.chr${i}.combined.vcf.gz

#merge samples back into one VCF
bcftools merge ${species_list[0]}.chr${i}.combined.vcf.gz ${species_list[1]}.chr${i}.combined.vcf.gz ${species_list[2]}.chr${i}.combined.vcf.gz ${species_list[3]}.chr${i}.combined.vcf.gz ${species_list[4]}.chr${i}.combined.vcf.gz ${species_list[5]}.chr${i}.combined.vcf.gz ${species_list[6]}.chr${i}.combined.vcf.gz ${species_list[7]}.chr${i}.combined.vcf.gz >Macaque_merged.chr${i}.filtered.vcf

bgzip -f Macaque_merged.chr${i}.filtered.vcf
tabix -f -p vcf Macaque_merged.chr${i}.filtered.vcf.gz
```

These final filtered VCF files have been uploaded to the Auburn University institutional respository [AUrora](https://aurora.auburn.edu/) at the following [DOI](http://dx.doi.org/10.35099/aurora-67). Any use of these files should cite this github repository, the Aurora DOI and our [BioRxiv preprint](https://www.biorxiv.org/content/10.1101/2020.05.18.102251v1.abstract) describing these data.  
