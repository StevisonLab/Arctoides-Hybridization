#!/bin/bash

#PBS -l nodes=1:ppn=13,walltime=20:00:00,mem=52gb
#PBS -W x=FLAGS:ADVRES:lss0021_lab.197924
#PBS -q general
#PBS -m abe
#PBS -d .

module load bcftools/1.11 
module load gffread/2
module load samtools/1.11
module load gnu-parallel/20181022
module load seqkit/0.8.1
module load R/4.0.3

# Define array for sample names

samples=(`cat samples/sample_list.txt`)

# Create reference genomes for variant samples using vcf-consensus 

parallel 'cat References_files/rheMac8.masked.fa | bcftools consensus -M N -H 1 -s {} samples/Macaque_merged.allSNPs.vcf.gz > References_files/{}_1.masked.fa' ::: ${samples[@]}
parallel 'cat References_files/rheMac8.masked.fa | bcftools consensus -M N -H 2 -s {} samples/Macaque_merged.allSNPs.vcf.gz > References_files/{}_2.masked.fa' ::: ${samples[@]}

# Index the variant genomes to allow for easier extraction of genes

parallel 'samtools faidx References_files/{}_1.masked.fa' ::: ${samples[@]}
parallel 'samtools faidx References_files/{}_2.masked.fa' ::: ${samples[@]}

# Extract the coding regions all of transcripts in a given annotation

parallel 'gffread -C -x {}_1.coding.masked.fa -g References_files/{}_1.masked.fa References_files/rheMac8.ensGene.gtf' ::: ${samples[@]}
parallel 'gffread -C -x {}_2.coding.masked.fa -g References_files/{}_2.masked.fa References_files/rheMac8.ensGene.gtf' ::: ${samples[@]}

# Generate list of transcripts in order that they appear in the FASTA files (FASTA should be replaced by a random particular FASTA file

grep ">" Malaya_2.coding.masked.fa | sed 's/>//' > Transcripts.txt

# Use python script to change gene transcript names to include gene and transcript ID
# This python script is very slow so should be used only once on a list of transcript names obtained from a CDS FASTA
# This will create a reference of gene-transcript names in the same order as they appear in the CDS FASTA file

python3 ReplaceFastaNames.py Transcripts.txt References_files/ensemblToGeneName.txt > GenesNTranscripts.txt

# Create R scripts that will quickly replace gene names in each CDS FASTA and then run the script
# Both steps are very quick and may be run outside this script if needed
# The R script itself should be modified to reflect haplotype number

for x in ${samples[@]}; do

sed "s/SAMPLE/${x}/g" ReplaceFastaNames_1.R > ${x}.ReplaceFastaNames_1.R
sed "s/SAMPLE/${x}/g" ReplaceFastaNames_2.R > ${x}.ReplaceFastaNames_2.R

Rscript ${x}.ReplaceFastaNames_1.R
Rscript ${x}.ReplaceFastaNames_2.R

done

# Use the gene-transcript list generated above to define an array

genes=(`cat GenesNTranscripts.txt`)
mkdir gene_sequences

# Extract gene-transcript sequences from gene FASTA files and combine homologous genes from all samples into gene-specific files

parallel 'seqkit grep -p {1} {2}_1.genes.masked.fa | sed "s/>{1}/>{2}_1/" >> gene_sequences/{1}.fa' ::: ${genes[@]} ::: ${samples[@]}
parallel 'seqkit grep -p {1} {2}_2.genes.masked.fa | sed "s/>{1}/>{2}_2/" >> gene_sequences/{1}.fa' ::: ${genes[@]} ::: ${samples[@]}

# Convert all the above files into phylip files and output error messages to a text file

for x in ${genes[@]}; do 
echo "${x}" >> F2P_Output.txt
./Fasta2Phylip.pl gene_sequences/${x}.fa gene_sequences/${x}.phy >> F2P_Output.txt
done

# Go to gene_sequences folder and reformat Phylip files generated above so they can be analyzed by PAML
# They will thus be "reborn", which I have named as such because this portion was scripted while listening to Kids See Ghosts

cd gene_sequences

files=(`ls *.phy`)

for y in ${files[@]}; do

        number=`grep -c "^[A-Z]" ${y}`
        number_minus_one=`expr ${number} - 1`
        samples=(`awk 'NR>1 {print $1}' ${y}`)
        sequences=(`awk 'NR>1 {print $2}' ${y}`)
        series=(`seq 0 ${number_minus_one}`)
        head -n 1 ${y} > reborn_${y}

        for x in ${series[@]}; do

                echo ${samples[${x}]} >> reborn_${y}
                echo ${sequences[${x}]} >> reborn_${y}

done
done
