# Genome Analysis Pipeline 

This file details the steps in genome analysis of raw reads. 

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
```

5. Merge alignments at the sample level into single bam file and sort/index.
