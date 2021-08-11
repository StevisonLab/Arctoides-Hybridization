#!/bin/bash

#PBS -l nodes=4:ppn=20,walltime=24:00:00:00,mem=100gb
#PBS -W x=FLAGS:ADVRES:lss0021_lab.197924
#PBS -q general
#PBS -m abe
#PBS -d .

# Load necessary modules

module load paml/4.9i
module load gnu-parallel/20181022

# Define phylip files for analysis and then enter results folder

cd ../gene_sequences
phylips=(`ls reborn*`)

cd ../paml_results/results

# Analyze each transcript in parallel and output the results into unique folder for that transcript

parallel '
prefix=`echo {} | sed "s/.phy//" | sed "s/reborn_//"` 
mkdir ${prefix}
sed "s/PHYLIP_FILE/{}/" ../codeml.ctl | sed "s/OUTPUT_FILE/${prefix}/" > ${prefix}/${prefix}.ctl
cd ${prefix}
codeml ${prefix}.ctl
cd ..' ::: ${phylips[@]}
