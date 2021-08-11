#!/bin/bash

# Define array of gene names based on unique ENSMMUG ID's, not HUGO genes names

Genes=(`awk '{print $4}' ../per_gene_stats/Gene.Topo.wSilenus.bed | sort -u | tail -n 5361`)

# Loop through these genes and define other parameters, including coordinates, transcripts, HUGO gene names, and topology proportions

for x in ${Genes[@]}
do

CHR=`grep "\s${x}\s" ../per_gene_stats/Gene.Topo.wSilenus.bed  | awk '{print $1}' | sort -u`
START=`grep "\s${x}\s" ../per_gene_stats/Gene.Topo.wSilenus.bed  | awk '{print $2}' | sort -u`
END=`grep "\s${x}\s" ../per_gene_stats/Gene.Topo.wSilenus.bed  | awk '{print $3}' | sort -u`
proportions=(`grep "\s${x}\s" ../per_gene_stats/Gene.Topo.wSilenus.bed | awk '{sumone+=$7; sumtwo+=$8; sumthree+=$9; sum=sumone+sumtwo+sumthree} END {OFS="\t" ; print $4,sumone/sum,sumtwo/sum,sumthree/sum}'`)
ENSMUT=(`awk -v TRANS=${x} '$2==TRANS {print $1}' ../References_files/ensmutToensmug.txt`)
GENE=`grep "\s${x}\s" ../per_gene_stats/Gene.Topo.wSilenus.bed  | awk '{print $5}' | sort -u`

# Loop through the transcripts and find their results in the correct respective result file

for a in ${ENSMUT[@]}
do

for z in "results" "results_2" "results_3" "results_4"
do
	output=(`ls ${z}/*_${a}/*mlc`)

# Loop through every mlc output file for each transcript

	for y in ${output[@]}
	do

# Define several variables. Namely first define the gene-transcript name by extracting from output header file
# Then extract the dN/dS and dS values for each topology of each branch
# These appear convoluted because paml output is convoluted. I verified these as the correct topologies and branches using FigTree
# That is to say you can compare the labels on your tree input file to labels provided by a PAML rst file
# Where the rst file contains a line saying "tree with node labels for Rod Page's TreeView"
# Page 16 of the PAML/4.9j manual refers to this and FigTree
# There are multiple posts on the PAML google discussion group elaborating on this, e.g. https://groups.google.com/g/pamlsoftware/c/pKkr0pAIPfs

		Transcript=`head -n 1 ${y} | awk '{print $8}' | sed "s/\/home\/npb0015\/arctoides_project\/gene_sequences\/reborn_//" | sed "s/.phy//"`
		Topo3_Arc_omega=`grep -A 60 "TREE #  1" ${y} | grep -m 1 "^\s\s29..30" | awk '{print $5}'`
		Topo2_Arc_omega=`grep -A 60 "TREE #  2" ${y} | grep -m 1 "^\s\s29..30" | awk '{print $5}'`
		Topo1_Arc_omega=`grep -A 60 "TREE #  3" ${y} | grep -m 1 "^\s\s28..39" | awk '{print $5}'`
		Topo3_Arc_dS=`grep -A 60 "TREE #  1" ${y} | grep -m 1 "^\s\s29..30" | awk '{print $7}'`
                Topo2_Arc_dS=`grep -A 60 "TREE #  2" ${y} | grep -m 1 "^\s\s29..30" | awk '{print $7}'`
                Topo1_Arc_dS=`grep -A 60 "TREE #  3" ${y} | grep -m 1 "^\s\s28..39" | awk '{print $7}'`
		Topo3_Sin_omega=`grep -A 60 "TREE #  1" ${y} | grep -m 1 "^\s\s28..39" | awk '{print $5}'`
                Topo2_Sin_omega=`grep -A 60 "TREE #  2" ${y} | grep -m 1 "^\s\s29..34" | awk '{print $5}'`
                Topo1_Sin_omega=`grep -A 60 "TREE #  3" ${y} | grep -m 1 "^\s\s29..35" | awk '{print $5}'`
                Topo3_Sin_dS=`grep -A 60 "TREE #  1" ${y} | grep -m 1 "^\s\s28..39" | awk '{print $7}'`
                Topo2_Sin_dS=`grep -A 60 "TREE #  2" ${y} | grep -m 1 "^\s\s29..34" | awk '{print $7}'`
                Topo1_Sin_dS=`grep -A 60 "TREE #  3" ${y} | grep -m 1 "^\s\s29..35" | awk '{print $7}'`
		Topo3_Fas_omega=`grep -A 60 "TREE #  1" ${y} | grep -m 1 "^\s\s29..34" | awk '{print $5}'`
                Topo2_Fas_omega=`grep -A 60 "TREE #  2" ${y} | grep -m 1 "^\s\s28..38" | awk '{print $5}'`
                Topo1_Fas_omega=`grep -A 60 "TREE #  3" ${y} | grep -m 1 "^\s\s29..30" | awk '{print $5}'`
                Topo3_Fas_dS=`grep -A 60 "TREE #  1" ${y} | grep -m 1 "^\s\s29..34" | awk '{print $7}'`
                Topo2_Fas_dS=`grep -A 60 "TREE #  2" ${y} | grep -m 1 "^\s\s28..38" | awk '{print $7}'`
                Topo1_Fas_dS=`grep -A 60 "TREE #  3" ${y} | grep -m 1 "^\s\s29..30" | awk '{print $7}'`

# Take weighted average of dN/dS values per topology based on proportion of topology value provided from Twisst
# Variable will be empty for later exclusion if dS is 0 or dN/dS is greater than 10

		Topo3_Arc_avg=`awk -v DS=${Topo3_Arc_dS} -v Topo3=${Topo3_Arc_omega} -v Proportion=${proportions[3]} 'BEGIN {if (DS > 0.0000 && Topo3 < 10) print Topo3*Proportion}'`
		Topo2_Arc_avg=`awk -v DS=${Topo2_Arc_dS} -v Topo2=${Topo2_Arc_omega} -v Proportion=${proportions[2]} 'BEGIN {if (DS > 0.0000 && Topo2 < 10) print Topo2*Proportion}'`
		Topo1_Arc_avg=`awk -v DS=${Topo1_Arc_dS} -v Topo1=${Topo1_Arc_omega} -v Proportion=${proportions[1]} 'BEGIN {if (DS > 0.0000 && Topo1 < 10) print Topo1*Proportion}'`
		Topo3_Sin_avg=`awk -v DS=${Topo3_Sin_dS} -v Topo3=${Topo3_Sin_omega} -v Proportion=${proportions[3]} 'BEGIN {if (DS > 0.0000 && Topo3 < 10) print Topo3*Proportion}'`
                Topo2_Sin_avg=`awk -v DS=${Topo2_Sin_dS} -v Topo2=${Topo2_Sin_omega} -v Proportion=${proportions[2]} 'BEGIN {if (DS > 0.0000 && Topo2 < 10) print Topo2*Proportion}'`
                Topo1_Sin_avg=`awk -v DS=${Topo1_Sin_dS} -v Topo1=${Topo1_Sin_omega} -v Proportion=${proportions[1]} 'BEGIN {if (DS > 0.0000 && Topo1 < 10) print Topo1*Proportion}'`
		Topo3_Fas_avg=`awk -v DS=${Topo3_Fas_dS} -v Topo3=${Topo3_Fas_omega} -v Proportion=${proportions[3]} 'BEGIN {if (DS > 0.0000 && Topo3 < 10) print Topo3*Proportion}'`
                Topo2_Fas_avg=`awk -v DS=${Topo2_Fas_dS} -v Topo2=${Topo2_Fas_omega} -v Proportion=${proportions[2]} 'BEGIN {if (DS > 0.0000 && Topo2 < 10) print Topo2*Proportion}'`
                Topo1_Fas_avg=`awk -v DS=${Topo1_Fas_dS} -v Topo1=${Topo1_Fas_omega} -v Proportion=${proportions[1]} 'BEGIN {if (DS > 0.0000 && Topo1 < 10) print Topo1*Proportion}'`

# Combine dN/dS values for all topology weights per branch but only if all topology weights per branch fulfilled dS and dN/dS criteria defined above
# If not, they will output as NA

		Transcript_Arc_omega=`awk -v Topo3_Avg=${Topo3_Arc_avg} -v Topo2_Avg=${Topo2_Arc_avg} -v Topo1_Avg=${Topo1_Arc_avg} 'BEGIN {if (Topo3_Avg > 0 && Topo2_Avg > 0 && Topo1_Avg > 0) print Topo3_Avg+Topo2_Avg+Topo1_Avg}'`
		Transcript_Sin_omega=`awk -v Topo3_Avg=${Topo3_Sin_avg} -v Topo2_Avg=${Topo2_Sin_avg} -v Topo1_Avg=${Topo1_Sin_avg} 'BEGIN {if (Topo3_Avg > 0 && Topo2_Avg > 0 && Topo1_Avg > 0) print Topo3_Avg+Topo2_Avg+Topo1_Avg}'`
		Transcript_Fas_omega=`awk -v Topo3_Avg=${Topo3_Fas_avg} -v Topo2_Avg=${Topo2_Fas_avg} -v Topo1_Avg=${Topo1_Fas_avg} 'BEGIN {if (Topo3_Avg > 0 && Topo2_Avg > 0 && Topo1_Avg > 0) print Topo3_Avg+Topo2_Avg+Topo1_Avg}'`

		if [ -z ${Transcript_Arc_omega} ]
		then
		Transcript_Arc_omega="NA"
		fi
		if [ -z ${Transcript_Sin_omega} ]
                then
                Transcript_Sin_omega="NA"
                fi
		if [ -z ${Transcript_Fas_omega} ]
                then
                Transcript_Fas_omega="NA"
                fi

		echo -e "${CHR}\t${START}\t${END}\t${x}\t${Transcript}\t${Transcript_Arc_omega}\t${Transcript_Sin_omega}\t${Transcript_Fas_omega}" >> dNdS_results.wSilenus.txt

done
done
done

# Sort out unique transcripts (to account for the fact that a single transcript could be present in multiple results files above)
# and then average uniquely for each gene excluding NA

ArcAvg=`grep "\s${x}\s" dNdS_results.wSilenus.txt | sort -k5,5 -u  | awk '{print $6}' | grep -v "NA" | awk '{Arc+=$1} END {print Arc/NR}'`
SinAvg=`grep "\s${x}\s" dNdS_results.wSilenus.txt | sort -k5,5 -u  | awk '{print $7}' | grep -v "NA" | awk '{Sin+=$1} END {print Sin/NR}'`
FasAvg=`grep "\s${x}\s" dNdS_results.wSilenus.txt | sort -k5,5 -u  | awk '{print $8}' | grep -v "NA" | awk '{Fas+=$1} END {print Fas/NR}'`

# If there are no values to average for any of these branches they should be excluded as NA

if [ -z ${ArcAvg} ]
then
ArcAvg="NA"
fi
if [ -z ${SinAvg} ]
then
SinAvg="NA"
fi
if [ -z ${FasAvg} ]
then
FasAvg="NA"
fi

# Finally output unique gene dN/dS values in bed format and then remove the temporary file generated above

grep "\s${x}\s" dNdS_results.wSilenus.txt | sort -k4,4 -u  | awk -v gene=${x} -v GENE=${GENE} -v Arc=${ArcAvg} -v Sin=${SinAvg} -v Fas=${FasAvg} '{OFS="\t" ; print $1,$2,$3,gene,GENE,Arc,Sin,Fas}' | sort -u >> dNdS.results.wSilenus.bed
done
rm dNdS_results.wSilenus.txt
