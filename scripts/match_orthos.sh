#! /bin/sh

rm mouse_list.orthologs.txt

#hpo="HPO_Taylor_merged_genes.csv"
mpo="MP_combined_genes.txt"
orth="Ensembl_gene_list_rhemac8_human_mouse.tsv"
#orth="human_rhesus_orthologs.clean.txt"

#echo $hpo, $orth

genes=(`awk -F, 'NR>1 {print $1}' $mpo`)

#echo "Mouse_gene\tRhesus_gene" >mouse_list.orthologs.txt

for i in ${genes[@]}
do
#i="ADAR"
#echo $i
match=(`awk '{print $4}' $orth | grep -n "^$i$" | awk -F: '{print $1}'`)
#echo ${match[0]}
rhe_gene=`awk -v ln=${match[0]} 'NR==ln {print $2}' $orth`

echo "$rhe_gene" >> mouse_list.orthologs.txt

done 

awk '$0!~/^$/' mouse_list.orthologs.txt >mpo_orthos.txt
