library(phylotools)

old_name=get.fasta.name("SAMPLE_2.coding.masked.fa")
new_name=read.table("GenesNTranscripts.txt")
Ref_Table=data.frame(old_name,new_name)

write.csv(Ref_Table, "SAMPLE_2_Ref_Table.csv")

rename.fasta(infile = "SAMPLE_2.coding.masked.fa", ref_table=Ref_Table, outfile = "SAMPLE_2.genes.masked.fa")
