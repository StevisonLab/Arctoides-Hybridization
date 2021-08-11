library(phylotools)

old_name=get.fasta.name("SAMPLE_1.coding.masked.fa")
new_name=read.table("GenesNTranscripts.txt")
Ref_Table=data.frame(old_name,new_name)

write.csv(Ref_Table, "SAMPLE_1_Ref_Table.csv")

rename.fasta(infile = "SAMPLE_1.coding.masked.fa", ref_table=Ref_Table, outfile = "SAMPLE_1.genes.masked.fa")
