# Candidate Gene Analysis

Below is code used to identify outliers and manage candidate gene lists. The R script is located in the scripts folder and is called `The_One_Script.R`

Candidate Gene lists are available as CSV files in the data directory:

- Genital_Gene_Candidate_List.csv (genes based on the literature)
- genes_for_HP_0000811.csv
- genes_for_HP_0010461.csv
- genes_for_HP_0010460.csv
- genes_for_HP_0000812.csv
- genes_for_HP_0001827.csv
- genes_for_HP_0012244.csv
- genes_for_HP_0003252.csv
- MP_0002210_genes.txt
- MP_0003936_genes.txt
- MP_0009198_genes.txt
- MP_0009208_genes.txt


Population Genetics Stats were read into R and intersected into a final dataset:

```R
div=read.csv(file = "combined.wSilenus.divergence.txt",header=T,stringsAsFactors = T,na.strings = "nan")
div2=read.csv(file = "combined.wSilenus.Arc_divergence.txt",header=T,stringsAsFactors = T,na.strings = "nan")
intro=read.csv(file="combined.wSilenus.genes.txt",header=T,stringsAsFactors = T,na.strings = "nan")
genes=read.table(file="Genital_Gene_Candidate_List.wRec.bed",header=F,stringsAsFactors = T,sep="\t", na.strings =".")
colnames(genes)=c("scaffold","start","end","ensGeneID","geneID","orientation","GC_class","RecRate")
dnds=read.table(file="../dNdS/dNdS.results.wSilenus.bed", header=F, stringsAsFactors = T,sep="\t")
colnames(dnds)=c("chr","start","end","ensGeneID","GeneID","Arc-dNdS","Sin-dNdS","Fas-dNdS")
#make column for sorting
div$conc=paste(div$scaffold,div$start,div$end,sep=".")
div2$conc=paste(div2$scaffold,div2$start,div2$end,sep=".")
intro$conc=paste(intro$scaffold,intro$start,intro$end,sep=".")
genes$conc=paste(genes$scaffold,genes$start,genes$end,sep=".")
#dnds$conc=paste(dnds$chr,dnds$start,dnds$end,sep=".")
colnames(dnds)=c("chr", "start","end","ensGeneID","geneID","Arc-dNdS","Sin-dNdS","Fas-dNdS")

#merge div and intro
stats1=merge(div,div2,by=c("conc","scaffold","start","end","mid","sites"),all=T)
stats2=merge(stats1,intro,by=c("conc","scaffold","start","end","mid","sites"),all=T)
stats3=merge(stats2,genes,by=c("conc","scaffold","start","end"),all=T)
merged=merge(stats3,dnds,by="ensGeneID",all=T)

keep=c("ensGeneID","scaffold","start.x","end.x","sites","pi_Arctoides.x","pi_Sinica","pi_Nemestrina",
       "pi_Fascicularis","dxy_Arctoides_Sinica","dxy_Arctoides_Nemestrina",   
       "dxy_Arctoides_Fascicularis","dxy_Sinica_Nemestrina","dxy_Sinica_Fascicularis",
       "dxy_Nemestrina_Fascicularis","Fst_Arctoides_Sinica","Fst_Arctoides_Nemestrina",
       "Fst_Arctoides_Fascicularis","Fst_Sinica_Nemestrina","Fst_Sinica_Fascicularis",
       "Fst_Nemestrina_Fascicularis","pi_Arctoides.y","pi_Macaque",
       "dxy_Arctoides_Macaque","Fst_Arctoides_Macaque","sitesUsed","ABBA","BABA",
       "D","fd","fdM","geneID.x","orientation","GC_class","RecRate","Arc-dNdS","Sin-dNdS","Fas-dNdS")

merged = merged[keep]
```

Next, a shell script was used to match the gene names from either human or mouse to the rheMac8 ensGeneID. The script is called `match_orthos.sh` and relies on a biomart output of all the rheMac8 geneIDs with the human/mouse ortholog gene name.

The outputs of the sets of orthologs were read into R and a new column was added to the merged dataframe to label genes based on whether they were in this list or not:

```R
#intersect updated GC list with merged dataset
hpo_orhtos=read.table(file="hpo_plusTaylor_orthos.txt",header=F)
mpo_orthos=read.table(file="mpo_orthos.txt",header=F)
orthos=rbind(hpo_orhtos,mpo_orthos)

merged$NEW_GC_class=ifelse(merged$ensGeneID %in% orthos$V1,"GCgene","nonGCgene")

#get list of GC genes
gc_genes=merged$ensGeneID[(merged$NEW_GC_class=="GCgene")]
write.csv(merged,file="Combined_perGene_popGen_Stats_updated.csv",quote = F,row.names = F)
```

Biomart output file read into R to add this info to gene lists to facilitate comparisons:

```R
#add info on which candidate gene list each gene is on using ensembl file
ensembl=read.table(file="Ensembl_gene_list_rhemac8_human_mouse.tsv", sep="\t",header=T,na.strings = "")

#intersect CG list with ensembl table
ensembl$CGgene=ifelse(ensembl$Gene.stable.ID.version %in% gc_genes,1,0)
my_data=subset(ensembl,ensembl$CGgene==1)

HPO811$rheMac8=ifelse(HPO811$GENE_SYMBOL %in% my_data$Human.gene.name,my_data$Gene.stable.ID.version[match(HPO811$GENE_SYMBOL,my_data$Human.gene.name)],NA)
HPO10461$rheMac8=ifelse(HPO10461$GENE_SYMBOL %in% my_data$Human.gene.name,my_data$Gene.stable.ID.version[match(HPO10461$GENE_SYMBOL,my_data$Human.gene.name)],NA)
HPO10460$rheMac8=ifelse(HPO10460$GENE_SYMBOL %in% my_data$Human.gene.name,my_data$Gene.stable.ID.version[match(HPO10460$GENE_SYMBOL,my_data$Human.gene.name)],NA)
HPO812$rheMac8=ifelse(HPO812$GENE_SYMBOL %in% my_data$Human.gene.name,my_data$Gene.stable.ID.version[match(HPO812$GENE_SYMBOL,my_data$Human.gene.name)],NA)
HPO1827$rheMac8=ifelse(HPO1827$GENE_SYMBOL %in% my_data$Human.gene.name,my_data$Gene.stable.ID.version[match(HPO1827$GENE_SYMBOL,my_data$Human.gene.name)],NA)
HPO12244$rheMac8=ifelse(HPO12244$GENE_SYMBOL %in% my_data$Human.gene.name,my_data$Gene.stable.ID.version[match(HPO12244$GENE_SYMBOL,my_data$Human.gene.name)],NA)
HPO3252$rheMac8=ifelse(HPO3252$GENE_SYMBOL %in% my_data$Human.gene.name,my_data$Gene.stable.ID.version[match(HPO3252$GENE_SYMBOL,my_data$Human.gene.name)],NA)
MP_0002210$rheMac8=ifelse(MP_0002210$V1 %in% my_data$Mouse.gene.name,my_data$Gene.stable.ID.version[match(MP_0002210$V1,my_data$Mouse.gene.name)],NA)
MP_0003936$rheMac8=ifelse(MP_0003936$V1 %in% my_data$Mouse.gene.name,my_data$Gene.stable.ID.version[match(MP_0003936$V1,my_data$Mouse.gene.name)],NA)
MP_0009198$rheMac8=ifelse(MP_0009198$V1 %in% my_data$Mouse.gene.name,my_data$Gene.stable.ID.version[match(MP_0009198$V1,my_data$Mouse.gene.name)],NA)
MP_0009208$rheMac8=ifelse(MP_0009208$V1 %in% my_data$Mouse.gene.name,my_data$Gene.stable.ID.version[match(MP_0009208$V1,my_data$Mouse.gene.name)],NA)
taylor$rheMac8=ifelse(taylor$Gene.Name %in% my_data$Human.gene.name,my_data$Gene.stable.ID.version[match(taylor$Gene.Name,my_data$Human.gene.name)],NA)
```

dN/dS outliers:

```R
#clean up dataframe
my_data=merged[c("ensGeneID","geneID.x","NEW_GC_class","Arc-dNdS","RecRate")]
my_data<- my_data[complete.cases(my_data),]
length(my_data$ensGeneID)
length(my_data$ensGeneID[my_data$NEW_GC_class=="GCgene"])

#dnds top 5% 
length=length(my_data$`Arc-dNdS`)
sorted=sort(my_data$`Arc-dNdS`)
cutoff=sorted[length-round(length*0.05,1)]
top_5_percent=subset(my_data,my_data$`Arc-dNdS`>=cutoff)
dnds_outliers=top_5_percent$ensGeneID[(top_5_percent$NEW_GC_class=="GCgene")]
```

all other outliers (pi used as example):

```R
#repeat for pi
#clean up dataframe
my_data=merged[c("ensGeneID","geneID.x","RecRate","pi_Arctoides.x","NEW_GC_class")]
my_data<- my_data[complete.cases(my_data),]

#get outliers
mod=lm(my_data$pi_Arctoides~my_data$RecRate)
summary(mod)

fitted <- mod$coefficients[2]*my_data$RecRate + mod$coefficients[1]
#cbind fitted to data
my_data <- cbind(my_data,fitted)

# calculate cooks d for all data
cooksd <- cooks.distance(mod)
#cbind cooksd to data
my_data <- cbind(my_data,cooksd)

#get outliers
pi_outliers=my_data$ensGeneID[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                                (my_data$pi_Arctoides < my_data$fitted)]

#intersect with gc_genes
pi_outliers2=intersect(pi_outliers,gc_genes)
```

Each set of outliers were intersected with the list of GC genes and combined at the end:

```R
#combine GC outlier lists
outliers=list(dnds_outliers,fdm_upper,fdm_lower,pi_outliers2,fst_arc_sin,dxy_arc_sin,
              fst_arc_fas,dxy_arc_fas)
outlier_table=merged[merged$ensGeneID %in% combined,c("ensGeneID","geneID.x","pi_Arctoides.x","fdM","RecRate","Arc-dNdS","Fst_Arctoides_Fascicularis","dxy_Arctoides_Fascicularis","Fst_Arctoides_Sinica","dxy_Arctoides_Sinica","NEW_GC_class")]

#add outlier columns
outlier_table$pi_outlier=ifelse(outlier_table$ensGeneID %in% pi_outliers2,1,0)
outlier_table$dnds_outlier=ifelse(outlier_table$ensGeneID %in% dnds_outliers,1,0)
outlier_table$fdm_upper_outlier=ifelse(outlier_table$ensGeneID %in% fdm_upper,1,0)
outlier_table$fdm_lower_outlier=ifelse(outlier_table$ensGeneID %in% fdm_lower,1,0)
outlier_table$fst_arc_sin_outlier=ifelse(outlier_table$ensGeneID %in% fst_arc_sin,1,0)
outlier_table$fst_arc_fas_outlier=ifelse(outlier_table$ensGeneID %in% fst_arc_fas,1,0)
outlier_table$dxy_arc_sin_outlier=ifelse(outlier_table$ensGeneID %in% dxy_arc_sin,1,0)
outlier_table$dxy_arc_fas_outlier=ifelse(outlier_table$ensGeneID %in% dxy_arc_fas,1,0)

#add columns to outlier table to specify which list each is on
outlier_table$HPO811_list=ifelse(outlier_table$ensGeneID %in% HPO811$rheMac8,1,0)
outlier_table$HPO10461_list=ifelse(outlier_table$ensGeneID %in% HPO10461$rheMac8,1,0)
outlier_table$HPO10460_list=ifelse(outlier_table$ensGeneID %in% HPO10460$rheMac8,1,0)
outlier_table$HPO812_list=ifelse(outlier_table$ensGeneID %in% HPO812$rheMac8,1,0)
outlier_table$HPO1827_list=ifelse(outlier_table$ensGeneID %in% HPO1827$rheMac8,1,0)
outlier_table$HPO12244_list=ifelse(outlier_table$ensGeneID %in% HPO12244$rheMac8,1,0)
outlier_table$HPO3252_list=ifelse(outlier_table$ensGeneID %in% HPO3252$rheMac8,1,0)
outlier_table$Lit_CG_list=ifelse(outlier_table$ensGeneID %in% taylor$rheMac8,1,0)
outlier_table$MP_0002210_list=ifelse(outlier_table$ensGeneID %in% MP_0002210$rheMac8,1,0)
outlier_table$MP_0003936_list=ifelse(outlier_table$ensGeneID %in% MP_0003936$rheMac8,1,0)
outlier_table$MP_0009198_list=ifelse(outlier_table$ensGeneID %in% MP_0009198$rheMac8,1,0)
outlier_table$MP_0009208_list=ifelse(outlier_table$ensGeneID %in% MP_0009208$rheMac8,1,0)

#add human/mouse gene name to the table
outlier_table$Mouse_gene_name=ifelse(outlier_table$geneID.x %in% my_data$Gene.name,my_data$Mouse.gene.name[match(outlier_table$geneID.x,my_data$Gene.name)],NA)
outlier_table$Human_gene_name=ifelse(outlier_table$geneID.x %in% my_data$Gene.name,my_data$Human.gene.name[match(outlier_table$geneID.x,my_data$Gene.name)],NA)
outlier_table$mouse_outlier=ifelse(outlier_table$MP_0003936_list==1,1,ifelse(outlier_table$MP_0009208_list==1,1,ifelse(outlier_table$MP_0002210_list==1,1,ifelse(outlier_table$MP_0009198_list==1,1,0))))
outlier_table$human_outlier=ifelse(outlier_table$HPO811_list==1,1,ifelse(outlier_table$HPO10461_list==1,1,ifelse(outlier_table$HPO10460_list==1,1,ifelse(outlier_table$HPO812_list==1,1,ifelse(outlier_table$HPO12244_list==1,1,ifelse(outlier_table$HPO1827_list==1,1,ifelse(outlier_table$HPO3252_list==1,1,ifelse(outlier_table$Lit_CG_list==1,1,0))))))))

outlier_table$Mouse_gene_name[outlier_table$mouse_outlier == 0] <- NA
outlier_table$Human_gene_name[outlier_table$human_outlier == 0] <- NA

write.csv(outlier_table,file="Outliers_Combined.csv",row.names = F,quote = F)

#intersect with Phenotype Ontology results and SNPeff results
outliers=read.csv(file="updated_candidate_gene_list/Outliers_Combined.csv",header=T,stringsAsFactors = T)
phen_terms=read.csv(file="updated_candidate_gene_list/Phenotype_Ontology_Results.csv",header=T,stringsAsFactors = T)
outliers_merged=merge(outliers,phen_terms,by=c("ensGeneID","geneID.x","NEW_GC_class","mouse_outlier","human_outlier"))

#now to intersect SNPeff data
SNPeff=read.csv(file="../SnpEff.transcripts.csv",header=T,stringsAsFactors = T)

#need to remove duplicated gene rows with multiple transcripts
#will use largest number of amino acids
df=arrange(SNPeff,SNPeff$GeneID,-SNPeff$aa_length)
SNPeff_trimmed=df[!duplicated(df$GeneId),]

#need to get ensembl gene ID in outlier table
outliers_merged$GeneId=ifelse(outliers_merged$ensGeneID %in% my_data$Gene.stable.ID.version,my_data$Gene.stable.ID[match(outliers_merged$ensGeneID,my_data$Gene.stable.ID.version)],NA)

#merge with outliers
outliers_final=merge(outliers_merged,SNPeff_trimmed,by="GeneId",all=T)
write.csv(outliers_final,file="Outliers_Combined_Final.csv",row.names = F,quote = F)
```

Script `SNPeff_plot.R` was run to create Figure 5.
