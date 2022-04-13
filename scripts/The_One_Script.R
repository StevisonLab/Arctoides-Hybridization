#candidate gene Analysis
library(ggplot2)
library(ggrepel)
library(interactions)
library(dplyr)

#read in all the gene lists:
#taylor's list
taylor=read.csv(file="Genital_Gene_Candidate_List.csv")

#HPO lists
HPO811=read.csv(file="updated_candidate_gene_list/genes_for_HP_0000811.csv")
HPO10461=read.csv(file="updated_candidate_gene_list/genes_for_HP_0010461.csv")
HPO10460=read.csv(file="updated_candidate_gene_list/genes_for_HP_0010460.csv")
HPO812=read.csv(file="updated_candidate_gene_list/genes_for_HP_0000812.csv")
HPO1827=read.csv(file="updated_candidate_gene_list/genes_for_HP_0001827.csv")
HPO12244=read.csv(file="updated_candidate_gene_list/genes_for_HP_0012244.csv")
HPO3252=read.csv(file="updated_candidate_gene_list/genes_for_HP_0003252.csv")

#repeat for mp lists
MP_0002210=read.table(file="updated_candidate_gene_list/MP_0002210_genes.txt")
MP_0003936=read.table(file="updated_candidate_gene_list/MP_0003936_genes.txt")
MP_0009198=read.table(file="updated_candidate_gene_list/MP_0009198_genes.txt")
MP_0009208=read.table(file="updated_candidate_gene_list/MP_0009208_genes.txt")

#compute data summary to add to plots
data_summary <- function (datum) {
  m = mean(datum)
  ymin = m-sd(datum)
  ymax = m+sd(datum)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

#read in pop gen data
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

#offline - run match_orthos.sh to get the rheMac8 EnsGeneID for each Human Gene Name

#intersect updated GC list with merged dataset
hpo_orhtos=read.table(file="updated_candidate_gene_list/hpo_plusTaylor_orthos.txt",header=F)
mpo_orthos=read.table(file="updated_candidate_gene_list/mpo_orthos.txt",header=F)
orthos=rbind(hpo_orhtos,mpo_orthos)

merged$NEW_GC_class=ifelse(merged$ensGeneID %in% orthos$V1,"GCgene","nonGCgene")

#get list of GC genes
gc_genes=merged$ensGeneID[(merged$NEW_GC_class=="GCgene")]
write.csv(merged,file="Combined_perGene_popGen_Stats_updated.csv",quote = F,row.names = F)


#add info on which candidate gene list each gene is on using ensembl file
ensembl=read.table(file="updated_candidate_gene_list/Ensembl_gene_list_rhemac8_human_mouse.tsv", sep="\t",header=T,na.strings = "")

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

#update Table 2
length(taylor$rheMac8[!is.na(taylor$rheMac8)])

length(HPO811$rheMac8[!is.na(HPO811$rheMac8)])
length(HPO10461$rheMac8[!is.na(HPO10461$rheMac8)])
length(HPO10460$rheMac8[!is.na(HPO10460$rheMac8)])
length(HPO812$rheMac8[!is.na(HPO812$rheMac8)])
length(HPO1827$rheMac8[!is.na(HPO1827$rheMac8)])
length(HPO12244$rheMac8[!is.na(HPO12244$rheMac8)])
length(HPO3252$rheMac8[!is.na(HPO3252$rheMac8)])

length(MP_0009198$rheMac8[!is.na(MP_0009198$rheMac8)])
length(MP_0009208$rheMac8[!is.na(MP_0009208$rheMac8)])
length(MP_0002210$rheMac8[!is.na(MP_0002210$rheMac8)])
length(MP_0003936$rheMac8[!is.na(MP_0003936$rheMac8)])



# dNdS outliers -----------------------------------------------------------
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
length(top_5_percent$ensGeneID)

top_5_percent$geneID.x[(top_5_percent$NEW_GC_class=="GCgene")]

dnds_outliers=top_5_percent$ensGeneID[(top_5_percent$NEW_GC_class=="GCgene")]

#write ranked list to a file             
temp=top_5_percent[,c("geneID.x","Arc-dNdS")]
#[1] PTPRC    ASPM     CD96     KIAA0586 DACT1    MAN2B2   HFM1     BTBD8    EXO1     PCNT     CSPP1    MYLK    
#[13] TOPAZ1   TLR2     ADAMTS3  SDCCAG8  CPLANE1  NPHP4    INPP5B   EYS      RP1L1 

mean(my_data$`Arc-dNdS`[my_data$NEW_GC_class=="GCgene"],na.rm=T)
#0.1596277
mean(my_data$`Arc-dNdS`[my_data$NEW_GC_class=="nonGCgene"],na.rm=T)
#0.1426137

#more values for Table 2 (background counts for dN/dS per list)
sum(my_data$ensGeneID %in% HPO811$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10461$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10460$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO812$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO1827$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO12244$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO3252$rheMac8,1,0)

sum(my_data$ensGeneID %in% MP_0009198$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0009208$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0002210$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0003936$rheMac8,1,0)
sum(my_data$ensGeneID %in% taylor$rheMac8)

# fdM outliers -----------------------------------------------------------
#repeat for fdm
#Note: when doing fdm, a cutoff of only 10 sites was used. Need to filter to get to cutoff of 100!
merged2=subset(merged,merged$sites>=100)

#clean up dataframe
my_data=merged2[c("ensGeneID","geneID.x","fdM","RecRate","NEW_GC_class")]
my_data<- my_data[complete.cases(my_data),]
length(my_data$ensGeneID)
length(my_data$ensGeneID[my_data$NEW_GC_class=="GCgene"])

#try to get outliers
mod=lm(my_data$fdM~my_data$RecRate)
summary(mod)

fitted <- mod$coefficients[2]*my_data$RecRate + mod$coefficients[1]
#cbind fitted to data
my_data <- cbind(my_data,fitted)

# calculate cooks d for all data
cooksd <- cooks.distance(mod)
#cbind cooksd to data
my_data <- cbind(my_data,cooksd)


#get upper outliers
nrow(my_data[(cooksd >=3*mean(cooksd, na.rm=T)) & # cooksD is high
               (my_data$fdM > my_data$fitted),]) # the value is an upper outlier

fdm_outliers=my_data$ensGeneID[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                                 (my_data$fdM > my_data$fitted)]

nrow(my_data[(my_data$NEW_GC_class=="nonGCgene"),])

#same as below, but neater
fdm_upper=intersect(fdm_outliers,gc_genes)

#outlier genes from GC list
my_data$geneID.x[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                   (my_data$fdM > my_data$fitted) &
                   (my_data$NEW_GC_class=="GCgene") ]

#get lower outliers
nrow(my_data[(cooksd >=3*mean(cooksd, na.rm=T)) & # cooksD is high
               (my_data$fdM < my_data$fitted),]) # the value is an upper outlier

fdm_outliers2=my_data$ensGeneID[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                                  (my_data$fdM < my_data$fitted)]

#same as below, but neater
fdm_lower=intersect(fdm_outliers2,gc_genes)

#outlier genes from GC list
my_data$geneID.x[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                   (my_data$fdM < my_data$fitted) &
                   (my_data$NEW_GC_class=="GCgene") ]

mean(my_data$fdM[my_data$NEW_GC_class=="GCgene"],na.rm=T)
mean(my_data$fdM[my_data$NEW_GC_class=="nonGCgene"],na.rm=T)


sum(my_data$ensGeneID %in% HPO811$rheMac8)
sum(my_data$ensGeneID %in% HPO10461$rheMac8)
sum(my_data$ensGeneID %in% HPO10460$rheMac8)
sum(my_data$ensGeneID %in% HPO812$rheMac8)
sum(my_data$ensGeneID %in% HPO1827$rheMac8)
sum(my_data$ensGeneID %in% HPO12244$rheMac8)
sum(my_data$ensGeneID %in% HPO3252$rheMac8)

sum(my_data$ensGeneID %in% MP_0009198$rheMac8)
sum(my_data$ensGeneID %in% MP_0009208$rheMac8)
sum(my_data$ensGeneID %in% MP_0002210$rheMac8)
sum(my_data$ensGeneID %in% MP_0003936$rheMac8)
sum(my_data$ensGeneID %in% taylor$rheMac8)

# pi outliers -----------------------------------------------------------
#repeat for pi
#clean up dataframe
my_data=merged[c("ensGeneID","geneID.x","RecRate","pi_Arctoides.x","NEW_GC_class")]
my_data<- my_data[complete.cases(my_data),]
length(my_data$ensGeneID)
length(my_data$ensGeneID[my_data$NEW_GC_class=="GCgene"])

#try to get outliers
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
nrow(my_data[(cooksd >=3*mean(cooksd, na.rm=T)) & # cooksD is high
               (my_data$pi_Arctoides < my_data$fitted),]) # the value is an upper outlier

pi_outliers=my_data$ensGeneID[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                                (my_data$pi_Arctoides < my_data$fitted)]

#same as below, but neater
pi_outliers2=intersect(pi_outliers,gc_genes)

#outlier genes from GC list
my_data$geneID.x[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                   (my_data$pi_Arctoides < my_data$fitted) &
                   (my_data$NEW_GC_class=="GCgene") ]

#CEP19  SLC9A3 GATA3  PTPN23 SNAI2  HMGA1

mean(my_data$pi_Arctoides[my_data$NEW_GC_class=="GCgene"],na.rm=T)
#0.0509636
mean(my_data$pi_Arctoides[my_data$NEW_GC_class=="nonGCgene"],na.rm=T)
#0.05472103

sum(my_data$ensGeneID %in% HPO811$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10461$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10460$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO812$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO1827$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO12244$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO3252$rheMac8,1,0)

sum(my_data$ensGeneID %in% MP_0009198$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0009208$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0002210$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0003936$rheMac8,1,0)


# FST outliers -----------------------------------------------------------
#fst outliers Arc-Sin
my_data=merged[c("ensGeneID","geneID.x","RecRate","Fst_Arctoides_Sinica","NEW_GC_class")]
my_data<- my_data[complete.cases(my_data),]
length(my_data$ensGeneID)
length(my_data$ensGeneID[my_data$NEW_GC_class=="GCgene"])

mod=lm(my_data$Fst_Arctoides_Sinica~my_data$RecRate)
summary(mod)

fitted <- mod$coefficients[2]*my_data$RecRate + mod$coefficients[1]
#cbind fitted to data
my_data <- cbind(my_data,fitted)

# calculate cooks d for all data
cooksd <- cooks.distance(mod)
#cbind cooksd to data
my_data <- cbind(my_data,cooksd)


#get outliers
nrow(my_data[(cooksd >=3*mean(cooksd, na.rm=T)) & # cooksD is high
               (my_data$Fst_Arctoides_Sinica > my_data$fitted),]) # the value is an upper outlier

fst_outliers_arc_sin=my_data$ensGeneID[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                                         (my_data$Fst_Arctoides_Sinica > my_data$fitted)]

#same as below, but neater
fst_arc_sin=intersect(fst_outliers_arc_sin,gc_genes)

#outlier genes from GC list
my_data$geneID.x[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                   (my_data$Fst_Arctoides_Sinica > my_data$fitted) &
                   (my_data$NEW_GC_class=="GCgene") ]

#[1] SLC35D1 MST1R   CEP19   TFAP2A  GATA3   HSPA4   ZCWPW1  PTPN23  CSPP1   CAMKV   UBE2D2  SKIDA1

mean(my_data$Fst_Arctoides_Sinica[my_data$NEW_GC_class=="GCgene"],na.rm=T)
#0.3647024
mean(my_data$Fst_Arctoides_Sinica[my_data$NEW_GC_class=="nonGCgene"],na.rm=T)
#0.3385975

sum(my_data$ensGeneID %in% HPO811$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10461$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10460$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO812$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO1827$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO12244$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO3252$rheMac8,1,0)

sum(my_data$ensGeneID %in% MP_0009198$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0009208$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0002210$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0003936$rheMac8,1,0)

#fst outliers Arc-Fas
my_data=merged[c("ensGeneID","geneID.x","RecRate","Fst_Arctoides_Fascicularis","NEW_GC_class")]
my_data<- my_data[complete.cases(my_data),]
length(my_data$ensGeneID)
length(my_data$ensGeneID[my_data$NEW_GC_class=="GCgene"])

mod=lm(my_data$Fst_Arctoides_Fascicularis~my_data$RecRate)
summary(mod)

fitted <- mod$coefficients[2]*my_data$RecRate + mod$coefficients[1]
#cbind fitted to data
my_data <- cbind(my_data,fitted)

# calculate cooks d for all data
cooksd <- cooks.distance(mod)
#cbind cooksd to data
my_data <- cbind(my_data,cooksd)


#get outliers
nrow(my_data[(cooksd >=3*mean(cooksd, na.rm=T)) & # cooksD is high
               (my_data$Fst_Arctoides_Fascicularis > my_data$fitted),]) # the value is an upper outlier

fst_outliers_arc_fas=my_data$ensGeneID[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                                         (my_data$Fst_Arctoides_Fascicularis > my_data$fitted)]

#same as below, but neater
fst_arc_fas=intersect(fst_outliers_arc_fas,gc_genes)

#outlier genes from GC list
my_data$geneID.x[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                   (my_data$Fst_Arctoides_Fascicularis > my_data$fitted) &
                   (my_data$NEW_GC_class=="GCgene") ]

#CEP19  AIRE   NR2F2  GATA3  TDRKH  PTPN23 SNAI2  HMGA1 

mean(my_data$Fst_Arctoides_Fascicularis[my_data$NEW_GC_class=="GCgene"],na.rm=T)
#0.4721669
mean(my_data$Fst_Arctoides_Fascicularis[my_data$NEW_GC_class=="nonGCgene"],na.rm=T)
#0.4521043

sum(my_data$ensGeneID %in% HPO811$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10461$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10460$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO812$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO1827$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO12244$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO3252$rheMac8,1,0)

sum(my_data$ensGeneID %in% MP_0009198$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0009208$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0002210$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0003936$rheMac8,1,0)


# dxy outliers -----------------------------------------------------------
#dxy outliers Arc-Sin
my_data=merged[c("ensGeneID","geneID.x","RecRate","dxy_Arctoides_Sinica","NEW_GC_class")]
my_data<- my_data[complete.cases(my_data),]
length(my_data$GeneID)
length(my_data$ensGeneID[my_data$NEW_GC_class=="GCgene"])

#dxy outliers
mod=lm(my_data$dxy_Arctoides_Sinica~my_data$RecRate)
summary(mod)

fitted <- mod$coefficients[2]*my_data$RecRate + mod$coefficients[1]
#cbind fitted to data
my_data <- cbind(my_data,fitted)

# calculate cooks d for all data
cooksd <- cooks.distance(mod)
#cbind cooksd to data
my_data <- cbind(my_data,cooksd)


#get outliers
nrow(my_data[(cooksd >=3*mean(cooksd, na.rm=T)) & # cooksD is high
               (my_data$dxy_Arctoides_Sinica > my_data$fitted),]) # the value is an upper outlier

dxy_outliers_arc_sin=my_data$ensGeneID[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                                         (my_data$dxy_Arctoides_Sinica > my_data$fitted)]

#same as below, but neater
dxy_arc_sin=intersect(dxy_outliers_arc_sin,gc_genes)

#outlier genes from GC list
my_data$geneID.x[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                   (my_data$dxy_Arctoides_Sinica > my_data$fitted) &
                   (my_data$NEW_GC_class=="GCgene") ]

#KDM3B   CYP21A2 DND1    HSPA4   DPY19L2 PSMB8   ZCWPW1  CSPP1   TACSTD2 HFE

mean(my_data$dxy_Arctoides_Sinica[my_data$NEW_GC_class=="GCgene"],na.rm=T)
#0.1479928
mean(my_data$dxy_Arctoides_Sinica[my_data$NEW_GC_class=="nonGCgene"],na.rm=T)
#0.1480222

sum(my_data$ensGeneID %in% HPO811$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10461$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10460$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO812$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO1827$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO12244$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO3252$rheMac8,1,0)

sum(my_data$ensGeneID %in% MP_0009198$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0009208$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0002210$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0003936$rheMac8,1,0)




#dxy outliers Arc-Fas
my_data=merged[c("ensGeneID","geneID.x","RecRate","dxy_Arctoides_Fascicularis","NEW_GC_class")]
my_data<- my_data[complete.cases(my_data),]
length(my_data$GeneID)
length(my_data$ensGeneID[my_data$NEW_GC_class=="GCgene"])

#dxy outliers 
mod=lm(my_data$dxy_Arctoides_Fascicularis~my_data$RecRate)
summary(mod)

fitted <- mod$coefficients[2]*my_data$RecRate + mod$coefficients[1]
#cbind fitted to data
my_data <- cbind(my_data,fitted)

# calculate cooks d for all data
cooksd <- cooks.distance(mod)
#cbind cooksd to data
my_data <- cbind(my_data,cooksd)


#get outliers
nrow(my_data[(cooksd >=3*mean(cooksd, na.rm=T)) & # cooksD is high
               (my_data$dxy_Arctoides_Fascicularis > my_data$fitted),]) # the value is an upper outlier

dxy_outliers_arc_fas=my_data$ensGeneID[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                                         (my_data$dxy_Arctoides_Fascicularis > my_data$fitted)]

#same as below, but neater
dxy_arc_fas=intersect(dxy_outliers_arc_fas,gc_genes)

#outlier genes from GC list
my_data$geneID.x[(cooksd >=3*mean(cooksd, na.rm=T)) & 
                   (my_data$dxy_Arctoides_Fascicularis > my_data$fitted) &
                   (my_data$NEW_GC_class=="GCgene") ]

#[1] STAG1    CEP19    AIRE     DND1     PIK3CA   TDRKH    SIN3A    TOGARAM1 PSMB8    PTPN23   TMED10   PRPF3   
#[13] CAMKV    CCR8     NODAL    MAB21L2

mean(my_data$dxy_Arctoides_Fascicularis[my_data$NEW_GC_class=="GCgene"],na.rm=T)
#0.1982903
mean(my_data$dxy_Arctoides_Fascicularis[my_data$NEW_GC_class=="nonGCgene"],na.rm=T)
#0.193532

sum(my_data$ensGeneID %in% HPO811$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10461$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO10460$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO812$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO1827$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO12244$rheMac8,1,0)
sum(my_data$ensGeneID %in% HPO3252$rheMac8,1,0)

sum(my_data$ensGeneID %in% MP_0009198$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0009208$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0002210$rheMac8,1,0)
sum(my_data$ensGeneID %in% MP_0003936$rheMac8,1,0)




#Summary Outliers

#combine GC outlier lists
outliers=list(dnds_outliers,fdm_upper,fdm_lower,pi_outliers2,fst_arc_sin,dxy_arc_sin,
              fst_arc_fas,dxy_arc_fas)

#make overlap table
sapply(outliers, function(x) sapply(outliers, function(y) sum(y %in% x)))

combined=unique(sort(c(as.character(dnds_outliers),fdm_upper,fdm_lower,pi_outliers2,fst_arc_sin,dxy_arc_sin,fst_arc_fas,dxy_arc_fas)))

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
outlier_table$mouse_outlier=ifelse(outlier_table$MP_0003936_list==1,1,
                                   ifelse(outlier_table$MP_0009208_list==1,1,
                                          ifelse(outlier_table$MP_0002210_list==1,1,
                                                 ifelse(outlier_table$MP_0009198_list==1,1,0))))
outlier_table$human_outlier=ifelse(outlier_table$HPO811_list==1,1,
                                   ifelse(outlier_table$HPO10461_list==1,1,
                                          ifelse(outlier_table$HPO10460_list==1,1,
                                                 ifelse(outlier_table$HPO812_list==1,1,
                                                        ifelse(outlier_table$HPO12244_list==1,1,
                                                               ifelse(outlier_table$HPO1827_list==1,1,
                                                                      ifelse(outlier_table$HPO3252_list==1,1,
                                                                             ifelse(outlier_table$Lit_CG_list==1,1,0))))))))

outlier_table$Mouse_gene_name[outlier_table$mouse_outlier == 0] <- NA
outlier_table$Human_gene_name[outlier_table$human_outlier == 0] <- NA

write.csv(outlier_table,file="Outliers_Combined.csv",row.names = F,quote = F)


#intersect with Phenotype Ontology results and SNPeff results
outliers=read.csv(file="Outliers_Combined.csv",header=T,stringsAsFactors = T)

phen_terms=read.csv(file="updated_candidate_gene_list/Phenotype_Ontology_Results.csv",header=T,stringsAsFactors = T)

outliers_merged=merge(outliers,phen_terms,by=c("ensGeneID","geneID.x","NEW_GC_class","mouse_outlier","human_outlier"))

#now to intersect SNPeff data
#SNPeff=read.csv(file="../SnpEff.transcripts.csv",header=T,stringsAsFactors = T)
SNPeff=read.table(file="../SnpEff/Outlier_Genes.lengths.stats.genes.txt",header=T,stringsAsFactors = T)

#need to remove duplicated gene rows with multiple transcripts
#will use largest number of amino acids
df=arrange(SNPeff,SNPeff$GeneID,-SNPeff$aa_length)
SNPeff_trimmed=df[!duplicated(df$GeneId),]

#need to get ensembl gene ID in outlier table
outliers_merged$GeneId=ifelse(outliers_merged$ensGeneID %in% ensembl$Gene.stable.ID.version,ensembl$Gene.stable.ID[match(outliers_merged$ensGeneID,ensembl$Gene.stable.ID.version)],NA)

#merge with outliers
outliers_final=merge(outliers_merged,SNPeff_trimmed,by.x="ensGeneID",by.y="GeneId",all=T)
write.csv(outliers_final,file="Outliers_Combined_Final.csv",row.names = F,quote = F)


#run script `SNPeff_plot.R`
