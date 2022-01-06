#SNPeff plot

SNPeff_outliers=subset(outliers_final,outliers_final$Lit_CG_list==1)
SNPeff_outliers=SNPeff_outliers[,c(3,44:66)]

SNPeff_outliers$Up=100*(SNPeff_outliers$variants_effect_upstream_gene_variant/SNPeff_outliers$nt_length)

SNPeff_outliers$Down=100*(SNPeff_outliers$variants_effect_downstream_gene_variant/SNPeff_outliers$nt_length)

SNPeff_outliers$Intron=100*(SNPeff_outliers$variants_effect_intron_variant/SNPeff_outliers$nt_length)

SNPeff_outliers$Splice=100*(((SNPeff_outliers$variants_effect_splice_acceptor_variant+SNPeff_outliers$variants_effect_splice_donor_variant+SNPeff_outliers$variants_effect_splice_region_variant))/SNPeff_outliers$nt_length)

SNPeff_outliers$ThreeUTR=100*(SNPeff_outliers$variants_effect_3_prime_UTR_variant/SNPeff_outliers$nt_length)

SNPeff_outliers$FiveUTR=100*(SNPeff_outliers$variants_effect_5_prime_UTR_variant/SNPeff_outliers$nt_length)

#impact
SNPeff_outliers$HIGH=100*(SNPeff_outliers$variants_impact_HIGH/SNPeff_outliers$nt_length)
SNPeff_outliers$MODERATE=100*(SNPeff_outliers$variants_impact_MODERATE/SNPeff_outliers$nt_length)
SNPeff_outliers$LOW=100*(SNPeff_outliers$variants_impact_LOW/SNPeff_outliers$nt_length)
SNPeff_outliers$MODIFIER=100*(SNPeff_outliers$variants_impact_MODIFIER/SNPeff_outliers$nt_length)

#subset
seo=SNPeff_outliers[,c("geneID.x","HIGH","MODERATE","LOW","MODIFIER")]
barplot(as.matrix(seo[,c(2:5)]),xlab="Percent of effects by impact",ylab="Percent",main="SNPeff Results",beside=TRUE,col=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628"))

legend("topleft",legend=seo$geneID.x,cex=0.3,col=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628"),fill=1:6)


#functional class
SNPeff_outliers$Missense=100*(SNPeff_outliers$variants_effect_missense_variant/SNPeff_outliers$aa_length)

SNPeff_outliers$Silent=100*(SNPeff_outliers$variants_effect_synonymous_variant/SNPeff_outliers$aa_length)

SNPeff_outliers$Nonsense=100*((SNPeff_outliers$variants_effect_start_lost+SNPeff_outliers$variants_effect_stop_gained+SNPeff_outliers$variants_effect_stop_retained_variant)/SNPeff_outliers$aa_length)


seo2=SNPeff_outliers[,c("geneID.x","Missense","Nonsense","Silent")]
gene_cols=c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628")

png(filename = "SNPeff_results.png",width=6.5,units="in",height=4,res=300)
par(mar = c(4.1, 4.1, 0.2, 4),                                  # Specify par parameters
    xpd = TRUE)
barplot(as.matrix(seo2[,c(2:4)]),xlab="SNP Functional class",ylab="Percent",beside=TRUE,col=gene_cols)

legend("topright",legend=seo2$geneID.x,cex=0.6,col=gene_cols,fill=gene_cols,bty="n",inset = c(- 0.15, 0))
dev.off()
