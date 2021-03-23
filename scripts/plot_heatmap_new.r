
combined=read.table(file="combined.50kb.txt",header=T,na.strings = "nan",stringsAsFactors = T,sep=",")

chr.labels=c(1:20, "X")

combined[combined$scaffold=="chr1","CHROM"]="1"
combined[combined$scaffold=="chr2","CHROM"]="2"
combined[combined$scaffold=="chr3","CHROM"]="3"
combined[combined$scaffold=="chr4","CHROM"]="4"
combined[combined$scaffold=="chr5","CHROM"]="5"
combined[combined$scaffold=="chr6","CHROM"]="6"
combined[combined$scaffold=="chr7","CHROM"]="7"
combined[combined$scaffold=="chr8","CHROM"]="8"
combined[combined$scaffold=="chr9","CHROM"]="9"
combined[combined$scaffold=="chr10","CHROM"]="10"
combined[combined$scaffold=="chr11","CHROM"]="11"
combined[combined$scaffold=="chr12","CHROM"]="12"
combined[combined$scaffold=="chr13","CHROM"]="13"
combined[combined$scaffold=="chr14","CHROM"]="14"
combined[combined$scaffold=="chr15","CHROM"]="15"
combined[combined$scaffold=="chr16","CHROM"]="16"
combined[combined$scaffold=="chr17","CHROM"]="17"
combined[combined$scaffold=="chr18","CHROM"]="18"
combined[combined$scaffold=="chr19","CHROM"]="19"
combined[combined$scaffold=="chr20","CHROM"]="20"
combined[combined$scaffold=="chrX","CHROM"]="21"
#combined[combined$scaffold=="chrM","CHROM"]="22"

combined$CHROM=as.numeric(combined$CHROM)

head(combined)

combined3=subset(combined,!is.na(combined$D))
combined4=subset(combined,!is.na(combined$fdM))
#combined4=subset(combined,!is.na(combined$fd))

#fd outliers
#top_fdM=subset(combined4,combined4$fdM<=-0.8 | combined4$fdM>=0.45)

#hist(combined2$fdM.2,xlab="fdM",main="Distribution of fdM values")
#top_fdM_sinica=subset(combined4,combined4$fdM.2<=-0.8)
#write.csv(top_fdM_sinica,file="top_fdM_regions_sinica.50kb.csv")


png(file="hist.50kb.png", height=5.5,width=7.5,res=300,units="in",pointsize=12)
par(mfrow=c(1,2))
hist(combined4$fdM,xlab="fdM",main="Distribution of fdM values")
hist(combined3$D,xlab="D",main="Distribution of D values")
dev.off()

#top_fdM_fascicularis=subset(combined4,combined4$fdM.2>=0.45)
#write.csv(top_fdM_fascicularis,file="top_fdM_regions_fascicularis.50kb.csv")

new=data.frame(combined4$CHROM,combined4$start,combined4$fdM)
library(reshape2)

x=dcast(new,combined4.start~combined4.CHROM)
mat=as.matrix(x[,-1])

#breaks_neg=seq(min(combined4$fdM)-0.01,0, length.out=11)
#breaks_pos=seq(0,max(combined4$fdM)+0.01,length.out=11)
#breaks=c(breaks_neg,breaks_pos[-1])
breaks_neg=seq(-1,0,length.out=11)
breaks_pos=seq(0,1,length.out=11)
breaks=c(breaks_neg,breaks_pos[-1])
breaks

#cols2=c("navy","mediumblue","royalblue1","steelblue1","lightgoldenrod","orange","orangered","red")
cols_red=rainbow(10,start=0, end=.167)
cols_blue=rainbow(10, start=0.51, end=.66)
cols2=c(cols_red,cols_blue)

ranges=vector(mode="character",length=length(breaks)-1)
for (i in 1:length(breaks)-1) { 
  ranges[i]=paste(round(breaks[i],2),",",round(breaks[i+1],2))
}

ranges

png(file="heatmap_mosaic.fdM.50kb.png", height=5.5,width=10,res=300,units="in",pointsize=12)
#par(mar=c(0,0,7,0))
#par(mfrow=c(2,1))
heatmap(mat,Rowv=NA,Colv=NA,na.rm=T,keep.dendro = FALSE,scale="none",
        col = cols2, breaks=breaks,
        labCol=chr.labels,labRow = FALSE, 
        add.expr = legend("topright",legend=rev(ranges),fill=rev(cols2),cex=0.65,bty="n",title="Sinica (-) vs. Fascicularis (+) ancestry"))
dev.off()


new2=data.frame(combined3$CHROM,combined3$start,combined3$D)
#library(reshape2)

y=dcast(new2,combined3.start~combined3.CHROM)
mat2=as.matrix(y[,-1])


png(file="heatmap_mosaic.D.50kb.png", height=5.5,width=7.5,res=300,units="in",pointsize=12)
#par(mar=c(0,0,7,0))
heatmap(mat2,Rowv=NA,Colv=NA,na.rm=T,keep.dendro = FALSE,scale="none",
        col = cols2, breaks=breaks,
        labCol=chr.labels,labRow = FALSE, 
        add.expr = legend("topright",legend=rev(ranges),fill=rev(cols2),cex=0.65,bty="n",title="Sinica (-) vs. Fascicularis (+) ancestry"))
dev.off()