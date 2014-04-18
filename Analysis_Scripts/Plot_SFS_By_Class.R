#   An R script to plot a site frequency spectrum for various partitions of
#   SNPs. 
#   The data files are in the format
#   SNPNAME MAF
#   ...
#   separated by tabs
#   Fill in the correct values for your files, and your desired output file

#   Read the MAFs of the various SNP classes
SFS.noncoding <- read.table("Noncoding_AFS.txt", header=F)
SFS.synonymous <- read.table("Synonymous_AFS.txt", header=F)  
SFS.tolerated <- read.table("Tolerated_AFS.txt", header=F)
SFS.deleterious <- read.table("Deleterious_AFS.txt", header=F)

#   Bin them up!
bins = seq(0.0, 1.0, by=0.2)
SFS.noncoding.bins <- cut(SFS.noncoding[,2], breaks=bins, include.lowest=TRUE)
SFS.synonymous.bins <- cut(SFS.synonymous[,2], breaks=bins, include.lowest=TRUE)
SFS.tolerated.bins <- cut(SFS.tolerated[,2], breaks=bins, include.lowest=TRUE)
SFS.deleterious.bins <- cut(SFS.deleterious[,2], breaks=bins, include.lowest=TRUE)

#   Count up the frequencies
SFS.noncoding.counts <- table(SFS.noncoding.bins)/length(SFS.noncoding[,2])
SFS.synonymous.counts <- table(SFS.synonymous.bins)/length(SFS.synonymous[,2])
SFS.tolerated.counts <- table(SFS.tolerated.bins)/length(SFS.tolerated[,2])
SFS.deleterious.counts <- table(SFS.deleterious.bins)/length(SFS.deleterious[,2])

#   Stick them together, in the order that we want them to appear
SFS.alldata <- as.data.frame(cbind(SFS.noncoding.counts, SFS.synonymous.counts, SFS.tolerated.counts, SFS.deleterious.counts))

#   Open a handle to a PDF file for the plot
pdf(file="DM_Barley_UnfoldedSFS_Classes.pdf", width=8, height=6,family="Helvetica",pointsize=16)
#   Make the barplot, we have to do a separate x-axis later
barplot(t(SFS.alldata) ,ylim=c(0,0.55),beside=TRUE, axisnames=F,xlab="Derived Allele Frequency",ylab="Proportion",col=c('black', 'grey', 'blue', 'red'))
#labels <- c("0-0.05","0.05-0.10","0.10-0.15","0.15-0.20","0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40","0.40-0.45","0.45-0.50")
#labels <- c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
labels <- c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0")

#   You'll have to play with these numbers to get the labels to line up with
#   the tick marks. As far as I know, there isn't a good way to figure this out
#   other than trial and error
#at <- c(3, 8, 13, 18, 23, 28, 33, 38, 43, 48)
at <- c(3, 8, 13, 18, 23)
axis(side=1, at=at, las=1,labels=labels,font=1,padj=1, cex.axis=0.8)
legend(inset=0,cex=1.1,"topright",c("Noncoding", "Synonymous", "Nonsynonymous Tolerated", "Nonsynonymous Deleterious"),fill=c('black', 'grey', 'blue', 'red'))
dev.off ()
