#   Read the MAFs of the various SNP classes
SFS <- read.table("Derived_SFS.txt", header=F)
#Sim <- read.table("Simulated_Derived_SFS.txt", header=F)

#   Bin them up!
#bins = seq(0.0, 0.5, by=0.05)
bins = seq(0.0, 1.0, by=0.2)
SFS.bins <- cut(SFS$V2, breaks=bins, include.lowest=TRUE)

#   Count up the frequencies
SFS.counts <- table(SFS.bins)/length(SFS$V2)
#Sim.counts <- Sim$V2/sum(Sim$V2)

#SFS.alldata <- as.data.frame(cbind(SFS.counts, Sim.counts))
#   Open a handle to a PDF file for the plot
pdf(file="DM_Barley_UnfoldedSFS.pdf", width=8, height=6,family="Helvetica",pointsize=16)
#   Make the barplot, we have to do a separate x-axis
barplot(t(SFS.counts) ,ylim=c(0,0.45),beside=TRUE, axisnames=F,xlab="Derived Allele Frequency",ylab="Proportion",col=grey(0.8))
#labels <- c("0-0.05","0.05-0.10","0.10-0.15","0.15-0.20","0.20-0.25","0.25-0.30","0.30-0.35","0.35-0.40","0.40-0.45","0.45-0.50")
#labels <- c("0-0.1","0.1-0.2","0.2-0.3","0.3-0.4","0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
labels <- c("0-0.2","0.2-0.4","0.4-0.6","0.6-0.8","0.8-1.0")

#at <- c(2, 5, 8, 11, 14)
at <- c(1.5, 3.5, 5.5, 7.5, 9.5)
axis(side=1, at=at, las=1,labels=labels,font=1,padj=1, cex.axis=0.8)
dev.off ()
