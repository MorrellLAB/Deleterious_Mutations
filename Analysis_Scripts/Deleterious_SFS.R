#	Plot the SFS for deleterious mutations, according to prediction method.
#	Will produce two figures:
#		1) Side-by-side SFS for each program (SIFT, PPH2, LRT)
#		2) Side-by-side SFS for variants in 1, 2, or all prediction methods.
#	Currently goes off of minor allele frequency, but will be modified to
#	go on derived allele frequency as soon as it is ready

#	Set the LRT significance threshold
lrt_sig <- 0.01/18000
#	Read in the table that JCF sent
snp_table <- read.csv("Barley_Exome_LRT.csv", header=T)

#	Get those deleterious by each method
sift <- snp_table[snp_table$SIFT == "DELETERIOUS", ]
pph <- snp_table[snp_table$PPH == "deleterious", ]
lrt <- snp_table[snp_table$MaskedP.value <= lrt_sig, ]

#	We can plot the SFS for each one
bins <- seq(0.0, 0.5, by = 0.1)
sfs.sift <- cut(sift$MAF, breaks=bins, include.lowest=TRUE)
sfs.pph <- cut(pph$MAF, breaks=bins, include.lowest=TRUE)
sfs.lrt <- cut(lrt$MAF, breaks=bins, include.lowest=TRUE)

sfs.sift <- table(sfs.sift)/length(sift$MAF)
sfs.pph <- table(sfs.pph)/length(pph$MAF)
sfs.lrt <- table(sfs.lrt)/length(lrt$MAF)

sfs.data <- as.data.frame(
	cbind(
		sfs.sift,
		sfs.pph,
		sfs.lrt
		)
	)

pdf(
	file="Deleterious_SFS_By_Program.pdf",
	width=8,
	height=6,
	family="Helvetica",
	pointsize=16
	)
plt <- barplot(
	t(sfs.data),
	ylim=c(0, 0.5),
	beside=TRUE,
	axisnames=F,
	xlab="Minor Allele Frequency",
	ylab="Proportion",
	col=c("red", "blue", "green")
	)
labels <- c(
	"[0, 0.1)",
	"[0.1, 0.2)",
	"[0.2, 0.3)",
	"[0.3, 0.4)",
	"[0.4, 0.5]"
	)
at <- apply(plt, 2, mean)
axis(
	 side=1,
	 at=at,
	 labels=labels,
	 font=1,
	 cex.axis=0.75
	 )
legend(
	"topright",
	c("SIFT", "PolyPhen2", "LRT"),
	fill=c("red", "blue", "green"),
	cex=1.0
	)
dev.off()

#	Now, build data frames for SNPs that are called deleterious by one, two
#	or three prediction methods.
#		One prediction method only
onlysift <- as.vector(snp_table[(snp_table$SIFT == "DELETERIOUS" & snp_table$PPH != "deleterious" & snp_table$MaskedP.value > lrt_sig), "MAF" ])
onlypph <- as.vector(snp_table[(snp_table$SIFT != "DELETERIOUS" & snp_table$PPH == "deleterious" & snp_table$MaskedP.value > lrt_sig), "MAF" ])
onlylrt <- (snp_table[(snp_table$SIFT != "DELETERIOUS" & snp_table$PPH != "deleterious" & snp_table$MaskedP.value <= lrt_sig), "MAF" ])
sfs.one <- c(
	onlysift,
	onlypph,
	onlylrt
	)
#		Two methods
sfs.two <- c(
	as.vector(sift[sift$SNPID %in% pph$SNPID, "MAF"]),
	as.vector(sift[sift$SNPID %in% lrt$SNPID, "MAF"]),
	as.vector(pph[pph$SNPID %in% lrt$SNPID, "MAF"])
	)
#		All three
sfs.three <- as.vector(
	sift[(sift$SNPID %in% pph$SNPID) & (sift$SNPID %in% lrt$SNPID), "MAF"]
	)

#	Bin them up
sfs.one.bin <- cut(sfs.one, breaks=bins, include.lowest=TRUE)
sfs.two.bin <- cut(sfs.two, breaks=bins, include.lowest=TRUE)
sfs.three.bin <- cut(sfs.three, breaks=bins, include.lowest=TRUE)

sfs.one.freq <- table(sfs.one.bin)/length(sfs.one)
sfs.two.freq <- table(sfs.two.bin)/length(sfs.two)
sfs.three.freq <- table(sfs.three.bin)/length(sfs.three)

sfs.pred.data <- as.data.frame(
	cbind(
		sfs.one.freq,
		sfs.two.freq,
		sfs.three.freq
		)
	)

pdf(
	file="Deleterious_SFS_By_Intersect.pdf",
	width=8,
	height=6,
	family="Helvetica",
	pointsize=16
	)
plt <- barplot(
	t(sfs.pred.data),
	ylim=c(0, 0.5),
	beside=TRUE,
	axisnames=F,
	xlab="Minor Allele Frequency",
	ylab="Proportion",
	col=c("red", "blue", "green")
	)
labels <- c(
	"[0, 0.1)",
	"[0.1, 0.2)",
	"[0.2, 0.3)",
	"[0.3, 0.4)",
	"[0.4, 0.5]"
	)
at <- apply(plt, 2, mean)
axis(
	 side=1,
	 at=at,
	 labels=labels,
	 font=1,
	 cex.axis=0.75
	 )
legend(
	"topright",
	c("One Method", "Two Methods", "Three Methods"),
	fill=c("red", "blue", "green"),
	cex=1.0
	)
dev.off()


