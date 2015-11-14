#	Plot the SFS for deleterious mutations, according to prediction method.
#	Will produce two figures:
#		1) Side-by-side SFS for each program (SIFT, PPH2, LRT)
#		2) Side-by-side SFS for variants in 1, 2, or all prediction methods.

#	Take arguments
args <- commandArgs(TRUE)
snp_table <- read.table(args[1], header=T)

#	Set the LRT significance threshold
#		We tested 59,277 codons
lrt_sig <- 0.05/59277

#	Get those deleterious by each method
sift <- snp_table[snp_table$SIFT == "DELETERIOUS", ]
pph <- snp_table[snp_table$PPH == "deleterious", ]
lrt <- snp_table[snp_table$MaskedP.value <= lrt_sig, ]

#	We can plot the SFS for each one
bins <- seq(0.0, 1.0, by = 0.1)
sfs.sift <- cut(sift$DAF, breaks=bins, include.lowest=TRUE)
sfs.pph <- cut(pph$DAF, breaks=bins, include.lowest=TRUE)
sfs.lrt <- cut(lrt$DAF, breaks=bins, include.lowest=TRUE)

sfs.sift <- table(sfs.sift)/sum(!is.na(sift$DAF))
sfs.pph <- table(sfs.pph)/sum(!is.na(pph$DAF))
sfs.lrt <- table(sfs.lrt)/sum(!is.na(lrt$DAF))

sfs.data <- as.data.frame(
	cbind(
		sfs.sift,
		sfs.pph,
		sfs.lrt
		)
	)

pdf(
	file="Deleterious_SFS_By_Program.pdf",
	width=12,
	height=6,
	family="Helvetica",
	pointsize=16
	)
plt <- barplot(
	t(sfs.data),
	ylim=c(0, 0.4),
	beside=TRUE,
	axisnames=F,
	xlab="Derived Allele Frequency",
	ylab="Proportion",
	col=c("red", "blue", "green")
	)
labels <- c(
	"[0, 0.1]",
	"(0.1, 0.2]",
	"(0.2, 0.3]",
	"(0.3, 0.4]",
	"(0.4, 0.5]",
	"(0.5, 0.6]",
	"(0.6, 0.7]",
	"(0.7, 0.8]",
	"(0.8, 0.9]",
	"(0.9, 1.0]"
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
onlysift <- as.vector(snp_table[(snp_table$SIFT == "DELETERIOUS" & snp_table$PPH != "deleterious" & snp_table$MaskedP.value > lrt_sig), "DAF" ])
onlypph <- as.vector(snp_table[(snp_table$SIFT != "DELETERIOUS" & snp_table$PPH == "deleterious" & snp_table$MaskedP.value > lrt_sig), "DAF" ])
onlylrt <- (snp_table[(snp_table$SIFT != "DELETERIOUS" & snp_table$PPH != "deleterious" & snp_table$MaskedP.value <= lrt_sig), "DAF" ])
sfs.one <- c(
	onlysift,
	onlypph,
	onlylrt
	)
#		Two methods
sfs.two <- c(
	as.vector(sift[sift$SNPID %in% pph$SNPID, "DAF"]),
	as.vector(sift[sift$SNPID %in% lrt$SNPID, "DAF"]),
	as.vector(pph[pph$SNPID %in% lrt$SNPID, "DAF"])
	)
#		All three
sfs.three <- as.vector(
	sift[(sift$SNPID %in% pph$SNPID) & (sift$SNPID %in% lrt$SNPID), "DAF"]
	)

#	Bin them up
sfs.one.bin <- cut(sfs.one, breaks=bins, include.lowest=TRUE)
sfs.two.bin <- cut(sfs.two, breaks=bins, include.lowest=TRUE)
sfs.three.bin <- cut(sfs.three, breaks=bins, include.lowest=TRUE)

sfs.one.freq <- table(sfs.one.bin)/sum(!is.na(sfs.one))
sfs.two.freq <- table(sfs.two.bin)/sum(!is.na(sfs.two))
sfs.three.freq <- table(sfs.three.bin)/sum(!is.na(sfs.three))

sfs.pred.data <- as.data.frame(
	cbind(
		sfs.one.freq,
		sfs.two.freq,
		sfs.three.freq
		)
	)

pdf(
	file="Deleterious_SFS_By_Intersect.pdf",
	width=12,
	height=6,
	family="Helvetica",
	pointsize=16
	)
plt <- barplot(
	t(sfs.pred.data),
	ylim=c(0, 0.4),
	beside=TRUE,
	axisnames=F,
	xlab="Derived Allele Frequency",
	ylab="Proportion",
	col=c("red", "blue", "green")
	)
labels <- c(
	"[0, 0.1]",
	"(0.1, 0.2]",
	"(0.2, 0.3]",
	"(0.3, 0.4]",
	"(0.4, 0.5]",
	"(0.5, 0.6]",
	"(0.6, 0.7]",
	"(0.7, 0.8]",
	"(0.8, 0.9]",
	"(0.9, 1.0]"
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
