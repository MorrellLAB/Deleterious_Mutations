#	Plot the unfolded site frequency spectrum for synonymous, tolerated, and
#	deleterious SNPs. Due to the nature of the ancestral state inference, all
#	SNPs are in coding regions.
#		Synonymous: Does not alter amino acid sequence
#		Tolerated: Alters AA sequence, not deleterious by at least 1 method
#		Deleterious: Alters AA sequence, deleterious by SIFT,PPH,LRT

args <- commandArgs(TRUE)

#	Read in the effect prediction table
effects <- read.table(args[1], header=T)

#	Set the LRT significance threshold
#		We tested 59,865 codons for barley
#lrt_sig <- 0.05/59865
#		And tested 64,087 for soy
lrt_sig <- 0.05/64087
#	We only want to consider codons with at least 10 sequences represented
minseq <- 10
#	And those with constraint < 1
maxconstraint <- 1

#	What is the DAF for synonymous SNPs
daf.syn <- effects[effects["Silent"] == "Yes", ]

#	If at least one method is tolerated, it gets put into the tolerated bin
nonsyn <- effects[effects["Silent"] == "No",]
tolerated <- (nonsyn["SIFT"] == "TOLERATED") | (nonsyn["PPH"] == "neutral") | (nonsyn["MaskedP.value"] >= lrt_sig & nonsyn["SeqCount"] >= 10)
daf.tol <- nonsyn[tolerated,]
#	If it is deleterious by all three methods, then it is inferred deleteirous
deleterious <- (nonsyn["SIFT"] == "DELETERIOUS") & (nonsyn["PPH"] == "deleterious") & (nonsyn["MaskedP.value"] < lrt_sig & nonsyn["SeqCount"] >= 10 & nonsyn["MaskedConstraint"] < maxconstraint & (nonsyn["Rn"] == 1 | nonsyn["An"] == 0))
daf.del <- nonsyn[deleterious,]

daf.syn.nomiss <- daf.syn[!is.na(daf.syn["DAF"]), "DAF"]
daf.tol.nomiss <- daf.tol[!is.na(daf.tol["DAF"]), "DAF"]
daf.del.nomiss <- daf.del[!is.na(daf.del["DAF"]), "DAF"]

length(daf.syn.nomiss)
length(daf.tol.nomiss)
length(daf.del.nomiss)

#	We can plot the SFS for each one
#		Again, for barley, we will bin by 10% classes
bins <- seq(0.0, 1.0, by = 0.1)
#	We can plot the SFS for each one
#		And 20% for soy.
#bins <- seq(0.0, 1.0, by = 0.2)
#	And then the SFS
sfs.syn <- cut(daf.syn.nomiss, breaks=bins, include.lowest=TRUE)
sfs.tol <- cut(daf.tol.nomiss, breaks=bins, include.lowest=TRUE)
sfs.del <- cut(daf.del.nomiss, breaks=bins, include.lowest=TRUE)

sfs.syn.prop <- table(sfs.syn)/length(sfs.syn)
sfs.tol.prop <- table(sfs.tol)/length(sfs.tol)
sfs.del.prop <- table(sfs.del)/length(sfs.del)
sfs.data <- cbind(
	sfs.syn.prop,
	sfs.tol.prop,
	sfs.del.prop)

#	And make the plot
pdf(
	file="SFS_by_Functional_Impact.pdf",
	width=6,
	height=6,
	family="Helvetica",
	pointsize=12
	)
plt <- barplot(
	t(sfs.data),
	ylim=c(0, 0.5),
	beside=TRUE,
	axisnames=F,
	xlab="Derived Allele Frequency",
	ylab="Proportion",
	col=c("black", "blue", "red")
	)
#	For barley, we bin up into 10% classes
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
#	For soy, we bin up into 20% classes, since our sample size is a lot
#	smaller
# labels <- c(
# 	"[0, 0.2]",
# 	"(0.2, 0.4]",
# 	"(0.4, 0.6]",
# 	"(0.6, 0.8]",
# 	"(0.8, 1.0]"
# 	)
at <- apply(plt, 2, mean)
axis(
	 side=1,
	 at=at,
	 labels=labels,
	 font=1,
	 cex.axis=0.75,
	 las=1
	 )
legend(
	"topright",
	c("Synonymous", "Tolerated", "Deleterious"),
	fill=c("black", "blue", "red"),
	cex=0.75
	)
dev.off()

