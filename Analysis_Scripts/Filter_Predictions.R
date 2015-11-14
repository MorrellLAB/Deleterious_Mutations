#	Filter the predictions to the list of deleterious variants as the intersect
#	of three prediction methods.

args <- commandArgs(TRUE)

#	Read in the effect prediction table
effects <- read.table(args[1], header=T)
#	Create the filename for output
outfile <- gsub(".txt", "_Deleterious.txt", args[1])

# Set the LRT significance threshold
#		First is to calcuate the significance threshold. We use a conservative
#		cutoff of 0.05/N_tests, which is 0.05/N_NS-SNPs
n_nsnps <- sum(effects["Silent"] == "No")
lrt_sig <- 0.05/n_nsnps
#		Then the minimum number of sequences we use for prediction. We start
#		with 10 sequences minimum
min_seq <- 10

#	Slice down the table to nonsynonymous only
nonsyn <- effects[effects["Silent"] == "No",]
#	If it is deleterious by all three methods, then it is inferred deleteirous
deleterious <- (nonsyn["SIFT"] == "DELETERIOUS") & (nonsyn["PPH"] == "deleterious") & (nonsyn["MaskedP.value"] <= lrt_sig & nonsyn["SeqCount"] >= min_seq)
deleterious_table <- nonsyn[deleterious,]

#	Then write it out
write.table(
	deleterious_table,
	file=outfile,
	sep="\t",
	quote=F,
	row.names=F,
	col.names=T
	)
