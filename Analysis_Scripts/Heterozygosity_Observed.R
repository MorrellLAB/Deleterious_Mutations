#   Plot the distribution of heterozygosity (as the number of individuals that
#   are heterozygous) for deleterious and non-deleterious SNPs.
#   Takes three arguments:
#       1) Deleterious SNP IDs
#       2) Nonsense SNP IDs
#       3) Sample alleles (from Sample_Allele_States.py)

#   Take arguments
args <- commandArgs(TRUE)

deleterious <- read.table(args[1], header=T)
deleterious <- deleterious$SNPID

nonsense <- read.table(args[2], header=T)
nonsense <- nonsense$SNPID

sample_alleles <- read.table(args[3], header=T)
sample_alleles <- sample_alleles[, c("SNPID", "Num_Ref", "Num_Het", "Num_Alt")]

#   Which are deleterious?
del_snps <- which(sample_alleles$SNPID %in% deleterious)
del_snps <- sample_alleles[del_snps, ]
non_snps <- which(sample_alleles$SNPID %in% nonsense)
non_snps <- sample_alleles[non_snps, ]
#   and tolerated?
tol_snps <- which(!(sample_alleles$SNPID %in% deleterious))
tol_snps <- sample_alleles[tol_snps, ]

#   Get the number of heterozygous accessions
del_Hobv <- as.numeric(del_snps$Num_Het)
non_Hobv <- as.numeric(non_snps$Num_Het)
tol_Hobv <- as.numeric(tol_snps$Num_Het)

#   Remove missing values
del_Hobv <- del_Hobv[!is.na(del_Hobv)]
non_Hobv <- non_Hobv[!is.na(non_Hobv)]
tol_Hobv <- tol_Hobv[!is.na(tol_Hobv)]

print(c(mean(del_Hobv), median(del_Hobv)))
print(c(mean(non_Hobv), median(non_Hobv)))
print(c(mean(tol_Hobv), median(tol_Hobv)))

#   Bin them up into side-by-side barplots like an SFS
bins <- seq(-1, 15, by=1)
del_Hobv <- cut(del_Hobv, breaks=bins, include.lowest=TRUE)
del_Hobv <- table(del_Hobv)/nrow(del_snps)
non_Hobv <- cut(non_Hobv, breaks=bins, include.lowest=TRUE)
non_Hobv <- table(non_Hobv)/nrow(non_snps)
tol_Hobv <- cut(tol_Hobv, breaks=bins, include.lowest=TRUE)
tol_Hobv <- table(tol_Hobv)/nrow(tol_snps)

ho_data <- cbind(tol_Hobv, del_Hobv, non_Hobv)

#   And plot some densities
pdf(file="Deleterious_Heterozygosity.pdf", 8, 8)
plt <- barplot(
    t(ho_data),
    beside=TRUE,
    ylim=c(0, 0.55),
    col=c("black", "red", "blue"),
    axisnames=FALSE,
    xlab="Number of Heterozygous Accessions",
    ylab="Proportion of SNPs",
    main="Distribution of Heterozygosity for Deleterious and Tolerated SNPs",
    )
at <- apply(plt, 2, mean)
axis(
     side=1,
     at=at,
     labels=as.character(seq(0,15,by=1)),
     font=1,
     cex.axis=1)
legend(
    "topright",
    c("Tolerated and Silent", "Deleterious", "Nonsense"),
    fill=c("black", "red", "blue")
    )
dev.off()
