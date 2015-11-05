#   Plot cM/Mb distance and proportion of deleterious variation over physical
#   distance in windows. This will only plot chromosome 1 of soybean, and will
#   use the physical and genetic map positions from the SoySNP6k (Lee et al.,
#   2015). Silent SNPs should not be considered.

#   To have plots with multiple axes
library(plotrix)
library(grDevices)
#   To run things in parallel
library(foreach)
library(doParallel)

#   Set up the parallel environment
cl <- makeCluster(8)
registerDoParallel(cl)

#   Define some parameters
#       Pericentromeric start and stop, in Mb
#       This is ugly, but oh well
pericentromeres <- data.frame(
    Chrom=c("Gm01", "Gm02", "Gm03", "Gm04", "Gm05", "Gm06", "Gm07", "Gm08",
            "Gm09", "Gm10", "Gm11", "Gm12", "Gm13", "Gm14", "Gm15", "Gm16",
            "Gm17", "Gm18", "Gm19", "Gm20"),
    Start=c(5700000, 15700000, 6300000, 10000000, 4400000, 18100000, 17800000,
            21800000, 7800000, 7600000, 17400000, 9100000, 8700000, 9600000,
            15700000, 7500000, 16100000, 7500000, 7200000, 3500000),
    End=c(47900000, 41800000, 33400000, 41200000, 31600000, 46300000, 35100000,
          41000000, 36700000, 37000000, 34600000, 32600000, 24000000, 44300000,
          48900000, 27600000, 38100000, 52800000, 34000000, 32500000)
    )

args <- commandArgs(TRUE)

genetic_map <- read.table(args[1], header=TRUE)
chroms <- unique(as.character(genetic_map$Chromosome))
#   Then read in the SNP table
snp_table <- read.table(args[2], header=TRUE)
for(k in chroms) {
    gm_curr <- genetic_map[genetic_map$Chromosome == k, ]
    keep <- grep(paste("^", k, sep=""), snp_table$Position, perl=TRUE, value=FALSE)
    snp_table_curr <- snp_table[keep, ]
    bp_pos <- strsplit(as.character(snp_table_curr$Position), ":")
    bp_pos <- as.numeric(sapply(bp_pos, "[[", 2))
    snp_table_curr <- data.frame(snp_table_curr, BP=bp_pos)
    #   And then, we want to build the list of postions to plot, and the respective
    #   recombination rates.
    plot_pos <- c()
    cMMb <- c()
    prop_del <- c()
    for(i in 1:(nrow(gm_curr)-1)) {
        j <- i + 1
        #   Generate the points to plot, which will be midpoints of intervals.
        #   Scale it all to Mb scale
        plot_pos[i] <- mean(gm_curr$V1Pos[i], gm_curr$V1Pos[j])/1000000
        Mb <- (gm_curr$V1Pos[j] - gm_curr$V1Pos[i]) / 1000000
        cM <- gm_curr$GeneticMapPos[j] - gm_curr$GeneticMapPos[i]
        cMMb[i] <- cM/Mb
        #   And then get the number of deleterious SNPs in this window
        in_window <- snp_table_curr$BP > as.numeric(gm_curr$V1Pos[i]) & snp_table_curr$BP < as.numeric(gm_curr$V1Pos[j])
        snps_in_window <- snp_table_curr[in_window, ]
        #   Which are deleterious?
        del <- snps_in_window$Deleterious == "Yes"
        del <- as.character(snp_table_curr$SNPID)[del]
        del <- unique(del)
        del <- length(del)
        #   And which are not?
        tol <- snps_in_window$Deleterious == "No"
        tol <- as.character(snp_table_curr$SNPID)[tol]
        tol <- unique(tol)
        tol <- length(tol)
        #   what is the proportion?
        prop_del[i] <- del/(del + tol)
    }
    #   make a plot
    pdf(file=paste("Recomb_DM_", k, ".pdf", sep=""), width=10, height=8)
    twoord.plot(
        lx=plot_pos,
        ly=cMMb,
        rx=plot_pos,
        ry=prop_del,
        main=paste("Proportion of Nonsynonymous SNPs Called\nDeleterious and Recombination Rate on ", k, sep=""),
        lcol="black",
        rcol="red",
        ylab="cM/Mb",
        rylab="D.S. Proportion",
        xlab="Physical Position (Mb)",
        lylim=c(-2, 20),
        rylim=c(-0.1, 1),
        type="p",
        lpch=19,
        rpch=19
        )
    abline(
        v=pericentromeres[pericentromeres$Chrom == k, 2]/1000000,
        lty=3,
        lwd=4,
        col="#aaaaaa"
        )
    abline(
        v=pericentromeres[pericentromeres$Chrom == k, 3]/1000000,
        lty=3,
        lwd=4,
        col="#aaaaaa"
        )
    dev.off()
}

#   And then, do a genome wide correlation
recomb <- c()
DM <- c()
for(k in 1:length(chroms)) {
    current_chr <- chroms[k]
    gm_curr <- genetic_map[genetic_map$Chromosome == current_chr, ]
    keep <- grep(paste("^", current_chr, sep=""), snp_table$Position, perl=TRUE, value=FALSE)
    snp_table_curr <- snp_table[keep, ]
    #   Now, we want to build a vector of basepair positions where the SNPs land
    bp_pos <- strsplit(as.character(snp_table_curr$Position), ":")
    bp_pos <- as.numeric(sapply(bp_pos, "[[", 2))
    snp_table_curr <- data.frame(snp_table_curr, BP=bp_pos)
    for(i in 1:(nrow(gm_curr)-1)) {
        j <- i + 1
        #   Generate the points to plot, which will be midpoints of intervals.
        #   Scale it all to Mb scale
        Mb <- (gm_curr$V1Pos[j] - gm_curr$V1Pos[i]) / 1000000
        cM <- gm_curr$GeneticMapPos[j] - gm_curr$GeneticMapPos[i]
        recomb <- c(recomb, (cM/Mb))
        #   And then get the number of deleterious SNPs in this window
        in_window <- snp_table_curr$BP > as.numeric(gm_curr$V1Pos[i]) & snp_table_curr$BP < as.numeric(gm_curr$V1Pos[j])
        snps_in_window <- snp_table_curr[in_window, ]
        #   Which are deleterious?
        del <- snps_in_window$Deleterious == "Yes"
        del <- as.character(snp_table_curr$SNPID)[del]
        del <- unique(del)
        del <- length(del)
        #   And which are not?
        tol <- snps_in_window$Deleterious == "No"
        tol <- as.character(snp_table_curr$SNPID)[tol]
        tol <- unique(tol)
        tol <- length(tol)
        #   what is the proportion?
        prop_del_win <- del/(del+tol)
        DM <- c(DM, prop_del_win)
    }
}

#   For the correlation, we do a little bit of pruning
rec_DM <- data.frame(
    recombination=recomb,
    prop_deleterious=DM
    )
#   First, remove rows that have negative recombination rate
#   These are markers that are in the wrong order
wrong_order <- rec_DM$recombination < 0
rec_DM <- rec_DM[!wrong_order, ]
#   Then, remove rows with recombination > 20, these are caused by other
#   markers being in the wrong order (HEURISTIC)
too_high <- rec_DM$recombination > 20
rec_DM <- rec_DM[!too_high, ]

#   And plot a regression
reg <- lm(rec_DM$prop_deleterious ~ rec_DM$recombination)
CI <- confint(reg)
reg_summary <- summary(reg)
rsq <- cor(rec_DM$recombination, rec_DM$prop_deleterious, method="pearson", use="pairwise.complete.obs")^2
#   And get a p_value
emp <- foreach(icount(1000000)) %dopar% {
    shuf_DM <- sample(rec_DM$prop_deleterious)
    shuf_recom <- sample(rec_DM$recombination)
    new_cor <- cor(shuf_DM, shuf_recom, method="pearson", use="pairwise.complete.obs")
    new_cor^2
}
emp <- unlist(emp)
alt <- sum(rsq <= emp)
pval <- alt/length(emp)
pdf(file="Cor_Resampling.pdf", 6, 6)
hist(emp, main="Bootstrap Distribution of Correlation Coefficients", xlab=expression(r^2), ylab="Frequency")
abline(v=rsq, col="red")
dev.off()
pdf(file="Recomb_DM_Corr.pdf", width=6, height=6)
plot(
    rec_DM$recombination,
    rec_DM$prop_deleterious,
    main="Proportion of Nonsyn. SNPs Called \nDeleterious Over Recombination Rate",
    xlab="cM/Mb",
    ylab="Proportion",
    type="p",
    pch=19,
    cex=0.3,
    col=rgb(0, 0, 0, 0.25)
    )
abline(reg, col="red", lty=1, lwd=1)
abline(coef=CI[, 1], col="darkgreen", lty=1, lwd=2)
abline(coef=CI[, 2], col="darkgreen", lty=1, lwd=2)
legend(
    "topright",
    bty="n",
    legend=c(
        paste("R2 = ", format(rsq, digits=3)),
        paste("P = ", format(pval, digits=3))
        )
    )
dev.off()
stopCluster(cl)
