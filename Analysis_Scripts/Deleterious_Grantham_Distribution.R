#   Take arguments
args <- commandArgs(TRUE)

#	Define the Grantham matrix. What a pain.
grantham <- matrix(
c(0, 112, 111, 126, 195, 91, 107, 60, 86, 94, 96, 106, 84, 113, 27, 99, 58, 148, 112, 64,
112, 0, 86, 96, 180, 43, 54, 125, 29, 97, 102, 26, 91, 97, 103, 110, 71, 101, 77, 96,
111, 86, 0, 23, 139, 46, 42, 80, 68, 149, 153, 94, 142, 158, 91, 46, 65, 174, 143, 133,
126, 96, 23, 0, 154, 61, 45, 94, 81, 168, 172, 101, 160, 177, 108, 65, 85, 181, 160, 152,
195, 180, 139, 154, 0, 154, 170, 159, 174, 198, 198, 202, 196, 205, 169, 112, 149, 215, 194, 192,
91, 43, 46, 61, 154, 0, 29, 87, 24, 109, 113, 53, 101, 116, 76, 68, 42, 130, 99, 96,
107, 54, 42, 45, 170, 29, 0, 98, 40, 134, 138, 56, 126, 140, 93, 80, 65, 152, 122, 121,
60, 125, 80, 94, 159, 87, 98, 0, 98, 135, 138, 127, 127, 153, 42, 56, 59, 184, 147, 109,
86, 29, 68, 81, 174, 24, 40, 98, 0, 94, 99, 32, 87, 100, 77, 89, 47, 115, 83, 84,
94, 97, 149, 168, 198, 109, 134, 135, 94, 0, 5, 102, 10, 21, 95, 142, 89, 61, 33, 29,
96, 102, 153, 172, 198, 113, 138, 138, 99, 5, 0, 107, 15, 22, 98, 145, 92, 61, 36, 32,
106, 26, 94, 101, 202, 53, 56, 127, 32, 102, 107, 0, 95, 102, 103, 121, 78, 110, 85, 97,
84, 91, 142, 160, 196, 101, 126, 127, 87, 10, 15, 95, 0, 28, 87, 135, 81, 67, 36, 21,
113, 97, 158, 177, 205, 116, 140, 153, 100, 21, 22, 102, 28, 0, 114, 155, 103, 40, 22, 50,
27, 103, 91, 108, 169, 76, 93, 42, 77, 95, 98, 103, 87, 114, 0, 74, 38, 147, 110, 68,
99, 110, 46, 65, 112, 68, 80, 56, 89, 142, 145, 121, 135, 155, 74, 0, 58, 177, 144, 124,
58, 71, 65, 85, 149, 42, 65, 59, 47, 89, 92, 78, 81, 103, 38, 58, 0, 128, 92, 69,
148, 101, 174, 181, 215, 130, 152, 184, 115, 61, 61, 110, 67, 40, 147, 177, 128, 0, 37, 88,
112, 77, 143, 160, 194, 99, 122, 147, 83, 33, 36, 85, 36, 22, 110, 144, 92, 37, 0, 55,
64, 96, 133, 152, 192, 96, 121, 109, 84, 29, 32, 97, 21, 50, 68, 124, 69, 88, 55, 0),
nrow=20,
ncol=20,
dimnames=list(
    c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"),
    c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
    )
)

#   A function to just pull values out of the Grantham matrix
GranthamScore <- function(x) {
    if(x["AA1"] == "*" | x["AA2"] == "*") {
        return(NA)
    }
    else {
        return(grantham[x["AA1"], x["AA2"]])
    }
}

#   Read the predictions table
preds <- read.table(args[1], header=T)

#   Figure out the LRT significance threshold. We will use 0.05 divided by the
#   number of nonsynonymous variants
lrt_sig <- 0.05/59277
lrt_min_seq <- 10

#   What are the SIFT deleterious SNPs?
sift_del <- preds[which(preds$SIFT == "DELETERIOUS"),]
#   And PPH?
pph_del <- preds[which(preds$PPH == "deleterious"),]
#   And LRT
lrt_del <- preds[which(preds$MaskedP.value <= lrt_sig & preds$SeqCount >= lrt_min_seq), ]
#   Finally, get the intersection of all these 
intersect_del <- intersect(
    intersect(sift_del$SNPID, pph_del$SNPID),
    lrt_del$SNPID)
intersect_del <- preds[preds$SNPID %in% intersect_del, ]

#   Grantham scores by prediction method. Be sure to remove missing data
sift_grantham <- apply(sift_del, 1, GranthamScore)
sift_grantham <- sift_grantham[!is.na(sift_grantham)]
pph_grantham <- apply(pph_del, 1, GranthamScore)
pph_grantham <- pph_grantham[!is.na(pph_grantham)]
lrt_grantham <- apply(lrt_del, 1, GranthamScore)
lrt_grantham <- lrt_grantham[!is.na(lrt_grantham)]
intersect_grantham <- apply(intersect_del, 1, GranthamScore)
intersect_grantham <- intersect_grantham[!is.na(intersect_grantham)]

#   And plot them
pdf(file="Grantham_By_Method.pdf", 8, 8)
plot(
    density(sift_grantham),
    col="red",
    lwd=2,
    lty=1,
    xlim=c(0, 250),
    ylim=c(0,0.02),
    xlab="Grantham Score",
    ylab="Density",
    main="Distribution of Grantham Scores"
    )
abline(
    v=mean(sift_grantham),
    col="red",
    lwd=1,
    lty=1
    )
lines(
    density(pph_grantham),
    col="blue",
    lwd=2,
    lty=2
    )
abline(
    v=mean(pph_grantham),
    col="blue",
    lwd=1,
    lty=2
)
lines(
    density(lrt_grantham),
    col="green",
    lwd=2,
    lty=3
    )
abline(
    v=mean(lrt_grantham),
    col="green",
    lwd=1,
    lty=3
)
lines(
    density(intersect_grantham),
    col="orange",
    lwd=2,
    lty=4
    )
abline(
    v=mean(intersect_grantham),
    col="orange",
    lwd=1,
    lty=4
)
legend(
    "topright",
    c("SIFT Deleterious", "PPH Deleterious", "LRT Deleterious", "Intersect"),
    col=c("red", "blue", "green", "orange"),
    lwd=2,
    lty=c(1, 2, 3, 4),
    cex=1
    )
dev.off()
