# required packages
require(gaston) # for many functions
## To install pcamethods
# source("https://bioconductor.org/biocLite.R")
# biocLite("pcaMethods")
require(pcaMethods) # for probablistic PCA (a package in Bioconductor)
require(BGLR) # for Bayesian regression
# read file
bm <- read.vcf("HDRA-G6-4-RDP1-RDP2-NIAS.AGCT_MAF005MIS001NR.vcf.gz")
bm
# read phenotypic data and line data
pheno <- read.csv("RiceDiversityPheno4GWASGS.csv", row.names = 1)
str(pheno)
head(bm@ped)
head(bm@snps)
# select samples with phenotypic data
bm.wp <- bm[bm@ped$id %in% rownames(pheno),]
bm.wp
# remove remove the monomorphic markers
bm.wp <- select.snps(bm.wp, maf > 0)
bm.wp
# extract the score matrix
gt.score <- as.matrix(bm.wp)
head(gt.score)[,1:10]
dim(gt.score)
# reorder the scores in the same order as the phenotypic data
gt.score <- gt.score[rownames(pheno),]
dim(gt.score)


# k-means clustering
tmp <- na.omit(t(gt.score))
dim(tmp)
km <- kmeans(t(tmp), centers = 5, nstart = 10, iter.max = 100)
grp <- as.factor(km$cluster)

# choose a trait: here "Seed.length.width.ratio" 
bm.wp@ped$pheno <- pheno[bm.wp@ped$id, "Seed.length.width.ratio"]
# remove lines with missing trait
bm.wom <- select.inds(bm.wp, !is.na(bm.wp@ped$pheno))
bm.wom
# keep marker a large enough MAF (>0.05) and low missing rate (callrate>0.9)
bm.wom <- select.snps(bm.wom, maf > 0.05)
bm.wom <- select.snps(bm.wom, callrate > 0.9)
bm.wom

# Compute K (=Genetic relationship matrix for the lines with observed/non missing trait)
K <- GRM(bm.wom)
dim(K)
# Compute fixed effects from k-means clustering
Xkm <- model.matrix(~grp[bm.wom@ped$id])
head(Xkm)

gwa <- association.test(bm.wom, method = "lmm", response = "quantitative", test = "lrt", eigenK = eigen(K), p = 4)
head(gwa)
gwa.wald <- association.test(bm.wom, method = "lmm", response = "quantitative", test = "wald", eigenK = eigen(K), p = 4)
head(gwa.wald)
gwa.score <- association.test(bm.wom, method = "lmm", response = "quantitative", test = "score", K = K, eigenK = eigen(K), p = 4)
head(gwa.score)

plot(-log10(gwa$p), -log10(gwa.wald$p))
abline(0,1)
plot(-log10(gwa$p), -log10(gwa.score$p))
abline(0,1)

manhattan(gwa, pch = 20)

p.adj <- p.adjust(gwa$p, method = "bonferroni")
gwa[p.adj < 0.05, ]

fdr <- p.adjust(gwa$p, method = "BH")
gwa[fdr < 0.05, ]
col <- rep("black", nrow(gwa))
col[gwa$chr %% 2 == 0] <- "gray50"
col[fdr < 0.05] <- "green"
manhattan(gwa, pch = 20, col = col)
col <- rep("black", nrow(gwa))
col[gwa$chr %% 2 == 0] <- "gray50"
col[p.adj < 0.05] <- "green"
manhattan(gwa, pch = 20, col = col)

gwa.km <- association.test(bm.wom, method = "lmm", response = "quantitative", test = "lrt", X = Xkm, eigenK = eigen(K))
p.adj.km <- p.adjust(gwa.km$p, method = "bonferroni")
col <- rep("black", nrow(gwa.km))
col[gwa.km$chr %% 2 == 0] <- "gray50"
col[p.adj.km < 0.05] <- "green"
manhattan(gwa.km, pch = 20, col = col)

gwa.km.p <- association.test(bm.wom, method = "lmm", response = "quantitative", test = "lrt", X = Xkm, eigenK = eigen(K), p=4)

gwa.no <- association.test(bm.wom, method = "lmm", response = "quantitative", test = "lrt", eigenK = eigen(K))
qqplot.pvalues(gwa$p, pch = 20, ylim = range(-log10(c(gwa$p, gwa.km$p, gwa.no$p))))
qqplot.pvalues(gwa.km$p, pch = 20, ylim = range(-log10(c(gwa$p, gwa.km$p, gwa.no$p))))
qqplot.pvalues(gwa.no$p, pch = 20, ylim = range(-log10(c(gwa$p, gwa.km$p, gwa.no$p))))
gwa.lm <- association.test(bm.wom, method = "lm", response = "quantitative", test = "wald")
qqplot.pvalues(gwa.lm$p, pch = 20)


selector <- gwa$chr == 5 & fdr < 0.05
gwa[selector,]
from <- 15882 - 1
to <- 15912 + 1

ld <- LD(bm.wom, c(from, to), measure = "r2")
LD.plot(ld, snp.positions = bm.wom@snps$pos[from:to])
LD.plot(ld, snp.positions = bm.wom@snps$pos[from:to], pdf.file = "ldplot.pdf")

bm.wp@ped$pheno <- pheno[bm.wp@ped$id, "Leaf.pubescence"]
bm.wom <- select.inds(bm.wp, !is.na(bm.wp@ped$pheno))
bm.wom <- select.snps(bm.wom, maf > 0.05)
bm.wom <- select.snps(bm.wom, callrate > 0.9)
bm.wom
K <- GRM(bm.wom)
gwa <- association.test(bm.wom, method = "lmm", response = "binary", test = "score", K = K, eigenK = eigen(K), p = 4)
manhattan(gwa, pch = 20)
qqplot.pvalues(gwa$p, pch = 20)
fdr <- p.adjust(gwa$p, method = "BH")
gwa[fdr < 0.05, ]

bm.thin <- LD.thin(bm.wp, threshold = 0.4)
bm.thin

bm.thin@ped$pheno <- pheno[bm.thin@ped$id, "Seed.length.width.ratio"]
bm.thin.wom <- select.inds(bm.thin, !is.na(bm.thin@ped$pheno))
bm.thin.wom <- select.snps(bm.thin.wom, maf > 0.05)
bm.thin.wom <- select.snps(bm.thin.wom, callrate > 0.9)
bm.thin.wom
X <- as.matrix(bm.thin.wom) - 1
for(i in 1:ncol(X)) {
	X[,i][is.na(X[,i])] <- mean(X[,i], na.rm = T)
}
Xkm <- model.matrix(~grp[bm.thin.wom@ped$id])
grm <- GRM(bm.thin.wom)
ETA <- list(list(X = X, model = "BayesB"), list(X = Xkm, model = "FIXED"), list(K = grm, model = "RKHS"))
y <- bm.thin.wom@ped$pheno
model <- BGLR(y = y, ETA = ETA, nIter = 4000, burnIn = 1000)
col <- rep("black", ncol(bm.thin.wom))
col[bm.thin.wom@snps$chr %% 2 == 0] <- "gray50"
plot(abs(model$ETA[[1]]$b), pch = 20, col = col)

