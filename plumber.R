#
# This is a Plumber API. In RStudio 1.2 or newer you can run the API by
# clicking the 'Run API' button above.
#
# In RStudio 1.1 or older, see the Plumber documentation for details
# on running the API.
#
# Find out more about building APIs with Plumber here:
#
#    https://www.rplumber.io/
#

library(plumber)


#* @apiTitle GWAS API
#* @apiDescription First test for GWAS API

# required packages
require(gaston) # for many functions

#####################################################################
# load data
# read file
bm <- read.vcf("data/HDRA-G6-4-RDP1-RDP2-NIAS.AGCT_MAF005MIS001NR.vcf.gz")
# read phenotypic data and line data
pheno <- read.csv("data/RiceDiversityPheno4GWASGS.csv", row.names = 1)
# select markers with traits and remove monomorphic markers
bm.wp <- bm[bm@ped$id %in% rownames(pheno),]
bm.wp <- select.snps(bm.wp, maf > 0)

# extract the score matrix and reorder the scores in the same order as the phenotypic data
gt.score <- as.matrix(bm.wp)
gt.score <- gt.score[rownames(pheno),]

# k-means clustering
tmp <- na.omit(t(gt.score))
km <- kmeans(t(tmp), centers = 5, nstart = 10, iter.max = 100)
grp <- as.factor(km$cluster)

## EXAMPLE
trait="Seed.length.width.ratio"
trait_type="quantitative"
test = "lrt"
fixed = 4
adj_method = "bonferroni"
######################################################################

#* Echo back the input
#* @param msg The message to echo
#* @get /echo
function(msg=""){
  list(msg = paste0("The message is: '", msg, "'"))
}

##### GWAS #####

#* Compute the GWAS object
#* @tag GWAS
#* @param test The testing method (lrt, Wald or score)
#* @param trait The trait to be analyzed
#* @param trait_type The trait type: quantitative or binary
#* @param fixed The option chosen for fixed effect (number of PC, or kmeans, or none)
#* @get /gwas
function(test, trait, trait_type, fixed){
  # put trait in bed matrix 
  bm.wp@ped$pheno <- pheno[bm.wp@ped$id, trait]
  # remove lines with missing trait
  bm.wom <- select.inds(bm.wp, !is.na(bm.wp@ped$pheno))
  # keep marker with a large enough MAF (>0.05) and low missing rate (callrate>0.9)
  bm.wom <- select.snps(bm.wom, maf > 0.05)
  bm.wom <- select.snps(bm.wom, callrate > 0.9)
  # Compute K (=Genetic relationship matrix for the lines with observed/non missing trait)
  K <- GRM(bm.wom)
  # Compute fixed effects from k-means clustering
  Xkm <- model.matrix(~grp[bm.wom@ped$id])
 
  if (is.numeric(fixed)){
    if (test!="score"){gwa=association.test(bm.wom, method = "lmm", response = trait_type, test = test, eigenK = eigen(K), p = fixed)
    }else{
      gwa=association.test(bm.wom, method = "lmm", response = trait_type, test = "score", K = K, eigenK = eigen(K), p = fixed)}
  }else if(fixed=="kmeans"){
    if (test!="score"){gwa=association.test(bm.wom, method = "lmm", response = trait_type, test = test, X=Xkm, eigenK = eigen(K))
    }else{
      gwa=association.test(bm.wom, method = "lmm", response = trait_type, test = "score", K = K, X=Xkm, eigenK = eigen(K))}
  }else if (is.null(fixed)){
    if (test!="score"){gwa=association.test(bm.wom, method = "lmm", response = trait_type, test = test, eigenK = eigen(K))
    }else{
      gwa=association.test(bm.wom, method = "lmm", response = trait_type, test = "score", K = K, eigenK = eigen(K))}
  }  
  gwa
  }

##### Plots #####

#* Manhattan plot
#* @tag ManhattanPlot
#* @param adj_method either bonferroni or FDR
#* @param gwa the GWAS object
#* @png
#* @get /manplot
function(gwa, adj_method){
  p.adj <- p.adjust(gwa$p, method = adj_method)
  col <- rep("black", nrow(gwa))
  col[gwa$chr %% 2 == 0] <- "gray50"
  col[p.adj < 0.05] <- "green"
  manhattan(gwa, pch = 20, col = col)
}

#* LD plot
#* @tag LDPlot
#* @param bm.wom
#* @param from
#* @param to
#* @png
#* @get /LDplot
function(bm.wom, from, to){
  ld <- LD(bm.wom, c(from, to), measure = "r2")
LD.plot(ld, snp.positions = bm.wom@snps$pos[from:to])
}

##### Tables output #####

#* Table of selected SNPs
#* @tag SNP Table
#* @get /datatable
#* @serializer htmlwidget
function(gwa, adj_method){
  p.adj <- p.adjust(gwa$p, method = adj_method)
  datatable(gwa[p.adj < 0.05, ])
}

sel=gt.score[,colnames(gt.score)%in% gwa[p.adj < 0.05, "id"]]