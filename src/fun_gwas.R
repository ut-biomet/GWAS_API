# Author: Julien Diot juliendiot@ut-biomet.org
# 2020 The University of Tokyo
#
# Description:
# Functions of the GWAS API related to data loading


#' perform GWAS
#'
#' @param data (markerData, phenoData, grMatrix)
#' @param trait
#' @param maf (0 < maf < 0.5)
#' @param callrate (0 < callrate <= 1)
#'
gwas <- function(data, trait, test, fixed, thresh_maf, thresh_callrate) {
  logger <- logger$new("r-gwas()")
	### GET DATA
  logger$log("aggregate data ...")
	bm <- data$markerData
	bm@ped$pheno <- data$phenoData[, trait]
	K <- data$grMatrix
	logger$log("aggregate data DONE")

	### FILTER SAMPLES
	# remove samples with missing phenotypic values
	logger$log("remove samples with missing phenotypic values ...")
	bm <- select.inds(bm, !is.na(pheno))
	K <- K[bm@ped$id, bm@ped$id]
	logger$log("remove samples with missing phenotypic values DONE")


	### FILTER SNPs
	# keep marker with a large enough MAF (>0.05)
	# and low missing rate (callrate>0.9)
	logger$log("filter SNPs ...")
	bm <- select.snps(bm, maf > thresh_maf)
	bm <- select.snps(bm, callrate > thresh_callrate)
	logger$log("filter SNPs DONE")

	### FIT MODEL
	logger$log("fit model ...")
   	if (test != "score") {
    	gwa <- association.test(
        	bm,
        	method = "lmm",
        	response = "quantitative",
        	test = test,
        	eigenK = eigen(K),
        	p = fixed)
    } else {
    	gwa <- association.test(
        	bm,
        	method = "lmm",
        	response = "quantitative",
        	test = "score",
        	K = K,
        	p = fixed)
    }
	logger$log("fit model DONE")

	logger$log("DONE, return output.")

	return(gwa)

}
