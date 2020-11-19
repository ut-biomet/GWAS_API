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
gwas <- function(data, trait, test, fixed, thresh.maf, thresh.callrate) {

	### GET DATA
  cat(as.character(Sys.time()), "-",
      " r-gwas(): aggregate data ... \n")
	bm <- data$markerData
	bm@ped$pheno <- data$phenoData[, trait]
	K <- data$grMatrix
	cat(as.character(Sys.time()), "-",
	    " r-gwas(): aggregate data DONE \n")

	### FILTER SAMPLES
	# remove samples with missing phenotypic values
	cat(as.character(Sys.time()), "-",
	    " r-gwas(): remove samples with missing phenotypic values ... \n")
	bm <- select.inds(bm, !is.na(pheno))
	K <- K[bm@ped$id, bm@ped$id]
	cat(as.character(Sys.time()), "-",
	    " r-gwas(): remove samples with missing phenotypic values DONE \n")


	### FILTER SNPs
	# keep marker with a large enough MAF (>0.05)
	# and low missing rate (callrate>0.9)
	cat(as.character(Sys.time()), "-",
	    " r-gwas(): filter SNPs ... \n")
	bm <- select.snps(bm, maf > thresh.maf)
	bm <- select.snps(bm, callrate > thresh.callrate)
	cat(as.character(Sys.time()), "-",
	    " r-gwas(): filter SNPs DONE \n")

	### FIT MODEL
	cat(as.character(Sys.time()), "-",
	    " r-gwas(): fit model ... \n")
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
	cat(as.character(Sys.time()), "-",
	    " r-gwas(): fit model DONE \n")

	cat(as.character(Sys.time()), "-",
	    " r-gwas(): DONE, return output.\n")

	return(gwa)

}
