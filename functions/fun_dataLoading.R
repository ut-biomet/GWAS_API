# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 The University of Tokyo
#
# Description:
# Functions of the GWAS API related to data loading


#' get markers data
#'
#' @param dataId
#'
#' @return bed.matrix
getMarkerData <- function(dataId) {
  cat(as.character(Sys.time()), "-",
      " r-getMarkerData(): Read data file ... \n")
  dta <- read.vcf(paste0("data/markers/", dataId, ".vcf.gz"), verbose = FALSE)
  cat(as.character(Sys.time()), "-",
      " r-getMarkerData(): Read data file DONE \n")
  cat(as.character(Sys.time()), "-",
      " r-getMarkerData(): DONE, return output.\n")
  dta
}

#' get phenotypic data
#'
#' @param dataId
#'
#' @return data.frame
getPhenoData <- function(dataId){
  cat(as.character(Sys.time()), "-",
      " r-getPhenoData(): Read data file ... \n")
  dta <- read.csv(paste0("data/pheno/", dataId,".csv"),
                    row.names = 1)
  cat(as.character(Sys.time()), "-",
      " r-getPhenoData(): Read data file DONE \n")


  cat(as.character(Sys.time()), "-",
      " r-getPhenoData(): DONE, return output.\n")
  dta
}



#' get data for GWAS model
#'
#' @param markerDataId
#' @param phenoDataId
#'
#' @return
loadData <- function(markerDataId, phenoDataId){

  cat(as.character(Sys.time()), "-",
      " r-loadData(): get marker data ...\n")
  mDta <- getMarkerData(markerDataId)
  cat(as.character(Sys.time()), "-",
      " r-loadData(): get marker data DONE\n")

  cat(as.character(Sys.time()), "-",
      " r-loadData(): get pheno data ...\n")
  pDta <- getPhenoData(phenoDataId)
  cat(as.character(Sys.time()), "-",
      " r-loadData(): get pheno data DONE\n")


  # select individuals with traits
  cat(as.character(Sys.time()), "-",
      " r-loadData(): extract individuals having phenotype for the corresponding trait ...\n")
  mDta <- select.inds(mDta, id %in% rownames(pDta))
  cat(as.character(Sys.time()), "-",
      " r-loadData(): extract individuals having phenotype for the corresponding trait DONE\n")


  # reorder phenotypic data with id in bed matrix
  cat(as.character(Sys.time()), "-",
  " r-loadData(): reorder matrix ...\n")
  pDta <- pDta[mDta@ped$id,]
  cat(as.character(Sys.time()), "-",
      " r-loadData(): reorder matrix DONE\n")


  # remove monomorphic markers
  cat(as.character(Sys.time()), "-",
      " r-loadData(): remove monomorphic markers ...\n")
  mDta <- select.snps(mDta, maf > 0)
  cat(as.character(Sys.time()), "-",
      " r-loadData(): remove monomorphic markers DONE\n")


  # calculate genetic relatinoal matrix
  cat(as.character(Sys.time()), "-",
      " r-loadData(): calculate genetic relatinoal matrix ...\n")
  grm <- GRM(mDta)
  cat(as.character(Sys.time()), "-",
      " r-loadData(): calculate genetic relatinoal matrix DONE\n")

  cat(as.character(Sys.time()), "-",
       " r-loadData(): DONE, return output.\n")

  list(markerData = mDta, phenoData = pDta, grMatrix = grm)
}

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
