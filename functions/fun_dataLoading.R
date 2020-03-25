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
  read.vcf(paste0("data/markers/", dataId, ".vcf.gz"), verbose = FALSE)
}

#' get phenotypic data
#'
#' @param dataId 
#' @param traits 
#'
#' @return data.frame
getPhenoData <- function(dataId, trait){
  # get data
  dta <- read.csv(paste0("data/pheno/", dataId,".csv"),
                    row.names = 1)
  dta <- as.data.frame(dta[, trait],
                       row.names = row.names(dta))
  names(dta) <- trait
  dta
}


#' get phenotypic data (type 2)
#'
#' @param dataId 
#'
#' @return data.frame
getPhenoData2 <- function(dataId){
  # get data
  dta <- read.csv(paste0("data/pheno/", dataId,".csv"),
                    row.names = 1)
  dta
}


#' get data for GWAS model
#'
#' @param markerDataId 
#' @param phenoDataId 
#' @param traits 
#'
#' @return 
loadData <- function(markerDataId, phenoDataId, trait){
  mDta <- getMarkerData(markerDataId)
  pDta <- getPhenoData(phenoDataId, trait)
  
  # select individuals with traits
  dta <- mDta[mDta@ped$id %in% rownames(pDta),]
  # reorder bed matrix like the phenotypic data
  dta@ped <- dta@ped[match(rownames(pDta), dta@ped$id),]
  
  # remove monomorphic markers
  dta <- select.snps(dta, maf > 0)
  
  # add phenotype to the bed matrix
  dta@ped$pheno <- pDta[dta@ped$id, 1]
  # remove lines with missing trait
  dta <- select.inds(dta, !is.na(dta@ped$pheno))
  
  dta
}


#' get data for GWAS model (type 2)
#'
#' @param markerDataId 
#' @param phenoDataId
#'
#' @return 
loadData2 <- function(markerDataId, phenoDataId){
  mDta <- getMarkerData(markerDataId)
  pDta <- getPhenoData2(phenoDataId)
  
  # select individuals with traits
  mDta <- select.inds(mDta, id %in% rownames(pDta))
  # reorder phenotypic data with id in bed matrix
  pDta <- pDta[mDta@ped$id,]
  
  # remove monomorphic markers
  mDta <- select.snps(mDta, maf > 0)
  
  # calculate genetic relatinoal matrix
  grm <- GRM(mDta)
  
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
	bm <- data$markerData
	bm@ped$pheno <- data$phenoData[, trait]
	K <- data$grMatrix
	
	### FILTER SAMPLES
	# remove samples with missing phenotypic values
	bm <- select.inds(bm, !is.na(pheno))
	K <- K[bm@ped$id, bm@ped$id]
	
	### FILTER SNPs
	# keep marker with a large enough MAF (>0.05) 
	# and low missing rate (callrate>0.9)
	bm <- select.snps(bm, maf > thresh.maf)
	bm <- select.snps(bm, callrate > thresh.callrate)
  
  
	### FIT MODEL
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
  
	
	return(gwa)

}
