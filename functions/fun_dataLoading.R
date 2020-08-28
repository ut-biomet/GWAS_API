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
