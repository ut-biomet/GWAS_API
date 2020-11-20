# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 The University of Tokyo
#
# Description:
# Functions of the GWAS API related to data loading


#' get markers data
#'
#' @param dtaS3Path url of the markers data file (.vcf.gz file)
#'
#' @return bed.matrix
getMarkerData <- function(dtaS3Path) {
  logger <- logger$new("r-getMarkerData()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "markers",
                        tmpdir = tempdir(),
                        fileext = ".vcf.gz")
  logger$log("Download markers file ...")
  download.file(dtaS3Path, localFile)

  logger$log("Read markers file ... ")
  dta <- read.vcf(localFile, verbose = FALSE)
  logger$log("Read markers file DONE ")
  logger$log("DONE, return output.")
  dta
}

#' get phenotypic data
#'
#' @param dtaS3Path url of the phenotypic data file (csv file)
#'
#' @return data.frame
getPhenoData <- function(dtaS3Path){
  logger <- logger$new("r-getPhenoData()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "pheno",
                        tmpdir = tempdir(),
                        fileext = ".csv")
  logger$log("Download phenotypic file ... ")
  download.file(dtaS3Path, localFile)

  logger$log("Read phenotypic file ... ")
  dta <- read.csv(localFile,
                  row.names = 1)
  logger$log("Read phenotypic file DONE ")


  logger$log("DONE, return output.")
  dta
}



#' get data for GWAS model
#'
#' @param markerDataId
#' @param phenoDataId
#'
#' @return
loadData <- function(markerS3Path, phenoS3Path){
  logger <- logger$new("r-loadData()")
  logger$log("get marker data ...")
  mDta <- getMarkerData(markerS3Path)
  logger$log("get marker data DONE")

  logger$log("get pheno data ...")
  pDta <- getPhenoData(phenoS3Path)
  logger$log("get pheno data DONE")


  # Remove from marker data individuals that are not in phenotypic data-set
  logger$log("Remove from marker data individuals that are not in phenotypic data-set ...")
  mDta <- select.inds(mDta, id %in% rownames(pDta))
  logger$log("Remove from marker data individuals that are not in phenotypic data-set DONE")


  # reorder phenotypic data with id in bed matrix
  logger$log("reorder matrix ...")
  pDta <- pDta[mDta@ped$id,]
  logger$log("reorder matrix DONE")


  # remove monomorphic markers
  logger$log("remove monomorphic markers ...")
  mDta <- select.snps(mDta, maf > 0)
  logger$log("remove monomorphic markers DONE")


  # calculate genetic relatinoal matrix
  logger$log("calculate genetic relatinoal matrix ...")
  grm <- GRM(mDta)
  logger$log("calculate genetic relatinoal matrix DONE")

  logger$log("DONE, return output.")

  list(markerData = mDta, phenoData = pDta, grMatrix = grm)
}


#' load a gwas model
#'
#' @param modelS3Path url of the model data file (rds file)
#' @return
loadModel <- function(modelS3Path){
  logger <- logger$new("r-loadModel()")
  logger$log("Create local temp file ... ")
  localFile <- tempfile(pattern = "downloadedModel",
                        tmpdir = tempdir(),
                        fileext = ".json")
  logger$log("Download model file ... ")
  download.file(modelS3Path, localFile)

  logger$log("Read model file ... ")
  gwa <- readLines(localFile)
  logger$log("Read model file DONE ")
  logger$log("Convert Json to data.frame ... ")
  gwa <- data.frame(fromJSON(gwa))
  logger$log("Convert Json to data.frame DONE ")

  logger$log("DONE, return output.")
  gwa
}
