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
  cat(as.character(Sys.time()), "-",
      " r-getMarkerData(): Create local temp file ... \n")
  localFile <- tempfile(pattern = "markers",
                        tmpdir = tempdir(),
                        fileext = ".vcf.gz")
  cat(as.character(Sys.time()), "-",
      " r-getMarkerData(): Download markers file ... \n")
  download.file(dtaS3Path, localFile)

  cat(as.character(Sys.time()), "-",
      " r-getMarkerData(): Read markers file ... \n")
  dta <- read.vcf(localFile, verbose = FALSE)
  cat(as.character(Sys.time()), "-",
      " r-getMarkerData(): Read markers file DONE \n")
  cat(as.character(Sys.time()), "-",
      " r-getMarkerData(): DONE, return output.\n")
  dta
}

#' get phenotypic data
#'
#' @param dtaS3Path url of the phenotypic data file (csv file)
#'
#' @return data.frame
getPhenoData <- function(dtaS3Path){
  cat(as.character(Sys.time()), "-",
      " r-getPhenoData(): Create local temp file ... \n")
  localFile <- tempfile(pattern = "pheno",
                        tmpdir = tempdir(),
                        fileext = ".csv")
  cat(as.character(Sys.time()), "-",
      " r-getPhenoData(): Download phenotypic file ... \n")
  download.file(dtaS3Path, localFile)

  cat(as.character(Sys.time()), "-",
      " r-getPhenoData(): Read phenotypic file ... \n")
  dta <- read.csv(localFile,
                  row.names = 1)
  cat(as.character(Sys.time()), "-",
      " r-getPhenoData(): Read phenotypic file DONE \n")


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
loadData <- function(markerS3Path, phenoS3Path){

  cat(as.character(Sys.time()), "-",
      " r-loadData(): get marker data ...\n")
  mDta <- getMarkerData(markerS3Path)
  cat(as.character(Sys.time()), "-",
      " r-loadData(): get marker data DONE\n")

  cat(as.character(Sys.time()), "-",
      " r-loadData(): get pheno data ...\n")
  pDta <- getPhenoData(phenoS3Path)
  cat(as.character(Sys.time()), "-",
      " r-loadData(): get pheno data DONE\n")


  # Remove from marker data individuals that are not in phenotypic data-set
  cat(as.character(Sys.time()), "-",
      " r-loadData(): Remove from marker data individuals that are not in phenotypic data-set ...\n")
  mDta <- select.inds(mDta, id %in% rownames(pDta))
  cat(as.character(Sys.time()), "-",
      " r-loadData(): Remove from marker data individuals that are not in phenotypic data-set DONE\n")


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


#' load a gwas model
#'
#' @param modelS3Path url of the model data file (rds file)
#' @return
loadModel <- function(modelS3Path){
  cat(as.character(Sys.time()), "-",
      " r-loadModel(): Create local temp file ... \n")
  localFile <- tempfile(pattern = "downloadedModel",
                        tmpdir = tempdir(),
                        fileext = ".json")
  cat(as.character(Sys.time()), "-",
      " r-loadModel(): Download model file ... \n")
  download.file(modelS3Path, localFile)

  cat(as.character(Sys.time()), "-",
      " r-loadModel(): Read model file ... \n")
  gwa <- readLines(localFile)
  cat(as.character(Sys.time()), "-",
      " r-loadModel(): Read model file DONE \n")
  cat(as.character(Sys.time()), "-",
      " r-loadModel(): Convert Json to data.frame ... \n")
  gwa <- data.frame(fromJSON(gwa))
  cat(as.character(Sys.time()), "-",
      " r-loadModel(): Convert Json to data.frame DONE \n")

  cat(as.character(Sys.time()), "-",
      " r-loadModel(): DONE, return output.\n")
  gwa
}
