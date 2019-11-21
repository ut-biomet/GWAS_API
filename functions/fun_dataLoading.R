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
getMarkerData <- function(dataId){
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
