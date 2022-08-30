# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Definition of the Progenies' blup expected values and variance calculation's
# # endpoints' parameters and handlers




# /progenyBlupCalc ----

progenyBlupCalc_params <- list(
  "geno_url" = list(
    desc = "url of the phased genotype file of the parents",
    type = "string",
    required = TRUE,
    isArray = FALSE),

  "crossTable_url" = list(
    desc = paste("url of the crossing table data file",
                 "(`csv` file of 2 or 3 columns).",
                 "It must contain the names of the variables",
                 "as its first line. The column 1 and 2 will",
                 "be interpreted as the parents ids.",
                 "The optional third column will be interpreted",
                 "as the offspring base name."),
    type = "string",
    required = TRUE,
    isArray = FALSE),

  "SNPcoord_url" = list(
    desc = paste("url of the SNPs coordinates",
                 "file (`csv` file). This `.csv` file should have",
                 "4 named columns:\n",
                 "- `chr`: chromosome name holding the SNP\n",
                 "- `physPos`: physical position of the SNP on",
                 "the chromosome\n",
                 "- `linkMapPos`: linkage map position of the SNP",
                 "on the chromosome in Morgan\n",
                 "- `SNPid`: ID of the SNP\n"),
    type = "string",
    required = TRUE,
    isArray = FALSE),

  "markerEffects_url" = list(
    desc = paste("url of the marker effects file",
                 "(`csv` file). This `.csv` file should",
                 "have 2 named columns:\\n",
                 "- `SNPid`: Marker id\\n",
                 "- `effects`: effect of the corresponding",
                 "marker\\n"),
    type = "string",
    required = TRUE,
    isArray = FALSE),

  "upload_url" = list(
    desc = "url of the PUT request for saving the results",
    type = "string",
    required = FALSE,
    isArray = FALSE)
)



progenyBlupCalc_handler <- function(res,
                                    geno_url,
                                    crossTable_url,
                                    SNPcoord_url,
                                    markerEffects_url,
                                    upload_url = NA) {
  logger <- logger$new("/progenyBlupCalc")

  # save call time.
  callTime <- Sys.time()

  # log inputs params
  inputParamsNames <- names(formals(rlang::current_fn()))
  inputParamsNames <- inputParamsNames[!inputParamsNames %in% c('res')]
  inputParams <- as.list(environment())[inputParamsNames]
  out <- list(
    inputParams = inputParams
  )
  logger$log("call with parameters:")
  logger$log(time = FALSE, context = FALSE,
             paste0(names(out$inputParams), ": ", out$inputParams,
                    collapse = '\n\t')
  )

  # calculate progenies BLUPs variance and expected values
  logger$log("Progenies blup calculation ...")
  progeniesBlups <- calc_progenyBlupEstimation(
    genoUrl = geno_url,
    crossTableUrl = crossTable_url,
    SNPcoordUrl = SNPcoord_url,
    markerEffectsUrl = markerEffects_url
  )
  logger$log("Progenies blup calculation DONE")
  progeniesBlupsId <- gsub("\\.", "-",
                  paste0("progenies-blup_",
                         as.numeric(callTime)))
  logger$log(time = FALSE, context = FALSE,
             "progeniesBlupsId: ", progeniesBlupsId)

  # save results information
  logger$log("Save results information ...")
  out$resultsInfo <- list(
    progeniesBlupsId = progeniesBlupsId,
    upload_url = upload_url,
    creationTime = callTime,
    genoUrl = geno_url,
    crossTableUrl = crossTable_url,
    SNPcoordUrl = SNPcoord_url,
    markerEffectsUrl = markerEffects_url,
    resultRobjectMD5 = digest::digest(progeniesBlups)
  )
  logger$log("Save results information DONE")



  ### UPLOAD results
  if (!is.na(upload_url)) {
    logger$log("Upload results ...:")

    logger$log("Save results in tmp dir...")
    localFile <- tempfile(fileext = '.json')
    save_dataFrame_as_json(progeniesBlups, localFile)
    logger$log("Save results in tmp dir DONE")


    logger$log("Make PUT request ...")
    putResult <- httr::PUT(url = upload_url,
                           body = httr::upload_file(localFile, type = ""))
    if (putResult$status_code != 200) {
      logger$log("Error, PUT request's satus code is different than 200: ",
                 putResult$status)
      res$status <- putResult$status_code
      out <- as_list(content(putResult))
      out$r_geno_tools_api_error <- "error PUT request didn't get status code 200"
      logger$log('Exit with error code ', putResult$status)
      logger$log("END")
      return(out)
    }
    logger$log("Upload results DONE:")
    logger$log(time = FALSE, context = FALSE,
               "putUrl: ", upload_url)

    ### RESPONSE
    logger$log("Create response ...")
    res$status <- 201 # status for good post response
    out$message <- "Progenies blup file created"
    out$progBlupId <- progBlupId
    logger$log("Create response DONE")
    logger$log("END")
    return(out)

  } else {
    logger$log("Create response ...")

    out <- jsonlite::toJSON(progeniesBlups,
                            dataframe = "rows",
                            pretty = T,
                            digits = NA,
                            na = 'string')
    res$status <- 200
    logger$log("Create response DONE")
    logger$log("END")
    return(out)
  }
}

