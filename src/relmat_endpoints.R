# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Definition of the relationship matrix endpoints' parameters and handlers





# /relmat-ped ----
relmatped_params <- list(
  "ped_url" = list(
    desc = "url of the pedigree file (.csv)",
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "upload_url" = list(
    desc = "url of the PUT request for saving the results",
    type = "string",
    required = FALSE,
    isArray = FALSE),
  "header" = list(
    desc = 'a logical value indicating whether the file contains the names of the variables as its first line. The default value is TRUE. In any cases, the column 1 will be interpreted as the individual id, column 2 as the first parent, column 3 as the second parent.',
    type = "string",
    required = FALSE,
    isArray = FALSE),
  "unknown_string" = list(
    desc = 'a character vector of strings which are to be interpreted as "unknown parent". By default: missing value in the file.',
    type = "string",
    required = TRUE,
    isArray = FALSE)
)

relmatped_handler <- function(res,
                              ped_url,
                              upload_url = NA,
                              header = TRUE,
                              unknown_string = ''){
  logger <- logger$new("/relmat-ped")
  # save call time.
  callTime <- Sys.time()
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


  ### relationship calculation
  logger$log("Calculate relationship ...")
  relmat <- calc_pedRelMAt(pedUrl = ped_url,
                           unknown_string = unknown_string,
                           header = header,
                           outFormat = 'json')
  localFile <- relmat$file
  logger$log("Calculate relationship DONE:")

  relMatId <- gsub("\\.", "-",
                  paste0("relationship-matrix_",
                         as.numeric(callTime)))
  logger$log(time = FALSE, context = FALSE,
             "relMatId: ", relMatId)



  # save results information
  logger$log("Save results information ...")
  out$modelInfo <- list(
    relMatId = relMatId,
    upload_url = upload_url,
    creationTime = callTime,
    ped_url = ped_url,
    resultRobjectMD5 = digest::digest(relmat$relmat)
  )
  logger$log("Save results information DONE")


  ### UPLOAD results
  if (!is.na(upload_url)) {
    logger$log("Upload results ...:")
    logger$log("Save results in tmp dir...")

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
    out$message <- "Relationship matrix created"
    out$relMatId <- relMatId
    logger$log("Create response DONE")
    logger$log("END")
    return(out)
  } else {
    logger$log("Create response ...")
    relmat <- list(relMat = as.data.frame(relmat$relMat),
                   metadata = relmat$metadata)
    out <- jsonlite::toJSON(relmat,
                            # dataframe = "rows",
                            # pretty = T,
                            digits = NA,
                            na = 'string')
    res$status <- 200
    logger$log("Create response DONE")
    logger$log("END")
    return(out)
  }
}







# /relmat-geno ----
relmatgeno_params <- list(
  "geno_url" = list(
    desc = "url of the genotype file (`.vcf` or `.vcf.gz`)",
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "upload_url" = list(
    desc = "url of the PUT request for saving the results",
    type = "string",
    required = FALSE,
    isArray = FALSE)
)

relmatgeno_handler <- function(res,
                               geno_url,
                               upload_url = NA){
  logger <- logger$new("/relmat-geno")
  # save call time.
  callTime <- Sys.time()
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


  ### relationship calculation
  logger$log("Calculate relationship ...")
  relmat <- calc_genoRelMAt(genoUrl = geno_url,
                           outFormat = 'json')
  localFile <- relmat$file
  logger$log("Calculate relationship DONE:")

  relMatId <- gsub("\\.", "-",
                  paste0("relationship-matrix_",
                         as.numeric(callTime)))
  logger$log(time = FALSE, context = FALSE,
             "relMatId: ", relMatId)



  # save results information
  logger$log("Save results information ...")
  out$modelInfo <- list(
    relMatId = relMatId,
    upload_url = upload_url,
    creationTime = callTime,
    geno_url = geno_url,
    resultRobjectMD5 = digest::digest(relmat$relmat)
  )
  logger$log("Save results information DONE")


  ### UPLOAD results
  if (!is.na(upload_url)) {
    logger$log("Upload results ...:")
    logger$log("Save results in tmp dir...")

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
    out$message <- "Relationship matrix created"
    out$relMatId <- relMatId
    logger$log("Create response DONE")
    logger$log("END")
    return(out)
  } else {
    logger$log("Create response ...")
    relmat <- list(relMat = as.data.frame(relmat$relMat),
                   metadata = relmat$metadata)
    out <- jsonlite::toJSON(relmat,
                            # dataframe = "rows",
                            # pretty = T,
                            digits = NA,
                            na = 'string')
    res$status <- 200
    logger$log("Create response DONE")
    logger$log("END")
    return(out)
  }
}





# /relmat-combined ----
relmatCombined_params <- list(
  "genoRelMat_url" = list(
    desc = "url of the genomic relationship matrix file. This url must end by either `.csv` or `.json`",
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "pedRelMat_url" = list(
    desc = "url of the pedigree relationship matrix file. This url must end by either `.csv` or `.json`",
    type = "string",
    required = FALSE,
    isArray = FALSE)
  )

relmatCombined_handler <- function(res,
                                   genoRelMat_url,
                                   pedRelMat_url,
                                   upload_url = NA){
  logger <- logger$new("/relmat-geno")
  # save call time.
  callTime <- Sys.time()
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


  ### relationship calculation
  logger$log("Calculate relationship ...")
  relmat <- calc_combinedRelMat(pedRelMatUrl = pedRelMat_url,
                                genoRelMatUrl = genoRelMat_url,
                                outFormat = 'json')
  localFile <- relmat$file
  logger$log("Calculate relationship DONE:")

  relMatId <- gsub("\\.", "-",
                  paste0("relationship-matrix_",
                         as.numeric(callTime)))
  logger$log(time = FALSE, context = FALSE,
             "relMatId: ", relMatId)



  # save results information
  logger$log("Save results information ...")
  out$modelInfo <- list(
    relMatId = relMatId,
    upload_url = upload_url,
    creationTime = callTime,
    pedRelMat_url = pedRelMat_url,
    genoRelMat_url = genoRelMat_url,
    resultRobjectMD5 = digest::digest(relmat$relmat)
  )
  logger$log("Save results information DONE")


  ### UPLOAD results
  if (!is.na(upload_url)) {
    logger$log("Upload results ...:")
    logger$log("Save results in tmp dir...")

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
    out$message <- "Relationship matrix created"
    out$relMatId <- relMatId
    logger$log("Create response DONE")
    logger$log("END")
    return(out)
  } else {
    logger$log("Create response ...")
    relmat <- list(relMat = as.data.frame(relmat$relMat),
                   metadata = relmat$metadata)
    out <- jsonlite::toJSON(relmat,
                            # dataframe = "rows",
                            # pretty = T,
                            digits = NA,
                            na = 'string')
    res$status <- 200
    logger$log("Create response DONE")
    logger$log("END")
    return(out)
  }
}
