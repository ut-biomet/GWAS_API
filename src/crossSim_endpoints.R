# Author: Julien Diot juliendiot@ut-biomet.org
# 2022 The University of Tokyo
#
# Description:
# Definition of the crossing simulation endpoints' parameters and handlers




# /crossing-sim ----
crossingSim_params <- list(
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
  "nCross" = list(
    desc = paste("Number of cross to simulate for each parent",
                 "pair defined in the crossing table."),
    type = "integer",
    required = FALSE,
    isArray = FALSE),
  "upload_url" = list(
    desc = "url of the PUT request for saving the simulated vcf file",
    type = "string",
    required = FALSE,
    isArray = FALSE)
)



crossingSim_handler <- function(res,
                                geno_url,
                                crossTable_url,
                                SNPcoord_url,
                                upload_url = NA,
                                nCross = 10) {
  logger <- Logger$new("/crossing-simulation")
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

  ### CHECK PARAMETERS
  # Convert to numeric
  nCross <- as.numeric(nCross)

  ### Crossing simulation
  logger$log("Crossing simulation ...")
  outFile <- tempfile(fileext = ".vcf.gz")
  crossingSimulation(genoUrl = geno_url,
                     crossTableUrl = crossTable_url,
                     SNPcoordUrl = SNPcoord_url,
                     nCross = nCross,
                     outFile = outFile)


  ### UPLOAD results
  if (!is.na(upload_url)) {
    logger$log("Upload simulated file ...:")
    logger$log("Save simulated vcf in tmp dir...")

    logger$log("Make PUT request ...")
    putResult <- httr::PUT(url = upload_url,
                           body = httr::upload_file(outFile, type = ""))
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
    logger$log("Upload simulated vcf DONE:")
    logger$log(time = FALSE, context = FALSE,
               "putUrl: ", upload_url)

    ### RESPONSE
    logger$log("Create response ...")
    res$status <- 201 # status for good post response
    out$message <- "Model created"
    logger$log("Create response DONE")
    logger$log("END")
    return(out)
  } else {
    logger$log("Return simulated vcf file.")
    res$status <- 201
    logger$log("Create response DONE")
    logger$log("END")
    return(include_file(outFile, res))
  }
}
