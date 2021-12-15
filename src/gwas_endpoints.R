# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Definition of the endpoints' parameters and handlers


# /gwas ----
gwas_params <- list(
  "geno_url" = list(
    desc = "url of the markers data file (.vcf.gz file)",
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "pheno_url" = list(
    desc = "url of the phenotypic data file (csv file)",
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "upload_url" = list(
    desc = "url of the PUT request for saving the model",
    type = "string",
    required = FALSE,
    isArray = FALSE),
  "trait" = list(
    desc = "name of the trait to analyze. Must be a column name of the phenotypic file.",
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "test" = list(
    desc = 'Which test to use. Either `"score"`, `"wald"` or `"lrt"`. For binary phenotypes, test = `"score"` is mandatory. For more information about this parameters see: https://www.rdocumentation.org/packages/gaston/versions/1.4.9/topics/association.test',
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "fixed" = list(
    desc = 'Number of Principal Components to include in the model with fixed effect (for test = `"wald"` or `"lrt"`). Default value is 0. For more information about this parameters see: https://www.rdocumentation.org/packages/gaston/versions/1.4.9/topics/association.test',
    type = "integer",
    required = FALSE,
    isArray = FALSE),
  "response" = list(
    desc = 'Either "quantitative" or "binary". Is the trait a quantitative or a binary phenotype? Default value is "quantitative".',
    type = "string",
    required = FALSE,
    isArray = FALSE),
  "thresh_maf" = list(
    desc = 'Threshold for filtering markers. Only markers with minor allele frequency > `thresh_maf` will be kept for the GWAS analysis. (default 0.05)',
    type = "number",
    required = FALSE,
    isArray = FALSE),
  "thresh_callrate" = list(
    desc = 'Threshold for filtering markers. Only markers with a callrate > `thresh_callrate` will be kept.',
    type = "number",
    required = FALSE,
    isArray = FALSE)
)

gwas_handler <- function(res,
                         geno_url,
                         pheno_url,
                         upload_url = NA,
                         trait,
                         test,
                         fixed = 0,
                         response = "quantitative",
                         thresh_maf = 0.05,
                         thresh_callrate = 0.9){

  logger <- logger$new("/gwas")
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
  options(warn = -1) # disable warnings. for "as.numeric" with character
  logger$log("Check parameters...")
  if (!is.na(as.numeric(fixed))) {
    fixed <- as.numeric(fixed)
  } else {
    logger$log('Error: "fixed" cannot be converted to numeric.')
    res$status <- 400 # bad request
    out$error <- '"fixed" should be a numeric value.'
    logger$log('Exit with error code 400')
    logger$log("END")
    return(out)
  }

  if (!is.na(as.numeric(thresh_maf))) {
    thresh_maf <- as.numeric(thresh_maf)
  } else {
    logger$log('Error: "thresh_maf" cannot be converted to numeric.')
    res$status <- 400 # bad request
    out$error <- '"thresh_maf" should be a numeric value.'
    logger$log('Exit with error code 400')
    logger$log("END")
    return(out)
  }

  if (!is.na(as.numeric(thresh_callrate))) {
    thresh_callrate <- as.numeric(thresh_callrate)
  } else {
    logger$log('Error: "thresh_callrate" cannot be converted to numeric.')
    res$status <- 400 # bad request
    out$error <- '"thresh_callrate" should be a numeric value.'
    logger$log('Exit with error code 400')
    logger$log("END")
    return(out)
  }
  logger$log("Convert numeric parameters DONE.")
  options(warn = 0) # enable warnings

  ### GWAS
  logger$log("Generate Gwas model...")
  gwas <- run_gwas(genoFile = NULL,
                   phenoFile = NULL,
                   genoUrl = geno_url,
                   phenoUrl = pheno_url,
                   trait = trait,
                   test = test,
                   fixed = fixed,
                   response = response,
                   thresh_maf = thresh_maf,
                   thresh_callrate = thresh_callrate)
  localFile <- gwas$file
  logger$log("Generate Gwas model DONE:")

  modelId <- gsub("\\.", "-",
                  paste0("GWAS-model_",
                         trait,"_",
                         test,"_",
                         fixed,"_",
                         response,"_",
                         thresh_maf,"_",
                         thresh_callrate,"_",
                         as.numeric(callTime)))
  logger$log(time = FALSE, context = FALSE,
             "modelId: ", modelId)



  # save model information

  logger$log("Save model information ...")
  out$modelInfo <- list(
    modelId = modelId,
    upload_url = upload_url,
    creationTime = callTime,
    geno_url = geno_url,
    pheno_url = pheno_url,
    trait = trait,
    test = test,
    fixed = as.character(fixed),
    response = response,
    thresh_maf = as.character(thresh_maf),
    thresh_callrate = as.character(thresh_callrate),
    modelRobjectMD5 = digest::digest(gwas$gwasRes)
  )
  logger$log("Save model information DONE")


  ### UPLOAD model
  if (!is.na(upload_url)) {
    logger$log("Upload model ...:")
    logger$log("Save model in tmp dir...")

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
    logger$log("Upload model DONE:")
    logger$log(time = FALSE, context = FALSE,
               "putUrl: ", upload_url)

    ### RESPONSE
    logger$log("Create response ...")
    res$status <- 201 # status for good post response
    out$message <- "Model created"
    out$modelId <- modelId
    logger$log("Create response DONE")
    logger$log("END")#,
    # redis = TRUE,
    # status = "DONE",
    # action_type = "GWAS")
    return(out)
  } else {
    logger$log("Create response ...")
    gwasList <- list(gwas = gwas$gwasRes,
                     metadata = gwas$metadata)
    out <- jsonlite::toJSON(gwasList,
                            complex = "list",
                            # pretty = T,
                            digits = NA)
    res$status <- 200
    logger$log("Create response DONE")
    logger$log("END")
    return(out)
  }
}


# /adjustedResults ----
adjustedResults_params <- list(
  "gwas_url" = list(
    desc = "url of the result file saved by `/gwas` (json file)",
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "adj_method" = list(
    desc = 'correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none". Default: "bonferroni". (see https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust for more details)',
    type = "string",
    required = FALSE,
    isArray = FALSE),
  "filter_pAdj" = list(
    desc = 'threshold to remove points with pAdj > filter_pAdj from the plot (default no filtering)',
    type = "number",
    required = FALSE,
    isArray = FALSE),
  "filter_nPoints" = list(
    desc = 'threshold to keep only the filter_nPoints with the lowest p-values for the plot.',
    type = "number",
    required = FALSE,
    isArray = FALSE),
  "filter_quant" = list(
    desc = 'threshold to keep only the filter_quant*100 % of the points with the lowest p-values for the plot (default no filtering)',
    type = "number",
    required = FALSE,
    isArray = FALSE)
)

adjustedResults_handler <- function(res,
                                    gwas_url,
                                    adj_method = "bonferroni",
                                    filter_pAdj = 1,
                                    filter_nPoints = 10^100,
                                    filter_quant = 1) {
  logger <- logger$new("/adjustedResults")
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
  logger$log("Convert numeric parameters...")
  logger$log("Convert numeric parameters DONE.")


  # Adjust p-values
  logger$log("Adjust p-values ...")
  adj_gwas <- run_resAdjustment(gwasFile = NULL,
                                gwasUrl = gwas_url,
                                adj_method = "bonferroni",
                                filter_pAdj = filter_pAdj,
                                filter_nPoints = filter_nPoints,
                                filter_quant = filter_quant)
  dta <- jsonlite::fromJSON(adj_gwas$gwasAdjusted)
  logger$log("Adjust p-values DONE")

  ### RESPONSE
  logger$log("Create response ... ")
  gwasList <- list(gwasAdjusted = dta,
                   metadata = adj_gwas$metadata)
  out <- jsonlite::toJSON(gwasList,
                          complex = "list",
                          # pretty = T,
                          digits = NA)
  res$status <- 200
  logger$log("Create response DONE ")
  logger$log("END")#,
  # redis = TRUE,
  # status = "DONE",
  # action_type = "GWAS_RESULTS")
  out
}
