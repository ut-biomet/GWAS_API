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
  out <- list(
    inputParams = list(
      geno_url = geno_url,
      pheno_url = pheno_url,
      upload_url = upload_url,
      trait = trait,
      test = test,
      fixed = as.character(fixed),
      thresh_maf = as.character(thresh_maf),
      thresh_callrate = as.character(thresh_callrate)
    )
  )
  # TODO ask Shuei if it is nessesary to specify that parameters used default values (for fixed, thresh_maf, thresh_callrate)

  logger$log("call with parameters:")
  logger$log(time = FALSE, context = FALSE,
             "geno_url: ", geno_url,"\n",
             "\t pheno_url: ", pheno_url, "\n",
             "\t trait: ", trait, "\n",
             "\t test: ", test, "\n",
             "\t fixed: ", fixed, "\n",
             "\t thresh_maf: ", thresh_maf, "\n",
             "\t thresh_callrate: ", thresh_callrate)


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
      out$GWAS_API_error <- "error PUT request didn't get status code 200"
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
  "thresh_p" = list(
    desc = "threshold for filtering non-significative markers. If not specify return all the values",
    type = "number",
    required = FALSE,
    isArray = FALSE)
)

adjustedResults_handler <- function(res, gwas_url, adj_method = "bonferroni", thresh_p = NA){
  logger <- logger$new("/adjustedResults")
  out <- list(
    inputParams = list(
      gwas_url = gwas_url,
      adj_method = adj_method,
      thresh_p = as.character(thresh_p)
    )
  )

  logger$log("call with parameters:")
  logger$log(time = FALSE, context = FALSE,
             "gwas_url: ", gwas_url,"\n",
             "\t adj_method: ", adj_method, "\n",
             "\t thresh_p: ", thresh_p)


  ### CHECK PARAMETERS
  # Convert to numeric
  logger$log("Convert numeric parameters...")
  if (!is.na(as.numeric(thresh_p)) | is.na(thresh_p)) {
    thresh_p <- as.numeric(thresh_p)
  } else {
    logger$log('Error: "thresh_p" cannot be converted to numeric.')
    res$status <- 400 # bad request
    out$error <- '"thresh_p" should be a numeric value.'
    return(out)
  }
  logger$log("Convert numeric parameters DONE.")


  # Adjust p-values
  logger$log("Adjust p-values ...")
  adj_gwas <- run_resAdjustment(gwasFile = NULL,
                                gwasUrl = gwas_url,
                                adj_method = "bonferroni")
  dta <- jsonlite::fromJSON(adj_gwas$gwasAdjusted)
  logger$log("Adjust p-values DONE")

  ### RESPONSE
  logger$log("Create response ... ")
  if (is.na(thresh_p)) {
    dta <- dta
  } else {
    dta <- dta[dta$p_adj <= thresh_p, ]
  }
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
