#
# This is a Plumber API. In RStudio 1.2 or newer you can run the API by
# clicking the 'Run API' button above.
#
# In RStudio 1.1 or older, see the Plumber documentation for details
# on running the API.
#
# Find out more about building APIs with Plumber here:
#
#    https://www.rplumber.io/
#


##### API description #####
#* @apiTitle GWAS API
#* @apiDescription REST API for GWAS analysis
#* @apiVersion 0.0.1
#* @apiContact @email juliendiot@ut-biomet.org
#* @apiTag Utils Endpoints for checking the API
#* @apiTag GWAS Endpoints related to gwas analysis
#* @apiTag Plots Endpoints related to plots drawn from a GWAS model

cat(as.character(Sys.time()), "-",
    "Start Plumber API using 'plumber v'", as.character(packageVersion("plumber")), "\n")

# required packages
library(plumber)
stopifnot("digest" %in% rownames(installed.packages()))  # for MD5 sum calculation
# library(DT)
# library(gaston) # for many functions
# library(httr) # make HTTP requests
# library(xml2) # manage xml format
# library(rjson) # manage json format
# library(R6) # R6 objects
# library(manhattanly) # manhattan plot using plotly
# library(redux) # manage redis
# stopifnot("git2r" %in% rownames(installed.packages()))


# load API's functions
sapply(list.files("GWAS-Engine/src/",
                  pattern = "\\.R$", # all files ending by ".R"
                  full.names = TRUE),
       source)

# create initialization logger
initLog <- logger$new("GWAS-API-INIT")

# # connect to redis
# initLog$log("Connexion to redis server ...")
# REDIS <<- connectRedis()
# initLog$log("Connexion to redis server DONE")

# send initialization message to redis
initLog$log("GWAS_API started")#,
            # redis = TRUE,
            # status = "DONE")

rm("initLog")
##################################### Filter ###################################

#* Log some information about the incoming requests
#* @filter logger
function(req){
  logger <- logger$new("GWAS-API-REQUESTS")
  logger$log(req$REQUEST_METHOD, req$PATH_INFO, "-",
             req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR)
  plumber::forward()
}




################################### Endpoints ##################################

##### Utils #####

#* Echo back the input
#* @tag Utils
#* @param msg The message to echo
#* @serializer unboxedJSON
#* @get /echo
function(msg=""){
  logger <- logger$new("/echo")
  logger$log("call with parameters:")
  logger$log("msg: ", msg,
             time = FALSE, context = FALSE)
  list(msg = paste0("The message is: '", msg, "'"))
}


#* give information about current API version
#* @tag Utils
#* @serializer unboxedJSON
#* @get /version
function(){
  logger <- logger$new("/version")
  logger$log("call with no parameters")
  out <- list()

  if (dir.exists(".git")) {
    logger$log("extract last commit informations")
    out$lastCommit <- as.data.frame(git2r::last_commit(), "data.frame")
  } else {
    logger$log("git repository not found")
    out$lastCommit <- NA
  }

  logger$log("calculate MD5 sum of all the '.R' files")
  apiRfiles <- dir(all.files = T, pattern = ".R$", recursive = T)
  allFP <- sapply(apiRfiles, function(f){digest::digest(file = f)})
  out$RfilesFingerPrint <- digest::digest(allFP)
  logger$log("END")
  return(out)
}



##### GWAS #####
#* Fit a GWAS model. This endpoint take Urls of geno and pheno data (and values of other GWAS parameters) and write an a json file to the give Url using a PUT request. It had been disign to work with amazon S3 services.
#* @tag GWAS
#* @param geno_url url of the markers data file (.vcf.gz file)
#* @param pheno_url url of the phenotypic data file (csv file)
#* @param upload_url url of the PUT request for saving the model
#* @param trait The trait name of the trait to analyzed, should be a column name of the phenotypic data.
#* @param test Which test to use. Either `"score"`,  `"wald"` or `"lrt"`. For binary phenotypes, test = `"score"` is mandatory. For more information about this parameters see: documentation of the function `association.test` of the R package `gaston`
#* @param fixed Number of Principal Components to include in the model with fixed effect (for test = `"wald"` or `"lrt"`). Default value is 0. For more information about this parameters see: documentation of the function `association.test` of the R package `gaston`.
#* @param response Character of length 1, Either "quantitative" or "binary". Is the trait a quantitative or a binary phenotype? Default value is "quantitative".
#* @param thresh_maf Threshold for filtering markers. Only markers with minor allele frequency > `thresh_maf` will be kept.
#* @param thresh_callrate Threshold for filtering markers. Only markers with a callrate > `thresh_callrate` will be kept.
#* @serializer unboxedJSON
#* @post /gwas
function(res,
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


##### adjusted Results #####

#* Adjusted results. This endpoint calculate the adjusted p-values of the gwas analysis and return all the results or only the significant adjusted p-value. The results are return in json format.
#* @tag GWAS
#* @param gwas_url url of the result file saved by `/gwas` (json file)
#* @param adj_method correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none" (see R function `p.adjust`'s documentation for more details)
#* @param thresh_p threshold for filtering non-significative markers. If not specify return all the values
#* @serializer unboxedJSON
#* @get /adjustedResults
function(res, gwas_url, adj_method, thresh_p = NA){
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

##### Plots #####

#* Draw a Manhattan plot. This endpoint return the html code of a plotly interactive graph.
#* @tag Plots
#* @param gwas_url url of the result file saved by `/gwas` (json file)
#* @param adj_method correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none" (see R function `p.adjust`'s documentation for more details)
#* @param thresh_p
#* @param chr names of the chromosomes to show separated using comma. Show all chromosomes if nothing is specified.
#* @serializer htmlwidget
#* @get /manplot
function(res, gwas_url, adj_method, thresh_p = 0.05, chr = NA){
  # # save call time.
  # callTime <- Sys.time()
  logger <- logger$new("/manplot")
  out <- list(
    inputParams = list(
      gwas_url = gwas_url,
      adj_method = adj_method,
      thresh_p = as.character(thresh_p),
      chr = chr
    )
  )

  logger$log("call with parameters:")
  logger$log(time = FALSE, context = FALSE,
    "gwas_url: ", gwas_url,"\n",
    "\t adj_method: ", adj_method, "\n",
    "\t thresh_p: ", thresh_p, "\n",
    "\t chr: ", chr)


  ### CHECK PARAMETERS
  # Convert to numeric
  logger$log("Convert numeric parameters...")
  if (!is.na(as.numeric(thresh_p))) {
    thresh_p <- as.numeric(thresh_p)
  } else {
    logger$log('Error: "thresh_p" cannot be converted to numeric.')
    res$status <- 400 # bad request
    out$error <- '"thresh_p" should be a numeric value.'
    logger$log('Exit with error code 400')
    logger$log("END")
    return(out)
  }
  logger$log("Convert numeric parameters DONE.")

  # CREATE PLOT
  logger$log("Create plot ...")
  p <- draw_manhattanPlot(gwasFile = NULL,
                          gwasUrl = gwas_url,
                          adj_method = adj_method,
                          thresh_p = thresh_p,
                          chr = chr)
  logger$log("Create plot DONE")

  # RESPONSE
  logger$log("Create response ... ")
  res$status <- 200 # status for good GET response
  logger$log("Create response DONE ")
  logger$log("END")#,
             # redis = TRUE,
             # status = "DONE",
             # action_type = "MANHANTTAN_PLOT")
  p
}



#* Draw a LD plot. This endpoint return a png image.
#* @tag Plots
#* @param geno_url url of the markers data file (.vcf.gz file)
#* @param from lower bound of the range of SNPs for which the LD is computed. (`from` must be lower than `to`)
#* @param to upper bound of the range of SNPs for which the LD is computed (the total number of SNP should be lower than 50)
#* @serializer png
#* @get /LDplot
function(res, geno_url, from, to){
  logger <- logger$new("/LDplot")
  out <- list(
    inputParams = list(
      geno_url = geno_url,
      from = from,
      to = to
    )
  )

  logger$log("call with parameters:")
  logger$log(time=FALSE, context=FALSE,
    "geno_url: ", geno_url,"\n",
    "\t from: ", from, "\n",
    "\t to: ", to)

  ### CHECK PARAMETERS
  # Convert to numeric
  logger$log("Convert numeric parameters...")
  if (!is.na(as.numeric(from))) {
    from <- as.numeric(from)
  } else {
    logger$log('Error: "from" cannot be converted to numeric.')
    res$status <- 400 # bad request
    out$error <- '"from" should be a numeric value.'
    logger$log('Exit with error code 400')
    logger$log("END")
    return(out)
  }

  logger$log("Convert numeric parameters...")
  if (!is.na(as.numeric(to))) {
    to <- as.numeric(to)
  } else {
    logger$log('Error: "to" cannot be converted to numeric.')
    res$status <- 400 # bad request
    out$error <- '"to" should be a numeric value.'
    logger$log('Exit with error code 400')
    logger$log("END")
    return(out)
  }
  logger$log("Convert numeric parameters DONE.")


  logger$log('Check "from" < "to"...')
  if (from >= to) {
    logger$log('Error: "from" greater than "to".')
    res$status <- 400 # bad request
    out$error <- '"from" should be inferior than "to".'
    logger$log('Exit with error code 400')
    logger$log("END")
    return(out)
  }
  logger$log('Check "from" < "to" DONE')


  logger$log('Check number of SNP < 50...')
  if (to - from > 50) {
    logger$log('Error: number of SNP is > 50.')
    res$status <- 400 # bad request
    out$error <- 'number of SNP should be < 50'
    logger$log('Exit with error code 400')
    logger$log("END")
    return(out)
  }
  logger$log('Check number of SNP < 50 DONE')

  logger$log('Draw LD plot ...')
  draw_ldPlot(genoFile = NULL,
              genoUrl = geno_url,
              from = from,
              to = to,
              outFile = NULL)
  logger$log('Draw LD plot DONE')

  logger$log("Create response ... ")
  res$status <- 200 # status for good GET response
  logger$log("Create response DONE ")
  logger$log("END")#,
             # redis = TRUE,
             # status = "DONE",
             # action_type = "LD_PLOT")
  # p
}
