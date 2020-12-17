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
#* @apiDescription API for GWAS models fitting and use
#* @apiVersion 0.0.1
#* @apiContact @email juliendiot@ut-biomet.org
#* @apiTag Utils Endpoints for checking the API
#* @apiTag Models Endpoints related to model managements
#* @apiTag Plots Endpoints related to plots drawn from a GWAS model
#* @apiTag Data Endpoints related to model data

cat(as.character(Sys.time()), "-",
    "Start Plumber API using 'plumber v'", as.character(packageVersion("plumber")), "\n")

# required packages
library(plumber)
library(digest) # for MD5 sum calculation
library(DT)
library(gaston) # for many functions
library(httr) # make HTTP requests
library(xml2) # manage xml format
library(rjson) # manage json format
library(R6) # R6 objects
library(manhattanly) # manhattan plot using plotly
library(redux) # manage redis
stopifnot("git2r" %in% rownames(installed.packages()))


# load API's functions
sapply(list.files("src",
                  pattern = "\\.R$", # all files ending by ".R"
                  full.names = TRUE),
       source)

# create initialization logger
initLog <- logger$new("GWAS-API-INIT")

# connect to redis
initLog$log("Connexion to redis server ...")
REDIS <<- connectRedis()
initLog$log("Connexion to redis server DONE")

# send initialization message to redis
initLog$log("GWAS_API started",
            redis = TRUE,
            status = "DONE")

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
#* @tag utils
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
#* @tag utils
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
  allFP <- sapply(apiRfiles, function(f){digest(file = f)})
  out$RfilesFingerPrint <- digest(allFP)
  logger$log("END")
  return(out)
}



##### GWAS #####
#* Fit a GWAS model (type 2)
#* @tag Models
#* @param geno_url url of the markers data file (.vcf.gz file)
#* @param pheno_url url of the phenotypic data file (csv file)
#* @param upload_url url of the put request for saving the model
#* @param trait The trait to be analyzed
#* @param test The testing method (lrt, Wald or score)
#* @param fixed The option chosen for fixed effect (number of PC, or none (0))
#* @param thresh_maf keep markers with a MAF > thresh_maf
#* @param thresh_callrate keep markers with a callrate > thresh_callrate
#* @serializer unboxedJSON
#* @post /gwas
function(res,
         geno_url,
         pheno_url,
         upload_url,
         trait,
         test,
         fixed = 0,
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
  logger$log("Convert numeric parameters...")
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


  ### GET DATA
  logger$log("Load data...")
  data <- loadData(geno_url, pheno_url)
  logger$log("Load data DONE.")

  ### GWAS
  logger$log("Generate Gwas model...")
  # calc model
  model <- gwas(data, trait, test, fixed, thresh_maf, thresh_callrate)
  modelId <- gsub("\\.", "-",
                  paste0("GWAS-model_",
                         "_", trait,"_",
                         as.numeric(callTime)))
  logger$log("Generate Gwas model DONE:")
  logger$log(time = FALSE, context = FALSE,
             "modelId: ", modelId)

  # Upload model
  logger$log("Upload model ...:")
  logger$log("Save model in tmp dir...")

  localFile <- tempfile(pattern = "GWAS-Results",
                        tmpdir = tempdir(),
                        fileext = ".json")

  writeLines(toJSON(model), con = localFile)

  logger$log("Make PUT request to AWS S3 ...")
  putResult <- PUT(url = upload_url,
                   body = upload_file(localFile, type = ""))
  if (putResult$status_code != 200) {
    logger$log("Error, PUT request's satus code is different than 200: ", putResult$status)
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

  # TODO:
  # see : https://docs.google.com/document/d/1gaTazFm_a6klD9krZPKGHc5x_D-SWL8K93M6dwsTiLE/edit
  # write models's information in a database
  # so that endpoints to check already fitted model can be created

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
    thresh_maf = as.character(thresh_maf),
    thresh_callrate = as.character(thresh_callrate),
    modelRobjectMD5 = digest(model),
    modelFileMD5 = digest(file = localFile)
  )
  logger$log("Save model information DONE")



  ### RESPONSE
  logger$log("Create response ...")
  res$status <- 201 # status for good post response
  out$message <- "Model created"
  out$modelId <- modelId
  logger$log("Create response DONE")
  logger$log("END",
             redis = TRUE,
             status = "DONE",
             action_type = "GWAS")
  out
}



##### Plots #####

#* Manhattan plot (type 2)
#* @tag Plots
#* @param modelS3Path url of the model data file (rds file)
#* @param adj_method either bonferroni or FDR
#* @param thresh_p
#* @param chr names of the chromosomes to show separated using comma. Show all chromosomes if nothing is specified.
#* @serializer htmlwidget
#* @get /manplot
function(res, modelS3Path, adj_method, thresh_p = 0.05, chr = NA){
  # # save call time.
  # callTime <- Sys.time()
  logger <- logger$new("/manplot")
  out <- list(
    inputParams = list(
      modelS3Path = modelS3Path,
      adj_method = adj_method,
      thresh_p = as.character(thresh_p),
      chr = chr
    )
  )

  logger$log("call with parameters:")
  logger$log(time = FALSE, context = FALSE,
    "modelS3Path: ", modelS3Path,"\n",
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

  # LOAD MODEL
  logger$log("load model ...")
  gwa <- loadModel(modelS3Path)
  logger$log("load model DONE.")

  ### CHECK PARAMETERS 2
  logger$log('Check "chr" parameter ...')
  if (!is.na(chr)) {
    # extract chr names
    chr <- strsplit(chr, ",")[[1]]
    if (!all(chr %in% as.character(unique(gwa$chr)))) {
      logger$log('Warning: "chr" is not valid, chromosomes\'names specified are not in the data.')
      logger$log(time = FALSE, context = FALSE,
                 '"chr" is:', chr, "\n",
                 '"gwa$chr" is: ', as.character(unique(gwa$chr)))
      logger$log(time = FALSE, context = FALSE,
                 'mismatch will be deleted')
      chr <- chr[chr %in% as.character(unique(gwa$chr))]
      if (length(chr) == 0) {
        chr <- NA
      }
    } else {
      logger$log(time = FALSE, context = FALSE,
                       '"chr" is good')
    }
  }
  logger$log('Check "chr" parameter DONE')



  # CREATE PLOT
  logger$log("Create plot ...")
  p <- manPlot(gwa = gwa,
               adj_method = adj_method,
               thresh_p = thresh_p,
               chr = chr)
  logger$log("Create plot DONE")
  logger$log("Create response ... ")
  res$status <- 200 # status for good GET response
  logger$log("Create response DONE ")
  logger$log("END",
             redis = TRUE,
             status = "DONE",
             action_type = "MANHANTTAN_PLOT")
  p
}



#* LD plot
#* @tag Plots
#* @param geno_url url of the markers data file (.vcf.gz file)
#* @param from (total number of SNP should be < 50)
#* @param to (total number of SNP should be < 50)
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


  ### GET DATA
  logger$log("Load data...")
  bm.wom <- getMarkerData(geno_url)
  logger$log("Load data DONE")

  # COMPUTE LD
  logger$log("Compute LD ...")
  ld <- LD(bm.wom, c(from, to), measure = "r2")
  logger$log("Compute LD DONE")

  logger$log("Create LD plot ...")
  p <- LD.plot(ld, snp.positions = bm.wom@snps$pos[from:to])
  logger$log("Create LD plot DONE")

  logger$log("Create response ... ")
  res$status <- 200 # status for good GET response
  logger$log("Create response DONE ")
  logger$log("END",
             redis = TRUE,
             status = "DONE",
             action_type = "LD_PLOT")
  p
}

##### Tables output #####

#* Table of selected SNPs
#* @tag Data
#* @param modelS3Path url of the model data file (rds file)
#* @param adj_method either bonferroni or FDR
#* @param thresh_p threshold for p values. If not specify return all values
#* @serializer unboxedJSON
#* @get /datatable
function(res, modelS3Path, adj_method, thresh_p = NA){
  logger <- logger$new("/datatable")
  out <- list(
    inputParams = list(
      modelS3Path = modelS3Path,
      adj_method = adj_method,
      thresh_p = as.character(thresh_p)
    )
  )

  logger$log("call with parameters:")
  logger$log(time = FALSE, context = FALSE,
    "modelS3Path: ", modelS3Path,"\n",
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


  # LOAD MODEL
  logger$log("Load model...")
  gwa <- loadModel(modelS3Path)
  logger$log("Load model DONE")

  # CREATE DATATABLE
  logger$log("Adjust p-values ...")
  p.adj <- p.adjust(gwa$p, method = adj_method)
  logger$log("Adjust p-values DONE")

  ### RESPONSE
  logger$log("Create response ... ")
  res$status <- 200 # status for good GET response
  if (is.na(thresh_p)) {
    out$data <- gwa
  } else {
    out$data <- gwa[p.adj < thresh_p, ]
  }

  # datatable(gwa[p.adj < thresh_p, ])
  logger$log("Create response DONE ")
  logger$log("END",
             redis = TRUE,
             status = "DONE",
             action_type = "GWAS_RESULTS")
  out

}
