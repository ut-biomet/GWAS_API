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


# required packages
library(plumber)
library(digest) # for MD5 sum calculation
library(DT)
library(gaston) # for many functions
library(httr) # make HTTP requests
library(xml2) # manage xml format
stopifnot("git2r" %in% rownames(installed.packages()))

# load API's functions
sapply(list.files("functions",
                  pattern = "\\.R$", # all files ending by ".R"
                  full.names = TRUE),
       source)

# create models folder
cat(as.character(Sys.time()), "-",
    "API initialization: check 'data/models' folder\n")
if (!dir.exists("data/models")) {
  perm = "0774"
  cat(as.character(Sys.time()), "-",
      "API initialization: 'data/models' not found, create it with permissions", perm, "\n")

  dir.create("data/models", mode = perm)
}



##################################### Filter ###################################

#* Log some information about the incoming requests
#* @filter logger
function(req){
  cat(as.character(Sys.time()), "-",
      req$REQUEST_METHOD, req$PATH_INFO, "-",
      req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR, "\n")
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
  cat(as.character(Sys.time()), "-",
      "/echo: call with parameters parameters:\n")
  cat("\t msg: ", msg,"\n")
  list(msg = paste0("The message is: '", msg, "'"))
}


#* give information about current API version
#* @tag utils
#* @serializer unboxedJSON
#* @get /version
function(){
  cat(as.character(Sys.time()), "-",
      "/version: call\n")
  out <- list()

  if (dir.exists(".git")) {
    cat(as.character(Sys.time()), "-",
        "/version: extract last commit informations\n")
    out$lastCommit <- as.data.frame(git2r::last_commit(), "data.frame")
  } else {
    cat(as.character(Sys.time()), "-",
        "/version: git repository not found\n")
    out$lastCommit <- NA
  }

  cat(as.character(Sys.time()), "-",
      "/version: calculate MD5 sum of all the '.R' files\n")
  apiRfiles <- dir(all.files = T, pattern = ".R$", recursive = T)
  allFP <- sapply(apiRfiles, function(f){digest(file = f)})
  out$RfilesFingerPrint <- digest(allFP)
  cat(as.character(Sys.time()), "-",
      "/version: END \n")
  return(out)
}



##### GWAS #####
#* Fit a GWAS model (type 2)
#* @tag Models
#* @param markerS3Path url of the markers data file (.vcf.gz file)
#* @param phenoS3Path url of the phenotypic data file (csv file)
#* @param modelS3Path url of the put request for saving the model
#* @param trait The trait to be analyzed
#* @param test The testing method (lrt, Wald or score)
#* @param fixed The option chosen for fixed effect (number of PC, or none (0))
#* @param tresh.maf keep markers with a MAF > tresh.maf
#* @param tresh.callrate keep markers with a callrate > tresh.callrate
#* @serializer unboxedJSON
#* @post /gwas
function(res,
         markerS3Path,
         phenoS3Path,
         modelS3Path,
         trait,
         test,
         fixed = 0,
         tresh.maf = 0.05,
         tresh.callrate = 0.9){
  # save call time.
  callTime <- Sys.time()
  out <- list(
    inputParams = list(
      markerS3Path = markerS3Path,
      phenoS3Path = phenoS3Path,
      modelS3Path = modelS3Path,
      trait = trait,
      test = test,
      fixed = as.character(fixed),
      tresh.maf = as.character(tresh.maf),
      tresh.callrate = as.character(tresh.callrate)
    )
  )
  # TODO ask Shuei if it is nessesary to specify that parameters used default values (for fixed, tresh.maf, tresh.callrate)

  cat(as.character(Sys.time()), "-",
      "/gwas: call with parameters parameters:\n")
  cat(
    "\t markerS3Path: ", markerS3Path,"\n",
    "\t phenoS3Path: ", phenoS3Path, "\n",
    "\t trait: ", trait, "\n",
    "\t test: ", test, "\n",
    "\t fixed: ", fixed, "\n",
    "\t tresh.maf: ", tresh.maf, "\n",
    "\t tresh.callrate: ", tresh.callrate, "\n"
  )


  ### CHECK PARAMETERS
  # Convert to numeric
  options(warn = -1) # disable warnings. for "as.numeric" with character
  cat(as.character(Sys.time()), "-",
      "/gwas: Convert numeric parameters...\n")
  if (!is.na(as.numeric(fixed))) {
    fixed <- as.numeric(fixed)
  } else {
    cat(as.character(Sys.time()), "-",
        '/gwas: Error: "fixed" cannot be converted to numeric.\n')
    res$status <- 400 # bad request
    out$error <- '"fixed" should be a numeric value.'
    cat(as.character(Sys.time()), "-",
        '/gwas: Exit with error code 400\n')
    cat(as.character(Sys.time()), "-",
        "/gwas: END \n")
    return(out)
  }

  if (!is.na(as.numeric(tresh.maf))) {
    tresh.maf <- as.numeric(tresh.maf)
  } else {
    cat(as.character(Sys.time()), "-",
        '/gwas: Error: "tresh.maf" cannot be converted to numeric.\n')
    res$status <- 400 # bad request
    out$error <- '"tresh.maf" should be a numeric value.'
    cat(as.character(Sys.time()), "-",
        '/gwas: Exit with error code 400\n')
    cat(as.character(Sys.time()), "-",
        "/gwas: END \n")
    return(out)
  }

  if (!is.na(as.numeric(tresh.callrate))) {
    tresh.callrate <- as.numeric(tresh.callrate)
  } else {
    cat(as.character(Sys.time()), "-",
        '/gwas: Error: "tresh.callrate" cannot be converted to numeric.\n')
    res$status <- 400 # bad request
    out$error <- '"tresh.callrate" should be a numeric value.'
    cat(as.character(Sys.time()), "-",
        '/gwas: Exit with error code 400\n')
    cat(as.character(Sys.time()), "-",
        "/gwas: END \n")
    return(out)
  }
  cat(as.character(Sys.time()), "-",
      "/gwas: Convert numeric parameters DONE.\n")
  options(warn = 0) # enable warnings


  ### GET DATA
  cat(as.character(Sys.time()), "-",
      "/gwas: Load data...\n")
  data <- loadData(markerS3Path, phenoS3Path)
  cat(as.character(Sys.time()), "-",
      "/gwas: Load data DONE.\n")

  ### GWAS
  cat(as.character(Sys.time()), "-",
      "/gwas: Generate Gwas model...\n")
  # calc model
  model <- gwas(data, trait, test, fixed, tresh.maf, tresh.callrate)
  modelId <- gsub("\\.", "-",
                  paste0("GWAS-model_",
                         "_", trait,"_",
                         as.numeric(callTime)))
  cat(as.character(Sys.time()), "-",
      "/gwas: Generate Gwas model DONE:\n")
  cat(paste0("\t modelId: ", modelId,"\n"))

  # Upload model
  cat(as.character(Sys.time()), "-",
  "/gwas: Upload model ...: \n")
  cat(as.character(Sys.time()), "-",
      "\t Save model in tmp dir... \n")
  localFile <- tempfile(pattern = "GWAS-Model",
                        tmpdir = tempdir(),
                        fileext = ".rds")
  saveRDS(model, file = localFile)
  cat(as.character(Sys.time()), "-",
      "\t Make PUT request to AWS S3 ... \n")
  putResult <- PUT(url = modelS3Path,
                   body = upload_file(localFile, type = ""))

  if (putResult$status_code != 200) {
    cat(as.character(Sys.time()), "-",
        "/gwas: Error, PUT request's satus code is different than 200: ", putResult$status,"\n")
    res$status <- putResult$status_code
    out <- as_list(content(putResult))
    out$GWAS_API_error <- "error PUT request didn't get status code 200"
    cat(as.character(Sys.time()), "-",
        '/gwas: Exit with error code ', putResult$status, '\n')
    cat(as.character(Sys.time()), "-",
        "/gwas: END \n")
    return(out)
  }
  cat(as.character(Sys.time()), "-",
      "/gwas: Upload model DONE: \n")
  cat(paste0("\t putUrl: ", modelS3Path,"\n"))

  # TODO:
  # see : https://docs.google.com/document/d/1gaTazFm_a6klD9krZPKGHc5x_D-SWL8K93M6dwsTiLE/edit
  # write models's information in a database
  # so that endpoints to check already fitted model can be created

  # save model information
  cat(as.character(Sys.time()), "-",
      "/gwas: Save model information ... \n")
  out$modelInfo <- list(
    modelId = modelId,
    putUrl = modelS3Path,
    creationTime = callTime,
    markerS3Path = markerS3Path,
    phenoS3Path = phenoS3Path,
    trait = trait,
    test = test,
    fixed = as.character(fixed),
    tresh.maf = as.character(tresh.maf),
    tresh.callrate = as.character(tresh.callrate),
    modelRobjectMD5 = digest(model),
    modelFileMD5 = digest(file = localFile)
  )
  cat(as.character(Sys.time()), "-",
      "/gwas: Save model information DONE \n")



  ### RESPONSE
  cat(as.character(Sys.time()), "-",
      "/gwas: Create response ... \n")
  res$status <- 201 # status for good post response
  out$message <- "Model created"
  out$modelId <- modelId
  cat(as.character(Sys.time()), "-",
      "/gwas: Create response DONE \n")
  cat(as.character(Sys.time()), "-",
      "/gwas: END \n")
  out
}



##### Plots #####

#* Manhattan plot (type 2)
#* @tag Plots
#* @param modelS3Path url of the model data file (rds file)
#* @param adj_method either bonferroni or FDR
#* @param thresh.p
#* @png
#* @get /manplot
function(res, modelS3Path, adj_method, thresh.p = 0.05){
  # # save call time.
  # callTime <- Sys.time()

  out <- list(
    inputParams = list(
      modelS3Path = modelS3Path,
      adj_method = adj_method,
      thresh.p = as.character(thresh.p)
    )
  )

  cat(as.character(Sys.time()), "-",
      "/manplot: call with parameters parameters:\n")
  cat(
    "\t modelS3Path: ", modelS3Path,"\n",
    "\t adj_method: ", adj_method, "\n",
    "\t thresh.p: ", thresh.p, "\n"
  )


  ### CHECK PARAMETERS
  # Convert to numeric
  cat(as.character(Sys.time()), "-",
      "/manplot: Convert numeric parameters...\n")
  if (!is.na(as.numeric(thresh.p))) {
    thresh.p <- as.numeric(thresh.p)
  } else {
    cat(as.character(Sys.time()), "-",
        '/manplot: Error: "thresh.p" cannot be converted to numeric.\n')
    res$status <- 400 # bad request
    out$error <- '"thresh.p" should be a numeric value.'
    cat(as.character(Sys.time()), "-",
        '/manplot: Exit with error code 400\n')
    cat(as.character(Sys.time()), "-",
        "/manplot: END \n")
    return(out)
  }
  cat(as.character(Sys.time()), "-",
      "/manplot: Convert numeric parameters DONE.\n")

  # LOAD MODEL
  cat(as.character(Sys.time()), "-",
      "/manplot: load model ...\n")
  gwa <- loadModel(modelS3Path)
  cat(as.character(Sys.time()), "-",
      "/manplot: load model DONE.\n")

  # CREATE PLOT
  cat(as.character(Sys.time()), "-",
      "/manplot: Adjust p-values ...\n")
  p.adj <- p.adjust(gwa$p, method = adj_method)
  cat(as.character(Sys.time()), "-",
      "/manplot: Adjust p-values DONE\n")

  col <- rep("black", nrow(gwa))
  col[gwa$chr %% 2 == 0] <- "gray50"
  col[p.adj < thresh.p] <- "green"
  cat(as.character(Sys.time()), "-",
      "/manplot: Create plot ...\n")
  p <- manhattan(gwa, pch = 20, col = col,
                 main = "TO DO trait name", # extract trait name from database
            sub = modelS3Path)
  cat(as.character(Sys.time()), "-",
      "/manplot: Create plot DONE\n")
  cat(as.character(Sys.time()), "-",
      "/manplot: Create response ... \n")
  res$status <- 200 # status for good GET response
  cat(as.character(Sys.time()), "-",
      "/manplot: Create response DONE \n")
  cat(as.character(Sys.time()), "-",
      "/manplot: END \n")
  p
}



#* LD plot
#* @tag Plots
#* @param markerS3Path url of the markers data file (.vcf.gz file)
#* @param from (total number of SNP should be < 50)
#* @param to (total number of SNP should be < 50)
#* @png
#* @get /LDplot
function(res, markerS3Path, from, to){

  out <- list(
    inputParams = list(
      markerS3Path = markerS3Path,
      from = from,
      to = to
    )
  )

  cat(as.character(Sys.time()), "-",
      "/LDplot: call with parameters parameters:\n")
  cat(
    "\t markerS3Path: ", markerS3Path,"\n",
    "\t from: ", from, "\n",
    "\t to: ", to, "\n"
  )

  ### CHECK PARAMETERS
  # Convert to numeric
  cat(as.character(Sys.time()), "-",
      "/LDplot: Convert numeric parameters...\n")
  if (!is.na(as.numeric(from))) {
    from <- as.numeric(from)
  } else {
    cat(as.character(Sys.time()), "-",
        '/LDplot: Error: "from" cannot be converted to numeric.\n')
    res$status <- 400 # bad request
    out$error <- '"from" should be a numeric value.'
    cat(as.character(Sys.time()), "-",
        '/LDplot: Exit with error code 400\n')
    cat(as.character(Sys.time()), "-",
        "/LDplot: END \n")
    return(out)
  }

  cat(as.character(Sys.time()), "-",
      "/LDplot: Convert numeric parameters...\n")
  if (!is.na(as.numeric(to))) {
    to <- as.numeric(to)
  } else {
    cat(as.character(Sys.time()), "-",
        '/LDplot: Error: "to" cannot be converted to numeric.\n')
    res$status <- 400 # bad request
    out$error <- '"to" should be a numeric value.'
    cat(as.character(Sys.time()), "-",
        '/LDplot: Exit with error code 400\n')
    cat(as.character(Sys.time()), "-",
        "/LDplot: END \n")
    return(out)
  }
  cat(as.character(Sys.time()), "-",
      "/LDplot: Convert numeric parameters DONE.\n")


  cat(as.character(Sys.time()), "-",
      '/LDplot: Check "from" < "to"...\n')
  if (from >= to) {
    cat(as.character(Sys.time()), "-",
        '/LDplot: Error: "from" greater than "to".\n')
    res$status <- 400 # bad request
    out$error <- '"from" should be inferior than "to".'
    cat(as.character(Sys.time()), "-",
        '/LDplot: Exit with error code 400\n')
    cat(as.character(Sys.time()), "-",
        "/LDplot: END \n")
    return(out)
  }
  cat(as.character(Sys.time()), "-",
      '/LDplot: Check "from" < "to" DONE\n')


  cat(as.character(Sys.time()), "-",
      '/LDplot: Check number of SNP < 50...\n')
  if (to - from > 50) {
    cat(as.character(Sys.time()), "-",
        '/LDplot: Error: number of SNP is > 50.\n')
    res$status <- 400 # bad request
    out$error <- 'number of SNP should be < 50'
    cat(as.character(Sys.time()), "-",
        '/LDplot: Exit with error code 400\n')
    cat(as.character(Sys.time()), "-",
        "/LDplot: END \n")
    return(out)
  }
  cat(as.character(Sys.time()), "-",
      '/LDplot: Check number of SNP < 50 DONE\n')


  ### GET DATA
  cat(as.character(Sys.time()), "-",
      "/LDplot: Load data...\n")
  bm.wom <- getMarkerData(markerS3Path)
  cat(as.character(Sys.time()), "-",
      "/LDplot: Load data DONE\n")

  # COMPUTE LD
  cat(as.character(Sys.time()), "-",
      "/LDplot: Compute LD ...\n")
  ld <- LD(bm.wom, c(from, to), measure = "r2")
  cat(as.character(Sys.time()), "-",
      "/LDplot: Compute LD DONE\n")

  cat(as.character(Sys.time()), "-",
      "/LDplot: Create LD plot ...\n")
  p <- LD.plot(ld, snp.positions = bm.wom@snps$pos[from:to])
  cat(as.character(Sys.time()), "-",
      "/LDplot: Create LD plot DONE\n")

  cat(as.character(Sys.time()), "-",
      "/LDplot: Create response ... \n")
  res$status <- 200 # status for good GET response
  cat(as.character(Sys.time()), "-",
      "/LDplot: Create response DONE \n")
  cat(as.character(Sys.time()), "-",
      "/LDplot: END \n")
  p
}

##### Tables output #####

#* Table of selected SNPs
#* @tag Data
#* @param modelS3Path url of the model data file (rds file)
#* @param adj_method either bonferroni or FDR
#* @param thresh.p threshold for p values. If not specify return all values
#* @serializer unboxedJSON
#* @get /datatable
function(res, modelS3Path, adj_method, thresh.p = NA){

  out <- list(
    inputParams = list(
      modelS3Path = modelS3Path,
      adj_method = adj_method,
      thresh.p = as.character(thresh.p)
    )
  )

  cat(as.character(Sys.time()), "-",
      "/datatable: call with parameters parameters:\n")
  cat(
    "\t modelS3Path: ", modelS3Path,"\n",
    "\t adj_method: ", adj_method, "\n",
    "\t thresh.p: ", thresh.p, "\n"
  )


  ### CHECK PARAMETERS
  # Convert to numeric
  cat(as.character(Sys.time()), "-",
      "/datatable: Convert numeric parameters...\n")
  if (!is.na(as.numeric(thresh.p)) | is.na(thresh.p)) {
    thresh.p <- as.numeric(thresh.p)
  } else {
    cat(as.character(Sys.time()), "-",
        '/datatable: Error: "thresh.p" cannot be converted to numeric.\n')
    res$status <- 400 # bad request
    out$error <- '"thresh.p" should be a numeric value.'
    return(out)
  }
  cat(as.character(Sys.time()), "-",
      "/datatable: Convert numeric parameters DONE.\n")


  # LOAD MODEL
  cat(as.character(Sys.time()), "-",
      "/datatable: Load model...\n")
  gwa <- loadModel(modelS3Path)
  cat(as.character(Sys.time()), "-",
      "/datatable: Load model DONE\n")

  # CREATE DATATABLE
  cat(as.character(Sys.time()), "-",
      "/datatable: Adjust p-values ...\n")
  p.adj <- p.adjust(gwa$p, method = adj_method)
  cat(as.character(Sys.time()), "-",
      "/datatable: Adjust p-values DONE\n")

  ### RESPONSE
  cat(as.character(Sys.time()), "-",
      "/datatable: Create response ... \n")
  res$status <- 200 # status for good GET response
  if (is.na(thresh.p)) {
    out$data <- gwa
  } else {
    out$data <- gwa[p.adj < thresh.p, ]
  }

  # datatable(gwa[p.adj < thresh.p, ])
  cat(as.character(Sys.time()), "-",
      "/datatable: Create response DONE \n")
  cat(as.character(Sys.time()), "-",
      "/datatable: END \n")
  out

}


