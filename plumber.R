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

# required packages
library(plumber)
library(digest)
library(DT)
library(gaston) # for many functions

#* @apiTitle GWAS API
#* @apiDescription First test for GWAS API

# api <- plumber::plumb('plumber.R')
# api$run(port = 8080, host = '0.0.0.0', swagger = TRUE)

# load API's functions
sapply(list.files("functions",
                  pattern = "\\.R$", # all files ending by ".R"
                  full.names = TRUE),
       source)

# create models folder
if (!dir.exists("data/models")) {
  dir.create("data/models", mode = "0664")
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

#* Echo back the input
#* @param msg The message to echo
#* @serializer unboxedJSON
#* @get /echo
function(msg=""){
  cat(as.character(Sys.time()), "-",
      "/echo: call with parameters parameters:\n")
  cat("\t msg: ", msg,"\n")
  list(msg = paste0("The message is: '", msg, "'"))
}






##### GWAS #####
#* Fit a GWAS model (type 2)
#* @tag Model fitting GWAS
#* @param markerDataId id of the genetic data
#* @param phenoDataId id of the phenotypic data
#* @param trait The trait to be analyzed
#* @param test The testing method (lrt, Wald or score)
#* @param fixed The option chosen for fixed effect (number of PC, or none (0))
#* @param tresh.maf keep markers with a MAF > tresh.maf
#* @param tresh.callrate keep markers with a callrate > tresh.callrate
#* @serializer unboxedJSON
#* @post /gwas
function(res,
         markerDataId,
         phenoDataId,
         trait,
         test,
         fixed = 0,
         tresh.maf = 0.05,
         tresh.callrate = 0.9){
  # save call time.
  callTime <- Sys.time()
  out <- list(
    inputParams = list(
      markerDataId = markerDataId,
      phenoDataId = phenoDataId,
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
    "\t markerDataId: ", markerDataId,"\n",
    "\t phenoDataId: ", phenoDataId, "\n",
    "\t trait: ", trait, "\n",
    "\t test: ", test, "\n",
    "\t fixed: ", fixed, "\n",
    "\t tresh.maf: ", tresh.maf, "\n",
    "\t tresh.callrate: ", tresh.callrate, "\n"
  )


  ### CHECK PARAMETERS
  # Convert to numeric
  cat(as.character(Sys.time()), "-",
      "/gwas: Convert numeric parameters...\n")
  if (!is.na(as.numeric(fixed))) {
    fixed <- as.numeric(fixed)
  } else {
    cat(as.character(Sys.time()), "-",
        '/gwas: Error: "fixed" cannot be converted to numeric.\n')
    res$status <- 400 # bad request
    out$error <- '"fixed" should be a numeric value.'
    return(out)
  }

  if (!is.na(as.numeric(tresh.maf))) {
    tresh.maf <- as.numeric(tresh.maf)
  } else {
    cat(as.character(Sys.time()), "-",
        '/gwas: Error: "tresh.maf" cannot be converted to numeric.\n')
    res$status <- 400 # bad request
    out$error <- '"tresh.maf" should be a numeric value.'
    return(out)
  }

  if (!is.na(as.numeric(tresh.callrate))) {
    tresh.callrate <- as.numeric(tresh.callrate)
  } else {
    cat(as.character(Sys.time()), "-",
        '/gwas: Error: "tresh.callrate" cannot be converted to numeric.\n')
    res$status <- 400 # bad request
    out$error <- '"tresh.callrate" should be a numeric value.'
    return(out)
  }
  cat(as.character(Sys.time()), "-",
      "/gwas: Convert numeric parameters DONE.\n")



  ### GET DATA
  cat(as.character(Sys.time()), "-",
      "/gwas: Load data...\n")
  data <- loadData2(markerDataId, phenoDataId)
  cat(as.character(Sys.time()), "-",
      "/gwas: Load data DONE.\n")

  ### GWAS
  cat(as.character(Sys.time()), "-",
      "/gwas: Generate Gwas model...\n")
  # calc model
  model <- gwas(data, trait, test, fixed, tresh.maf, tresh.callrate)
  modelId <- gsub("\\.", "-",
                  paste0("GWAS_",
                         markerDataId, "_", phenoDataId, "_", trait,"_",
                         as.numeric(callTime)))
  cat(as.character(Sys.time()), "-",
      "/gwas: Generate Gwas model DONE:\n")
  cat(paste0("\t modelId: ", modelId,"\n"))

  # save model
  cat(as.character(Sys.time()), "-",
      "/gwas: Save model ... \n")
  modPath <- paste0("data/models/", modelId, ".rds")
  saveRDS(model, file = modPath)
  cat(as.character(Sys.time()), "-",
      "/gwas: Save model DONE: \n")
  cat(paste0("\t path: ", modPath,"\n"))

  # TODO:
  # see : https://docs.google.com/document/d/1gaTazFm_a6klD9krZPKGHc5x_D-SWL8K93M6dwsTiLE/edit
  # cat some information ID of the model,
  # write models's information in a database
  # so that endpoints to check already fitted model can be created
  #
  # Send all the information back so that, listenfield can write these in their database
  #
  # fp <- digest(file = modPath) # model's file finger print (hash)



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
#* @tag ManhattanPlot
#* @param modelId GWAS model id
#* @param adj_method either bonferroni or FDR
#* @param thresh.p
#* @png
#* @get /manplot
function(res, modelId, adj_method, thresh.p = 0.05){
  # # save call time.
  # callTime <- Sys.time()

  out <- list(
    inputParams = list(
      modelId = modelId,
      adj_method = adj_method,
      thresh.p = as.character(thresh.p)
    )
  )

  cat(as.character(Sys.time()), "-",
      "/manplot: call with parameters parameters:\n")
  cat(
    "\t modelId: ", modelId,"\n",
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
    return(out)
  }
  cat(as.character(Sys.time()), "-",
      "/manplot: Convert numeric parameters DONE.\n")

  # LOAD MODEL
  cat(as.character(Sys.time()), "-",
      "/manplot: load model ...\n")
  gwa <- readRDS(paste0("data/models/", modelId, ".rds"))
  # can also be done in an external function: gwa <- loadModel(modelId)
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
            sub = modelId)
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
#* @tag LDPlot
#* @param markerDataId
#* @param from (total number of SNP should be < 50)
#* @param to (total number of SNP should be < 50)
#* @png
#* @get /LDplot
function(res, markerDataId, from, to){

  out <- list(
    inputParams = list(
      markerDataId = markerDataId,
      from = from,
      to = to
    )
  )

  cat(as.character(Sys.time()), "-",
      "/LDplot: call with parameters parameters:\n")
  cat(
    "\t markerDataId: ", markerDataId,"\n",
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
    return(out)
  }
  cat(as.character(Sys.time()), "-",
      '/LDplot: Check number of SNP < 50 DONE\n')


  ### GET DATA
  cat(as.character(Sys.time()), "-",
      "/LDplot: Load data...\n")
  bm.wom <- getMarkerData(markerDataId)
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
#* @tag SNP Table
#* @param modelId GWAS model id
#* @param adj_method either bonferroni or FDR
#* @param thresh.p
#* @serializer unboxedJSON
#* @get /datatable
function(res, modelId, adj_method, thresh.p = 0.05){

  out <- list(
    inputParams = list(
      modelId = modelId,
      adj_method = adj_method,
      thresh.p = as.character(thresh.p)
    )
  )

  cat(as.character(Sys.time()), "-",
      "/datatable: call with parameters parameters:\n")
  cat(
    "\t modelId: ", modelId,"\n",
    "\t adj_method: ", adj_method, "\n",
    "\t thresh.p: ", thresh.p, "\n"
  )


  ### CHECK PARAMETERS
  # Convert to numeric
  cat(as.character(Sys.time()), "-",
      "/datatable: Convert numeric parameters...\n")
  if (!is.na(as.numeric(thresh.p))) {
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
  gwa <- readRDS(paste0("data/models/", modelId, ".rds"))
  # can also be done in an external function: gwa <- loadModel(modelId)
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
  out$data <- gwa[p.adj < thresh.p, ]
  # datatable(gwa[p.adj < thresh.p, ])
  cat(as.character(Sys.time()), "-",
      "/datatable: Create response DONE \n")
  cat(as.character(Sys.time()), "-",
      "/datatable: END \n")
  out

}

# sel=gt.score[,colnames(gt.score)%in% gwa[p.adj < thresh.p, "id"]]
