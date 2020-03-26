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
         fixed,
         tresh.maf = 0.05,
         tresh.callrate = 0.9){
  
  ### CHECK PARAMETERS
  if (!is.na(as.numeric(fixed))) {
    fixed <- as.numeric(fixed)
  } else {
  	fixed <- 0
  }
  
  ### GET DATA
  data <- loadData2(markerDataId, phenoDataId)
  
  ### GWAS 
  # calc model
  creatTime <- Sys.time()
  model <- gwas(data, trait, test, fixed, tresh.maf, tresh.callrate)
  modelId <- gsub("\\.", "-",
                  paste0("GWAS_",
                         markerDataId, "_", phenoDataId, "_", trait,"_",
                         as.numeric(creatTime)))
  
  # save model
  modPath <- paste0("data/models/", modelId, ".rds")
  saveRDS(model, file = modPath)

  # TODO:
  # write models's information in a database
  # so that endpoints to check already fitted model can be created
  #
  # fp <- digest(file = modPath) # model's file finger print (hash)

  res$status <- 201 # status for good post response
  list(
    message = "Model created",
    modelId = modelId
  )
  
}


##### Plots #####

#* Manhattan plot (type 2)
#* @tag ManhattanPlot
#* @param modelId GWAS model id
#* @param adj_method either bonferroni or FDR
#* @param thresh.p
#* @png
#* @get /manplot
function(modelId, adj_method, thresh.p){
  
  # LOAD MODEL
  gwa <- readRDS(paste0("data/models/", modelId, ".rds"))
  # can also be done in an external function: gwa <- loadModel(modelId) 

  # CREATE PLOT
  p.adj <- p.adjust(gwa$p, method = adj_method)
  col <- rep("black", nrow(gwa))
  col[gwa$chr %% 2 == 0] <- "gray50"
  col[p.adj < thresh.p] <- "green"
  manhattan(gwa, pch = 20, col = col,
            main = "TO DO trait name", # extract trait name from database
            sub = modelId)
}



#* LD plot
#* @tag LDPlot
#* @param markerDataId
#* @param from (total number of SNP should be < 50)
#* @param to (total number of SNP should be < 50)
#* @png
#* @get /LDplot
function(markerDataId, from, to){
  
  from <- as.numeric(from)
  to <- as.numeric(to)
  
  if (to - from > 50) {
    return("ERROR: number of SNP should be < 50.")
  }
  
  ### GET DATA 
  bm.wom <- getMarkerData(markerDataId)
  
  
  # COMPUTE LD
  ld <- LD(bm.wom, c(from, to), measure = "r2")
  LD.plot(ld, snp.positions = bm.wom@snps$pos[from:to])
}

##### Tables output #####

#* Table of selected SNPs
#* @tag SNP Table
#* @param modelId GWAS model id
#* @param adj_method either bonferroni or FDR
#* @param thresh.p
#* @serializer unboxedJSON
#* @get /datatable
function(modelId, adj_method, thresh.p){
  
  # LOAD MODEL
  gwa <- readRDS(paste0("data/models/", modelId, ".rds"))
  # can also be done in an external function: gwa <- loadModel(modelId)  
  
  # CREATE DATATABLE
  p.adj <- p.adjust(gwa$p, method = adj_method)
  # datatable(gwa[p.adj < thresh.p, ])
  gwa[p.adj < thresh.p, ]
}

# sel=gt.score[,colnames(gt.score)%in% gwa[p.adj < thresh.p, "id"]]
