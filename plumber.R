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





##################################### Filter ###################################

#* Log some information about the incoming requests
#* @filter logger
function(req){
  cat(as.character(Sys.time()), "-",
      req$REQUEST_METHOD, req$PATH_INFO, "-",
      req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR, "\n")
  plumber::forward()
}




######################################################################

#* Echo back the input
#* @param msg The message to echo
#* @serializer unboxedJSON
#* @get /echo
function(msg=""){
  list(msg = paste0("The message is: '", msg, "'"))
}






##### GWAS #####
 
#* Fit a GWAS model
#* @tag Model fitting GWAS
#* @param markerDataId id of the genetic data
#* @param phenoDataId id of the phenotypic data
#* @param test The testing method (lrt, Wald or score)
#* @param trait The trait to be analyzed
#* @param trait_type The trait type: quantitative or binary
#* @param fixed The option chosen for fixed effect (number of PC, or kmeans, or none)
#* @serializer unboxedJSON
#* @post /gwas
function(markerDataId, phenoDataId, test, trait, trait_type, fixed){
  
  ### CHECK PARAMETERS 
  if (!is.na(as.numeric(fixed))) {
    fixed <- as.numeric(fixed)
  }
  
  ### GET DATA
  bm.wom <- loadData(markerDataId, phenoDataId, trait)
  
  ### CLEAN DATA
  # keep marker with a large enough MAF (>0.05) and low missing rate (callrate>0.9)
  bm.wom <- select.snps(bm.wom, maf > 0.05)
  bm.wom <- select.snps(bm.wom, callrate > 0.9)
  
  
  ### ????
  # Compute K (=Genetic relationship matrix for the lines with observed/non missing trait)
  K <- GRM(bm.wom)
  
  ### FIT MODEL
  if (is.numeric(fixed)) {
    if (test != "score") {
      gwa <- association.test(
        bm.wom,
        method = "lmm",
        response = trait_type,
        test = test,
        eigenK = eigen(K),
        p = fixed)
    } else {
      gwa <- association.test(
        bm.wom,
        method = "lmm",
        response = trait_type,
        test = "score",
        K = K,
        eigenK = eigen(K),
        p = fixed
      )
    }
  } else if (fixed == "kmeans") {
    
    # extract the score matrix 
    gt.score <- as.matrix(bm.wom)
    
    # k-means clustering
    tmp <- na.omit(t(gt.score))
    km <- kmeans(t(tmp), centers = 5, nstart = 10, iter.max = 100)
    grp <- as.factor(km$cluster)
    
    # Compute fixed effects from k-means clustering
    Xkm <- model.matrix(~grp[bm.wom@ped$id])
    
    if (test != "score") {
      gwa <- association.test(
        bm.wom,
        method = "lmm",
        response = trait_type,
        test = test,
        X = Xkm,
        eigenK = eigen(K)
      )
    } else{
      gwa <- association.test(
        bm.wom,
        method = "lmm",
        response = trait_type,
        test = "score",
        K = K,
        X = Xkm,
        eigenK = eigen(K)
      )
    }
  } else if (is.null(fixed)) {
    if (test != "score") {
      gwa <- association.test(
        bm.wom,
        method = "lmm",
        response = trait_type,
        test = test,
        eigenK = eigen(K)
      )
    } else{
      gwa <- association.test(
        bm.wom,
        method = "lmm",
        response = trait_type,
        test = "score",
        K = K,
        eigenK = eigen(K)
      )
    }
  }
  
  # SAVE MODEL
  creatTime <- Sys.time()
  modName <- paste0("GWAS_",
                    markerDataId, "_", phenoDataId, "_", trait,"_",
                    as.numeric(creatTime))
  modName <- gsub("\\.", "-", modName)
  modPath <- paste0("data/models/", modName, ".Rdata")
  
  if (!dir.exists("data/models")) {
    dir.create("data/models", mode = "0664")
  }
  save(gwa, file = modPath)

  # TODO:
  # write models's information in a database
  # so that endpoints to check already fitted model can be created
  #
  # fp <- digest(file = modPath) # model's file finger print (hash)
  # 
  #
  
  list(
    message = "Model created !",
    modelId = modName
  )
}

##### Plots #####

#* Manhattan plot
#* @tag ManhattanPlot
#* @param adj_method either bonferroni or FDR
#* @param modelId GWAS model id
#* @png
#* @get /manplot
function(modelId, adj_method){
  
  # LOAD MODEL
  load(paste0("data/models/", modelId, ".Rdata"))
  # can also be done in an external function: gwa <- loadModel(modelId) 
  
  
  # CREATE PLOT
  p.adj <- p.adjust(gwa$p, method = adj_method)
  col <- rep("black", nrow(gwa))
  col[gwa$chr %% 2 == 0] <- "gray50"
  col[p.adj < 0.05] <- "green"
  manhattan(gwa, pch = 20, col = col)
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
#* @param adj_method either bonferroni or FDR
#* @param modelId GWAS model id
#* @get /datatable
function(modelId, adj_method){
  # LOAD MODEL
  load(paste0("data/models/", modelId, ".Rdata"))
  # can also be done in an external function: gwa <- loadModel(modelId) 
  
  # CREATE DATATABLE
  p.adj <- p.adjust(gwa$p, method = adj_method)
  # datatable(gwa[p.adj < 0.05, ])
  gwa[p.adj < 0.05, ]
}

# sel=gt.score[,colnames(gt.score)%in% gwa[p.adj < 0.05, "id"]]
