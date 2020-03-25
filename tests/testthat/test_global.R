# Author: Julien Diot juliendiot@ut-biomet.org
# 2020 The University of Tokyo
#
# Description:
# Test file for GWAS_API

# init ----
# create tempoary folder
dir.create("tmp")


# Test GET /echo ----
test_that("GET /echo", {
  
  # creat path and request
  path <- paste0(host,"/echo")
  query <- list(msg = "hello world")
  
  
  # send request
  resp <- GET(path,
              query = query)
  
  # test status
  expect_equal(resp$status_code, 200)
  
  # test response content
  respContent <- content(resp)
  expect_equal(respContent$msg, paste0("The message is: \'", query$msg, "\'"))
})





# Test POST /gwas ----
test_that("POST /gwas", {
  
  # create path and request
  path <- paste0(host,"/gwas")
  query <- list(
    fixed = 4,
    trait_type = "quantitative",
    trait = "Seed.length.width.ratio",
    test = "lrt", # (lrt, Wald or score)
    phenoDataId = "testPhenoData01",
    markerDataId = "testMarkerData01"
  )
  
  
  # send request
  resp <- POST(path,
               query = query)
  
  # test status
  expect_equal(resp$status_code, 201)
  
  # test response content
  # test message
  respContent <- content(resp)
  expect_equal(respContent$message, "Model created !")
  
  # test model id
  expect_match(respContent$modelId, "^GWAS_") # start by "GWAS_"
  expect_match(respContent$modelId,
               chartr(".", "-", query$phenoDataId)) # contain phenoDataId
  expect_match(respContent$modelId,
               chartr(".", "-", query$markerDataId)) # contain markerDataId
  expect_match(respContent$modelId,
               chartr(".", "-", query$trait)) # contain trait
  
  # save model id for future tests
  saveRDS(respContent$modelId, "tmp/modelID.rds")
  
})



# Test GET /manplot ----
test_that("GET /manplot", {
    
  # creat path and request
  path <- paste0(host,"/manplot")
  query <- list(
    modelId = readRDS("tmp/modelID.rds"),
    adj_method = "bonferroni"
  )
  
  
  # send request
  resp <- GET(path,
              query = query)
  
  # test status
  expect_equal(resp$status_code, 200)
  
  # test response content
  expect_equal(resp$headers$`content-type`, "image/png")
})



# Test GET /LDplot ----
test_that("GET /LDplot", {
    
  # creat path and request
  path <- paste0(host,"/LDplot")
  query <- list(
    markerDataId = "testMarkerData01",
    from = 1,
    to = 11
    )
  
  
  # send request
  resp <- GET(path,
              query = query)
  
  # test status
  expect_equal(resp$status_code, 200)
  
  # test response content
  expect_equal(resp$headers$`content-type`, "image/png")
})



# Test GET /datatable ----
test_that("GET /datatable", {
    
  # creat path and request
  path <- paste0(host,"/datatable")
  query <- list(
    modelId = readRDS("tmp/modelID.rds"),
    adj_method = "bonferroni"
  )
  
  
  # send request
  resp <- GET(path,
              query = query)
  
  # test status
  expect_equal(resp$status_code, 200)
  
  # test response content
  respContent <- content(resp)
  expect_equal(respContent, data.frame())
})



# end ----
unlink("tmp", recursive = TRUE, force = TRUE)