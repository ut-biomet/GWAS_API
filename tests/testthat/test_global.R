# Author: Julien Diot juliendiot@ut-biomet.org
# 2020 The University of Tokyo
#
# Description:
# Test file for r-geno-tools-api

# init ----

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
    geno_url = paste0(dtaPref, "/geno/testMarkerData01.vcf.gz"),
    pheno_url = paste0(dtaPref, "/pheno/testPhenoData01.csv"),
    # upload_url = NA,
    trait = "Flowering.time.at.Arkansas",
    test = "score",
    fixed = 0,
    response = "quantitative",
    thresh_maf = 0.05,
    thresh_callrate = 0.9
  )


  # send request
  resp <- POST(path,
               query = query)

  # test status
  expect_equal(resp$status_code, 200)

  # test response content
  res <- jsonlite::fromJSON(content(resp))
  expect_is(res, "list")
  expect_equal(names(res), c("gwas", "metadata"))
  expect_is(res$metadata, "list")
  expect_equal(names(res$metadata),
               c("genoFP",
                 "phenoFP",
                 "trait",
                 "test",
                 "fixed",
                 "response",
                 "thresh_maf",
                 "thresh_callrate",
                 "date"))
  expect_is(fromJSON(res$gwas), "data.frame")
})

# Test GET /adjustedResults ----
test_that("GET /adjustedResults", {
  # create path and request
  path <- paste0(host,"/adjustedResults")
  query <- list(
    gwas_url = paste0(dtaPref, "/results/gwasResult.json"),
    adj_method = "bonferroni",
    filter_pAdj = 1,
    filter_nPoints = 3000,
    filter_quant = 1)

  # send request
  resp <- GET(path, query = query)
  # test status
  expect_equal(resp$status_code, 200)
  # test response content
  res <- jsonlite::fromJSON(content(resp))
  expect_is(res, "list")
  expect_equal(names(res),
               c("gwasAdjusted", "metadata"))
  expect_equal(names(res$metadata),
               c("genoFP",
                 "phenoFP",
                 "trait",
                 "test",
                 "fixed",
                 "response",
                 "thresh_maf",
                 "thresh_callrate",
                 "date",
                 "adj_method"))
  expect_is(res$gwasAdjusted, "data.frame")
})




# Test GET /manplot ----
manplotExt <- c("", ".html", ".png")

for (ext in manplotExt) {
  test_that(paste0("GET /manplot", ext), {

    # create path and request
    path <- paste0(host, paste0("/manplot", ext))
    query <- list(
      gwas_url = paste0(dtaPref, "/results/gwasResult.json"),
      adj_method = "bonferroni",
      thresh_p = 0.05
      # chr = NA
    )


    # send request
    resp <- GET(path,
                query = query)

    # test status
    expect_equal(resp$status_code, 200)

    # test response content
    if (identical(ext, ".png")) {
      expect_equal(resp$headers$`content-type`, "image/png")
    } else {
      expect_equal(resp$headers$`content-type`, "text/html; charset=UTF-8")
    }
  })
}




# Test GET /LDplot ----
test_that("GET /LDplot", {

  # creat path and request
  path <- paste0(host,"/LDplot")
  query <- list(
    geno_url = paste0(dtaPref, "/geno/testMarkerData01.vcf.gz"),
    from = "1",
    to = "20"
    )


  # send request
  resp <- GET(path,
              query = query)

  # test status
  expect_equal(resp$status_code, 200)

  # test response content
  expect_equal(resp$headers$`content-type`, "image/png")
})

# Test GET /pedNetwork ----
test_that("GET /pedNetwork", {

  # creat path and request
  path <- paste0(host,"/pedNetwork")
  query <- list(
    ped_url = paste0(dtaPref, "/pedigree/testPedData_char.csv")
  )


  # send request
  resp <- GET(path,
              query = query)

  # test status
  expect_equal(resp$status_code, 200)

  # test response content
  expect_equal(resp$headers$`content-type`, "text/html; charset=UTF-8")
})



# Test GET /relmat-heatmap ----
manplotExt <- c(".html", ".png")

for (ext in manplotExt) {
  test_that(paste0("GET /relmat-heatmap", ext), {

    # create path and request
    path <- paste0(host, paste0("/relmat-heatmap", ext))
    query <- list(
      relmat_url =  paste0(dtaPref, "/results/pedigreeRelationship.json")
    )


    # send request
    resp <- GET(path,
                query = query)

    # test status
    expect_equal(resp$status_code, 200)

    # test response content
    if (identical(ext, ".png")) {
      expect_equal(resp$headers$`content-type`, "image/png")
    } else {
      expect_equal(resp$headers$`content-type`, "text/html; charset=UTF-8")
    }
  })
}


# Test POST /relmat-ped ----
test_that("POST /relmat-ped", {

  # create path and request
  path <- paste0(host,"/relmat-ped")
  query <- list(
    ped_url = paste0(dtaPref, "/pedigree/testPedData_char.csv")
  )


  # send request
  resp <- POST(path,
               query = query)

  # test status
  expect_equal(resp$status_code, 200)

  # test response content
  res <- jsonlite::fromJSON(content(resp))
  expect_is(res, "list")
  expect_equal(names(res), c("relMat", "metadata"))
  expect_is(res$metadata, "list")
  expect_equal(names(res$metadata),
               c("info", "date", "nInds", "pedFP"))
  expect_is(res$relMat, "data.frame")
})

# Test POST /relmat-geno ----
test_that("POST /relmat-geno", {

  # create path and request
  path <- paste0(host,"/relmat-geno")
  query <- list(
    geno_url = paste0(dtaPref, "/geno/breedGame_geno.vcf.gz")
  )

  # send request
  resp <- POST(path,
               query = query)

  # test status
  expect_equal(resp$status_code, 200)

  # test response content
  res <- jsonlite::fromJSON(content(resp))
  expect_is(res, "list")
  expect_equal(names(res), c("relMat", "metadata"))
  expect_is(res$metadata, "list")
  expect_equal(names(res$metadata),
               c("info", "date", "nInds", "genoFP"))
  expect_is(res$relMat, "data.frame")
})




# Test POST /crossing-simulation ----
test_that("POST /crossing-simulation", {

  # create path and request
  path <- paste0(host,"/crossing-simulation")

  query <- list(
    geno_url = paste0(dtaPref, '/geno/breedGame_phasedGeno.vcf.gz'),
    crossTable_url = paste0(dtaPref, '/crossingTable/breedGame_crossTable.csv'),
    SNPcoord_url = paste0(dtaPref, '/SNPcoordinates/breedingGame_SNPcoord.csv'),
    chromosomeInfo_url = paste0(dtaPref, '/chromosomesInformation/breedingGame_chrInfo.csv'),
    nCross = 10
  )

  # send request
  resp <- POST(path,
               query = query)

  # test status
  expect_equal(resp$status_code, 201)

  # test response content
  resFile <- tempfile(fileext = '.vcf')
  writeBin(content(resp), resFile)
  source('../../r-geno-tools-engine/src/readWriteData.R')
  source('../../r-geno-tools-engine/src/utils.R')
  expect_error({
    readGenoData(file = resFile)
  }, NA)
})
