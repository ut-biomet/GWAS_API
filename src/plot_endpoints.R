# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Definition of the endpoints' parameters and handlers

# /manplot ----
manplot_params <- list(
  "gwas_url" = list(
    desc = "url of the result file saved by `/gwas` (json file)",
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "adj_method" = list(
    desc = 'p-value correction method: "holm", "hochberg", "bonferroni", "BH", "BY", "fdr", "none". Default: "bonferroni". (see https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust for more details)',
    type = "string",
    required = FALSE,
    isArray = FALSE),
  "thresh_p" = list(
    desc = "p value significant threshold (default 0.05)",
    type = "number",
    required = FALSE,
    isArray = FALSE),
  "chr" = list(
    desc = "names of the chromosomes to show separated using comma. Show all chromosomes if nothing is specified.",
    type = "string",
    required = FALSE,
    isArray = FALSE),
  "filter_pAdj" = list(
    desc = 'threshold to remove points with pAdj > filter_pAdj from the plot (default no filtering)',
    type = "number",
    required = FALSE,
    isArray = FALSE),
  "filter_nPoints" = list(
    desc = 'threshold to keep only the filter_nPoints with the lowest p-values for the plot (default 3000 points)',
    type = "number",
    required = FALSE,
    isArray = FALSE),
  "filter_quant" = list(
    desc = 'threshold to keep only the filter_quant*100 % of the points with the lowest p-values for the plot (default no filtering)',
    type = "number",
    required = FALSE,
    isArray = FALSE)
)

manplot_handler <- function(res,
                            gwas_url,
                            adj_method = "bonferroni",
                            thresh_p = 0.05,
                            chr = NA,
                            filter_pAdj = 1,
                            filter_nPoints = 3000,
                            filter_quant = 1){
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
  if (!is.na(as.numeric(filter_quant))) {
    filter_quant <- as.numeric(filter_quant)
  } else {
    logger$log('Error: "filter_quant" cannot be converted to numeric.')
    res$status <- 400 # bad request
    out$error <- '"filter_quant" should be a numeric value.'
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
                          chr = chr,
                          filter_pAdj = filter_pAdj,
                          filter_nPoints = filter_nPoints,
                          filter_quant = filter_quant)
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



# /LDplot ----
LDplot_params <- list(
  "geno_url" = list(
    desc = "url of the markers data file (.vcf.gz file)",
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "from" = list(
    desc = 'lower bound of the range of SNPs for which the LD is computed (`from` must be lower than `to`)',
    type = "integer",
    required = TRUE,
    isArray = FALSE),
  "to" = list(
    desc = "upper bound of the range of SNPs for which the LD is computed (the total number of SNP should be lower than 50)",
    type = "integer",
    required = TRUE,
    isArray = FALSE)
)

LDplot_handler <- function(res, geno_url, from, to){
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
