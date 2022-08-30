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
    desc = 'threshold to keep only the filter_nPoints with the lowest p-values for the plot.',
    type = "number",
    required = FALSE,
    isArray = FALSE),
  "filter_quant" = list(
    desc = 'threshold to keep only the filter_quant*100 % of the points with the lowest p-values for the plot (default no filtering)',
    type = "number",
    required = FALSE,
    isArray = FALSE)
)

create_manplot_handler <- function(interactive){

  function(res,
           gwas_url,
           adj_method = "bonferroni",
           thresh_p = 0.05,
           chr = NA,
           filter_pAdj = 1,
           filter_nPoints = ifelse(interactive, 3000, 10^100),
           filter_quant = 1){
    # # save call time.
    # callTime <- Sys.time()
    if (interactive) {
      logger <- logger$new("/manplot.html")
    } else {
      logger <- logger$new("/manplot.png")
    }

  inputParamsNames <- names(formals(rlang::current_fn()))
  inputParamsNames <- inputParamsNames[!inputParamsNames %in% c('res')]
  inputParams <- as.list(environment())[inputParamsNames]

  out <- list(
    inputParams = inputParams
  )

  logger$log("call with parameters:")
  logger$log(time = FALSE, context = FALSE,
             paste0(names(out$inputParams), ": ", out$inputParams,
                    collapse = '\n\t')
  )


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
    if (!is.na(as.numeric(filter_pAdj))) {
      filter_pAdj <- as.numeric(filter_pAdj)
    } else {
      logger$log('Error: "filter_pAdj" cannot be converted to numeric.')
      res$status <- 400 # bad request
      out$error <- '"filter_pAdj" should be a numeric value.'
      logger$log('Exit with error code 400')
      logger$log("END")
      return(out)
    }
    if (!is.na(as.numeric(filter_nPoints))) {
      filter_nPoints <- as.numeric(filter_nPoints)
    } else {
      logger$log('Error: "filter_nPoints" cannot be converted to numeric.')
      res$status <- 400 # bad request
      out$error <- '"filter_nPoints" should be a numeric value.'
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
                            interactive = interactive,
                            filter_pAdj = filter_pAdj,
                            filter_nPoints = filter_nPoints,
                            filter_quant = filter_quant,
                            outFile = NULL)
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

  inputParamsNames <- names(formals(rlang::current_fn()))
  inputParamsNames <- inputParamsNames[!inputParamsNames %in% c('res')]
  inputParams <- as.list(environment())[inputParamsNames]

  out <- list(
    inputParams = inputParams
  )

  logger$log("call with parameters:")
  logger$log(time = FALSE, context = FALSE,
             paste0(names(out$inputParams), ": ", out$inputParams,
                    collapse = '\n\t')
  )


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


# /pedNetwork ----
pedNetwork_params <- list(
  "ped_url" = list(
    desc = "url of the pedigree file (.csv)",
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "header" = list(
    desc = 'a logical value indicating whether the file contains the names of the variables as its first line. The default value is TRUE. In any cases, the column 1 will be interpreted as the individual id, column 2 as the first parent, column 3 as the second parent.',
    type = "string",
    required = FALSE,
    isArray = FALSE),
  "unknown_string" = list(
    desc = 'a character vector of strings which are to be interpreted as "unknown parent". By default: missing value in the file.',
    type = "string",
    required = TRUE,
    isArray = FALSE)
)

pedNetwork_handler <- function(res, ped_url, header = TRUE, unknown_string = ''){
  logger <- logger$new("/pedNetwork")

  inputParamsNames <- names(formals(rlang::current_fn()))
  inputParamsNames <- inputParamsNames[!inputParamsNames %in% c('res')]
  inputParams <- as.list(environment())[inputParamsNames]

  out <- list(
    inputParams = inputParams
  )

  logger$log("call with parameters:")
  logger$log(time = FALSE, context = FALSE,
             paste0(names(out$inputParams), ": ", out$inputParams,
                    collapse = '\n\t')
  )


  logger$log('Draw pedigree network plot ...')
  p <- draw_pedNetwork(pedUrl = ped_url,
                       unknown_string = unknown_string,
                       header = header,
                       outFile = NULL)
  logger$log('Draw pedigree network plot DONE')

  logger$log("Create response ... ")
  res$status <- 200 # status for good GET response
  logger$log("Create response DONE ")
  logger$log("END")

  p
}




# /relmat-heatmap ----
relmatHeatmap_params <- list(
  "relmat_url" = list(
    desc = "url of the relationship matrix file (.json)",
    type = "string",
    required = TRUE,
    isArray = FALSE)
)


create_relmatHeatmap_handler <- function(interactive){

  function(res, relmat_url){
    logger <- logger$new("/relmat-heatmap")

    inputParamsNames <- names(formals(rlang::current_fn()))
    inputParamsNames <- inputParamsNames[!inputParamsNames %in% c('res')]
    inputParams <- as.list(environment())[inputParamsNames]

    out <- list(
      inputParams = inputParams
    )

    logger$log("call with parameters:")
    logger$log(time = FALSE, context = FALSE,
               paste0(names(out$inputParams), ": ", out$inputParams,
                      collapse = '\n\t')
    )


    logger$log('Draw relationship heatmap plot ...')
    p <- draw_relHeatmap(relMatUrl = relmat_url,
                         format = 'json',
                         interactive = interactive,
                         outFile = NULL)
    logger$log('Draw relationship heatmap plot DONE')

    logger$log("Create response ... ")
    res$status <- 200 # status for good GET response
    logger$log("Create response DONE ")
    logger$log("END")

    p
  }
}



# /progenyBlup-plot ----
progenyBlupPlot_params <- list(
  "progenyBlup_url" = list(
    desc = "url of the relationship matrix file (.json)",
    type = "string",
    required = TRUE,
    isArray = FALSE),
  "sorting" = list(
    desc =  paste(
      "method to sort the individuals (X axis) can be:\\n",
      '- "asc": sort the BLUP expected value in ascending order',
      "(from left to right)\\n",
      '- "dec": sort the BLUP expected value in decreasing order',
      "(from left to right)\\n",
      "- any other value will sort the individuals in alphabetical order",
      "(from left to right)\\n",
      "Default is : Alphabetical"),
    type = "string",
    required = FALSE,
    isArray = FALSE
  )
)


progenyBlupPlot_handler <- function(res,
                                    progenyBlup_url,
                                    sorting = "alphabetical"){

  logger <- logger$new("/progenyBlup-plot")

  inputParamsNames <- names(formals(rlang::current_fn()))
  inputParamsNames <- inputParamsNames[!inputParamsNames %in% c('res')]
  inputParams <- as.list(environment())[inputParamsNames]

  out <- list(
    inputParams = inputParams
  )

  logger$log("call with parameters:")
  logger$log(time = FALSE, context = FALSE,
             paste0(names(out$inputParams), ": ", out$inputParams,
                    collapse = '\n\t')
  )


  logger$log('Draw progenies blups plot ...')
  p <- draw_progBlupsPlot(progEstimUrl = progenyBlup_url,
                          sorting = sorting,
                          outFile = NULL)
  logger$log('Draw relationship heatmap plot DONE')

  logger$log("Create response ... ")
  res$status <- 200 # status for good GET response
  logger$log("Create response DONE ")
  logger$log("END")

  p

}
