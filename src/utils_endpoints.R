# Author: Julien Diot juliendiot@ut-biomet.org
# 2021 The University of Tokyo
#
# Description:
# Definition of the endpoints' parameters and handlers


# /echo ----
echo_params <- list(
  "msg" = list(
    desc = "The message to echo",
    type = "string",
    required = FALSE,
    isArray = FALSE)
)

echo_handler <- function(msg=""){
  logger <- logger$new("/echo")
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

    list(msg = paste0("The message is: '", msg, "'"))
}

# /version ----
version_params <- list()
version_handler <- function(){
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

  logger$log("calculate MD5 sum of all '.R' files of the api")
  apiRfiles <- c(
    dir(".", recursive = FALSE,
        full.names = TRUE, all.files = TRUE, pattern = ".R$"),
    dir("r-geno-tools-engine/src/", recursive = TRUE,
        full.names = TRUE, all.files = TRUE, pattern = ".R$"),
    dir("src/", recursive = TRUE,
        full.names = TRUE, all.files = TRUE, pattern = ".R$")
  )
  allFP <- sapply(apiRfiles, function(f){digest::digest(file = f)})
  out$RfilesFingerPrint <- digest::digest(allFP)
  logger$log("END")
  return(out)
}
