
#' Connect to redis
#'
#' @param ... see ?redux::hiredis
#'
#' @return redis_api object or NULL
connectRedis <- function(...){
  logger <- logger$new("r-connectRedis()")
  logger$log("Connection to redis server ...")

  redis <- tryCatch({
    redis <- redux::hiredis(...)
    logger$log("Connection to redis succeed ...")
    logger$log("Connection to redis server DONE")
    redis
  }, error = function(err) {
    logger$log("Connection to redis failed. Error message is:\n\t", err$message)
    logger$log("Connection to redis server FAILED")
    return(NULL)
  })
  redis
}


#' Publish to redis channel
#'
#' @param action_type
#' @param status
#' @param message
#' @param channel
#'
#' @examples
redPub <- function(action_type, status, message, channel){
  logger <- logger$new("r-redPub()")

  if (is.null(REDIS)) {
    logger$log(time = FALSE, context = FALSE,
               "No connexion to redis, message not sent")
    return(NULL)
  }
  message = toJSON(list(
    service = 'GWAS_API',
    action_type = action_type,
    status = status,
    message = message
  ))
  n <- REDIS$PUBLISH(channel, message)
  logger$log(time = FALSE, context = FALSE,
             "Redis message received by", n," clients.")
}





#' Subscribe to redis channel
#'
#' @param channel
#'
#' @details for doc see:  https://cran.r-project.org/web/packages/redux/vignettes/redux.html
#'
#' @examples
redSub <- function(channel){
  REDIS$subscribe(channel,
                  transform = function(x) {
                    print("message received:")
                    print(x$value)
                    x$value
                  },
                  terminate = function(x) {
                    if (identical(x, "STOP")) {
                      print("STOP message received ---")
                      return(TRUE)
                    } else FALSE
                  },
                  n = Inf)
}
