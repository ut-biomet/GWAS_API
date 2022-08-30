library(httr)
library(jsonlite)
library(testthat)

# Initialisation:

# setup global variables for tests

ip <- "127.0.0.1"
port <- 8080
host <- paste0("http://", ip, ":", port)

docker <- FALSE
if (docker) {
    dtaPref <- "file:///r-geno-tools-api/r-geno-tools-engine/data"
} else {
    dtaPref <- paste0("file://", getwd(), "/r-geno-tools-engine/data")
}

# Check if the API is available
testthat::expect_error(GET(host), NA)


# run tests
test_dir("tests/testthat", stop_on_failure = TRUE, stop_on_warning = TRUE)
