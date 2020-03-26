library(httr)
library(jsonlite)
library(testthat)

# Initialisation:

# setup global variables for tests

ip <- "127.0.0.1"
# ip <- "52.193.46.252" # SIP AWS server
port <- 8080
host <- paste0("http://", ip, ":", port)

# Check if the API is available
testthat::expect_error(GET(host), NA)


# run tests
test_dir("tests/testthat", stop_on_failure = TRUE, stop_on_warning = TRUE)
