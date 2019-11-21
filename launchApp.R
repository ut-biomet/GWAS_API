#!/usr/bin/Rscript
arg <- commandArgs(trailingOnly = TRUE)

if (length(arg) == 0) {
    port <- 8080
} else if (length(arg) == 1) {
    port <- as.numeric(arg)
    if (is.na(port)) stop("provide only the port number eg: Rscript launchApp.R 8080")
} else stop("Too many arguments provided.")

.libPaths("/home/userplumberapi/R/x86_64-pc-linux-gnu-library/3.4")

api <- plumber::plumb('/home/userplumberapi/plumberApi/GWAS_API/plumber.R')
api$run(port = port,
        host = '0.0.0.0',
        swagger = TRUE)
