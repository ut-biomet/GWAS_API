#!/usr/bin/Rscript
arg <- commandArgs()


# extract working directory :
pat <- "(?<=(--file=)).*(?=\\/[a-zA-Z0-9\\_\\-]*\\.R$)"
wd <- regmatches(arg[4], gregexpr(pat, arg[4], perl = TRUE))[[1]][1]
setwd(wd)
print(paste0("set working directory to: ", wd))
cat("\n")

# run API
api <- plumber::plumb("plumber.R")
api$run(port = 8000,
        host = '0.0.0.0',
        swagger = TRUE)
