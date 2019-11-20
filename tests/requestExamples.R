# Author: Julien Diot juliendiot@ut-biomet.org
# 2019 The University of Tokyo
#
# Description:
# Some request example to try the GWAS api

#### PACKAGES ####
library(httr)
library(magick)
library(DT)


# 1 -- Fit a new model:
fixed = 4
trait_type = "quantitative"
trait = "Seed.length.width.ratio"
test = "lrt" # (lrt, Wald or score)
phenoDataId <- "testPhenoData01"
markerDataId <- "testMarkerData01"



reqFitModel <- POST(paste0("http://127.0.0.1:6831/gwas?",
                           "fixed=", fixed,
                           "&trait_type=", trait_type,
                           "&trait=", trait,
                           "&test=", test,
                           "&phenoDataId=", phenoDataId,
                           "&markerDataId=", markerDataId))
print(reqFitModel)
print(content(reqFitModel))

modId <- content(reqFitModel)[[1]]$modelId





# 2 -- ManhattanPlot
adj_method = "bonferroni"

reqManhattan <- GET(paste0("http://127.0.0.1:6831/manplot?",
                            "modelId=", modId,
                            "&adj_method=", adj_method))
print(reqManhattan)

image_read(paste0("http://127.0.0.1:6831/manplot?",
                  "modelId=", modId,
                  "&adj_method=", adj_method))





# 3 -- LDplot
from <- 1
to <- 11

reqLDplot <- GET(paste0("http://127.0.0.1:6831/LDplot?",
                        "&markerDataId=", markerDataId,
                        "&from=", from,
                        "&to=", to))
print(reqLDplot)

image_read(paste0("http://127.0.0.1:6831/LDplot?",
                  "&markerDataId=", markerDataId,
                  "&from=", from,
                  "&to=", to))





# 4 -- DataTable
reqDT <- GET(paste0("http://127.0.0.1:6831/datatable?",
                    "modelId=", modId,
                    "&adj_method=", adj_method))
print(reqDT)

