#
# This is a Plumber API. In RStudio 1.2 or newer you can run the API by
# clicking the 'Run API' button above.
#
# In RStudio 1.1 or older, see the Plumber documentation for details
# on running the API.
#
# Find out more about building APIs with Plumber here:
#
#    https://www.rplumber.io/
#


# Initialisation ----
cat(as.character(Sys.time()), "-",
    "Start Plumber API using 'plumber v'",
    as.character(packageVersion("plumber")), "\n")

cat(as.character(Sys.time()), "-",
    "Load `plumber` package","\n")
library(plumber)

# load gwas-engine functions
cat(as.character(Sys.time()), "-",
    "Load `gwas-engine`","\n")
sapply(list.files("GWAS-Engine/src/",
                  pattern = "\\.R$", # all files ending by ".R"
                  full.names = TRUE),
       source)

# load gwas-api functions
cat(as.character(Sys.time()), "-",
    "Load `gwas-api`'s functions","\n")
sapply(list.files("src/",
                  pattern = "\\.R$", # all files ending by ".R"
                  full.names = TRUE),
       source)

# create initialization logger
initLog <- logger$new("GWAS-API-INIT")


# Define the default png serializer for the images
my_serializer_png <- serializer_png(width = 40,
                                    height = 30,
                                    units = 'cm',
                                    res = 177,
                                    pointsize = 20)

# create new plumber router
initLog$log("create new router")
gwasApi <- pr()





# Set api description ----
initLog$log("Set api description")
gwasApi$setApiSpec(
  function(spec) {
    spec$info <- list(
      title = "GWAS API",
      description = "REST API for GWAS analysis",
      # termsOfService = "",
      contact = list(name = "Laboratory of Biometry and Bioinformatics, Hiroyoshi Iwata",
                     email = "iwata@ut-biomet.org"),
      license = list(name = "MIT",
                     url = "https://opensource.org/licenses/MIT"),
      version = "0.0.1"
    )
    spec$tags <- list(
      list(name = "Utils",
           description = "Endpoints for checking the API"),
      list(name = "GWAS",
           description = "Endpoints related to gwas analysis"),
      list(name = "Plots",
           description = "Endpoints related to plots drawn from a GWAS model"))
    spec
  }
)


# Set filters ----
initLog$log("Set api filters")
filterLogger <- logger$new("GWAS-API-REQUESTS")
# Log some information about the incoming requests
gwasApi <- gwasApi %>%
  pr_filter("logger",
            function(req){
              logger <- filterLogger
              logger$log(req$REQUEST_METHOD, req$PATH_INFO, "-",
                         req$HTTP_USER_AGENT, "@", req$REMOTE_ADDR)
              plumber::forward()
            })

# Redirect request sent to `/manplot` to `/manplot.html`
gwasApi <- gwasApi %>%
  pr_filter("redirect manplot",
            function(req){
              if (identical(req$PATH_INFO, "/manplot")) {
                logger <- filterLogger
                logger$log("Request to /manplot detected, redirect to /manplot.html")
                req$PATH_INFO <- "/manplot.html"
              }
              plumber::forward()
            })


# Set endpoints ----
initLog$log("Set api endpoints")

## Utils endpoints----
initLog$log("Set `/echo`")
### /echo ----
gwasApi <- gwasApi %>% pr_get(
  path = "/echo",
  tags = "Utils",
  comments = "Echo back the input",
  params = echo_params,
  handler = echo_handler,
  serializer = serializer_unboxed_json(),
)


### /version ----
initLog$log("Set `/version`")
gwasApi <- gwasApi %>% pr_get(
  path = "/version",
  tags = "Utils",
  comments = "Give information about current API version",
  params = version_params,
  handler = version_handler,
  serializer = serializer_unboxed_json(),
)





## GWAS endpoints ----
### /gwas ----
initLog$log("Set `/gwas`")
gwasApi <- gwasApi %>% pr_post(
  path = "/gwas",
  tags = "GWAS",
  comments = "Fit a GWAS model. This endpoint take Urls of geno and pheno data (and values of other GWAS parameters) and write an a json file to the give Url using a PUT request. It had been disign to work with amazon S3 services.",
  params = gwas_params,
  handler = gwas_handler,
  serializer = serializer_unboxed_json(),
)



### /adjustedResults ----
initLog$log("Set `/adjustedResults`")
gwasApi <- gwasApi %>% pr_get(
  path = "/adjustedResults",
  tags = "GWAS",
  comments = "Adjusted results. This endpoint calculate the adjusted p-values of the gwas analysis and return all the results or only the significant adjusted p-value. The results are return in json format.",
  params = adjustedResults_params,
  handler = adjustedResults_handler,
  serializer = serializer_unboxed_json(),
)



## Plot endpoints ----
### /manplot ----
initLog$log("Set `/manplot`")
gwasApi <- gwasApi %>% pr_get(
  path = "/manplot",
  tags = "Plots",
  comments = "DEPRECATED. Please use `/manplot.html` or `/manplot.png`. All request sent to this endpoint will be redirect to `/manplot.html` for backward compatibility.",
  params = manplot_params,
  handler = function(){
    "You should not be able to see that, this endpoint have been deprecated."
  },
  serializer = serializer_htmlwidget(),
)

### /manplot.html ----
initLog$log("Set `/manplot.html`")
gwasApi <- gwasApi %>% pr_get(
  path = "/manplot.html",
  tags = "Plots",
  comments = "Draw a Manhattan plot. This endpoint return the html code of a plotly interactive graph. By default only the 3000 points with the lowest p-values are display on the graph.",
  params = manplot_params,
  handler = create_manplot_handler(interactive = TRUE),
  serializer = serializer_htmlwidget(),
)


### /manplot.png ----
initLog$log("Set `/manplot.png`")
gwasApi <- gwasApi %>% pr_get(
  path = "/manplot.png",
  tags = "Plots",
  comments = "Draw a Manhattan plot. This endpoint return png Image of the graph. By default all the points are display on the graph.",
  params = manplot_params,
  handler = create_manplot_handler(interactive = FALSE),
  serializer = my_serializer_png,
)


### /LDplot ----
initLog$log("Set `/LDplot`")
gwasApi <- gwasApi %>% pr_get(
  path = "/LDplot",
  tags = "Plots",
  comments = "Draw a LD plot. This endpoint return a png image.",
  params = LDplot_params,
  handler = LDplot_handler,
  serializer = my_serializer_png,
)


# Mark deprecated endpoints ----


initLog$log("Set deprecated endpoints")
gwasApi$setApiSpec(
  utils::modifyList(
    gwasApi$getApiSpec(),
    list(
      paths = list(
        `/manplot` = list(
          get = list(
            deprecated = TRUE)
        )
      )
    )
  )
)

