# GWAS_API

This is a REST api wrapper around "GWAS-Engine". For more detailed documentation, please check [GWAS-Engine's readme](https://github.com/ut-biomet/GWAS-Engine).



## API documentation

To generate the [`openapi.json`](./openapi.json) file one can run in `R`:

```R
source("gwas_api.R")
spec <- utils::modifyList(list(servers = list(list(url = ""))),
                               gwasApi$getApiSpec())
writeLines(jsonlite::toJSON(spec,
                            complex = "list",
                            auto_unbox = TRUE,
                            pretty = T,
                            digits = NA,
                            na = 'string'),
                            con = "openapi.json")
```

## Build and deploy

To this project use docker to build and deploy the API.

### Build

To build an image of the api, one can just run this command in the app folder:

```sh
docker build -t utbiomet/gwasapi ./
```

### Deploy

```sh
docker push utbiomet/gwasapi:latest
```

### Debug

```sh
docker run -d --rm --name gwasapi -p 8080:8080 utbiomet/gwasapi
```

```sh
docker exec -u 0 -it gwasapi bash
```

## Tests

To test the API, first launch a local instance of the API on port `8181`. This can be done through several ways, for example in `R` run:

```R
source('gwas_api.R');
gwasApi$run(port = 8080, host = '0.0.0.0')
```

Or you can run the Docker container (mapping the API on the port `8181`):

```sh
docker run -d --rm --name gwasapi -p 8080:8080 utbiomet/gwasapi
```

Then run the script `tests/testthat.R`: 
```sh
Rscript tests/testthat.R
```

For manually testing, you can use these urls when the api is running in docker:

```
file:///GWAS_API/GWAS-Engine/data/results/gwasResult.json
file:///GWAS_API/GWAS-Engine/data/geno/testMarkerData01.vcf.gz
file:///GWAS_API/GWAS-Engine/data/pheno/testPhenoData01.csv
```

For trait this one is in the pheno file: `Flowering.time.at.Arkansas`
