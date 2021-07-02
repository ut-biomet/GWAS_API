# GWAS_API

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
library(plumber)
plumber::plumb('plumber.R')$run(port = 8181)
```

Or you can run the Docker container (mapping the API on the port `8181`):

```sh
docker run -d --rm --name gwasapi -p 8181:8080 utbiomet/gwasapi
```

Then run the script `tests/testthat.R`:

```sh
Rscript tests/testthat.R
```
