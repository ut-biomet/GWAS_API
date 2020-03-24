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
docker run -d --rm --name gwasapi -p 8080:8080 utbiomet/gwasapi
```
