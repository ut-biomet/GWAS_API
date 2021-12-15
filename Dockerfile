FROM rstudio/r-base:4.1-focal
LABEL maintainer="Julien Diot <juliendiot@ut-biomet.org>"



RUN apt-get update --allow-releaseinfo-change && apt-get install -y \
    libsodium-dev \
    libssl-dev \
    gcc-8-base \
    libcurl4 \
    libcurl4-openssl-dev \
    libxml2-dev \
    libz-dev \
    curl \
    pandoc \
    libhiredis-dev \
    libpng-dev \
    python3

# install dependencies
# install renv package
ENV RENV_VERSION 0.14.0
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR /r-geno-tools-api

# install deps using renv
COPY renv.lock renv.lock
RUN R -e 'renv::init()'


# get application code
COPY . .

EXPOSE 8080

ENTRYPOINT []
CMD R -e "source('r-geno-tools-api.R'); genoApi\$run(port = 8080, host = '0.0.0.0')"
