# FROM rstudio/plumber
FROM rstudio/r-base:4.1-focal
LABEL maintainer="Julien Diot <juliendiot@ut-biomet.org>"

WORKDIR /GWAS_API

RUN apt-get update --allow-releaseinfo-change && apt-get install -y \
    libpng-dev \
    libsodium-dev \
    libssl-dev \
    gcc-8-base \
    libcurl4 \
    libcurl4-openssl-dev \
    libxml2-dev \
    libz-dev \
    curl \
    pandoc \
    libhiredis-dev

# install dependencies
# install renv package
ENV RENV_VERSION 0.14.0
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"
# install deps using renv
COPY renv.lock renv.lock
RUN R -e 'renv::init()'


# get application code
COPY . .


EXPOSE 8080

ENTRYPOINT []
CMD R -e "api <- plumber::plumb('plumber.R'); api\$run(port = 8080, host = '0.0.0.0')"
