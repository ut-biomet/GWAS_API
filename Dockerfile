FROM rstudio/plumber
LABEL maintainer="Julien Diot <juliendiot@ut-biomet.org>"

EXPOSE 8080

RUN useradd plumber

RUN apt-get update --allow-releaseinfo-change && apt-get install -y \
    gcc-8-base \
    libcurl4 \
    libcurl4-openssl-dev \
    libxml2-dev \
    curl \
    pandoc \
    libhiredis-dev


# set and check R package repository, repo will be saved in the ~/.Rprofile file
RUN SNAPSHOT="https://cran.microsoft.com/snapshot/2021-07-01/" && \
    if [ $(curl -s -o /dev/null -w "%{http_code}" $SNAPSHOT) = "200" ] ; then echo "OK: R pkg repository server accessible." ; else echo "ERROR: R pkg repository server not accessible. status = $(curl -s -o /dev/null -w "%{http_code}" $SNAPSHOT)" ; fi &&\
    touch /.Rprofile && \
    echo "options(repos = c(CRAN = '$SNAPSHOT'))" >> ~/.Rprofile

# get application code
COPY ./ /GWAS_API
# using git (app repo must be public)
# RUN git clone --depth=1 https://github.com/ut-biomet/GWAS_API.git

RUN chown plumber.plumber -R /GWAS_API
RUN chmod -R 774 /GWAS_API


# install dependencies
RUN Rscript /GWAS_API/installDeps.R
RUN Rscript /GWAS_API/GWAS-Engine/installDeps.R
RUN rm ~/.Rprofile


USER plumber

ENTRYPOINT []
CMD R -e "setwd('/GWAS_API'); api <- plumber::plumb('/GWAS_API/plumber.R'); api\$run(port = 8080, host = '0.0.0.0')"
