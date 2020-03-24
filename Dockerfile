FROM trestletech/plumber
LABEL maintainer="Julien Diot <juliendiot@ut-biomet.org>"

EXPOSE 8080

RUN useradd plumber

RUN apt-get update --allow-releaseinfo-change && apt-get install -y \
    libcurl4 \
    libcurl4-openssl-dev \
    libxml2-dev \
    curl \
    pandoc

# set and check R package repository
RUN SNAPSHOT="https://cran.microsoft.com/snapshot/2020-03-17/" && \
    if [ $(curl -s -o /dev/null -w "%{http_code}" $SNAPSHOT) = "200" ] ; then echo "OK: R pkg repository server accessible." ; else echo "ERROR: R pkg repository server not accessible. status = $(curl -s -o /dev/null -w "%{http_code}" $SNAPSHOT)" ; fi &&\
    touch /.Rprofile && \
    echo "options(repos = c(CRAN = '$SNAPSHOT'))" >> ~/.Rprofile

RUN R -e "install.packages(c('BGLR', 'digest', 'DT', 'gaston', 'httr', 'magick'))"
RUN rm ~/.Rprofile



# get application code
COPY ./ /GWAS_API

# using git (app repo must be public)
# RUN git clone --depth=1 https://github.com/ut-biomet/GWAS_API.git

RUN chown plumber.plumber -R /GWAS_API
RUN chmod -R 774 /GWAS_API



USER plumber

ENTRYPOINT []
CMD R -e "setwd('/GWAS_API'); api <- plumber::plumb('/GWAS_API/plumber.R'); api\$run(port = 8080, host = '0.0.0.0', swagger = TRUE)"
