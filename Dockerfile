FROM rocker/r-base
LABEL maintainer="Julien Diot <juliendiot@ut-biomet.org>"

RUN apt-get update --allow-releaseinfo-change && apt-get install -y \
    libcurl4 \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    pandoc

RUN R -e "install.packages(c('BGLR', 'digest', 'DT', 'gaston', 'httr', 'magick', 'plumber'))"



EXPOSE 8000


# V1: the running argument must be the plumber.R file
# pay attention to the working directory of the application
# ENTRYPOINT ["R", "-e", "print(commandArgs());pr <- plumber::plumb(commandArgs()[4]);print(regmatches(commandArgs()[4], gregexpr('.*(?=\\\\/[a-zA-Z0-9\\\\_\\\\-]*\\\\.R$)',commandArgs()[4], perl = TRUE))[[1]][1]);setwd(regmatches(commandArgs()[4], gregexpr('.*(?=\\\\/[a-zA-Z0-9\\\\_\\\\-]*\\\\.R$)',commandArgs()[4], perl = TRUE))[[1]][1]);pr$run(host='0.0.0.0', port=8000, swagger = TRUE)"]
# CMD ["/app/plumber.R"]

# V2: the running argument must be the launchApp.R file
# the good app's working directory must be mannage in the launchApp.R file
ENTRYPOINT ["Rscript"]

CMD ["/app/launchApp.R"]
