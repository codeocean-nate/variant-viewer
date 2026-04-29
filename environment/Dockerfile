# hash:sha256:977c934e665975b9bfd769c844dd9364f3d80755d3ac73cc7e1d06bfe3ca9e2c
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/bioconductor-320:RELEASE_3_20_NoGPU

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bcftools \
        build-essential \
        chromium-browser \
        libbz2-dev \
        libcairo2-dev \
        libcurl4-openssl-dev \
        libjpeg-dev \
        liblzma-dev \
        libpng-dev \
        libssl-dev \
        libxml2-dev \
        openjdk-17-jre-headless \
        pandoc \
        samtools \
        tabix \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("shiny")' \
    && Rscript -e 'remotes::install_version("shinydashboard")' \
    && Rscript -e 'remotes::install_version("shinyFiles")' \
    && Rscript -e 'remotes::install_version("shinyjs")' \
    && Rscript -e 'remotes::install_version("DT")' \
    && Rscript -e 'remotes::install_version("plotly")' \
    && Rscript -e 'remotes::install_version("arrow")' \
    && Rscript -e 'remotes::install_version("digest")' \
    && Rscript -e 'remotes::install_version("rmarkdown")' \
    && Rscript -e 'remotes::install_version("webshot2")' \
    && Rscript -e 'remotes::install_version("chromote")' \
    && Rscript -e 'remotes::install_version("remotes")' \
    && Rscript -e 'remotes::install_version("dplyr")'

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
