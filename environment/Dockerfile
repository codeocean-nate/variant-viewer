# hash:sha256:3124d2dc67c488f9973081c89a5795d2358d4e03505cf351359fb7512baa47bf
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/bioconductor-320:RELEASE_3_20_NoGPU

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bcftools=1.19-1build2 \
        build-essential=12.10ubuntu1 \
        chromium-browser=2:1snap1-0ubuntu2 \
        libbz2-dev=1.0.8-5.1build0.1 \
        libcairo2-dev=1.18.0-3build1 \
        libcurl4-openssl-dev=8.5.0-2ubuntu10.9 \
        libjpeg-dev=8c-2ubuntu11 \
        liblzma-dev=5.6.1+really5.4.5-1ubuntu0.2 \
        libpng-dev=1.6.43-5ubuntu0.6 \
        libssl-dev=3.0.13-0ubuntu3.9 \
        libxml2-dev=2.9.14+dfsg-1.3ubuntu3.7 \
        openjdk-17-jre-headless=17.0.18+8-1~24.04.1 \
        pandoc=3.8.2.1-1 \
        samtools=1.19.2-1build2 \
        shiny-server=1.5.20.1002-2build1 \
        tabix=1.19+ds-1.1build3 \
        zlib1g-dev=1:1.3.dfsg-3.1ubuntu2.1 \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("shiny", "1.13.0")' \
    && Rscript -e 'remotes::install_version("shinydashboard", "0.7.3")' \
    && Rscript -e 'remotes::install_version("shinyjs", "2.1.1")' \
    && Rscript -e 'remotes::install_version("DT", "0.34.0")' \
    && Rscript -e 'remotes::install_version("plotly", "4.12.0")' \
    && Rscript -e 'remotes::install_version("arrow", "23.0.1.2")' \
    && Rscript -e 'remotes::install_version("digest", "0.6.39")' \
    && Rscript -e 'remotes::install_version("rmarkdown", "2.31")' \
    && Rscript -e 'remotes::install_version("webshot2", "0.1.2")' \
    && Rscript -e 'remotes::install_version("chromote", "0.5.1")' \
    && Rscript -e 'remotes::install_version("remotes", "2.5.0")' \
    && Rscript -e 'remotes::install_version("dplyr", "1.2.1")' \
    && Rscript -e 'remotes::install_version("jsonlite", "2.0.0")'

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
