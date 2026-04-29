# hash:sha256:63d2297e7154dcb4f386f54d711227df1164fe4acc692df08426cd066bbc6671
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/bioconductor-320:RELEASE_3_20_NoGPU

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bcftools \
        chromium-browser \
        openjdk-17-jre-headless \
        pandoc \
        samtools \
        shiny-server \
        tabix \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'remotes::install_version("shiny")' \
    && Rscript -e 'remotes::install_version("shinyFiles")' \
    && Rscript -e 'remotes::install_version("shinyjs")' \
    && Rscript -e 'remotes::install_version("webshot2")' \
    && Rscript -e 'remotes::install_version("chromote")' \
    && Rscript -e 'remotes::install_version("arrow")' \
    && Rscript -e 'remotes::install_version("dplyr")' \
    && Rscript -e 'remotes::install_version("tidyr")' \
    && Rscript -e 'remotes::install_version("DT")' \
    && Rscript -e 'remotes::install_version("plotly")' \
    && Rscript -e 'remotes::install_version("qqman")' \
    && Rscript -e 'remotes::install_version("digest")' \
    && Rscript -e 'remotes::install_version("rmarkdown")'

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
