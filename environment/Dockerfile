# hash:sha256:509782dc505c56a64bc3014904f709da298ed31bc75e037c661c754be0934d8e
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/bioconductor-320:RELEASE_3_20_NoGPU

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bcftools=1.19-1build2 \
        chromium-browser=2:1snap1-0ubuntu2 \
        openjdk-17-jre-headless=17.0.18+8-1~24.04.1 \
        pandoc=3.8.2.1-1 \
        samtools=1.19.2-1build2 \
        shiny-server=1.5.20.1002-2build1 \
        tabix=1.19+ds-1.1build3 \
    && rm -rf /var/lib/apt/lists/*

COPY postInstall /
RUN --mount=type=secret,id=secrets . /run/secrets/secrets \
    && /postInstall
