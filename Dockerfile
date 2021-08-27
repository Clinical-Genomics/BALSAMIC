FROM continuumio/miniconda3:4.9.2-alpine

LABEL base_image="continuumio/miniconda3:4.9.2-alpine"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.maintainer="Hassan Foroughi hassan dot foroughi at scilifelab dot se" 
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"
LABEL about.version="8.0.1"

ENV LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
ENV PATH="/opt/conda/bin/:${PATH}"
RUN apk add --no-cache bash gcc git wget libc-dev zlib-dev linux-headers

ARG CONTAINER_NAME
ARG WORK_DIR=project

# Generate a custom locale
RUN wget -q -O /etc/apk/keys/sgerrand.rsa.pub https://alpine-pkgs.sgerrand.com/sgerrand.rsa.pub && \
    wget https://github.com/sgerrand/alpine-pkg-glibc/releases/download/2.33-r0/glibc-2.33-r0.apk && \
    apk add glibc-2.33-r0.apk && \
    wget https://github.com/sgerrand/alpine-pkg-glibc/releases/download/2.33-r0/glibc-bin-2.33-r0.apk && \
    wget https://github.com/sgerrand/alpine-pkg-glibc/releases/download/2.33-r0/glibc-i18n-2.33-r0.apk && \
    apk add glibc-bin-2.33-r0.apk glibc-i18n-2.33-r0.apk && \
    /usr/glibc-compat/bin/localedef -i ${LANG%%.*} -f ${LANG##*.} ${LANG} && \
    echo "export LANG=${LANG}" > /etc/profile.d/locale.sh && \
    rm \
        /etc/apk/keys/sgerrand.rsa.pub \
        glibc-*

# Copy all project files
COPY . /${WORK_DIR}

RUN cd /${WORK_DIR}/BALSAMIC/containers/${CONTAINER_NAME}/ && /bin/sh ${CONTAINER_NAME}.sh ${CONTAINER_NAME}
RUN if [ "$CONTAINER_NAME" = "balsamic" ]; then cd /project && pip install --no-cache-dir . ; fi

# Clean work environment
RUN rm -rf /project && conda clean --all --yes
