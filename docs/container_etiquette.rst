===================
Container etiquette
===================

BALSAMIC uses singularity containers to perform the bioinformatics analysis. These containers are built using Docker and pushed to Docker Hub.
For more details on building containers using docker, please refer to the official docker documentation: https://docs.docker.com/

**Structure of Docker recipe**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    FROM <CONTAINER>:<VERSION>

    LABEL base.image="<CONTAINER>:<VERSION>"
    LABEL maintainer="Clinical Genomics"
    LABEL about.contact="support@clinicalgenomics.se"
    LABEL software="<NAME_OF_THE_MAIN_SOFTWARE>"
    LABEL software.version="<VERSION_OF_THE_MAIN_SOFTWARE>"
    LABEL about.summary="<DESCRIPTION_OF_THE_MAIN_SOFTWARE>"
    LABEL about.home="<URL_OF_THE_MAIN_SOFTWARE>""
    LABEL about.documentation="<DOCS_URL_OF_THE_MAIN_SOFTWARE>"
    LABEL about.license="MIT License (MIT)"

    RUN apt-get update && apt-get -y upgrade && \
        apt-get -y install --no-install-recommends && \
        <SOFTWARE_1 SOFTWARE_2> && \
        apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

    RUN ....

    USER    ubuntu
    WORKDIR /home/ubuntu
    CMD ["/bin/bash"]
