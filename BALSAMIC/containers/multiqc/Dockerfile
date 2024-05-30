FROM python:3.11.3-slim

LABEL base.image="python:3.11.3-slim"
LABEL about.home="https://github.com/Clinical-Genomics/BALSAMIC"
LABEL about.documentation="https://balsamic.readthedocs.io/"
LABEL about.license="MIT License (MIT)"
LABEL about.description="Bioinformatic analysis pipeline for somatic mutations in cancer"

RUN python -m pip install --upgrade --no-cache-dir pip && \
    pip install multiqc==1.22.1

RUN adduser --disabled-password --gecos '' ubuntu && \
chsh -s /bin/bash && \
mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]
