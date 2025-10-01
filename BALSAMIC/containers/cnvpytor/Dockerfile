# ---- builder ----
FROM python:3.11-slim AS builder

ENV WORK_DIR=/opt
ENV VENV=${WORK_DIR}/venv
ENV PATH="${VENV}/bin:${PATH}"

RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get -y install --no-install-recommends git && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN python -m pip install --upgrade --no-cache-dir pip && \
    python -m venv "${VENV}"

# constraints file
RUN printf '%s\n' \
    'numpy==1.23.5' \
    'seaborn==0.12.2' \
    > "${WORK_DIR}/constraints.txt"

RUN pip install --no-cache-dir -c "${WORK_DIR}/constraints.txt" seaborn==0.12.2 numpy==1.23.5

# Clone CNVpytor
WORKDIR ${WORK_DIR}
RUN git clone https://github.com/abyzovlab/CNVpytor.git
WORKDIR ${WORK_DIR}/CNVpytor
RUN git checkout v1.3.1 && \
    pip install --no-cache-dir -c "${WORK_DIR}/constraints.txt" .

# Optional: download CNVpytor resources during build
RUN cnvpytor -download

# sanity check
RUN python - <<'PY'
import numpy, sys
print("Python:", sys.version.split()[0])
print("NumPy :", numpy.__version__)
assert numpy.__version__ == "1.23.5", "NumPy version drifted!"
PY


# ---- final ----
FROM python:3.11-slim

LABEL base.image="python:3.11-slim" \
      maintainer="Clinical Genomics" \
      about.contact="support@clinicalgenomics.se" \
      software="cnvpytor" \
      software.version="1.3.1" \
      about.summary="A tool for CNV discovery and genotyping from depth-of-coverage by mapped reads" \
      about.home="https://github.com/abyzovlab/CNVpytor" \
      about.documentation="https://github.com/abyzovlab/CNVpytor/wiki" \
      about.license="MIT License (MIT)"

RUN apt-get update && \
    apt-get -y upgrade && \
    apt-get -y install --no-install-recommends tabix && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV WORK_DIR=/opt
ENV VENV=${WORK_DIR}/venv
ENV PATH="${VENV}/bin:${PATH}"

COPY --from=builder ${WORK_DIR} ${WORK_DIR}

# Set user env
RUN adduser --disabled-password --gecos '' ubuntu && \
    chsh -s /bin/bash ubuntu && \
    mkdir -p /home/ubuntu

USER ubuntu
WORKDIR /home/ubuntu
CMD ["/bin/bash"]