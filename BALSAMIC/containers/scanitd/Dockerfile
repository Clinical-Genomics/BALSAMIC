FROM python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# toolchain helpers
RUN python -m pip install --upgrade pip setuptools wheel cython

# 👇 Pin Click to <8.2 before installing ScanITD
RUN pip install "click==8.1.7"

# Install ScanITD
RUN pip install "scanitd==0.9.1"

# (optional) drop to non-root at runtime
RUN useradd -m -u 1000 scanitd
USER scanitd
WORKDIR /home/scanitd

ENTRYPOINT ["scanitd"]
CMD ["--help"]