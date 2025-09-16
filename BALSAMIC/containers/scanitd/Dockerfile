FROM python:3.12-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1

RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --upgrade pip setuptools wheel cython

# Install as root → CLI ends up in /usr/local/bin (always in PATH)
RUN pip install "scanitd==0.9.1"

# Drop to non-root for runtime
RUN useradd -m -u 1000 scanitd
USER scanitd
WORKDIR /home/scanitd

ENTRYPOINT ["scanitd"]
CMD ["--help"]