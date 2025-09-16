FROM python:3.12-slim

ENV PIP_NO_CACHE_DIR=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1

# Install ScanITD from PyPI (0.9.1 is latest as of Nov 20, 2024)
# Supports Python 3.9–3.13 and provides the `scanitd` CLI.
# https://pypi.org/project/scanitd/
RUN pip install "scanitd==0.9.1"

# Default help
ENTRYPOINT ["scanitd"]
CMD ["--help"]
