#!/bin/bash

# Tumor-Normal Panel
balsamic config case \
  -t tests/test_data/fastq/S1_R_1.fastq.gz \
  -n tests/test_data/fastq/S2_R_1.fastq.gz  \
  --case-id TN_panel \
  --analysis-dir run_tests/ \
  -r reference/GRCh37/reference.json  \
  --singularity BALSAMIC/containers/BALSAMIC_latest.sif \
  -p tests/test_data/references/panel/panel.bed \
  --output-config balsamic_config.json

# Tumor-only Panel
balsamic config case \
  -t tests/test_data/fastq/S1_R_1.fastq.gz \
  --case-id T_panel \
  --analysis-dir run_tests/ \
  -r reference/GRCh37/reference.json  \
  --singularity BALSAMIC/containers/BALSAMIC_latest.sif \
  -p tests/test_data/references/panel/panel.bed \
  --output-config balsamic_config.json

# Tumor-Normal WGS
balsamic config case \
  -t tests/test_data/fastq/S1_R_1.fastq.gz \
  -n tests/test_data/fastq/S2_R_1.fastq.gz  \
  --case-id TN_wgs \
  --analysis-dir run_tests/ \
  -r reference/GRCh37/reference.json  \
  --singularity BALSAMIC/containers/BALSAMIC_latest.sif \
  --output-config balsamic_config.json

# Tumor-only WGS
balsamic config case \
  -t tests/test_data/fastq/S1_R_1.fastq.gz \
  --case-id T_wgs \
  --analysis-dir run_tests/ \
  -r reference/GRCh37/reference.json  \
  --singularity BALSAMIC/containers/BALSAMIC_latest.sif \
  --output-config balsamic_config.json
