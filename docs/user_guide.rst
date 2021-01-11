==============
Short tutorial
==============

Here a short toturial is provided for BALSAMIC (**version** = 6.0.4). 

.. contents::

Step 1. generate a reference
----------------------------


First reference files must be downloaded. Let's assume BALSAMIC is installed and available at `D_BALSAMIC-base_5.0.0`,
and a COSMIC key is generated via: https://cancer.sanger.ac.uk/cosmic/help/file_download 

The following commands will create and download reference directory at `./BALSAMIC_reference` (change this path if you
want it to be created in another location):

::

  # Given:
  # 1. COSMIC key is in variable $COSMIC_KEY
  # 2. reference directory is in variable: $reference_path
  # 3. BALSAMIC container file path is in $BALSAMIC_container

  balsamic init reference \
    --cosmic-key ${COSMIC_KEY} \
    --outdir ${reference_path} \
    --genome-version hg19 \
    --singularity ${BALSAMIC_container} \
    --quiet \
    --run-analysis
  

A `json` file with reference specifications is created at: `${reference_path}/BALSAMIC_version/hg19/reference.json`

Step 2. Running a test sample
-----------------------------
Now a config file for a test run must be created. Let's use the test data in `tests` directory:

::

  balsamic config case \
    --tumor tests/test_data/fastq/S1_R_1.fastq.gz \
    --normal tests/test_data/fastq/S2_R_1.fastq.gz \
    --case-id demo_run_balsamic \
    --analysis-dir demo/ \
    --panel-bed tests/test_data/references/panel/panel.bed \
    --reference-config ${reference_path}/BALSAMIC_version/hg19/reference.json \
    --singularity ${BALSAMIC_container} \
    --output-config demo_run_balsamic.json 

Notes:

- If you want to test tumor_only mode, remove the `--normal tests/test_data/fastq/S2_R_1.fastq.gz` line.
- `--output-config demo_run_balsamic.json` is also optional

Let's try a dry run and see everything is in place:

::

  balsamic run analysis --sample-config demo/demo_run_balsamic/demo_run_balsamic.json

Command above should exit a similar output as below:

::

  Job counts:
  count jobs
  1 BaseRecalibrator
  1 CollectAlignmentSummaryMetrics
  1 CollectHsMetrics
  1 CollectInsertSizeMetrics
  1 IndelRealigner
  1 MarkDuplicates
  1 RealignerTargetCreator
  1 all
  1 bwa_mem
  1 cnvkit_single
  1 fastp
  1 fastqc
  13  haplotypecaller
  1 haplotypecaller_merge
  1 manta_germline
  1 manta_tumor_only
  1 mergeBam_tumor
  1 mergeBam_tumor_gatk
  1 multiqc
  1 mutect2_merge
  13  mutect2_tumor_only
  1 sambamba_exon_depth
  1 sambamba_panel_depth
  1 samtools_sort_index
  1 somatic_snv_indel_vcf_merge
  1 split_bed_by_chrom
  1 strelka_germline
  1 vardict_merge
  13  vardict_tumor_only
  7 vep
  72
  This was a dry-run (flag -n). The order of jobs does not reflect the order of execution.
 
And now run balsamic through SLURM. Make sure you set your SLURM project account using `--account` if your local
settings require it:

::

  balsamic run analysis --sample-config demo/demo_run_balsamic/demo_run_balsamic.json \
    --profile slurm --qos low --account development --run-analysis

And now run balsamic through QSUB. Make sure you set your QSUB project account using `--account` if your local
settings require it: 

::

  balsamic run analysis --sample-config demo/demo_run_balsamic/demo_run_balsamic.json \
    --profile qsub --qos low --account development --run-analysis


And running workflow without submitting jobs. Set number of cores by passing an argument to snakemake as seen below:

::

  balsamic run analysis --sample-config demo/demo_run_balsamic/demo_run_balsamic.json \
    --run-mode local --snakemake-opt "--cores 8" --run-analysis
