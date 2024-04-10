==============
Short tutorial
==============

Here a short tutorial is provided for BALSAMIC (**version** = 15.0.0).

Regarding fastq-inputs
---------------------

Previous versions of BALSAMIC only accepted one fastq-pair per sample, which required concatenation of fastq-pairs if multiple existed.

The current version BALSAMIC takes ``--fastq-path`` instead which is a path to a directory containing ALL fastq-files you want to include in the analysis, for tumor and normal (if it exists for the analysis).

**NOTE**: The fastq-files in ``--fastq-path`` need to contain the names from ``--tumor-sample-name [sample_name]`` and ``--normal-sample-name [sample_name]`` as a sub-string in the fastq-names to correctly assign them to their respective sample.

Running a test sample
---------------------

Example config demo case:

::

  balsamic config case \
    --analysis-dir demo/
    --balsamic-cache ~/balsamic_cache
    --fastq-path tests/test_data/fastq/
    --case-id demo_run_balsamic
    --gender female
    --analysis-workflow balsamic
    --genome-version hg19
    --tumor-sample-name S1
    --panel-bed tests/test_data/references/panel/panel.bed





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
 
And now run balsamic through SLURM. Make sure you set your SLURM project account using ``--account`` if your local
settings require it:

::

  balsamic run analysis --sample-config demo/demo_run_balsamic/demo_run_balsamic.json \
    --profile slurm --qos low --account development --run-analysis

And now run balsamic through QSUB. Make sure you set your QSUB project account using ``--account`` if your local
settings require it: 

::

  balsamic run analysis --sample-config demo/demo_run_balsamic/demo_run_balsamic.json \
    --profile qsub --qos low --account development --run-analysis


And running workflow without submitting jobs. Set number of cores by passing an argument to snakemake as seen below:

::

  balsamic run analysis --sample-config demo/demo_run_balsamic/demo_run_balsamic.json \
    --run-mode local --snakemake-opt "--cores 8" --run-analysis
