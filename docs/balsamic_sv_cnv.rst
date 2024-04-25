************************************
Structural and Copy Number variants
************************************

Depending on the sequencing type, BALSAMIC is currently running the following structural and copy number variant callers:


.. list-table:: SV CNV callers
   :widths: 25 25 25 25 25
   :header-rows: 1

   * - Variant caller
     - Sequencing type
     - Analysis type
     - Somatic/Germline
     - Variant type
   * - AscatNgs
     - WGS
     - tumor-normal
     - somatic
     - CNV
   * - CNVkit
     - TGA, WES
     - tumor-normal, tumor-only
     - somatic
     - CNV
   * - Delly
     - TGA, WES, WGS
     - tumor-normal, tumor-only
     - somatic
     - SV, CNV
   * - Manta
     - TGA, WES, WGS
     - tumor-normal, tumor-only
     - somatic, germline
     - SV
   * - TIDDIT
     - WGS
     - tumor-normal, tumor-only
     - somatic
     - SV
   * - CNVpytor
     - WGS
     - tumor-only
     - somatic
     - CNV
   * - igh_dux4 (see note below)
     - WGS
     - tumor-normal, tumor-only
     - somatic
     - SV

Further details about a specific caller can be found in the links for the repositories containing the documentation for SV and CNV callers along with the links for the articles are listed in `bioinfo softwares <https://balsamic.readthedocs.io/en/latest/bioinfo_softwares.html>`_.

Note that igh_dux4 is not a variant caller itself. This is a custom script that uses samtools to detect read pairs supporting IGH::DUX4 rearrangements. In short, the command identifies discordant reads mapping to the IGH region and to either DUX4 or its homologous DUX4-like regions (see references for details). The inclusion of this feature aims to alleviate the failure of callers to detect this rearrangement. It is important to note, however, that the reported breakpoints are fixed to the IGH and DUX4 coordinates and are, therefore, imprecise and uncertain. Therefore, we advise caution when interpreting this information.


It is mandatory to provide the gender of the sample from BALSAMIC version >= 10.0.0 For CNV analysis.



**Pre-merge Filtrations**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


The copy number variants, identified using ascatNgs and `dellycnv`, are converted to deletion and duplications before they are merged using `SVDB` with `--bnd_distance = 5000` (distance between end points for the variants from different callers) and  `--overlap = 0.80` (percentage for overlapping bases for the variants from different callers).

Tumor and normal calls in `TIDDIT` are merged using `SVDB` with `--bnd_distance 500` and `--overlap = 0.80`.
Using a custom made script "filter_SVs.py", soft-filters are added to the calls based on the presence of the variant in the normal, with the goal of retaining only somatic variants as PASS.

Manta calls are filtered using bcftools to only keep variants that have evidence from 3 or more reads.

.. list-table:: SV filters
   :widths: 25 25 40
   :header-rows: 1

   * - Variant caller
     - Filter added
     - Filter expression
   * - TIDDIT
     - high_normal_af_fraction
     - (AF_N_MAX / AF_T_MAX) > 0.25
   * - TIDDIT
     - max_normal_allele_frequency
     - AF_N_MAX > 0.25
   * - TIDDIT
     - normal_variant
     - AF_T_MAX == 0 and ctg_t == False
   * - TIDDIT
     - in_normal
     - ctg_n == True and AF_N_MAX == 0 and AF_T_MAX <= 0.25
   * - Manta
     - low_pr_sr_count
     - SUM(FORMAT/PR[0:1]+FORMAT/SR[0:1]) < 4.0
   * - igh_dux4
     - samtools_igh_dux4
     - DV < 1


Further information regarding the TIDDIT tumor normal filtration: As translocation variants are represented by 2 BNDs in the VCF which allows for mixed assignment of soft-filters, a requirement for assigning soft-filters to translocations is that neither BND is PASS.


**Post-merge Filtrations**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`SVDB` prioritizes the merging of variants from SV and CNV callers to fetch position and genotype information,  in the following order:

.. list-table:: SVDB merge caller priority order
   :widths: 25 25 25 25
   :header-rows: 1

   * - TGA, WES
        tumor-normal
     - TGA, WES
        tumor-only
     - WGS
        tumor-normal
     - WGS
        tumor-only
   * - | 1. manta
       | 2. dellysv
       | 3. cnvkit
       | 4. dellycnv
     - | 1. manta
       | 2. dellysv
       | 3. cnvkit
       | 4. dellycnv
     - | 1. manta
       | 2. dellysv
       | 3. ascat
       | 4. dellycnv
       | 5. tiddit
       | 6. igh_dux4
     - | 1. manta
       | 2. dellysv
       | 3. dellycnv
       | 4. tiddit
       | 5. cnvpytor
       | 6. igh_dux4


The merged `SNV.somatic.<CASE_ID>.svdb.vcf.gz` file retains all the information for the variants from the caller in which the variants are identified, which are then annotated using `ensembl-vep`.
The SweGen and frequencies and the frequency of observed structural variants from clinical normal samples are annotated using `SVDB`.

The following filter applies for both tumor-normal and tumor-only samples in addition to caller specific filters.

*SWEGENAF*: SweGen Allele Frequency

::

    SWEGENAF <= 0.02  (or) SWEGENAF == "."

*Frq*: Frequency of observation of the variants from normal `Clinical` samples

::

    Frq <= 0.02  (or) Frq == "."

The variants scored as `PASS` are included in the final vcf file (`SNV.somatic.<CASE_ID>.svdb.<research/clinical>.filtered.pass.vcf.gz`).

The following command can be used to fetch the variants identified by a specific caller from merged structural and copy number variants.

::

  zgrep -E "#|<Caller>" <*.svdb.vcf.gz>


**Using GENS for WGS**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

GENS is a visualization tool similar to IGV, originally developed in Clinical Genomics Lund, and primarily for visualizing genomic copy number profiles from WGS samples.

To visualise the GENS-formatted files from BALSAMIC you need to have GENS installed, and to do this you can follow the instructions on the Clinical-Genomics-Lund GENS-repository:

- `Clinical-Genomics-Lund-GENS`_

.. _Clinical-Genomics-Lund-GENS: https://github.com/Clinical-Genomics-Lund/gens

Two files per sample are uploaded to GENS, one file with allele-frequencies from SNV & InDel germline-calls (BAF file) which can be used to aid the interpretation of the CN-profile, and one file with the Log2 copy number ratios normalized against a PON. Instructions for how to generate this PON using the BALSAMIC PON workflow can be found here:

`Generate GENS PON <https://balsamic.readthedocs.io/en/latest/balsamic_pon.html>`_.

There are three required arguments for creating the input files for GENS:
1. Genome interval file produced by GATK ``PreprocessIntervals`` (see instructions in GENS PON creation)
2. A gender specific PON (see instructions in GENS PON creation)
3. A population database VCF with variant positions to be reported in the BAF file.

We created the file mentioned in **3** using the file ``gnomad.genomes.r2.1.1.sites`` filtered with bcftools to only keep variants with an AF above 0.05.

.. code-block::

    bcftools view -i AF>=0.05 -Oz

To config BALSAMIC to run with GENS activated you supply these files like this:

::

  balsamic config case \
    --case-id <CASE_ID>
    --balsamic-cache </path/reference_cache/>
    --analysis-dir </path/analysis/>
    --fastq-path </path/fastq/>
    --gender <[male/female]>
    --analysis-workflow balsamic
    --genome-version hg19
    --tumor-sample-name <TUMOR_NAME>
    --genome-interval </path/genome_interval>
    --gens-coverage-pon </path/pon_file>
    --gnomad-min-af5 </path/population_vcf.vcf.gz>


**Genome Reference Files**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**How to generate genome reference files for ascatNGS**

Detailed information is available from `ascatNGS <https://github.com/cancerit/ascatNgs>`_ documentation

The file *SnpGcCorrections.tsv* prepared from the 1000 genome SNP panel.

**GC correction file:**

First step is to download the 1000 genome snp file and convert it from .vcf to .tsv. The detailed procedure to for this step is available from `ascatNGS-reference-files <https://github.com/cancerit/ascatNgs/wiki/Human-reference-files-from-1000-genomes-VCFs>`_ (Human reference files from 1000 genomes VCFs)

.. code-block::

    export TG_DATA=ftp://ftp.ensembl.org/pub/grch37/release-83/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz


Followed by:

.. code-block::

    curl -sSL $TG_DATA | zgrep -F 'E_Multiple_observations' | grep -F 'TSA=SNV' |\
    perl -ane 'next if($F[0] !~ m/^\d+$/ && $F[0] !~ m/^[XY]$/);\
    next if($F[0] eq $l_c && $F[1]-1000 < $l_p); $F[7]=~m/MAF=([^;]+)/;\
    next if($1 < 0.05); printf "%s\t%s\t%d\n", $F[2],$F[0],$F[1];\
    $l_c=$F[0]; $l_p=$F[1];' > SnpPositions_GRCh37_1000g.tsv


--or--

.. code-block::

    curl -sSL $TG_DATA | zgrep -F 'E_Multiple_observations' | grep -F 'TSA=SNV' |\
    perl -ane 'next if($F[0] !~ m/^\d+$/ && $F[0] !~ m/^[XY]$/); $F[7]=~m/MAF=([^;]+)/;\
    next if($1 < 0.05); next if($F[0] eq $l_c && $F[1]-1000 < $l_p);\
    printf "%s\t%s\t%d\n", $F[2],$F[0],$F[1]; $l_c=$F[0]; $l_p=$F[1];'\
    > SnpPositions_GRCh37_1000g.tsv

Second step is to use *SnpPositions.tsv* file and generate *SnpGcCorrections.tsv* file, more details see `ascatNGS-convert-snppositions <https://github.com/cancerit/ascatNgs/wiki/Convert-SnpPositions.tsv-to-SnpGcCorrections.tsv>`_

.. code-block::

    ascatSnpPanelGcCorrections.pl genome.fa SnpPositions.tsv > SnpGcCorrections.tsv

