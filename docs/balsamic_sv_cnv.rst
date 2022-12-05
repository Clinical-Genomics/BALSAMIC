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

Further details about a specific caller can be found in the links for the repositories containing the documentation for SV and CNV callers along with the links for the articles are listed in `bioinfo softwares <https://balsamic.readthedocs.io/en/latest/bioinfo_softwares.html>`_.

It mandatory to provide the gender of the sample from BALSAMIC version >= 10.0.0 For CNV analysis.

The copy number variants, identified using ascatNgs and `dellycnv`, are converted to deletion and duplications before they are merged using `SVDB` with `--bnd_distance = 5000` (distance between end points for the variants from different callers) and  `--overlap = 0.80` (percentage for overlapping bases for the variants from different callers). `SVDB` prioritizes the merging of variants from SV and CNV callers to fetch position and genotype information,  in the following order:

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
     - | 1. manta
       | 2. dellysv
       | 3. dellycnv
       | 4. tiddit
       | 5. cnvpytor



The merged `SNV.somatic.<CASE_ID>.svdb.vcf.gz` file retains all the information for the variants from the caller in which the variants are identified, which are then annotated using `ensembl-vep`.
The SweGen and frequencies and the frequency of observed structural variants from clinical normal samples are annotated using `SVDB`.

The following filter applies for both tumor-normal and tumor-only samples in addition to caller specific filters.

*SWEGENAF*: SweGen Allele Frequency

::

    SWEGENAF <= 0.02  (or) SWEGENAF == "."

*Frq*: Frequency of observation of the variants from normal `Clinical` samples

::

    Frq <= 0.02  (or) Frq == "."

The variants scored as `PASS` or `MaxDepth` are included in the final vcf file (`SNV.somatic.<CASE_ID>.svdb.<research/clinical>.filtered.pass.vcf.gz`).

The following command can be used to fetch the variants identified by a specific caller from merged structural and copy number variants.

::

  zgrep -E "#|<Caller>" <*.svdb.vcf.gz>



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

**Attention:**
**BALSAMIC >= v11.0.0 removes unmapped reads from the bam and cram files for all the workflows.**
