************************************
Structural and Copy Number Variants
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

Further details about a specific caller can be found in the links for the repositories containing the documentation for SV and CNV callers along with the links for the articles are listed in `bioinfo softwares <https://github.com/Clinical-Genomics/BALSAMIC/blob/master/docs/bioinfo_softwares.rst>`_.

The copy number variants, identified using ascatNgs and `dellycnv`, are converted to deletion and duplications before they are merged using `SVDB` with `--bnd_distance = 5000` (distance between the end points for a variant from different callers) and  `--overlap = 0.80` (percentage for overlapping bases between variants from diffent callers). `SVDB` prioritizes the merging of variants from SV and CNV callers to fetch position and genotype information,  in the following order:

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


The merged `*.svdb.vcf.gz` file retains all the information for the variants from the caller in which the variants are identified, which are then annotated using `ensembl-vep`.

The following command can be used to fetch the variants identified by a specific caller from merged structural and copy number variants.

::

  zgrep -E "#|<Caller>" <*.svdb.vcf.gz>