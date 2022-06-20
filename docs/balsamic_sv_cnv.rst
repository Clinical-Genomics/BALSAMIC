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
     - Panel, WES
     - tumor-normal, tumor-only
     - somatic
     - CNV
   * - Delly
     - Panel, WES, WGS
     - tumor-normal, tumor-only
     - somatic
     - SV, CNV
   * - Manta
     - Panel, WES, WGS
     - tumor-normal, tumor-only
     - somatic, germline
     - SV
   * - TIDDIT
     - WGS
     - tumor-normal, tumor-only
     - somatic
     - SV

The links for the repositories containing the documentation for SV and CNV callers along with the links for the articles are listed in `bioinfo softwares <https://github.com/Clinical-Genomics/BALSAMIC/blob/master/docs/bioinfo_softwares.rst>`_.

The copy number variants, identified using ascatNgs and `dellycnv`, are converted to deletion and duplications before they are merged using `SVDB` with `--bnd_distance = 5000` and  `--overlap = 0.80`. `SVDB` priorityzes the merging of variants from SV and CNV callers to fetch position and genotype information,  in the following order:

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
   * - manta, dellysv, cnvkit, dellycnv
     - manta, dellysv, cnvkit, dellycnv
     - manta, dellysv, ascat, dellycnv,tiddit
     - manta, dellysv, dellycnv, tiddit


The merged `*.svdb.vcf.gz` file retains all the information for the variants from the caller in which the variants are identified, which are then annotated using `ensembl-vep`.

The following command can be used to fetch the variants identified by specific caller from merge structural and copy number.

::

  zgrep -E "#|<Caller>" <*.svdb.vcf.gz>