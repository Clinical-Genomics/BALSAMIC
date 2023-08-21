**********************
Annotation resources
**********************

BALSAMIC annotates somatic single nucleotide variants (SNVs) using ``ensembl-vep`` and ``vcfanno``. Somatic structural variants (SVs), somatic copy-number variants (CNVs) and germline single nucleotide variants are annotated using only ``ensembl-vep``. All SVs and CNVs are merged using ``SVDB`` before annotating for `Target Genome Analysis (TGA)` or `Whole Genome Sequencing (WGS)` analyses.

`BALSAMIC` adds the following annotation from `gnomAD` database using ``vcfanno``.

.. list-table:: gnomAD
   :widths: 50 50
   :header-rows: 1

   * - VCF tag
     - description
   * - GNOMADAF_popmax
     - maximum allele frequency across populations
   * - GNOMADAF
     - fraction of the reads supporting the alternate allele, allelic frequency

`BALSAMIC` adds the following annotation from `ClinVar` database using ``vcfanno``.

.. list-table:: ClinVar
   :widths: 50 50
   :header-rows: 1

   * - VCF tag
     - description
   * - CLNACC
     - Variant Accession and Versions
   * - CLNREVSTAT
     - ClinVar review status for the Variation ID
   * - CLNSIG
     - Clinical significance for this single variant
   * - CLNVCSO
     - Sequence Ontology id for variant type
   * - CLNVC
     - Variant type
   * - ORIGIN
     - Allele origin

The values for `ORIGIN` are described below:

.. list-table:: ORIGIN
   :widths: 25 25
   :header-rows: 1

   * - Value
     - Annotation
   * - 0
     - unknown
   * - 1
     - germline
   * - 2
     - somatic
   * - 4
     - inherited
   * - 8
     - paternal
   * - 16
     - maternal
   * - 32
     - *de-novo*
   * - 64
     - biparental
   * - 128
     - uniparental
   * - 256
     - not-tested
   * - 512
     - tested-inconclusive
   * - 1073741824
     - other

`BALSAMIC` uses `ensembl-vep` to add the following annotation from `COSMIC` database.

.. list-table:: COSMIC
   :widths: 50 50
   :header-rows: 1

   * - VCF tag
     - description
   * - COSMIC_CDS
     - CDS annotation
   * - COSMIC_GENE
     - gene name
   * - COSMIC_STRAND
     - strand
   * - COSMIC_CNT
     - number of samples with this mutation in the `COSMIC` database
   * - COSMIC_AA
     - peptide annotation

`BALSAMIC` adds the following annotation for SNVs from `CADD` database using ``vcfanno``.

.. list-table:: CADD
   :widths: 50 50
   :header-rows: 1

   * - VCF tag
     - description
   * - CADD
     - Combined Annotation Dependent Depletion

Where relevant, `BALSAMIC` uses `ensembl-vep` to annotate somatic and germline SNVs and somatic SVs/CNVs from `1000genomes (phase3)`, `ClinVar`, `ESP, HGMD-PUBLIC`, `dbSNP`, `gencode`, `gnomAD`, `polyphen`, `refseq`, and `sift` databases.
The following annotations are added by `ensembl-vep`.

VEP has a setting for the maximum size of a structural variant that it will annotate, currently this is set to the size of the size of chromosome 1 (in hg19) (`--max_sv_size 249250621`).


.. list-table:: ensembl-vep
   :widths: 10 60
   :header-rows: 1

   * - Annotation
     - description
   * - Allele
     - the variant allele used to calculate the consequence
   * - Gene
     - Ensembl stable ID of affected gene
   * - Feature
     - Ensembl stable ID of feature
   * - Feature type
     - type of feature. Currently one of Transcript, RegulatoryFeature, MotifFeature.
   * - Consequence
     - consequence type of this variant
   * - Position in cDNA
     - relative position of base pair in cDNA sequence
   * - Position in CDS
     - relative position of base pair in coding sequence
   * - Position in protein
     - relative position of amino acid in protein
   * - Amino acid change
     - only given if the variant affects the protein-coding sequence
   * - Codon change
     - the alternative codons with the variant base in upper case
   * - Co-located variation
     - identifier of any existing variants
   * - VARIANT_CLASS
     - Sequence Ontology variant class
   * - SYMBOL
     - the gene symbol
   * - SYMBOL_SOURCE
     - the source of the gene symbol
   * - STRAND
     - the DNA strand (1 or -1) on which the transcript/feature lies
   * - ENSP
     - the Ensembl protein identifier of the affected transcript
   * - FLAGS
     - | transcript quality flags:
       | cds_start_NF: CDS 5' incomplete
       | cds_end_NF: CDS 3' incomplete
   * - SWISSPROT
     - Best match UniProtKB/Swiss-Prot accession of protein product
   * - TREMBL
     - Best match UniProtKB/TrEMBL accession of protein product
   * - UNIPARC
     - Best match UniParc accession of protein product
   * - HGVSc
     - the HGVS coding sequence name
   * - HGVSp
     - the HGVS protein sequence name
   * - HGVSg
     - the HGVS genomic sequence name
   * - HGVS_OFFSET
     - Indicates by how many bases the HGVS notations for this variant have been shifted
   * - SIFT
     - the SIFT prediction and/or score, with both given as prediction(score)
   * - PolyPhen
     - the PolyPhen prediction and/or score
   * - MOTIF_NAME
     - The source and identifier of a transcription factor binding profile aligned at this position
   * - MOTIF_POS
     - The relative position of the variation in the aligned TFBP
   * - HIGH_INF_POS
     - A flag indicating if the variant falls in a high information position of a transcription factor binding profile (TFBP)
   * - MOTIF_SCORE_CHANGE
     - The difference in motif score of the reference and variant sequences for the TFBP
   * - CANONICAL
     - a flag indicating if the transcript is denoted as the canonical transcript for this gene
   * - CCDS
     - the CCDS identifer for this transcript, where applicable
   * - INTRON
     - the intron number (out of total number)
   * - EXON
     - the exon number (out of total number)
   * - DOMAINS
     - the source and identifer of any overlapping protein domains
   * - DISTANCE
     - Shortest distance from variant to transcript
   * - AF
     - Frequency of existing variant in 1000 Genomes
   * - AFR_AF
     - Frequency of existing variant in 1000 Genomes combined African population
   * - AMR_AF
     - Frequency of existing variant in 1000 Genomes combined American population
   * - EUR_AF
     - Frequency of existing variant in 1000 Genomes combined European population
   * - EAS_AF
     - Frequency of existing variant in 1000 Genomes combined East Asian population
   * - SAS_AF
     - Frequency of existing variant in 1000 Genomes combined South Asian population
   * - AA_AF
     - Frequency of existing variant in NHLBI-ESP African American population
   * - EA_AF
     - Frequency of existing variant in NHLBI-ESP European American population
   * - gnomAD_AF
     - Frequency of existing variant in gnomAD exomes combined population
   * - gnomAD_AFR_AF
     - Frequency of existing variant in gnomAD exomes African/American population
   * - gnomAD_AMR_AF
     - Frequency of existing variant in gnomAD exomes American population
   * - gnomAD_ASJ_AF
     - Frequency of existing variant in gnomAD exomes Ashkenazi Jewish population
   * - gnomAD_EAS_AF
     - Frequency of existing variant in gnomAD exomes East Asian population
   * - gnomAD_FIN_AF
     - Frequency of existing variant in gnomAD exomes Finnish population
   * - gnomAD_NFE_AF
     - Frequency of existing variant in gnomAD exomes Non-Finnish European population
   * - gnomAD_OTH_AF
     - Frequency of existing variant in gnomAD exomes combined other combined populations
   * - gnomAD_SAS_AF
     - Frequency of existing variant in gnomAD exomes South Asian population
   * - MAX_AF
     - Maximum observed allele frequency in 1000 Genomes, ESP and gnomAD
   * - MAX_AF_POPS
     - Populations in which maximum allele frequency was observed
   * - CLIN_SIG
     - ClinVar clinical significance of the dbSNP variant
   * - BIOTYPE
     - Biotype of transcript or regulatory feature
   * - APPRIS
     - Annotates alternatively spliced transcripts as primary or alternate based on a range of computational methods. NB: not available for GRCh37
   * - TSL
     - Transcript support level. NB: not available for GRCh37
   * - PUBMED
     - Pubmed ID(s) of publications that cite existing variant
   * - SOMATIC
     - Somatic status of existing variant(s); multiple values correspond to multiple values in the Existing_variation field
   * - PHENO
     - Indicates if existing variant is associated with a phenotype, disease or trait; multiple values correspond to multiple values in the Existing_variation field
   * - GENE_PHENO
     - Indicates if overlapped gene is associated with a phenotype, disease or trait
   * - BAM_EDIT
     - Indicates success or failure of edit using BAM file
   * - GIVEN_REF
     - Reference allele from input
   * - REFSEQ_MATCH
     - | the RefSeq transcript match status; contains a number of flags indicating whether this RefSeq transcript matches the underlying reference sequence and/or an Ensembl transcript (more information):

       - rseq_3p_mismatch: signifies a mismatch between the RefSeq transcript and the underlying primary genome assembly sequence. Specifically, there is a mismatch in the 3' UTR of the RefSeq model with respect to the primary genome assembly (e.g. GRCh37/GRCh38).
       - rseq_5p_mismatch: signifies a mismatch between the RefSeq transcript and the underlying primary genome assembly sequence. Specifically, there is a mismatch in the 5' UTR of the RefSeq model with respect to the primary genome assembly.
       - rseq_cds_mismatch: signifies a mismatch between the RefSeq transcript and the underlying primary genome assembly sequence. Specifically, there is a mismatch in the CDS of the RefSeq model with respect to the primary genome assembly.
       - rseq_ens_match_cds: signifies that for the RefSeq transcript there is an overlapping Ensembl model that is identical across the CDS region only. A CDS match is defined as follows: the CDS and peptide sequences are identical and the genomic coordinates of every translatable exon match. Useful related attributes are: rseq_ens_match_wt and rseq_ens_no_match.
       - rseq_ens_match_wt: signifies that for the RefSeq transcript there is an overlapping Ensembl model that is identical across the whole transcript. A whole transcript match is defined as follows: 1) In the case that both models are coding, the transcript, CDS and peptide sequences are all identical and the genomic coordinates of every exon match. 2) In the case that both transcripts are non-coding the transcript sequences and the genomic coordinates of every exon are identical. No comparison is made between a coding and a non-coding transcript. Useful related attributes are: rseq_ens_match_cds and rseq_ens_no_match.
       - rseq_ens_no_match: signifies that for the RefSeq transcript there is no overlapping Ensembl model that is identical across either the whole transcript or the CDS. This is caused by differences between the transcript, CDS or peptide sequences or between the exon genomic coordinates. Useful related attributes are: rseq_ens_match_wt and rseq_ens_match_cds.
       - rseq_mrna_match: signifies an exact match between the RefSeq transcript and the underlying primary genome assembly sequence (based on a match between the transcript stable id and an accession in the RefSeq mRNA file). An exact match occurs when the underlying genomic sequence of the model can be perfectly aligned to the mRNA sequence post polyA clipping.
       - rseq_mrna_nonmatch: signifies a non-match between the RefSeq transcript and the underlying primary genome assembly sequence. A non-match is deemed to have occurred if the underlying genomic sequence does not have a perfect alignment to the mRNA sequence post polyA clipping. It can also signify that no comparison was possible as the model stable id may not have had a corresponding entry in the RefSeq mRNA file (sometimes happens when accessions are retired or changed). When a non-match occurs one or several of the following transcript attributes will also be present to provide more detail on the nature of the non-match: rseq_5p_mismatch, rseq_cds_mismatch, rseq_3p_mismatch, rseq_nctran_mismatch, rseq_no_comparison
       - rseq_nctran_mismatch: signifies a mismatch between the RefSeq transcript and the underlying primary genome assembly sequence. This is a comparison between the entire underlying genomic sequence of the RefSeq model to the mRNA in the case of RefSeq models that are non-coding.
       - rseq_no_comparison: signifies that no alignment was carried out between the underlying primary genome assembly sequence and a corresponding RefSeq mRNA. The reason for this is generally that no corresponding, unversioned accession was found in the RefSeq mRNA file for the transcript stable id. This sometimes happens when accessions are retired or replaced. A second possibility is that the sequences were too long and problematic to align (though this is rare).
   * - CHECK_REF
     - Reports variants where the input reference does not match the expected reference
   * - HGNC_ID
     - A unique ID provided by the HGNC for each gene with an approved symbol
   * - MANE
     - indicating if the transcript is the MANE Select or MANE Plus Clinical transcript for the gene.
   * - miRNA
     - Reports where the variant lies in the miRNA secondary structure.


`BALSAMIC` adds the following annotation from `Swegen` database using ``vcfanno`` for SNVs and SVDB for SVs.

.. list-table:: Swegen SNV
   :widths: 50 150 50
   :header-rows: 1

   * - VCF tag
     - description
     - variant type
   * - SWEGENAF
     - allele frequency from 1000 Swedish genomes project
     - SNV, SV
   * - SWEGENAAC_Hom
     - allele counts of homozygous variants
     - SNV
   * - SWEGENAAC_Het
     - allele counts of heterozygous variants
     - SNV
   * - SWEGENAAC_Hemi
     - allele counts of hemizygous variants
     - SNV
   * - swegen_obs
     - allele count
     - SV

`BALSAMIC` adds the following annotation from database of normal `Clinical` samples using ``vcfanno`` for SNVs and SVDB for SVs.

.. list-table:: Clinical Normal samples SNV
   :widths: 50 150 50
   :header-rows: 1

   * - VCF tag
     - description
     - variant type
   * - Frq
     - Frequency of observation of the variants from normal `Clinical` samples
     - SNV, SV
   * - Obs
     - allele counts of the variant in normal `Clinical` samples
     - SNV
   * - Hom
     - allele counts of the homozygous variant in normal `Clinical` samples
     - SNV
   * - clin_obs
     - allele counts
     - SV

