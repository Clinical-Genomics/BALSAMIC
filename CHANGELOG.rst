[X.X.X]
-------

Added:
^^^^^^
* Fastq concatenation https://github.com/Clinical-Genomics/BALSAMIC/pull/1069
* `CADD` SNV references https://github.com/Clinical-Genomics/BALSAMIC/pull/1126
* `CADD` SNV annotation https://github.com/Clinical-Genomics/BALSAMIC/pull/1150
* Samtools stats, flagstat, idxstat to WGS workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/1176
* Functionality for dynamically assigning fastq-info to sample dict in config from input fastq-dir https://github.com/Clinical-Genomics/BALSAMIC/pull/1176
* Annotate SNVs with cancer germline SNV observations from loqusDB https://github.com/Clinical-Genomics/BALSAMIC/pull/1178
* Annotate SNVs with somatic SNV observations from loqusDB https://github.com/Clinical-Genomics/BALSAMIC/pull/1187
* Tests for Annotation with Cancer germline, somatic and clinical observations, and swegen frequencies https://github/Clinical-Genomics/BALSAMIC/pull/1190
* Annotate SVs with somatic SV observations from loqusDB https://github.com/Clinical-Genomics/BALSAMIC/pull/1194
* Support singularity bind paths with different destination directories https://github/Clinical-Genomics/BALSAMIC/pull/1211
* Added `--rerun-trigger mtime` option to Snakemake command https://github.com/Clinical-Genomics/BALSAMIC/pull/1217
* `CADD` container https://github.com/Clinical-Genomics/BALSAMIC/pull/1222
* `Container ettiquette` to ReadtheDocs https://github.com/Clinical-Genomics/BALSAMIC/pull/1132
* `htslib` (samtools, bcftools tabix) container https://github.com/Clinical-Genomics/BALSAMIC/pull/1234
* Release version support for cache generation https://github.com/Clinical-Genomics/BALSAMIC/pull/1231
* `CADD` scores for INDELs https://github.com/Clinical-Genomics/BALSAMIC/pull/1238
* `CADD` reference to tests https://githuc.com/Clinical-Genomics/BALSAMIC/pull/1241
* Add cache version option to config case https://github.com/Clinical-Genomics/BALSAMIC/pull/1244
* `cnvkit` container https://github.com/Clinical-Genomics/BALSAMIC/pull/1252
* `PureCN` container https://github.com/Clinical-Genomics/BALSAMIC/pull/1255
* `GATK` container https://github.com/Clinical-Genomics/BALSAMIC/pull/1266
* Resolved FASTQ paths to sample dictionary (balsamic logging) https://github.com/Clinical-Genomics/BALSAMIC/pull/1275
* Picard HsMetrics and CollectGcBiasMetrics for WGS https://github.com/Clinical-Genomics/BALSAMIC/pull/1288
* Command-line arguments and rules for creation of GENS files https://github.com/Clinical-Genomics/BALSAMIC/pull/1279


Changed:
^^^^^^^^
* Changed CN header field in cnvpytor in cnvpytor_tumor_only to be Float instead of Integer https://github.com/Clinical-Genomics/BALSAMIC/pull/1182
* Changed samples in case_config.json from being a dict to a list of dicts  https://github.com/Clinical-Genomics/BALSAMIC/pull/1176
* Updated snakemake version to 7.25.0 https://github.com/Clinical-Genomics/BALSAMIC/pull/1099
* Updated cryptography version to 41.0.1 https://github.com/Clinical-Genomics/BALSAMIC/pull/1173
* Refactor bam and fastq inputs in snakemake to call pydantic model functions https://github.com/Clinical-Genomics/BALSAMIC/pull/1176
* Standardised alignment workflows to WGS-workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/1176
* Implemented parallel trimming and alignment in all workflows per lane https://github.com/Clinical-Genomics/BALSAMIC/pull/1176
* All bam-QC tools take the final dedup.realign bamfile as input https://github.com/Clinical-Genomics/BALSAMIC/pull/1176
* Validation of pydantic models done both during config and run https://github.com/Clinical-Genomics/BALSAMIC/pull/1176
* Refactored fastp rules, and changed order of UMI-trimming and quality trimming https://github.com/Clinical-Genomics/BALSAMIC/pull/1176
* Fix pydantic version (<2.0) https://github.com/Clinical-Genomics/BALSAMIC/pull/1191
* Refactor constants https://github.com/Clinical-Genomics/BALSAMIC/pull/1174
* Move models to their own folder https://github.com/Clinical-Genomics/BALSAMIC/pull/1176
* Balsamic init workflow refactoring https://github.com/Clinical-Genomics/BALSAMIC/pull/1188
* Updated cryptography version to 41.0.2 https://github.com/Clinical-Genomics/BALSAMIC/pull/1205
* Refactor snakemake executable command generation https://github/Clinical-Genomics/BALSAMIC/pull/1211
* Updated Python version to 3.11 and its dependencies https://github.com/Clinical-Genomics/BALSAMIC/pull/1216
* Tools versions in doc https:/github.com/Clinical-Genomics/BALSAMIC/pull/1239
* Reuse common Balsamic CLI options https://github.com/Clinical-Genomics/BALSAMIC/pull/1242
* Update `reference.json` file to use relative paths https://github.com/Clinical-Genomics/BALSAMIC/pull/1251
* Update pydantic to v2 while maintaining support for v1 models https://github.com/Clinical-Genomics/BALSAMIC/pull/1253
* `PCT_PF_READS_IMPROPER_PAIRS` QC threshold lowered to 5% https://github.com/Clinical-Genomics/BALSAMIC/issues/1265
* Migrate Metrics models to pydantic v2 https://github.com/Clinical-Genomics/BALSAMIC/pull/1270
* Migrate Snakemake models to pydantic v2 https://github.com/Clinical-Genomics/BALSAMIC/pull/1268
* Migrate Cache models to pydantic v2 https://github.com/Clinical-Genomics/BALSAMIC/pull/1274
* Made BALSAMIC compatible with multiple PON creation workflows https://github.com/Clinical-Genomics/BALSAMIC/pull/1279


Fixed:
^^^^^^
* vcf2cytosure container https://github.com/Clinical-Genomics/BALSAMIC/pull/1159
* Link external fastqs to case folder & create case directory https://github.com/Clinical-Genomics/BALSAMIC/pull/1195
* vcf2cytosure container missing constants https://github.com/Clinical-Genomics/BALSAMIC/pull/1198
* Bash commands in vep_somatic_clinical_snv https://github.com/Clinical-Genomics/BALSAMIC/pull/1200
* Fix SVDB annotation intermediate rule https://github.com/Clinical-Genomics/BALSAMIC/pull/1218
* Broken documentation links https://github.com/Clinical-Genomics/BALSAMIC/pull/1226
* Updated contributors in main README https://github.com/Clinical-Genomics/BALSAMIC/pull/1237
* CNVpytor container https://github.com/Clinical-Genomics/BALSAMIC/pull/1246
* Restored balsamic container in UMI concatenation rule https://github.com/Clinical-Genomics/BALSAMIC/pull/1261
* CNVpytor container, fixing numpy version https://github.com/Clinical-Genomics/BALSAMIC/pull/1273
* QC workflow store https://github.com/Clinical-Genomics/BALSAMIC/pull/1295

Removed:
^^^^^^^^
* Config folder https://github.com/Clinical-Genomics/BALSAMIC/pull/1175
* Quality trimming of fastqs for UMI workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/1176
* Balsamic container https://github.com/Clinical-Genomics/BALSAMIC/pull/1230
* Plugin CLI https://github.com/Clinical-Genomics/BALSAMIC/pull/1245
* Realignment step for TGA workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/1272
* Archived/outdated workflows and scripts https://github.com/Clinical-Genomics/BALSAMIC/pull/1296

[12.0.2]
--------

Fixed:
^^^^^^
* Missing `Number` in VCF header for SVs https://github.com/Clinical-Genomics/BALSAMIC/pull/1203

Changed:
^^^^^^^^
* Fix cyvcf2 to version 0.30.22 https://github.com/Clinical-Genomics/BALSAMIC/pull/1206
* Fix pydantic version (<2.0) https://github.com/Clinical-Genomics/BALSAMIC/pull/1206
* Update varcall-cnvkit container versions https://github.com/Clinical-Genomics/BALSAMIC/pull/1207

[12.0.1]
--------

Added:
^^^^^^
* WGS QC criteria for `PCT_PF_READS_IMPROPER_PAIRS` (condition: <= 0.1) https://github.com/Clinical-Genomics/BALSAMIC/pull/1164

Fixed:
^^^^^^
* Logged version of Delly (changing it to v1.0.3)  https://github.com/Clinical-Genomics/BALSAMIC/pull/1170

[12.0.0]
--------

Added:
^^^^^^
* PIP specific missing tools to config https://github.com/Clinical-Genomics/BALSAMIC/pull/1096
* Filtering script to remove normal variants from TIDDIT https://github.com/Clinical-Genomics/BALSAMIC/pull/1120
* Store TMB files in HK https://github.com/Clinical-Genomics/BALSAMIC/pull/1144

Changed:
^^^^^^^^
* Fixed all conda container dependencies https://github.com/Clinical-Genomics/BALSAMIC/pull/1096
* Changed --max_sv_size in VEP params to the size of chr1 for hg19 https://github.com/Clinical-Genomics/BALSAMIC/pull/1124
* Increased time-limit for sambamba_exon_depth and picard_markduplicates to 6 hours https://github.com/Clinical-Genomics/BALSAMIC/pull/1143
* Update cosmicdb to v97 https://github.com/Clinical-Genomics/BALSAMIC/pull/1147
* Updated read the docs with the changes relevant to mention https://github.com/Clinical-Genomics/BALSAMIC/pull/1153

Fixed:
^^^^^^
* Update cryptography version (39.0.1) due to security alert https://github.com/Clinical-Genomics/BALSAMIC/pull/1087
* Bump cryptography to v40.0.2 and gsutil to v5.23 https://github.com/Clinical-Genomics/BALSAMIC/pull/1154
* Pytest file saved in balsamic directory https://github.com/Clinical-Genomics/BALSAMIC/pull/1093
* Fix varcall_py3 container bcftools dependency error https://github.com/Clinical-Genomics/BALSAMIC/pull/1097
* AscatNgs container https://github.com/Clinical-Genomics/BALSAMIC/pull/1155

[11.2.0]
--------

Fixed:
^^^^^^
* Number of variants are increased with triallelic_site https://github.com/Clinical-Genomics/BALSAMIC/pull/1089

[11.1.0]
--------

Added:
^^^^^^
* Added somalier integration and relatedness check: https://github.com/Clinical-Genomics/BALSAMIC/pull/1017
* Cluster resources for CNVPytor tumor only https://github.com/Clinical-Genomics/BALSAMIC/pull/1083

Changed:
^^^^^^^^
* Parallelize download of reference files https://github.com/Clinical-Genomics/BALSAMIC/pull/1065
* Parallelize download of container images https://github.com/Clinical-Genomics/BALSAMIC/pull/1068

Fixed:
^^^^^^
* triallelic_site in quality filter for SNV https://github.com/Clinical-Genomics/BALSAMIC/pull/1052
* Compression of SNV, research and clinical, VCF files https://github.com/Clinical-Genomics/BALSAMIC/pull/1060
* `test_write_json` failing locally https://github.com/Clinical-Genomics/BALSAMIC/pull/1063
* Container build and push via github actions by setting buildx `provenance` flag to false https://github.com/Clinical-Genomics/BALSAMIC/pull/1071
* Added buildx to the submodule workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/1072
* Change user in somalier container to defaultuser https://github.com/Clinical-Genomics/BALSAMIC/pull/1080
* Reference files for hg38 https://github.com/Clinical-Genomics/BALSAMIC/pull/1081

[11.0.2]
--------

Changed:
^^^^^^^^
* Code owners https://github.com/Clinical-Genomics/BALSAMIC/pull/1050

Fixed:
^^^^^^
* MaxDepth in quality filter for SV https://github.com/Clinical-Genomics/BALSAMIC/pull/1051

[11.0.1]
--------

Fixed:
^^^^^^
* Incorrect raw `TNscope` VCF delivered https://github.com/Clinical-Genomics/BALSAMIC/pull/1042

[11.0.0]
--------

Added:
^^^^^^
* Use of PON reference, if exists for CNVkit tumor-normal analysis https://github.com/Clinical-Genomics/BALSAMIC/pull/982
* Added PON version to CLI and config.json https://github.com/Clinical-Genomics/BALSAMIC/pull/983
* `cnvpytor` to varcallpy3 container https://github.com/Clinical-Genomics/BALSAMIC/pull/991
* `cnvpytor` for tumor only workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/994
* R packages to cnvkit container https://github.com/Clinical-Genomics/BALSAMIC/pull/996
* Missing R packages to cnvkit container https://github.com/Clinical-Genomics/BALSAMIC/pull/997
* add rlang to cnvkit container https://github.com/Clinical-Genomics/BALSAMIC/pull/998
* AnnotSV and bedtools to annotate container https://github.com/Clinical-Genomics/BALSAMIC/pull/1005
* cosmicdb to TNscope for tumor only and tumor normal workflows https://github.com/Clinical-Genomics/BALSAMIC/pull/1006
* `loqusDB` dump files to the config through the balsamic config case CLI https://github.com/Clinical-Genomics/BALSAMIC/pull/992
* Pre-annotation quality filters for SNVs annd added `research` to output files https://github.com/Clinical-Genomics/BALSAMIC/pull/1007
* Annotation of snv_clinical_observations for somatic snv https://github.com/Clinical-Genomics/BALSAMIC/pull/1012
* Annotation of sv_clinical_observations  for somatic sv and SV CNV filter rules https://github.com/Clinical-Genomics/BALSAMIC/pull/1013
* Swegen SNV and SV frequency database for WGS https://github.com/Clinical-Genomics/BALSAMIC/pull/1014
* triallelic_sites and variants with MaxDepth to the VCFs https://github.com/Clinical-Genomics/BALSAMIC/pull/1021
* Clinical VCF for TGA workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/1024
* CNVpytor plots into the CNV PDF report https://github.com/Clinical-Genomics/BALSAMIC/pull/1023
* Research and clinical housekeeper tags https://github.com/Clinical-Genomics/BALSAMIC/pull/1023
* Cluster configuration for rules https://github.com/Clinical-Genomics/BALSAMIC/pull/1028
* Variant filteration using loqusDB and Swegen annotations https://github.com/Clinical-Genomics/BALSAMIC/pull/1029
* Annotation resources to readsthedocs https://github.com/Clinical-Genomics/BALSAMIC/pull/1031
* Delly CNV rules for TGA workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/103
* cnvpytor container and removed cnvpytor from varcallpy3 https://github.com/Clinical-Genomics/BALSAMIC/pull/1037

Changed:
^^^^^^^^
* Added version number to the PON reference filename (`.cnn`) https://github.com/Clinical-Genomics/BALSAMIC/pull/982
* Update `TIDDIT` to v3.3.0, `SVDB` to v2.6.4, `delly` to v1.1.3, `vcf2cytosure` to v0.8 https://github.com/Clinical-Genomics/BALSAMIC/pull/987
* toml config file for vcfanno https://github.com/Clinical-Genomics/BALSAMIC/pull/1012
* Split `vep_germline` rule into `tumor` and `normal` https://github.com/Clinical-Genomics/BALSAMIC/pull/1018
* Extract number of variants from clinical files https://github.com/Clinical-Genomics/BALSAMIC/pull/1022

Fixed:
^^^^^^
* Reverted `pandas` version (from `1.3.5` to `1.1.5`) https://github.com/Clinical-Genomics/BALSAMIC/pull/1018
* Mate in realigned bam file https://github.com/Clinical-Genomics/BALSAMIC/pull/1019
* samtools command in merge bam and names in toml for vcfanno https://github.com/Clinical-Genomics/BALSAMIC/pull/1020
* If statement in `vep_somatic_clinical_snv` rule https://github.com/Clinical-Genomics/BALSAMIC/pull/1022
* Invalid flag second of pair validation error https://github.com/Clinical-Genomics/BALSAMIC/pull/1025
* Invalid flag second of pair validation error using picardtools https://github.com/Clinical-Genomics/BALSAMIC/pull/1027
* Samtools command for mergetype tumor https://github.com/Clinical-Genomics/BALSAMIC/pull/1030
* `varcall_py3` container building https://github.com/Clinical-Genomics/BALSAMIC/pull/1036
* Picard and fastp commands params and cluster config for umi workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/1032
* Set channels in `varcall_py3` container https://github.com/Clinical-Genomics/BALSAMIC/pull/1035
* Delly command for tumor-normal analysis https://github.com/Clinical-Genomics/BALSAMIC/pull/1039
* tabix command in bcftools_quality_filter_TNscope_umi_tumor_only rule https://github.com/Clinical-Genomics/BALSAMIC/pull/1040

Removed:
^^^^^^^^
* case ID from the PON `.cnn` output file https://github.com/Clinical-Genomics/BALSAMIC/pull/983
* `TNhaplotyper` for paired WGS analysis https://github.com/Clinical-Genomics/BALSAMIC/pull/988
* `TNhaplotyper` for tumor only WGS analysis https://github.com/Clinical-Genomics/BALSAMIC/pull/1006
* `TNhaplotyper` for TGS https://github.com/Clinical-Genomics/BALSAMIC/pull/1022

[10.0.5]
--------

Changed:
^^^^^^^^
* Update `vcf2cytosure` version to v0.8 https://github.com/Clinical-Genomics/BALSAMIC/pull/1010
* Update GitHub action images to `ubuntu-20.04` https://github.com/Clinical-Genomics/BALSAMIC/pull/1010
* Update GitHub actions to their latest versions https://github.com/Clinical-Genomics/BALSAMIC/pull/1010

[10.0.4]
---------

Fixed:
^^^^^^
* Increase `sambamba_exon_depth` rule run time https://github.com/Clinical-Genomics/BALSAMIC/pull/1001

[10.0.3]
---------
Fixed:
^^^^^^

* Input VCF files for cnvkit rules, cnvkit command and container https://github.com/Clinical-Genomics/BALSAMIC/pull/995

[10.0.2]
---------

Fixed:
^^^^^^

* TIDDIT delivery rule names (undo rule name changes made in Balsamic 10.0.1) https://github.com/Clinical-Genomics/BALSAMIC/pull/977
* BALSAMIC readthedocs CLI documentation generation  https://github.com/Clinical-Genomics/BALSAMIC/issues/965

[10.0.1]
---------

Fixed:
^^^^^^

* Command and condition for TIDDIT and fixed ReadtheDocs https://github.com/Clinical-Genomics/BALSAMIC/pull/973
* ReadtheDocs and updated the header https://github.com/Clinical-Genomics/BALSAMIC/pull/973


Changed:
^^^^^^^^

* Time allocation in cluster configuration for SV rules https://github.com/Clinical-Genomics/BALSAMIC/pull/973



[10.0.0]
---------

Added:
^^^^^^

* New option `analysis-workflow` to balsamic config case CLI https://github.com/Clinical-Genomics/BALSAMIC/pull/932
* New python script to edit INFO tags in `vardict` and `tnscope_umi` VCF files https://github.com/Clinical-Genomics/BALSAMIC/pull/948
* Added `cyvcf2` and `click` tools to the `varcallpy3` container https://github.com/Clinical-Genomics/BALSAMIC/pull/948
* Delly TIDDIT and vcf2cytosure for WGS https://github.com/Clinical-Genomics/BALSAMIC/pull/947
* `Delly` `TIDDIT` `vcf2cytosure` and method to process SVs and CNVs for WGS https://github.com/Clinical-Genomics/BALSAMIC/pull/947
* SV and CNV analysis and `TIDDIT` to balsamic ReadtheDocs https://github.com/Clinical-Genomics/BALSAMIC/pull/951
* Gender to `config.json` https://github.com/Clinical-Genomics/BALSAMIC/pull/955
* Provided gender as input for `vcf2cyosure` https://github.com/Clinical-Genomics/BALSAMIC/pull/955
* SV CNV doc to balsamic READTHEDOCS https://github.com/Clinical-Genomics/BALSAMIC/pull/960
* Germline normal SNV VCF file header renaming to be compatible with genotype uploads https://github.com/Clinical-Genomics/BALSAMIC/issues/882
* Add tabix and gzip to vcf2cytosure container https://github.com/Clinical-Genomics/BALSAMIC/pull/969

Changed:
^^^^^^^^

* UMI-workflow for panel cases to be run only with `balsamic-umi` flag https://github.com/Clinical-Genomics/BALSAMIC/issues/896
* Update `codecov` action version to @v2 https://github.com/Clinical-Genomics/BALSAMIC/pull/941
* QC-workflow for panel cases to be run only with `balsamic-qc` https://github.com/Clinical-Genomics/BALSAMIC/pull/942
* `get_snakefile` function takes the argument `analysis_workflow` to trigger the QC workflow when necessary https://github.com/Clinical-Genomics/BALSAMIC/pull/942
* `bcftools_counts` input depending on `analysis_workflow` https://github.com/Clinical-Genomics/BALSAMIC/pull/942
* UMI output filename `TNscope_umi` is changed to `tnscope_umi` https://github.com/Clinical-Genomics/BALSAMIC/pull/948
* Update `delly` to v1.0.3 https://github.com/Clinical-Genomics/BALSAMIC/pull/950
* Update versions of `delly` in ReadtheDocs https://github.com/Clinical-Genomics/BALSAMIC/pull/951
* Provided gender as input for `ascat` and `cnvkit` https://github.com/Clinical-Genomics/BALSAMIC/pull/955
* Update QC criteria for panel and wgs analysis according to https://github.com/Clinical-Genomics/project-planning/issues/338#issuecomment-1132643330. https://github.com/Clinical-Genomics/BALSAMIC/pull/952
* For uploads to scout, increasing the number of variants failing threshold from 10000 to 50000 https://github.com/Clinical-Genomics/BALSAMIC/pull/952

Fixed:
^^^^^^

* GENOME_VERSION set to the different genome_version options and replaced with config["reference"]["genome_version"] https://github.com/Clinical-Genomics/BALSAMIC/pull/942
* `run_validate.sh` script https://github.com/Clinical-Genomics/BALSAMIC/pull/952
* Somatic SV tumor normal rules https://github.com/Clinical-Genomics/BALSAMIC/pull/959
* Missing `genderChr` flag for `ascat_tumor_normal` rule https://github.com/Clinical-Genomics/BALSAMIC/pull/963
* Command in vcf2cytosure rule and updated ReadtheDocs https://github.com/Clinical-Genomics/BALSAMIC/pull/966
* Missing name `analysis_dir` in QC.smk https://github.com/Clinical-Genomics/BALSAMIC/pull/970
* Remove `sample_type` wildcard from the `vcfheader_rename_germline` rule and change genotype file name https://github.com/Clinical-Genomics/BALSAMIC/pull/971

Removed
^^^^^^^

* Removed `qc_panel` config in favor of standard config https://github.com/Clinical-Genomics/BALSAMIC/pull/942
* Removed cli `--analysis_type` for `balsamic report deliver` command and `balsamic run analysis` https://github.com/Clinical-Genomics/BALSAMIC/pull/942
* Removed `analysis_type`: `qc_panel` and replace the trigger for QC workflow by `analysis_workflow`: `balsamic-qc` https://github.com/Clinical-Genomics/BALSAMIC/pull/942
* Outdated balsamic report files (`balsamic_report.html` & `balsamic_report.md`) https://github.com/Clinical-Genomics/BALSAMIC/pull/952

[9.0.1]
-------

Fixed:
^^^^^^

* Revert `csvkit` tool in align_qc container https://github.com/Clinical-Genomics/BALSAMIC/pull/928
* Automatic version update for balsamic methods https://github.com/Clinical-Genomics/BALSAMIC/pull/930

[9.0.0]
--------

Added:
^^^^^^

* Snakemake workflow to create canfam3 reference https://github.com/Clinical-Genomics/BALSAMIC/pull/843
* Call umi variants using TNscope in bed defined regions https://github.com/Clinical-Genomics/BALSAMIC/issues/821
* UMI duplication metrics to report in multiqc_picard_dups.json https://github.com/Clinical-Genomics/BALSAMIC/issues/844
* Option to use PON reference in cnv calling for TGA tumor-only cases https://github.com/Clinical-Genomics/BALSAMIC/pull/851
* QC default validation conditions (for not defined capture kits) https://github.com/Clinical-Genomics/BALSAMIC/pull/855
* SVdb to the varcall_py36 container https://github.com/Clinical-Genomics/BALSAMIC/pull/872
* SVdb to WGS workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/873
* Docker container for vcf2cytosure https://github.com/Clinical-Genomics/BALSAMIC/pull/869
* Snakemake rule for creating `.cgh` files from `CNVkit` outputs https://github.com/Clinical-Genomics/BALSAMIC/pull/880
* SVdb to TGA workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/879
* SVdb merge SV and CNV https://github.com/Clinical-Genomics/BALSAMIC/pull/886
* Readthedocs for BALSAMIC method descriptions https://github.com/Clinical-Genomics/BALSAMIC/pull/906
* Readthedocs for BALSAMIC variant filters for WGS somatic callers https://github.com/Clinical-Genomics/BALSAMIC/pull/906
* bcftools counts to varcall filter rules https://github.com/Clinical-Genomics/BALSAMIC/pull/899
* Additional WGS metrics to be stored in ``<case>_metrics_deliverables.yaml`` https://github.com/Clinical-Genomics/BALSAMIC/pull/907
* ascatNGS copynumber file https://github.com/Clinical-Genomics/BALSAMIC/pull/914
* ReadtheDocs for BALSAMIC annotation resources https://github.com/Clinical-Genomics/BALSAMIC/pull/916
* Delly CNV for tumor only workflow https://github.com/Clinical-Genomics/BALSAMIC/pull/923
* Delly CNV Read-depth profiles for tumor only workflows https://github.com/Clinical-Genomics/BALSAMIC/pull/924
* New metric to be extracted and validated: ``NUMBER_OF_SITES`` (``bcftools`` counts) https://github.com/Clinical-Genomics/BALSAMIC/pull/925

Changed:
^^^^^^^^

* Merge QC metric extraction workflows https://github.com/Clinical-Genomics/BALSAMIC/pull/833
* Changed the base-image for balsamic container to 4.10.3-alpine https://github.com/Clinical-Genomics/BALSAMIC/pull/869
* Updated SVdb to 2.6.0 https://github.com/Clinical-Genomics/BALSAMIC/pull/901
* Upgrade black to 22.3.0
* For UMI workflow, post filter `gnomad_pop_freq` value is changed from `0.005` to `0.02` https://github.com/Clinical-Genomics/BALSAMIC/pull/919
* updated delly to 0.9.1 https://github.com/Clinical-Genomics/BALSAMIC/pull/920
* container base_image (align_qc, annotate, coverage_qc, varcall_cnvkit, varcall_py36) to 4.10.3-alpine https://github.com/Clinical-Genomics/BALSAMIC/pull/921
* update container (align_qc, annotate, coverage_qc, varcall_cnvkit,varcall_py36) bioinfo tool versions  https://github.com/Clinical-Genomics/BALSAMIC/pull/921
* update tool versions (align_qc, annotate, coverage_qc, varcall_cnvkit) in methods and softwares docs https://github.com/Clinical-Genomics/BALSAMIC/pull/921
* Updated the list of files to be stored and delivered https://github.com/Clinical-Genomics/BALSAMIC/pull/915
* Moved ``collect_custom_qc_metrics`` rule from ``multiqc.rule`` https://github.com/Clinical-Genomics/BALSAMIC/pull/925

Fixed:
^^^^^^
* Automate balsamic version for readthedocs install page https://github.com/Clinical-Genomics/BALSAMIC/pull/888
* ``collect_qc_metrics.py`` failing for WGS cases with empty ``capture_kit`` argument https://github.com/Clinical-Genomics/BALSAMIC/pull/850
* QC metric validation for different panel bed version https://github.com/Clinical-Genomics/BALSAMIC/pull/855
* Fixed development version of ``fpdf2`` to ``2.4.6`` https://github.com/Clinical-Genomics/BALSAMIC/issues/878
* Added missing svdb index file https://github.com/Clinical-Genomics/BALSAMIC/issues/848

Removed
^^^^^^^

* ``--qc-metrics/--no-qc-metrics`` flag from the ``balsamic report deliver`` command https://github.com/Clinical-Genomics/BALSAMIC/pull/833
* Unused pon option for SNV calling with TNhaplotyper tumor-only https://github.com/Clinical-Genomics/BALSAMIC/pull/851
* SV and CNV callers from annotation and filtering https://github.com/Clinical-Genomics/BALSAMIC/pull/889
* vcfanno and COSMIC from SV annotation https://github.com/Clinical-Genomics/BALSAMIC/pull/891
* Removed `MSK_impact` and `MSK_impact_noStrelka` json files from config https://github.com/Clinical-Genomics/BALSAMIC/pull/903
* Cleanup of `strelka`, `pindel` , `mutect2` variables from BALSAMIC https://github.com/Clinical-Genomics/BALSAMIC/pull/903
* bcftools_stats from vep https://github.com/Clinical-Genomics/BALSAMIC/issues/898
* QC delivery report workflow (generating the ``<case>_qc_report.html`` file) https://github.com/Clinical-Genomics/BALSAMIC/issues/878
* ``--sample-id-map`` and ``--case-id-map`` flags from the ``balsamic report deliver`` command https://github.com/Clinical-Genomics/BALSAMIC/issues/878
* Removed `gatk_haplotypecaller` for reporting panel germline variants https://github.com/Clinical-Genomics/BALSAMIC/issues/918

[8.2.10]
--------

Added:
^^^^^^
* `libopenblas=0.3.20` dependency to annotate container for fixing bcftools #909

Fixes:
^^^^^^

* bcftools version locked at `1.10` #909

Changed:
^^^^^^^^
* base image of balsamic container to `4.10.3-alphine` #909
* Replaced annotate container tests with new code #909

Removed:
^^^^^^^^
* Removed failed `vcf2cytosure` installation from annotate container #909

[8.2.9]
-------

Added:
^^^^^^

* Added slurm qos tag `express` #885
* Included more text about UMI-workflow variant calling settings to the readthedocs #888
* Extend QCModel to include `n_base_limit` which outputs in config json `QC` dict

Fixes:
^^^^^^
* Automate balsamic version for readthedocs install page #888

Changed:
^^^^^^^^
* Upgrade black to 22.3.0
* fastp default setting of `n_base_limit` is changed to `50` from `5`

[8.2.8]
--------

Added:
^^^^^^
* Added the readthedocs page for BALSAMIC variant-calling filters #867
* Project requirements (setup.py) to build the docs #874
* Generate cram from umi-consensus called bam files #865

Changed:
^^^^^^^^
* Updated the bioinfo tools version numbers in BALSAMIC readthedocs #867
* Sphinx version fixed to <0.18 #874
* Sphinx GitHub action triggers only on master branch PRs
* VAF filter for reporting somatic variants (Vardict) is minimised to 0.7% from 1% #876

Fixes:
^^^^^^
* cyvcf2 mock import for READTHEDOCS environment #874

[8.2.7]
-------
Fixes:
^^^^^^
* Fixes fastqc timeout issues for wgs cases #861
* Fix cluster configuration for vep and vcfanno #857

[8.2.6]
-------

Fixes:
^^^^^^

* Set right qos in scheduler command #856

[8.2.5]
-------

* balsamic.sif container installation during cache generation #841

Fixed:
^^^^^^

* Execution of `create_pdf` python script inside the balsamic container #841

[8.2.4]
-------

Added:
^^^^^^

* ``--hgvsg`` annotation to VEP #830
* ``ascatNgs`` PDF delivery (plots & statistics) #828

[8.2.3]
-------
Fixed:
^^^^^^

* Add default for gender if ``purecn`` captures dual gender values #824

Changed:
^^^^^^^^
* Updated ``purecn`` and its dependencies to latest versions

[8.2.2]
-------
Added:
^^^^^^

* ``ascatNGS`` tumor normal delivery #810

Changed:
^^^^^^^^
* QC metrics delivery tag #820
* Refactor tmb rule that contains redundant line #817

[8.2.1]
-------

Fixed:
^^^^^^

* ``cnvkit`` gender comparison operator bug #819

[8.2.0]
-------

Added:
^^^^^^

* Added various basic filters to all variant callers irregardless of their delivery status #750
* BALSAMIC container #728
* BALSAMIC reference generation via cluster submission for both reference and container #686
* Container specific tests #770
* BALSAMIC quality control metrics extraction and validation #754
* Delly is added as a submodule and removed from rest of the conda environments #787
* Store research VCFs for all filtered and annotated VCF files
* Added `.,PASS` to all structural variant filter rules to resolve the issues with missing calls in filtered file
* Handling of QC metrics validation errors #783
* Github Action workflow that builds the docs using Sphinx #809
* Zenodo integration to create citable link #813
* Panel BED specific QC conditions #800
* Metric extraction to a YAML file for Vogue #802

Changed:
^^^^^^^^

* refactored main workflow with more readible organization #614
* refactored conda envs within container to be on base and container definition is uncoupled #759
* renamed umi output file names to fix issue with picard HSmetrics #804
* locked requirements for graphviz io 0.16 #811
* QC metric validation is performed across all metrics of each of the samples #800

Removed:
^^^^^^^^

* The option of running umiworkflow independently with balsamic command-line option "-a umi"
* Removed source activate from reference and pon workflows #764

Fixed:
^^^^^^

* Pip installation failure inside balsamic container #758
* Fixed issue #768 with missing ``vep_install`` command in container
* Fixed issue #765 with correct input bam files for SV rules
* Continuation of CNVkit even if ``PURECN`` fails and fix ``PureCN`` conda paths #774 #775
* Locked version for ``cryptography`` package
* Bumped version for ``bcftools`` in cnvkit container
* Fixed issues #776 and #777 with correct install paths for gatk and manta
* Fixed issue #782 for missing AF in the vcf INFO field
* Fixed issues #748 #749 with correct sample names
* Fixed issue #767 for ascatngs hardcoded values
* Fixed missing output option in bcftools filters for tnhaplotyper #793
* Fixed issue #795 with increasing resources for vep and filter SV prior to vep
* Building ``wheel`` for ``cryptography`` bug inside BALSAMIC container #801
* Fixed badget for docker container master and develop status
* ReadtheDocs building failure due to dependencies, fixed by locking versions #773
* Dev requirements installation for Sphinx docs (Github Action) #812
* Changed path for main Dockerfile version in ``.bumpversion.cfg``

[8.1.0]
-------

Added:
^^^^^^

* Workflow to check PR tiltes to make easier to tell PR intents #724
* ``bcftools stats``  to calculate Ti/Tv for all post annotate germline and somatic calls #93
* Added reference download date to ``reference.json`` #726
* ``ascatngs`` hg38 references to constants #683
* Added ClinVar as a source to download and to be annotated with VCFAnno #737

Changed:
^^^^^^^^

* Updated docs for git FAQs #731
* Rename panel of normal filename Clinical-Genomics/cgp-cancer-cnvcall#10


Fixed:
^^^^^^

* Fixed bug with using varcall_py36 container with VarDict #739
* Fixed a bug with VEP module in MultiQC by excluding #746
* Fixed a bug with ``bcftools stats`` results failing in MultiQC #744

[8.0.2]
-------

Fixed:
^^^^^^

* Fixed breaking shell command for VEP annotation rules #734

[8.0.1]
-------

Fixed:
^^^^^^

* Fixed context for Dockerfile for release content #720

[8.0.0]
-------

Added:
^^^^^^

* ``samtools`` flagstats and stats to workflow and MultiQC
* ``delly v0.8.7`` somatic SV caller #644
* ``delly`` containter #644
* ``bcftools v1.12`` to ``delly`` container #644
* ``tabix v0.2.6`` to ``delly`` container #644
* Passed SV calls from Manta to clinical delivery
* An extra filter to VarDict tumor-normal to remove variants with STATUS=Germline, all other will still be around
* Added ``vcf2cytosure`` to annotate container
* ``git`` to the container definition
* prepare_delly_exclusion rule
* Installation of ``PureCN`` rpackage in ``cnvkit`` container
* Calculate tumor-purity and ploidy using ``PureCN`` for ``cnvkit`` call
* ``ascatngs`` as a submodule #672
* GitHub action to build and test ``ascatngs`` container
* Reference section to ``docs/FAQ.rst``
* ``ascatngs`` download references from reference_file repository #672
* ``delly`` tumor only rule #644
* ``ascatngs`` download container #672
* Documentation update on setting sentieon env variables in ``bashrc``
* ``ascatngs`` tumor normal rule for wgs cases #672
* Individual rules (i.e. ngs filters) for cnv and sv callers. Only Manta will be delivered and added to the list of output files. #708
* Added "targeted" and "wgs" tags to variant callers to provide another layer of separation. #708
* ``manta`` convert inversion #709
* Sentieon version to bioinformatic tool version parsing #685
* added ``CITATION.cff`` to cite BALSAMIC


Changed:
^^^^^^^^

* Upgrade to latest sentieon version 202010.02
* New name ``MarkDuplicates`` to ``picard_markduplicates`` in ``bwa_mem`` rule and ``cluster.json``
* New name rule ``GATK_contest`` to ``gatk_contest``
* Avoid running pytest github actions workflow on ``docs/**`` and ``CHANGELOG.rst`` changes
* Updated ``snakemake`` to ``v6.5.3`` #501
* Update ``GNOMAD`` URL
* Split Tumor-only ``cnvkit batch`` into individual commands
* Improved TMB calculation issue #51
* Generalized ascat, delly, and manta result in workflow. #708
* Generalized workflow to eliminate duplicate entries and code. #708
* Split Tumor-Normal ``cnvkit batch`` into individual commands
* Moved params that are used in multiple rules to constants #711
* Changed the way conda and non-conda bioinfo tools version are parsed
* Python code formatter changed from Black to YAPF #619


Fixed:
^^^^^^

* post-processing of the umi consensus in handling BI tags
* vcf-filtered-clinical tag files will have all variants including PASS
* Refactor snakemake ``annotate`` rules according to snakemake etiquette #636
* Refactor snakemake ``align`` rules according to snakemake etiquette #636
* Refactor snakemake ``fastqc`` ``vep`` contest and ``mosdepth`` rules according to ``snakemake`` etiquette #636
* Order of columns in QC and coverage report issue #601
* ``delly`` not showing in workflow at runtime #644
* ``ascatngs`` documentation links in ``FAQs`` #672
* ``varcall_py36`` container build and push #703
* Wrong spacing in reference json issue #704
* Refactor snakemake ``quality control`` rules according to snakemake etiquette #636

Removed:
^^^^^^^^

* Cleaned up unused container definitions and conda environment files
* Remove cnvkit calling for WGS cases
* Removed the install.sh script

[7.2.5]
-------

Changed:
^^^^^^^^

* Updated COSMIC path to use version 94

[7.2.5]
-------

Changed:
^^^^^^^^

* Updated path for gnomad and 1000genomes to a working path from Google Storage

[7.2.4]
-------

Changed:
^^^^^^^^

* Updated sentieon util sort in umi to use Sentieon 20201002 version

[7.2.3]
-------

Fixed:
^^^^^^

* Fixed memory issue with vcfanno in vep_somatic rule fixes #661

[7.2.2]
-------

Fixed:
^^^^^^

* An error with Sentieon for better management of memory fixes #621

[7.2.1]
-------

Changed:
^^^^^^^^

* Rename Github actions to reflect their content

[7.2.0]
-------

Added:
^^^^^^

* Changelog reminder workflow to Github
* Snakemake workflow for created PON reference
* Balsamic cli config command(pon) for creating json for PON analysis
* tumor lod option for passing tnscope-umi final variants
* Git guide to make balsamic release in FAQ docs

Changed:
^^^^^^^^

* Expanded multiqc result search dir to whole analysis dir
* Simple test for docker container

Fixed:
^^^^^^

* Correctly version bump for Dockerfile

Removed:
^^^^^^^^

* Removed unused Dockerfile releases
* Removed redundant genome version from ``reference.json``

[7.1.10]
--------

Fixed:
^^^^^^

* Bug in ``ngs_filter`` rule set for tumor-only WGS
* Missing delivery of tumor only WGS filter

[7.1.9]
-------


Changed:
^^^^^^^^

* only pass variants are not part of delivery anymore
* delivery tag file ids are properly matched with sample_name
* tabix updated to 0.2.6
* fastp updated to 0.20.1
* samtools updated to 1.12
* bedtools updated to 2.30.0

Removed:
^^^^^^^^

* sentieon-dedup rule from delivery
* Removed all pre filter pass from delivery


[7.1.8]
-------

Fixed:
^^^^^^

* Target coverage (Picard HsMetrics) for UMI files is now correctly calculated.

Changed:
^^^^^^^^


* TNscope calculated AF values are fetched and written to AFtable.txt.

[7.1.7]
-------

Added:
^^^^^^

* ngs_filter_tnscope is also part of deliveries now

Changed:
^^^^^^^^

* rankscore is now a research tag instead of clinical
* Some typo and fixes in the coverage and constant metrics
* Delivery process is more verbose

Fixed:
^^^^^^

* CNVKit output is now properly imported in the deliveries and workflow

[7.1.6]
-------

Fixed:
^^^^^^

* CSS style for qc coverage report is changed to landscape

[7.1.5]
-------

Changed:
^^^^^^^^

* update download url for 1000genome WGS sites from ftp to http

[7.1.4]
-------

Changed:
^^^^^^^^

* bump picard to version 2.25.0

[7.1.3]
-------

Fixed:
^^^^^^

* ``assets`` path is now added to bind path

[7.1.2]
-------

Fixed:
^^^^^^

* umi_workflow config json is set as true for panel and wgs as false.
* Rename umiconsensus bam file headers from {samplenames} to TUMOR/NORMAL.
* Documentation autobuild on RTFD


[7.1.1]
-------

Fixed:
^^^^^^

* Moved all requirements to setup.py, and added all package_data there. Clean up unused files.

[7.1.0]
-------

Removed
^^^^^^^

* ``tnsnv`` removed from WGS analysis, both tumor-only and tumor-normal
* GATK-BaseRecalibrator is removed from all workflows

Fixed
^^^^^

* Fixed issue 577 with missing ``tumor.merged.bam`` and ``normal.merged.bam``
* Issue 448 with lingering tmp_dir. It is not deleted after analysis is properly finished.

Changed
^^^^^^^

* All variant calling rules use proper ``tumor.merged.bam`` or ``normal.merged.bam`` as inputs

[7.0.2]
-------

Added
^^^^^

* Updated docs with FAQ for UMI workflow

Fixed
^^^^^

* fix job scheduling bug for benchmarking
* rankscore's output is now a proper vcf.gz file
* Manta rules now properly make a sample_name file


[7.0.1]
-------

Added
^^^^^

* github action workflow to autobuild release containers


[7.0.0]
-------

Added
^^^^^

* ``balsamic init`` to download reference and related containers done in PRs #464 #538
* ``balsamic config case`` now only take a cache path instead of container and reference #538
* UMI workflow added to main workflow in series of PRs #469 #477 #483 #498 #503 #514 #517
* DRAGEN for WGS applications in PR #488
* A framework for QC check PR #401
* ``--quiet``` option for ``run analysis`` PR #491
* Benchmark SLURM jobs after the analysis is finished PR #534
* One container per conda environment (i.e. decouple containers) PR #511 #525 #522
* ``--disable-variant-caller`` command for ``report deliver`` PR #439
* Added genmod and rankscore in series of two PRs #531 and #533
* Variant filtering to Tumor-Normal in PR #534
* Split SNV/InDels and SVs from TNScope variant caller PR #540
* WGS Tumor only variant filters added in PR #548

Changed
^^^^^^^

* Update Manta to 1.6.0 PR #470
* Update FastQC to 0.11.9 PR #532
* Update BCFTools to 1.11 PR #537
* Update Samtools to 1.11 PR #537
* Increase resources and runtime for various workflows in PRs #482
* Python package dependenicies versions fixed in PR #480
* QoL changes to workflow in series of PR #471
* Series of documentation updates in PRs #489 #553
* QoL changes to scheduler script PR #491
* QoL changes to how temporary directories are handlded PR #516
* TNScope model apply rule merged with TNScope variant calling for tumor-normal in WGS #540
* Decoupled ``fastp`` rule into two rules to make it possible to use it for UMI runs #570


Fixed
^^^^^

* A bug in Manta variant calling rules that didn't name samples properly to TUMOR/NORMAL in the VCF file #572


[6.1.2]
-------

Changed
^^^^^^^
* Changed hk delivery tag for coverage-qc-report


[6.1.1]
-------

Fixed
^^^^^

* No UMI trimming for WGS applications #486
* Fixed a bug where BALSAMIC was checking for sacct/jobid file in local mode PR #497
* ``readlink`` command in ``vep_germline``, ``vep_somatic``, ``split_bed``, and ``GATK_popVCF`` #533
* Fix various bugs for memory handling of Picardtools and its executable in PR #534
* Fixed various issues with ``gsutils`` in PR #550

Removed
^^^^^^^

* ``gatk-register`` command removed from installing GATK PR #496

[6.1.1]
-------

* Fixed a bug with missing QC templates after ``pip install``


[6.1.0]
-------

Added
^^^^^
* CLI option to expand report generation for TGA and WES runs. Please see ``balsamic report deliver --help``
* BALSAMIC now generates a custom HTML report for TGA and WES cases.


[6.0.4]
-------

Changed
^^^^^^^

* Reduces MQ cutoff from 50 to 40 to only remove obvious artifacts PR #535
* Reduces AF cutoff from 0.02 to 0.01 PR #535

[6.0.3]
-------

Added
^^^^^

* ``config case`` subcommand now has ``--tumor-sample-name`` and ``--normal-sample-name``

Fixed
^^^^^

* Manta resource allocation is now properly set PR #523
* VarDict resource allocation in cluster.json increased (both core and time allocation) PR #523
* minimum memory request for GATK mutect2 and haplotypecaller is removed and max memory increased PR #523

[6.0.2]
-------

Added
^^^^^

* Document for Snakemake rule grammar PR #489


Fixed
^^^^^

* removed ``gatk3-register`` command from Dockerfile(s) PR #508


[6.0.1]
-------

Added
^^^^^
* A secondary path for latest jobids submitted to cluster (slurm and qsub) PR #465

[6.0.0]
-------

Added
^^^^^
* UMI workflow using Sentieon tools. Analysis run available via `balsamic run analysis --help` command. PR #359
* VCFutils to create VCF from flat text file. This is for internal purpose to generate validation VCF. PR #349
* Download option for hg38 (not validated) PR #407
* Option to disable variant callers for WES runs. PR #417

Fixed
^^^^^
* Missing cyvcf2 dependency, and changed conda environment for base environment PR #413
* Missing numpy dependency PR #426

Changed
^^^^^^^
* COSMIC db for hg19 updated to v90 PR #407
* Fastp trimming is now a two-pass trimming and adapter trimming is always enabled. This might affect coverage slightly PR #422
* All containers start with a clean environment #425
* All Sentieon environment variables are now added to config when workflow executes #425
* Branching model will be changed to gitflow

[5.1.0]
-------

Fixed
^^^^^
* Vardict-java version fixed. This is due to bad dependency and releases available on conda. Anaconda is not yet update with vardict 1.8, but vardict-java 1.8 is there. This causes various random breaks with Vardict's TSV output. #403

Changed
^^^^^^^
* Refactored Docker files a bit, preparation for decoupling #403

Removed
^^^^^^^
* In preparation for GATK4, IndelRealigner is removed #404


[5.0.1]
-------

Added
^^^^^
* Temp directory for various rules and workflow wide temp directory #396

Changed
^^^^^^^
* Refactored tags for housekeeper delivery to make them unique #395
* Increased core requirements for mutect2 #396
* GATK3.8 related utils run via jar file instead of gatk3 #396


[5.0.0]
-------

Added
^^^^^
* Config.json and DAG draph included in Housekeeper report #372
* New output names added to cnvkit_single and cnvkit_paired #372
* New output names added to vep.rule #372
* Delivery option to CLI and what to delivery with delivery params in rules that are needed to be delivered #376
* Reference data model with validation #371
* Added container path to install script #388

Changed
^^^^^^^
* Delivery file format simplified #376
* VEP rules have "all" and "pass" as output #376
* Downloaded reference structure changed #371
* genome/refseq.flat renamed to genome/refGene.flat #371
* reverted CNVKit to version 0.9.4 #390

Fixed
^^^^^
* Missing pygments to requirements.txt to fix travis CI #364
* Wildcard resolve for deliveries of vep_germline #374
* Missing index file from deliverables #383
* Ambiguous deliveries in vep_somatic and ngs_filters #387
* Updated documentation to match with installation #391

Removed
^^^^^^^
* Temp files removed from list of outputs in vep.rule #372
* samtools.rule and merged it with bwa_mem #375


[4.5.0]
-------

Added
^^^^^
* Models to build config case JSON. The models and descriptions of their contents can now be found
  in BALSAMIC/utils/models.py
* Added analysis_type to `report deliver` command
* Added report and delivery capability to Alignment workflow
* run_validate.sh now has -d to handle path to analysis_dir (for internal use only) #361

Changed
^^^^^^^

* Fastq files are no longer being copied as part of creation of the case config file.
  A symlink is now created at the destination path instead
* Config structure is no longer contained in a collestion of JSON files.
  The config models are now built using Pydantic and are contained in BALSAMIC/utils/models.py

Removed
^^^^^^^

* Removed command line option "--fastq-prefix" from config case command
* Removed command line option "--config-path" from config case command.
  The config is now always saved with default name "case_id.json"
* Removed command line option "--overwrite-config" from config-case command
  The command is now always executed with "--overwrite-config True" behavior

Refactored
^^^^^^^^^^

* Refactored BALSAMIC/commands/config/case.py:
  Utility functions are moved to BALSAMIC/utils/cli.py
  Models for config fields can be found at BALSAMIC/utils/models.py
  Context aborts and logging now contained in pilot function
  Tests created to support new architecture
* Reduce analysis directory's storage

Fixed
^^^^^
* Report generation warnings supressed by adding workdirectory
* Missing tag name for germline annotated calls #356
* Bind path is not added as None if analysis type is wgs #357
* Changes vardict to vardict-java #361


[4.4.0]
-------

Added
^^^^^

* pydantic to validate various models namely variant caller filters

Changed
^^^^^^^

* Variant caller filters moved into pydantic
* Install script and setup.py
* refactored install script with more log output and added a conda env suffix option
* refactored docker container and decoupled various parts of the workflow


[4.3.0]
-------


Added
^^^^^

* Added cram files for targeted sequencing runs fixes #286
* Added `mosdepth` to calculate coverage for whole exome and targeted sequencing
* Filter models added for tumor-only mode
* Enabling adapter trim enables pe adapter trim option for fastp
* Annotate germline variant calls
* Baitset name to picard hsmetrics

Deprecated
^^^^^^^^^^

* Sambamba coverage and rules will be deprecated

Fixed
^^^^^

* Fixed latest tag in install script
* Fixed lack of naming final annotated VCF TUMOR/NORMAL


Changed
^^^^^^^

* Increased run time for various slurm jobs fixes #314
* Enabled SV calls for VarDict tumor-only
* Updated `ensembl-vep` to v100.2

[4.2.4]
-------


Fixed
^^^^^

* Fixed sort issue with bedfiles after 100 slop


[4.2.3]
-------

Added
^^^^^


* Added Docker container definition for release and bumpversion

Changed
^^^^^^^


* Quality of life change to rtfd docs

Fixed
^^^^^


* Fix Docker container with faulty git checkout

[4.2.2]
-------

Added
^^^^^


* Add "SENTIEON_TMPDIR" to wgs workflow

[4.2.1]
-------

Changed
^^^^^^^


* Add docker container pull for correct version of install script

[4.2.0]
-------

Added
^^^^^


* CNV output as VCF
* Vep output for PASSed variants
* Report command with status and delivery subcommands

Changed
^^^^^^^


* Bed files are slopped 100bp for variant calling fix #262
* Disable vcfmerge
* Picard markduplicate output moved from log to output
* Vep upgraded to 99.1
* Removed SVs from vardict
* Refactored delivery plugins to produce a file with list of output files from workflow
* Updated snakemake to 5.13

Fixed
^^^^^


* Fixed a bug where threads were not sent properly to rules

Removed
^^^^^^^


* Removed coverage annotation from mutect2
* Removed source deactivate from rules to suppress conda warning
* Removed ``plugins delivery`` subcommand
* Removed annotation for germline caller results

[4.1.0]
-------

Added
^^^^^


* VEP now also produces a tab delimited file
* CNVkit rules output genemetrics and gene break file
* Added reference genome to be able to calculate AT/CG dropouts by Picard
* coverage plot plugin part of issue #75
* callable regions for CNV calling of tumor-only

Changed
^^^^^^^


* Increased time for indel realigner and base recalib rules
* decoupled vep stat from vep main rule
* changed qsub command to match UGE
* scout plugin updated

Fixed
^^^^^


* WGS qc rules - updated with correct options
  (picard - CollectMultipleMetrics, sentieon - CoverageMetrics)
* Log warning if WES workflow cannot find SENTIEON* env variables
* Fixes issue with cnvkit and WGS samples #268
* Fix #267 coverage issue with long deletions in vardict

[4.0.1] - 2019-11-08
--------------------

Added
^^^^^


* dependencies for workflow report
* sentieon variant callers germline and somatic for wes cases

Changed
^^^^^^^


* housekeeper file path changed from basename to absolute
* scout template for sample location changed from delivery_report to scout
* rule names added to benchmark files

[4.0.0] - 2019-11-04
--------------------

SGE qsub support release

Added
^^^^^


* ``install.sh`` now also downloads latest container
* Docker image for balsamic as part of ci
* Support for qsub alongside with slurm on ``run analysis --profile``

Changed
^^^^^^^


* Documentation updated
* Test fastq data and test panel bed file with real but dummy data

[3.3.1] - 2019-10-28
--------------------

Fixed
^^^^^


* Various links for reference genome is updated with working URL
* Config reference command now print correct output file

[3.3.0] - 2019-10-24
--------------------

somatic vcfmerge release

Added
^^^^^


* QC metrics for WGS workflow
* refGene.txt download to reference.json and reference workflow
* A new conda environment within container
* A new base container built via Docker (centos7:miniconda3_4_6_14)
* VCFmerge package as VCF merge rule (https://github.com/hassanfa/VCFmerge)
* A container for develop branch
* Benchmark rules to variant callers

Changed
^^^^^^^


* SLURM resource allocation for various variancalling rules optimized
* mergetype rule updated and only accepts one single tumor instead of multiple

[3.2.3] - 2019-10-24
--------------------

Fixed
^^^^^


* Removed unused output files from cnvkit which caused to fail on targetted analysis

[3.2.2] - 2019-10-23
--------------------

Fixed
^^^^^


* Removed target file from cnvkit batch

[3.2.1] - 2019-10-23
--------------------

Fixed
^^^^^


* CNVkit single missing reference file added

[3.2.0] - 2019-10-11
--------------------

Adds:
^^^^^


* CNVkit to WGS workflow
* get_thread for runs

Changed:
^^^^^^^^


* Optimized resources for SLURM jobs

Removed:
^^^^^^^^


* Removed hsmetrics for non-mark duplicate bam files

[3.1.4] - 2019-10-08
--------------------

Fixed
^^^^^


* Fixes a bug where missing capture kit bed file error for WGS cases

[3.1.3] - 2019-10-07
--------------------

Fixed
^^^^^


* benchmark path bug issue #221

[3.1.2] - 2019-10-07
--------------------

Fixed
^^^^^


* libreadline.so.6 symlinking and proper centos version for container

[3.1.1] - 2019-10-03
--------------------

Fixed
^^^^^


* Proper tag retrieval for release
  ### Changed
* BALSAMIC container change to latest and version added to help line

[3.1.0] - 2019-10-03
--------------------

TL;DR:


* QoL changes to WGS workflow
* Simplified installation by moving all tools to a container

Added
^^^^^


* Benchmarking using psutil
* ML variant calling for WGS
* ``--singularity`` option to ``config case`` and ``config reference``

Fixed
^^^^^


* Fixed a bug with boolean values in analysis.json

Changed
^^^^^^^


* ``install.sh`` simplified and will be depricated
* Singularity container updated
* Common somatic and germline variant callers are put in single file
* Variant calling workflow and analysis config files merged together

Removed
^^^^^^^


* ``balsamic install`` is removed
* Conda environments for py36 and py27 are removed

[3.0.1] - 2019-09-11
--------------------

Fixed
^^^^^


* Permissions on ``analysis/qc`` dir are 777 now

[3.0.0] - 2019-09-05
--------------------

This is major release.
TL;DR:


* Major changes to CLI. See documentation for updates.
* New additions to reference generation and reference config file generation and complete overhaul
* Major changes to reposityory structure, conda environments.

Added
^^^^^


* Creating and downloading reference files: ``balsamic config reference`` and ``balsamic run reference``
* Container definitions for install and running BALSAMIC
* Bunch of tests, setup coveralls and travis.
* Added Mutliqc, fastp to rule utilities
* Create Housekeeper and Scout files after analysis completes
* Added Sentieon tumor-normal and tumor only workflows
* Added trimming option while creating workflow
* Added multiple tumor sample QC analysis
* Added pindle for indel variant calling
* Added Analysis finish file in the analysis directory

Fixed
^^^^^


* Multiple fixes to snakemake rules

Changed
^^^^^^^


* Running analysis through: ``balsamic run analysis``
* Cluster account and email info added to ``balsamic run analysis``
* ``umi`` workflow through ``--umi`` tag. [workflow still in evaluation]
* ``sample-id`` replaced by ``case-id``
* Plan to remove FastQC as well

Removed
^^^^^^^


* ``balsamic config report`` and ``balsamic report``
* ``sample.config`` and ``reference.json`` from config directory
* Removed cutadapt from workflows

[2.9.8] - 2019-01-01
--------------------

Fixed
^^^^^


* picard hsmetrics now has 50000 cov max
* cnvkit single wildcard resolve bug fixed

[2.9.7] - 2019-02-28
--------------------

Fixed
^^^^^


* Various fixes to umi_single mode
* analysis_finish file does not block reruns anymore
* Added missing single_umi to analysis workflow cli

Changed
^^^^^^^


* vardict in single mode has lower AF threshold filter (0.005 -> 0.001)

[2.9.6] - 2019-02-25
--------------------

Fixed
^^^^^


* Reference to issue #141, fix for 3 other workflows
* CNVkit rule update for refflat file

[2.9.5] - 2019-02-25
--------------------

Added
^^^^^


* An analysis finish file is generated with date and time inside (%Y-%M-%d T%T %:z)

[2.9.4] - 2019-02-13
--------------------

Fixed
^^^^^


* picard version update to 2.18.11 github.com/hassanfa/picard

[2.9.3] - 2019-02-12
--------------------

Fixed
^^^^^


* Mutect single mode table generation fix
* Vardict single mode MVL annotation fix

[2.9.2] - 2019-02-04
--------------------

Added
^^^^^


* CNVkit single sample mode now in workflow
* MVL list from cheng et al. 2015 moved to assets

[2.9.1] - 2019-01-22
--------------------

Added
^^^^^


* Simple table for somatic variant callers for single sample mode added

Fixed
^^^^^


* Fixes an issue with conda that unset variables threw an error issue #141

[2.9.0] - 2019-01-04
--------------------

Changed
^^^^^^^


* Readme structure and example
* Mutect2's single sample output is similar to paired now
* cli path structure update

Added
^^^^^


* test data and sample inputs
* A dag PDF will be generated when config is made
* umi specific variant calling

[2.8.1] - 2018-11-28
--------------------

Fixed
^^^^^


* VEP's perl module errors
* CoverageRep.R now properly takes protein_coding transcatipts only

[2.8.0] - 2018-11-23
--------------------

UMI single sample align and QC

Added
^^^^^


* Added rules and workflows for UMI analysis: QC and alignment

[2.7.4] - 2018-11-23
--------------------

Germline single sample

Added
^^^^^


* Germline single sample addition
  ### Changed
* Minor fixes to some rules to make them compatible with tumor mode

[2.7.3] - 2018-11-20
--------------------

Fixed
^^^^^


* Various bugs with DAG to keep popvcf and splitbed depending on merge bam file
* install script script fixed and help added

[2.7.2] - 2018-11-15
--------------------

Changed
^^^^^^^


* Vardict, Strelka, and Manta separated from GATK best practice pipeline

[2.7.1] - 2018-11-13
--------------------

Fixed
^^^^^


* minro bugs with strelka_germline and freebayes merge
  ### Changed
* removed ERC from haplotypecaller

[2.7.0] - 2018-11-08
--------------------

Germline patch

Added
^^^^^


* Germline caller tested and added to the paired analysis workflow: Freebayes, HaplotypeCaller, Strelka, Manta

Changed
^^^^^^^


* Analysis config files updated
* Output directory structure changed
* vep rule is now a single rule
* Bunch of rule names updated and shortened, specifically in Picard and GATK
* Variant caller rules are all updated and changed
* output vcf file names are now more sensible: {SNV,SV}.{somatic,germline}.sampleId.variantCaller.vcf.gz
* Job limit increased to 300

Removed
^^^^^^^


* removed bcftools.rule for var id annotation

Changed
^^^^^^^

Fixed
^^^^^

[2.6.3] - 2018-11-01
--------------------

Changed
^^^^^^^


* Ugly and godforsaken ``runSbatch.py`` is now dumping sacct files with job IDs. Yikes!

[2.6.2] - 2018-10-31
--------------------

Fixed
^^^^^


* added ``--fastq-prefix`` option for ``config sample`` to set fastq prefix name. Linking is not changed.

[2.6.1] - 2018-10-29
--------------------

Fixed
^^^^^


* patched a bug for copying results for strelka and manta which was introduced in ``2.5.0``

[2.5.0] - 2018-10-22
--------------------

Changed
^^^^^^^


* ``variant_panel`` changed to ``capture_kit``
* sample config file takes balsamic version
* bioinfo tool config moved bioinfotool to cli_utils from ``config report``

Added
^^^^^


* bioinfo tool versions is now added to analysis config file

[2.4.0] - 2018-10-22
--------------------

Changed
^^^^^^^


* ``balsamic run`` has 3 stop points: paired variant calling, single mode variant calling, and QC/Alignment mode.
* ``balsamic run [OPTIONS] -S ...`` is depricated, but it supersedes ``analysis_type`` mode if provided.

[2.3.3] - 2018-10-22
--------------------

Added
^^^^^


* CSV output for variants in each variant caller based on variant filters
* DAG image of workflow
  ### Changed
* Input for variant filter has a default value
* ``delivery_report`` is no created during config generation
* Variant reporter R script cmd updated in ``balsamic report``

[2.3.2] - 2018-10-19
--------------------

Changed
^^^^^^^


* Fastq files are now always linked to ``fastq`` directory within the analysis directory

Added
^^^^^


* ``balsamic config sample`` now accepts individual files and paths. See README for usage.

[2.3.1] - 2018-09-25
--------------------

Added
^^^^^


* CollectHSmetric now run twice for before and after markduplicate

[2.3.0] - 2018-09-25
--------------------

Changed
^^^^^^^


* Sample config file now includes a list of chromosomes in the panel bed file

Fixed
^^^^^


* Non-matching chrom won't break the splitbed rule anymore
* collectqc rules now properly parse tab delimited metric files

[2.2.0] - 2018-09-11
--------------------

Added
^^^^^


* Coverage plot to report
* target coverage file to report json
* post-cutadapt fastqc to collectqc
* A header to report pdf
* list of bioinfo tools used in the analysis added to report
  ### Changed
* VariantRep.R now accepts multiple inputs for each parameter (see help)
* AF values for MSKIMPACT config
  ### Fixed
* Output figure for coverageplot is now fully square :-)

[2.1.0] - 2018-09-11
--------------------

Added
^^^^^


* normalized coverage plot script
* fastq file IO check for config creation
* added qos option to ``balsamic run``
  ### Fixed
* Sambamba depth coverage parameters
* bug with picard markduplicate flag

[2.0.2] - 2018-09-11
--------------------

Added
^^^^^


* Added qos option for setting qos to run jobs with a default value of low

[2.0.1] - 2018-09-10
--------------------

Fixed
^^^^^


* Fixed package dependencies with vep and installation

[2.0.0] - 2018-09-05
--------------------

Variant reporter patch and cli update

Added
^^^^^


* Added ``balsamic config sample`` and ``balsamic config report`` to generate run analysis and reporting config
* Added ``VariantRep.R`` script to information from merged variant table: variant summry, TMB, and much more
* Added a workflow for single sample mode alignment and QC only
* Added QC skimming script to qccollect to generate nicely formatted information from picard
  ### Changed
* Change to CLI for running and creating config
* Major overhaul to coverage report script. It's now simpler and more readable!
  ### Fixed
* Fixed sambamba depth to include mapping quality
* Markduplicate now is now by default on marking mode, and will NOT remove duplicates
* Minor formatting and script beautification happened

[1.13.1] - 2018-08-17
---------------------

Fixed
^^^^^


* fixed a typo in MSKMVL config
* fixed a bug in strelka_simple for correct column orders

[1.13.0] - 2018-08-10
---------------------

Added
^^^^^


* rule for all three variant callers for paired analysis now generate a simple VCF file
* rule for all three variant callers for paired analysis to convert VCF into table format
* MVL config file and MVL annotation to VCF calls for SNV/INDEL callers
* CALLER annotation added to SNV/INDEL callers
* exome specific option for strelka paired
* create_config subcommand is now more granular, it accepts all enteries from sample.json as commandline arguments
* Added tabQuery to the assets as a tool to query the tabulated output of summarized VCF
* Added MQ annotation field to Mutect2 output see #67
  ### Changed
* Leaner VCF output from mutect2 with coverage and MQ annotation according to #64
* variant ids are now updated from simple VCF file
  ### Fixed
* Fixed a bug with sambamba depth coverage reporting wrong exon and panel coverage see #68
* The json output is now properly formatted using yapf
* Strelka rule doesn't filter out PASS variants anymore fixes issue #63

[1.12.0] - 2018-07-06
---------------------

Coverage report patch

Added
^^^^^


* Added a new script to retrieve coverage report for a list of gene(s) and transcripts(s)
* Added sambamba exon depth rule for coverage report
* Added a new entry in reference json for exon bed file, this file generated using: https://github.com/hassanfa/GFFtoolkit
  ### Changed
* sambamba_depth rule changed to sambama_panel_depth
* sambamba depth now has fix-mate-overlaps parameter enabled
* sambamba string filter changed to ``unmapped or mate\_is\_unmapped) and not duplicate and not failed\_quality\_control``.
* sambamba depth for both panel and exon work on picard flag (rmdup or mrkdup).
  ### Fixed
* Fixed sambamba panel depth rule for redundant coverage parameter

[1.11.0] - 2018-07-05
---------------------

create config patch for single and paired mode

Changed
^^^^^^^


* create_config is now accepting a paired|single mode instead of analysis json template (see help for changes). It is
  not backward compatible
  ### Added
* analysis_{paired single}.json for creating config. Analysis.json is now obsolete.
  ### Fixed
* A bug with writing output for analysis config, and creating the path if it doesn't exist.
* A bug with manta rule to correctly set output files in config.
* A bug that strelka was still included in sample analysis.

[1.10.0] - 2018-06-07
---------------------

Added
^^^^^


* Markduplicate flag to analysis config

[1.9.0] - 2018-06-04
--------------------

Added
^^^^^


* Single mode for vardict, manta, and mutect.
* merge type for tumor only
  ### Changed
* Single mode variant calling now has all variant calling rules
  ### Fixed
* run_analaysis now accepts workflows for testing pyrposes

[1.8.0] - 2018-06-01
--------------------

Changed
^^^^^^^


* picard create bed interval rule moved into collect hsmetric
* split bed is dependent on bam merge rule
* vardict env now has specific build rather than URL download (conda doesn't support URLs anymore)
  ### Fixed
* new logs and scripts dirs are not re-created if they are empty

[1.7.0] - 2018-05-31
--------------------

Added
^^^^^


* A source altered picard to generated more quality metrics output is added to installation and rules

[1.6.0] - 2018-05-30
--------------------

Added
^^^^^


* report subcommand for generating a pdf report from a json input file
* Added fastqc after removing adapter
  ### Changed
* Markduplicate now has both REMOVE and MARK (rmdup vs mrkdup)
* CollectHSMetrics now has more steps on PCT_TARGET_BASES

[1.5.0] - 2018-05-28
--------------------

Changed
^^^^^^^


* New log and script directories are now created for each re-run
  ### Fixed
* Picardtools' memory issue addressed for large samples

[1.4.0] - 2018-05-18
--------------------

Added
^^^^^


* single sample analysis mode
* alignment and insert size metrics are added to the workflow
  ### Changed
* collectqc and contest have their own rule for paired (tumor vs normal) and single (tumor only) sample.

[1.3.0] - 2018-05-13
--------------------

Added
^^^^^


* bed file for panel analysis is now mandatory to create analaysis config

[1.2.3] - 2018-05-13
--------------------

Changed
^^^^^^^


* vep execution path
* working directory for snakemake

[1.2.2] - 2018-05-04
--------------------

Added
^^^^^


* sbatch submitter and cluster config now has an mail field
  ### Changed
* ``create_config`` now only requires sample and output json. The rest are optional

[1.2.0] - 2018-05-02
--------------------

Added
^^^^^


* snakefile and cluster config in run analysis are now optional with a default value

[1.1.2] - 2018-04-27
--------------------

Fixed
^^^^^


* vardict installation was failing without conda-forge channel
* gatk installation was failing without correct jar file

[1.1.1] - 2018-04-27
--------------------

Fixed
^^^^^


* gatk-register tmp directory

[1.1.0] - 2018-04-26
--------------------

Added
^^^^^


* create config sub command added as a new feature to create input config file
* templates to generate a config file for analysis added
* code style template for YAPF input created. see: https://github.com/google/yapf
* vt conda env added

Changed
^^^^^^^


* install script changed to create an output config
* README updated with usage

Fixed
^^^^^


* fastq location for analysis config is now fixed
* lambda rules removed from cutadapt and fastq

[1.0.3-rc2] - 2018-04-18
------------------------

Added
^^^^^


* Added sbatch submitter to handle it outside snakemake
  ### Changed
* sample config file structure changed
* coding styles updated

[1.0.2-rc2] - 2018-04-17
------------------------

Added
^^^^^


* Added vt environment
  ### Fixed
* conda envs are now have D prefix instead of P (develop vs production)
* install_conda subcommand now accepts a proper conda prefix

[1.0.1-rc2] - 2018-04-16
------------------------

Fixed
^^^^^


* snakemake rules are now externally linked

[1.0.0-rc2] - 2018-04-16
------------------------

Added
^^^^^


* run_analysis subcommand
* Mutational Signature R script with CLI
* unittest to install_conda
* a method to semi-dynamically retrieve suitable conda env for each rule

Fixed
^^^^^


* install.sh updated with gatk and proper log output
* conda environments updated
* vardict now has its own environment and it should not raise anymore errors

[1.0.0-rc1] - 2018-04-05
------------------------

Added
^^^^^


* install.sh to install balsamic
* balsamic barebone cli
* subcommand to install required environments
* README.md updated with basic installation instructions

Fixed
^^^^^


* conda environment yaml files
