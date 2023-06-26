===================
Method description
===================

Target Genome Analysis
~~~~~~~~~~~~~~~~~~~~~~

BALSAMIC :superscript:`1` (**version** = 12.0.1) was used to analyze the data from raw FASTQ files.
We first quality controlled FASTQ files using FastQC v0.11.9 :superscript:`2`.
Adapter sequences and low-quality bases were trimmed using fastp v0.23.2 :superscript:`3`.
Trimmed reads were mapped to the reference genome hg19 using BWA MEM v0.7.17 :superscript:`4`.
The resulted SAM files were converted to BAM files and sorted using samtools v1.15.1 :superscript:`5`.
Duplicated reads were marked using Picard tools MarkDuplicate v2.27.1 :superscript:`6`. The unmapped reads are removed
and promptly quality controlled using CollectHsMetrics, CollectInsertSizeMetrics and CollectAlignmentSummaryMetrics functionalities.
Results of the quality controlled steps were summarized by MultiQC v1.12 :superscript:`7`.
Small somatic mutations (SNVs and INDELs) were called for each sample using VarDict v2019.06.04 :superscript:`8`.
Apart from the Vardict filters to report the variants, the called-variants were also further second filtered using the criteria
(*MQ >= 40, DP >= 100, VD >= 5, Minimum AF >= 0.007, Maximum AF < 1, GNOMADAF_popmax <= 0.005, swegen AF < 0.01*).
Only those variants that fulfilled the filtering criteria and scored as `PASS` in the VCF file were reported.
Structural variants (SV) were called using Manta v1.6.0 :superscript:`9` and Dellyv1.0.3 :superscript:`10`.
Copy number variations (CNV) were called using CNVkit v0.9.9 :superscript:`11`.
The variant calls from CNVkit, Manta and Delly were merged using SVDB v2.8.1 :superscript:`12`.
The clinical set of SNV and SV is also annotated and filtered against loqusDB curated frequency of observed variants (frequency < 0.01) from normal samples.
All variants were annotated using Ensembl VEP v104.3 :superscript:`13`. We used vcfanno v0.3.3 :superscript:`14`
to annotate somatic variants for their population allele frequency from gnomAD v2.1.1 :superscript:`18`, CADD v1.6 :superscript:`24`, SweGen :superscript:`22` and frequency of observed variants in normal samples.

Whole Genome Analysis
~~~~~~~~~~~~~~~~~~~~~

BALSAMIC :superscript:`1` (**version** = 12.0.1) was used to analyze the data from raw FASTQ files.
We first quality controlled FASTQ files using FastQC v0.11.9 :superscript:`2`.
Adapter sequences and low-quality bases were trimmed using fastp v0.23.2 :superscript:`3`.
Trimmed reads were mapped to the reference genome hg19 using sentieon-tools :superscript:`15`.
The resulted SAM files were converted to BAM files and sorted using samtools v1.15.1 :superscript:`5`.
Duplicated reads were marked using Picard tools MarkDuplicate v2.27.1 :superscript:`6`. The unmapped reads are removed
and promptly quality controlled using CollectMultipleMetrics and CollectWgsMetrics functionalities.
Results of the quality controlled steps were summarized by MultiQC v1.12 :superscript:`7`.
Small somatic mutations (SNVs and INDELs) were called for each sample using Sentieon TNscope :superscript:`16`.
The called-variants were also further second filtered using the criteria (DP(tumor,normal) >= 10; AD(tumor) >= 3; AF(tumor) >= 0.05, Maximum AF(tumor < 1;  GNOMADAF_popmax <= 0.001; normalized base quality scores >= 20, read_counts of alt,ref alle > 0).
Structural variants were called using Manta v1.6.0 :superscript:`9`, Delly v1.0.3 :superscript:`10` and TIDDIT v3.3.2 :superscript:`12`.
Copy number variations (CNV) were called using ascatNgs v4.5.0 :superscript:`17` (tumor-normal), Delly v1.0.3 :superscript:`10` and CNVpytor v1.2.1 :superscript:`22` (tumor-only) and converted from CNV to deletions (DEL) and duplications (DUP).
The structural variant (SV) calls from Manta, Delly, TIDDIT ascatNgs (tumor-normal) and CNVpytor (tumor-only) were merged using SVDB v2.8.1 :superscript:`12`
The clinical set of SNV and SV is also annotated and filtered against loqusDB curated frequency of observed variants (frequency < 0.01) from normal samples.
All variants were annotated using Ensembl VEP v104.3 :superscript:`13`. We used vcfanno v0.3.3 :superscript:`14`
to annotate somatic single nucleotide variants for their population allele frequency from gnomAD v2.1.1 :superscript:`18`, CADD v1.6 :superscript:`24`, SweGen :superscript:`22`  and frequency of observed variants in normal samples.

UMI Data Analysis
~~~~~~~~~~~~~~~~~~~~~

BALSAMIC :superscript:`1` (**version** = 12.0.1) was used to analyze the data from raw FASTQ files.
We first quality controlled FASTQ files using FastQC v0.11.9 :superscript:`2`.
Adapter sequences and low-quality bases were trimmed using fastp v0.23.2 :superscript:`3`.
UMI tag extraction and consensus generation were performed using Sentieon tools v202010.02 :superscript:`15`.
The alignment of UMI extracted and consensus called reads to the human reference genome (hg19) was done by bwa-mem and
samtools using Sentieon utils. Consensus reads were filtered based on the number of minimum reads supporting each UMI tag group.
We applied a criteria filter of minimum reads `3,1,1`. It means that at least three UMI tag groups should be ideally considered from both DNA strands,
where a minimum of at least one UMI tag group should exist in each single-stranded consensus read.
The filtered consensus reads were quality controlled using Picard CollectHsMetrics v2.27.1 :superscript:`5`. Results of the quality controlled steps were summarized by MultiQC v1.12 :superscript:`6`.
For each sample, somatic mutations were called using Sentieon TNscope :superscript:`16`, with non-default parameters for passing the final list of variants
(--min_tumor_allele_frac 0.0005, --filter_t_alt_frac 0.0005, --min_init_tumor_lod 0.5, min_tumor_lod 4, --max_error_per_read 5  --pcr_indel_model NONE, GNOMADAF_popmax <= 0.02).
The clinical set of SNV and SV is also annotated and filtered against loqusDB curated frequency of observed variants (frequency < 0.01) from normal samples.
All variants were annotated using Ensembl VEP v104.3 :superscript:`7`. We used vcfanno v0.3.3 :superscript:`8` to annotate somatic variants for their population allele frequency from gnomAD v2.1.1 :superscript:`18`, CADD v1.6 :superscript:`24`, SweGen :superscript:`22` and frequency of observed variants in normal samples.
For exact parameters used for each software, please refer to  https://github.com/Clinical-Genomics/BALSAMIC.
We used three commercially available products from SeraCare [Material numbers: 0710-067110 :superscript:`19`, 0710-067211 :superscript:`20`, 0710-067312 :superscript:`21`] for validating the efficiency of the UMI workflow in identifying 14 mutation sites at known allelic frequencies.


**References**
~~~~~~~~~~~~~~~~

1. Foroughi-Asl, H., Jeggari, A., Maqbool, K., Ivanchuk, V., Elhami, K., & Wirta, V. BALSAMIC: Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer (Version v8.2.10) [Computer software]. https://github.com/Clinical-Genomics/BALSAMIC
2. Babraham Bioinformatics - FastQC A Quality Control tool for High Throughput Sequence Data. Accessed June 22, 2020. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
3. Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics. 2018;34(17):i884-i890. https://doi.org/10.1093/bioinformatics/bty560
4. Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997v2 [q-bio.GN]
5. Li H, Handsaker B, Wysoker A, Fennell T, Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup (2009) The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9. https://doi.org/10.1093/bioinformatics/btp352
6. Picard Tools - By Broad Institute. Accessed June 22, 2020. https://broadinstitute.github.io/picard/
7. Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-3048. https://doi.org/10.1093/bioinformatics/btw354
8. Lai Z, Markovets A, Ahdesmaki M, Chapman B, Hofmann O, McEwen R, Johnson J, Dougherty B, Barrett JC, and Dry JR. VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. Nucleic Acids Res. 2016. https://doi.org/10.1093/nar/gkw227
9. Chen, X. et al. (2016) Manta: rapid detection of structural variants and indels for germline and cancer sequencing applications. Bioinformatics, 32, 1220-1222. https://doi.org/10.1093/bioinformatics/btv710
10. Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stuetz, Vladimir Benes, Jan O. Korbel. DELLY: structural variant discovery by integrated paired-end and split-read analysis. Bioinformatics. 2012 Sep 15;28(18):i333-i339. https://doi.org/10.1093/bioinformatics/bts378
11. Talevich, E, Shain, A.H, Botton, T, & Bastian, B.C. CNVkit: Genome-wide copy number detection and visualization from targeted sequencing. PLOS Computational Biology. 2016, 12(4):e1004873. https://doi.org/10.1371/journal.pcbi.1004873
12. Jesper Eisfeldt et.al. TIDDIT, an efficient and comprehensive structural variant caller for massive parallel sequencing data. F1000 research. 2017. https://doi.org/10.12688/f1000research.11168.2
13. McLaren W, Gil L, Hunt SE, et al. The Ensembl Variant Effect Predictor. Genome Biology. 2016;17(1):122. https://doi.org/10.1186/s13059-016-0974-4
14. Pedersen BS, Layer RM, Quinlan AR. Vcfanno: fast, flexible annotation of genetic variants. Genome Biology. 2016;17(1):118. https://doi.org/10.1186/s13059-016-0973-5
15. Donald Freed, Rafael Aldana, Jessica A. Weber, Jeremy S. Edwards. The Sentieon Genomics Tools - A fast and accurate solution to variant calling from next-generation sequence data. Bioinformatics. 2016, Volume 32,Issue 8. https://doi.org/10.1093/bioinformatics/btv710
16. Donald Freed, Renke Pan, Rafael Aldana. TNscope: Accurate Detection of Somatic Mutations with Haplotype-based Variant Candidate Detection and Machine Learning Filtering. bioRvix. https://doi.org/10.1101/250647
17. Keiran MR, Peter VL, David CW, David J, Andrew M, Adam PB , Jon WT, Patrick T, Serena Nik-Zainal, Peter J C. ascatNgs: Identifying Somatically Acquired Copy-Number Alterations from Whole-Genome Sequencing Data. Curr Protoc Bioinformatics. 2016. https://doi.org/10.1002/cpbi.17
18. Karczewski, K.J., Francioli, L.C., Tiao, G. et al. The mutational constraint spectrum quantified from variation in 141,456 humans. Nature 581, 434–443. 2020. https://doi.org/10.1038/s41586-020-2308-7
19. https://www.seracare.com/Seraseq-ctDNA-Complete-Reference-Material-AF1-0710-0671/
20. https://www.seracare.com/Seraseq-ctDNA-Complete-Reference-Material-AF05-0710-0672/
21. https://www.seracare.com/Seraseq-ctDNA-Complete-Reference-Material-AF01-0710-0673/
22. Ameur, A., Dahlberg, J., Olason, P. et al. SweGen: a whole-genome data resource of genetic variability in a cross-section of the Swedish population. Eur J Hum Genet 25, 1253–1260 (2017). https://doi.org/10.1038/ejhg.2017.130
23. Milovan Suvakov, Arijit Panda, Colin Diesh, Ian Holmes, Alexej Abyzov, CNVpytor: a tool for copy number variation detection and analysis from read depth and allele imbalance in whole-genome sequencing, GigaScience, Volume 10, Issue 11, November 2021, giab074, https://doi.org/10.1093/gigascience/giab074
24. Rentzsch P., Witten D., Cooper G.M., Shendure J., Kircher M. CADD: predicting the deleteriousness of variants throughout the human genome. Nucleic Acids Res. 2018. https://doi.org/10.1093/nar/gky1016. PubMed PMID: 30371827.