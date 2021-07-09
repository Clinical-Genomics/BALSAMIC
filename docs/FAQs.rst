=================================
Frequently Asked Questions (FAQs)
=================================

**BALSAMIC**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



**UMIworkflow**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**What are UMIs**

- Unique Molecular Identifiers (UMIs) are short random nucleotide sequences (3-20 bases) that are ligated to the ends of DNA fragments prior to sequencing to greatly reduce the impact of PCR duplicates and sequencing errors on the variant calling process.

.. figure:: images/UMI.png

    Figure1: Design of UMI adapters in the library preparation. Ref_ 

.. _Ref: https://plone.bcgsc.ca/services/solseq/duplex-umi-documents/idt_analysisguideline_varcall-umis-dupseqadapters/

__ Ref_


**How is the UMIworkflow implemented**

- CG's UMIworkflow is implemented using the commercial software Sentieon. The Sentieon tools provide functionality for extracting UMI tags from fastq reads and performing barcode-aware consensus generation. The workflow is as described:

.. figure:: images/UMIworkflow.png

    Figure2: UMI workflow steps.

**How is the UMI structure defined**

Our pair-end sequencing read length is about 151 bp and the UMI structure is defined as`3M2S146T, 3M2S146T` where `3M` represents 3 UMI bases, `2S` represents 2 skipped bases,  `146T` represents 146 bases in the read.

**Are there any differences in the UMI read extraction if the read structure is defined as `3M2S146T, 3M2S146T` or `3M2S+T, 3M2S+T`?**

In theory, this should be the same if the read length is always 151bp. But the recommendation is to use `3M2S+T, 3M2S+T` so that UMIworkflow can handle any unexpected input data.

**How does the `umi extract` tool handle sequencing adapters?  Do the input reads always need to be adapter removed fastq reads**

The presence of 5' adapter sequences can cause issues for the Sentieon `umi extract` tool, as the extract tool will not correctly identify the UMI sequence. If 5' adapter contamination is found in the data, before processing with the `umi extract` tool, these adapter sequences needed to be removed with a third-party trimming tool. 
3' adapter contamination is much more common and can occur when the insert size is shorter than the sequence read length. The Sentieon `umi consensus` tool will correctly identify and handle 3' adapter/barcode contamination during consensus read creation.

**How does Sentieon `umi consensus` tool handles paired-end reads**

The `umi consensus` tool will merge overlapping read pairs when it can, but it is not possible for reads with an insert size greater than 2x the read length as there is some unknown intervening sequence. In this case, `umi consensus` will output a consensus read pair where each consensus read in the pair is constructed separately, while other reads in the dataset are collapsed/merged to single-end reads.

.. figure:: images/sentieon_consensus.jpg

    Figure3: Figure taken from Sentieon document. 

**Purpose of consensus-filtering step in the UMIworkflow**

Mainly to reduce the calling of false-positive variants. Consensus filtering is based on the setting of minimum raw reads (MinR) supporting each UMI group.  By default, `MinR` is set as 3,1,1, meaning that the minimum number of raw reads in both strands should be greater than 1 and the sum of both strands is greater than 3.   The default `3,1,1` is a good starting point at lower coverages. This setting can be further adjusted accordingly at higher coverages or if finding false-positive calls due to consensus reads with little read support.

**How is the performance of other variant callers for analysing UMI datasets**
UMI workflow is validated with two datasets (SeraCare and HapMap). The Vardict failed to call the true reference variants while the TNscope performed better. A more detailed analysis is summarized here_. 

.. _here: https://drive.google.com/file/d/1Y1kNPE5u9VvykjmNhG4RydVMUyezbqh5/view?usp=sharing

We are still investigating other UMI-aware variant callers and maybe in the future, if something works better, additional varcallers will be added to the UMIworkflow.

**Git Related Questions**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**How to make a new release of balsamic**

Our release model looks like in the figure below:

.. figure:: images/git_releasemodel.png

*Requirements:* bumpversion

`pip install bumpversion`

If changes are hotfixes/patches, make a release from master.
Else make a release from develop.

If decided to make a release version from `develop`, do the following to release a balsamic new version:

1. `git checkout master && git pull`
2. `git checkout develop`
3. `git merge master`
4.  Fix all possible conflicts and `git push`
5. `git checkout -B release_X.X.X`

Change X.X.X version of the release in the CHANGELOG.rst

Make a pull request to master at this point. After pull request is approved and merge it into master:

1. `git checkout master`
2. `git pull`
3. `bumpversion --verbose [major/minor/patch]`
4. `git push`
5. `git push --tags`


**References**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**How to generate reference files for ascatNGS**

Detailed information is available at the ascatNGS_

.. _ascatNGS: https://github.com/cancerit/ascatNgs

Briefly, ascatNGS needs gender loci file if gender information for the input sample is not available. The second file is `SnpGcCorrections.tsv`, which is prepared from the 1000 genome SNP panel.

1. **Gender loci file:**
  
GRCh37d5_Y.loci contains the following contents:

.. line-block::
    Y	4546684
    Y	2934912
    Y	4550107
    Y	4549638


2. **GC correction file:**
 
First step is to download the 1000 genome snp file and convert it from .vcf to .tsv. The detailed procedure to for this step is provided here_

.. here_: https://github.com/cancerit/ascatNgs/wiki/Human-reference-files-from-1000-genomes-VCFs

.. code:: console

    export TG_DATA=ftp://ftp.ensembl.org/pub/grch37/release-83/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz


Followed by:

.. code:: console

    curl -sSL $TG_DATA | zgrep -F 'E_Multiple_observations' | grep -F 'TSA=SNV' |\
    perl -ane 'next if($F[0] !~ m/^\d+$/ && $F[0] !~ m/^[XY]$/);\
    next if($F[0] eq $l_c && $F[1]-1000 < $l_p); $F[7]=~m/MAF=([^;]+)/;\
    next if($1 < 0.05); printf "%s\t%s\t%d\n", $F[2],$F[0],$F[1];\
    $l_c=$F[0]; $l_p=$F[1];' > SnpPositions_GRCh37_1000g.tsv


--or--

.. code:: console

    curl -sSL $TG_DATA | zgrep -F 'E_Multiple_observations' | grep -F 'TSA=SNV' |\
    perl -ane 'next if($F[0] !~ m/^\d+$/ && $F[0] !~ m/^[XY]$/); $F[7]=~m/MAF=([^;]+)/;\
    next if($1 < 0.05); next if($F[0] eq $l_c && $F[1]-1000 < $l_p);\
    printf "%s\t%s\t%d\n", $F[2],$F[0],$F[1]; $l_c=$F[0]; $l_p=$F[1];'\
    > SnpPositions_GRCh37_1000g.tsv

Second step is to use `SnpPositions.tsv` file and generate `SnpGcCorrections.tsv` file as described here_

.. here_: https://github.com/cancerit/ascatNgs/wiki/Convert-SnpPositions.tsv-to-SnpGcCorrections.tsv

.. code:: console

    ascatSnpPanelGcCorrections.pl genome.fa SnpPositions.tsv > SnpGcCorrections.tsv

