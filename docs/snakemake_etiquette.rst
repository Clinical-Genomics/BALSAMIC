========================
Snakemake Etiquette
========================

The bioinformatics core analysis in BALSAMIC is defined by set of rules written as a Snakemakefile (`*smk`). Each rule decomposes the workflow into small steps by specifying how to create sets of output files from sets of input files. Using `{wildcards}` Snakemake can automatically determine the dependencies between the rules by matching file names. The following guidelines describe the general conventions for naming and order of the rules, while writing a Snakemake file in BALSAMIC.


**Structure of Snakemake rules**
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    rule <program>_<function>_<optional_tag>_<optional_tag>: 
        input:
            <named_input_1> = ...,
            <named_input_2> = ...,
        output:
            <named_output_1> = ...,
        benchmark:
            Path(benchmark_dir, "<rule_name>_<{sample}/{case_name}>.tsv").as_posix()
        singularity:
            singularity_image
        params:
            <named_param_1> = ...,
            <named_param_1> = ...,
        threads:
            get_threads(cluster_config, '<rule_name>')
        message:
            ("Align fastq files with bwa-mem to reference genome and sort using samtools for sample: {sample}"
            "<second line is allowed to cover more description>")
        shell:
             """
        <first_command> <options>;
        
        <second_command> <options>;

        <a_long_command> <--option-1> <value_1> \
        <--option-2> <value_2> \
        <--option-3> <value_3>;
            """


**rulename**: Rule name briefly should outline the program and functions utilized inside the rule. The word length shouldnt exceed more than 3 or 4 words. Make sure rule names are updated within ``config/cluster.json`` and it is all lowercase. Examples: `picard_collecthsmetrics_umi`, `bcftools_query_calculateAFtable_umi`

**input**: 

**output**:

**benchmark**: Benchmark name is prefixed with rule name and suffixed with '.tsv' file extension.

**singularity**: Make sure the singularity image does contain a conda environment with required bioinformatics tools.

**params**: If the defined parameter is a threshold or globally used constant; add it to ``utils/constants.py``

**threads**: Make sure for each rule, the correct number of threads are assigned in ``config/cluster.json``.

**message**: A short message describing the function of rule. Add any relevant wildcard to message to make it readable and understandable.

**shell (run)**: Code inside the `shell/run` command should be left idented. Shell lines no longer than 100 characters. Break the long commands with ``\`` and followed by a new line. 

`Example:`

::

    java -jar -Djava.io.tmpdir=${{tmpdir}} -Xms8G -Xmx16G $CONDA_PREFIX/share/picard.jar \
    MarkDuplicates {input.named_input_1} {output.named_output_1}

    java -jar \
    -Djava.io.tmpdir=${{tmpdir}} \
    -Xms8G -Xmx16G \
    $CONDA_PREFIX/share/picard.jar \
    MarkDuplicates \
    {input.named_input_1} \
    {output.named_output_1};`

