========================
Snakemake Etiquette
========================

The bioinformatics core analysis in BALSAMIC is defined by set of rules written as a Snakemake rules (``*.rule``) and Snakemake workflow as (``*.smk``). Main ``balsamic.smk`` workflow uses these rules to create sets of output files from sets of input files. Using ``{wildcards}`` Snakemake can automatically determine the dependencies between the rules by matching file names. The following guidelines describe the general conventions for naming and order of the rules, while writing a Snakemake file in BALSAMIC. For further description of how Snakemake works, please refer to Snakemake official documentation: https://snakemake.readthedocs.io/


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

**Descriptions**
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**rulename**: Rule name briefly should outline the program and functions utilized inside the rule. Each word is seperated by a underscore ``_``. First word is the bioinformatic tool or script's name. The following words describe subcommand within that bioinformatic tool and then followed by workflow specific description. The word length shouldn't exceed more than 3 or 4 words. Make sure rule names are updated within ``config/cluster.json`` and it is all lowercase. Examples: ``picard_collecthsmetrics_umi``, ``bcftools_query_calculateaftable_umi``

**input**: It is strongly recommended to set input and output files as named. Refrain from introducing new wildcards as much as possible.

**output**: This should follow the same instructions as ``input``.

**benchmark**: Benchmark name is prefixed with rule name and suffixed with '.tsv' file extension.

**singularity**: Make sure the singularity image does contain a Conda environment with required bioinformatics tools. Do not use this field if ``run`` is used instead of ``shell``.

**params**: If the defined parameter is a threshold or globally used constant; add it to ``utils/constants.py``. Respective class models need to be updated in ``utils/models.py``. 

**threads**: Make sure for each rule, the correct number of threads are assigned in ``config/cluster.json``. Otherwise it will be assigned default values from ``config/cluster.json`` . If there is no need for multithreading, this field can be removed from rule.

**message**: A short message describing the function of rule. Add any relevant wildcard to message to make it readable and understandable. It is also recommended to use ``params`` to build a more descriptive ``message``

**shell (run)**: Code inside the `shell/run` command should be left indented. Shell lines no longer than 100 characters. Break the long commands with ``\`` and followed by a new line. Avoid having long Python code within ``run``, instead add it to ``utils/`` as a Python script and import the function.

Example:

::

java -jar \
-Djava.io.tmpdir=${{tmpdir}} \
-Xms8G -Xmx16G \
$CONDA_PREFIX/share/picard.jar \
MarkDuplicates \
{input.named_input_1} \
{output.named_output_1};


Example for external python scripts that can be saved as modules in ``utils/*.py`` and can use them as definitions in rules as:

:: 

  from BALSAMIC.utils.workflowscripts import get_densityplot
  get_densityplot(input.named_input1, params.named_params_1, output.named_output1 )

Similarly ``awk`` or ``R`` external scripts can be saved in ``assets/scripts/*awk`` and can be invoked using `get_script_path` as: 

::
  
  params: 
      consensusfilter_script = get_script_path("FilterDuplexUMIconsensus.awk")
  shell:
       """
  samtools view -h {input} | \
  awk -v MinR={params.minreads} \
  -v OFS=\'\\t\' -f {params.consensusfilter_script} | \
  samtools view -bh - > {output}
       """

**References**
~~~~~~~~~~~~~~~

1. https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html
2. https://snakemake.readthedocs.io/en/stable/snakefiles/writing_snakefiles.html
