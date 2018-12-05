========
BALSAMIC
========

Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer
(**version** = 2.8.1)

* [Installation](doc/install.rst)
* [CLI Guide](doc/cli.rst)
* [Run BALSAMIC](doc/user_guide.rst)

Introduction
======

BALSAMIC is basically a wrapper for its core workflow manager. The goal is to have a package with well defined cli to
make it predictable for user to run somatic calling regaradless of the workflow manger at its core. Right now, BALSAMIC
is using Snakemake as its core, but the goal is to make easily extensible to use other workflow managers such as
Nextflow.
