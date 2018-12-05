========
BALSAMIC
========

Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer
(**version** = 2.8.1)


* `Installation <doc/install.rst>`_
* `CLI Guide <doc/cli.rst>`_
* `Run BALSAMIC <doc/user_guide.rst>`_

Introduction
======

BALSAMIC is basically a wrapper for its core workflow manager. The goal is to have a package with well defined cli to
make it predictable for user to run somatic calling regaradless of the workflow manger at its core. Right now, BALSAMIC
is using Snakemake as its core, but the goal is to make easily extensible to use other workflow managers such as
Nextflow.

Essentially, one can run the sample using workflows available within this package and standard Snakemake cli.
