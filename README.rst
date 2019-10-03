.. image:: https://travis-ci.org/Clinical-Genomics/BALSAMIC.svg?branch=master
    :target: https://travis-ci.org/Clinical-Genomics/BALSAMIC

.. image:: https://coveralls.io/repos/github/Clinical-Genomics/BALSAMIC/badge.svg?branch=master 
    :target: https://coveralls.io/github/Clinical-Genomics/BALSAMIC 

.. image:: https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg
   :target: https://singularity-hub.org/collections/3005

========
BALSAMIC
========

Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer
(**version** = 3.0.1)


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

Development and branching model
======

BALSAMIC is using a development structure similar to GitHub Flow: https://guides.github.com/introduction/flow/ , where a development branch is maintained and features branch are merged into development branch. For deployment, a pull request is created from development branch into master. It is then validated, verified, and a bumpversion justified. Finally, the code owner(s) will approve the merge and merge it into the master branch and release.
