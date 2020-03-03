.. raw:: html

    <p align="center">
        <a href="https://github.com/Clinical-Genomics/BALSAMIC">
            <img  width=480 src="../BALSAMIC/assets/balsamic_logo.png">
        </a>
        <h3 align="center">Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer (v 4.1.0)</h3>
        <h3 align="center">FastQ to Annotated VCF</h3>
    </p>

.. image:: https://travis-ci.org/Clinical-Genomics/BALSAMIC.svg?branch=master
    :target: https://travis-ci.org/Clinical-Genomics/BALSAMIC
    :align: right

.. image:: https://coveralls.io/repos/github/Clinical-Genomics/BALSAMIC/badge.svg?branch=master 
    :target: https://coveralls.io/github/Clinical-Genomics/BALSAMIC 
    :align: center


Documentation
=======

BALSAMIC is basically a wrapper for its core workflow manager. The goal is to have a package with well defined cli to
make it predictable for user to run somatic calling regaradless of the workflow manger at its core. BALSAMIC
is using Snakemake as its core.

Essentially, one can run the sample using workflows available within this package and standard Snakemake cli given that
there is a proper config file created.

**Development and branching model**

BALSAMIC is using a development structure similar to GitHub Flow: https://guides.github.com/introduction/flow/ , where
features branch are merged into master branch. Releases will be managed from master branch. It is then validated,
verified, and a bumpversion justified.


* Logo by: Mikael A.
