.. raw:: html

    <p align="center">
        <a href="https://github.com/Clinical-Genomics/BALSAMIC">
            <img  width=480 src="https://raw.githubusercontent.com/Clinical-Genomics/BALSAMIC/master/BALSAMIC/assets/balsamic_logo.png">
        </a>
        <h3 align="center">Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer (v 4.1.0)</h3>
        <h3 align="center">FastQ to Annotated VCF</h3>
    </p>




   
BALSAMIC is basically a wrapper for its core workflow manager. The goal is to have a package with well defined cli to
make it reproducible for user to run somatic calling regaradless of the workflow manger at its core. Right now, BALSAMIC
is using Snakemake as its core. So one can run the sample using workflows available within this package and standard
Snakemake cli given that there is a proper config file created.


.. list-table:: 
   :widths: 20 50
   :header-rows: 0
   :stub-columns: 1

   * - Source code
     - https://github.com/Clinical-Genomics/BALSAMIC
   * - Version
     - 4.1.0
   * - Author
     - Hassan Foroughi Asl
   * - Development model
     - Github Flow, release branch
   * - Build status
     - |travis_status_badge|_
   * - Code coverage
     - |code_cov_badge|_
   * - Dependencies
     - |snakemake_badge|, |singularity_badge|
   * - Contributors
     - Sarath Murugan


.. |code_cov_badge| image:: https://coveralls.io/repos/github/Clinical-Genomics/BALSAMIC/badge.svg?branch=master 
.. _code_cov_badge: https://coveralls.io/github/Clinical-Genomics/BALSAMIC

.. |travis_status_badge| image:: https://travis-ci.org/Clinical-Genomics/BALSAMIC.svg?branch=master
.. _travis_status_badge: https://travis-ci.org/Clinical-Genomics/BALSAMIC

.. |snakemake_badge| image:: https://img.shields.io/badge/snakemake-%E2%89%A55.12.3-brightgreen.svg 

.. |singularity_badge| image:: https://img.shields.io/badge/singularity-%E2%89%A53.1.1-brightgreen.svg 
