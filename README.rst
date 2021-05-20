.. raw:: html

    <p align="center">
        <a href="https://github.com/Clinical-Genomics/BALSAMIC">
            <img  width=480 src="https://raw.githubusercontent.com/Clinical-Genomics/BALSAMIC/master/BALSAMIC/assets/balsamic_logo.png">
        </a>
        <h3 align="center">Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer (v 7.2.2)</h3>
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
     - |latest_tag|
   * - Author
     - Hassan Foroughi Asl
   * - Development model
     - Gitflow
   * - Build status
     - |test_status_badge|
   * - Container latest release status
     - |docker_latest_release_status|
   * - Container master status 
     - |docker_latest_build_status|
   * - Code coverage
     - |code_cov_badge|_
   * - Documentation
     - |rtfd_badge|_
   * - Dependencies
     - |snakemake_badge| |singularity_badge|
   * - Contributors
     - @ashwini06 , @mropat , @imsarath , @keyvanelhami


.. |code_cov_badge| image:: https://codecov.io/gh/Clinical-Genomics/BALSAMIC/branch/develop/graph/badge.svg?token=qP68U3PNwV 
.. _code_cov_badge: https://codecov.io/gh/Clinical-Genomics/BALSAMIC

.. |latest_tag| image:: https://img.shields.io/github/v/tag/clinical-genomics/BALSAMIC

.. |test_status_badge| image:: https://github.com/Clinical-Genomics/BALSAMIC/actions/workflows/pytest_and_coveralls.yml/badge.svg

.. |docker_latest_build_status| image:: https://github.com/Clinical-Genomics/BALSAMIC/actions/workflows/docker_build_push_master.yml/badge.svg 

.. |docker_latest_release_status| image:: https://github.com/Clinical-Genomics/BALSAMIC/actions/workflows/docker_build_push_release.yml/badge.svg?tag=v7.2.2 
  
.. |snakemake_badge| image:: https://img.shields.io/badge/snakemake-%E2%89%A55.12.3-brightgreen.svg 

.. |singularity_badge| image:: https://img.shields.io/badge/singularity-%E2%89%A53.1.1-brightgreen.svg

.. |rtfd_badge| image:: https://readthedocs.org/projects/balsamic/badge/?version=latest&style=flat
.. _rtfd_badge: https://balsamic.readthedocs.io/en/latest

