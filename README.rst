.. raw:: html

    <p align="center">
        <a href="https://github.com/Clinical-Genomics/BALSAMIC">
            <img  width=480 src="https://raw.githubusercontent.com/Clinical-Genomics/BALSAMIC/master/BALSAMIC/assets/balsamic_logo.png">
        </a>
        <h3 align="center">Bioinformatic Analysis pipeLine for SomAtic MutatIons in Cancer (v 4.2.4)</h3>
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
     - 4.2.4
   * - Author
     - Hassan Foroughi Asl
   * - Development model
     - Github Flow, release branch
   * - Build status
     - |travis_status_badge|_
   * - Container latest release status
     - |docker_latest_release_status|_
   * - Container build type
     - |docker_build_type|_
   * - Container master status 
     - |docker_latest_build_status|_
   * - Code coverage
     - |code_cov_badge|_
   * - Documentation
     - |rtfd_badge|_
   * - Dependencies
     - |snakemake_badge| |singularity_badge|
   * - Contributors
     - @imsarath , @keyvanelhami


.. |code_cov_badge| image:: https://coveralls.io/repos/github/Clinical-Genomics/BALSAMIC/badge.svg?branch=master 
.. _code_cov_badge: https://coveralls.io/github/Clinical-Genomics/BALSAMIC

.. |travis_status_badge| image:: https://travis-ci.org/Clinical-Genomics/BALSAMIC.svg?branch=master
.. _travis_status_badge: https://travis-ci.org/Clinical-Genomics/BALSAMIC

.. |docker_latest_build_status| image:: https://img.shields.io/docker/cloud/build/hassanf/balsamic
.. _docker_latest_build_status: https://hub.docker.com/r/hassanf/balsamic

.. |docker_latest_release_status| image:: https://img.shields.io/docker/v/hassanf/balsamic?sort=semver 
.. _docker_latest_release_status: https://hub.docker.com/r/hassanf/balsamic/tags 
  
.. |docker_build_type| image:: https://img.shields.io/docker/cloud/automated/hassanf/balsamic
.. _docker_build_type: https://hub.docker.com/r/hassanf/balsamic

.. |snakemake_badge| image:: https://img.shields.io/badge/snakemake-%E2%89%A55.12.3-brightgreen.svg 

.. |singularity_badge| image:: https://img.shields.io/badge/singularity-%E2%89%A53.1.1-brightgreen.svg

.. |rtfd_badge| image:: https://readthedocs.org/projects/balsamic/badge/?version=latest&style=flat
.. _rtfd_badge: https://balsamic.readthedocs.io/en/latest

