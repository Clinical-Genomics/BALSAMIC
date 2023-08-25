=======================
Documentation Guidline
=======================

BALSAMIC uses Sphinx to build the documentation, see the official documentation of Sphinx: https://www.sphinx-doc.org/en/master/index.html


Following steps explains how to build documents locally.

Create a conda environment:

.. code-block::

   conda create -n balsamic_doc -c bioconda -c conda-forge python=3.11 pip pygraphviz
   conda activate balsamic_doc

Install Sphinx and extensions:

.. code-block::

   cd /path/to/BALSAMIC
   python -m pip install --upgrade --upgrade-strategy eager --no-cache-dir .
   cd docs
   pip install -r requirements.txt -r ../requirements-dev.txt

Build docs:

.. code-block::

   sphinx-build -T -E -b html -d _build/doctrees-readthedocs -D language=en . _build/html

View docs (\ ``open`` or similar command from your OS):

.. code-block::

   open _build/html/index.html
