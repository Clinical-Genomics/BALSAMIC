=========
Build Doc
=========

Following steps explains how to build documents locally.

Create a conda environment:

.. code-block::

   conda create -n balsamic_doc -c bioconda -c conda-forge python=3.6 pip
   conda activate balsamic_doc

Install Sphinx and extensions:

.. code-block::

   python -m pip install --upgrade --upgrade-strategy eager --no-cache-dir .
   cd docs
   pip install -r requirements.txt -r ../requirements-dev.txt

Build docs:

.. code-block::

   sphinx-build -T -E -b html -d _build/doctrees-readthedocs -D language=en . _build/html

View docs (\ ``open`` or similar command from your OS):

.. code-block::

   open _build/html/index.html
