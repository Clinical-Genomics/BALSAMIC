conda env update -n base --file ${1}.yaml --prune
pip install --no-cache-dir cnvkit==0.9.9 biopython==1.79
