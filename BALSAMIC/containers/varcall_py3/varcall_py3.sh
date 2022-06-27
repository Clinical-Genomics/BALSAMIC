conda env update -n base --file ${1}.yaml --prune
pip install --no-cache-dir cyvcf2==0.30.15
