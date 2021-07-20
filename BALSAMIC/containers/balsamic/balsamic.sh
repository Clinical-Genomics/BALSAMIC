conda env update --name base --file ${1}.yaml --prune
pip install --no-cache-dir .
