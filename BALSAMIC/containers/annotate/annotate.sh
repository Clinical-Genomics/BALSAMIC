conda env create -n ${1} --file ${1}.yaml
source activate ${1}
pip install --no-cache-dir genmod==3.7.4
