conda env create -n ${1} --file ${1}.yaml
source activate ${1}
pip install --no-cache-dir cnvkit==0.9.4 biopython==1.76
