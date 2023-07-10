conda env update -n base --file "${1}".yaml
pip install --no-cache-dir fpdf2==2.4.6
pip install --no-cache-dir zipp==3.16.0
