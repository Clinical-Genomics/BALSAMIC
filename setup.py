import os
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []

# Load the package's __version__.py module as a dictionary.
about = {}
with open(os.path.join(here, "BALSAMIC", "__version__.py")) as f:
    exec(f.read(), about)

setup(
    name="BALSAMIC",
    version=about["__version__"],
    url="https://github.com/Clinical-Genomics/BALSAMIC",
    author="Hassan Foroughi Asl",
    author_email='hassan.foroughi@scilifelab.se',
    install_requires=[
        "click>=7.1.2",
        "colorclass>=2.2.0",
        "coloredlogs>=14.0",
        "cyvcf2<0.10.0",
        "graphviz>=0.14",
        "gsutil==4.50",
        "jinja2>=2.11.2",
        "matplotlib>=3.3.0",
        "networkx>=2.4",
        "numpy>=1.19.2",
        "pandas>1.1.0",
        "psutil>=5.7.0",
        "pydantic>=1.5.1",
        "pygments>=2.6.1",
        "pyyaml>=5.3.1",
        "six>=1.12.0",
        "snakemake==5.13.0",
        "yapf>=0.30.0",
        "h5py>=3.1.0",
        "PyPDF2>=1.26.0",
        "markdown==3.3.3",
        "cryptography<3.4",
    ],
    packages=find_packages(),
    package_data={
        "": [
            "*.toml", "*.json", "*.R", "*.model", "*.yaml", "*.sh", "*.rule",
            "*.smk", "*.awk", "*.html", "*.md"
        ],
    },
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': ['balsamic=BALSAMIC.commands.base:cli'],
    },
)
