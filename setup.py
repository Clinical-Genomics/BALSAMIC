import os
from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

# Load the package's __version__.py module as a dictionary.
about = {}
with open(os.path.join(here, "BALSAMIC", "__version__.py")) as f:
    exec(f.read(), about)

# Requirements
requirements = [
    "click==8.1.3",
    "colorclass>=2.2.0",
    "coloredlogs>=14.0",
    "graphviz==0.16",
    "gsutil==5.23",
    "jinja2>=2.11.2",
    "matplotlib==3.5.2",
    "networkx>=2.6.3",
    "numpy>=1.21.6",
    "pandas>=1.1.5",
    "psutil>=5.7.0",
    "pydantic>=1.9.0",
    "pygments>=2.6.1",
    "pyyaml>=5.3.1",
    "six>=1.12.0",
    "snakemake==6.5.3",
    "yapf>=0.30.0",
    "h5py>=3.6.0",
    "PyPDF2>=1.26.0",
    "markdown==3.3.3",
    "cryptography==40.0.2",
    "tabulate==0.8.10",
    "toml==0.10.2",
]

# The C libraries required to build numpy are not available on RTD
if not os.getenv("READTHEDOCS"):
    requirements.extend(["cyvcf2==0.30.22"])

setup(
    name="BALSAMIC",
    version=about["__version__"],
    url="https://github.com/Clinical-Genomics/BALSAMIC",
    author="Hassan Foroughi Asl",
    author_email="hassan.foroughi@scilifelab.se",
    install_requires=requirements,
    packages=find_packages(),
    package_data={
        "": [
            "*.toml",
            "*.json",
            "*.R",
            "*.model",
            "*.yaml",
            "*.sh",
            "*.rule",
            "*.smk",
            "*.awk",
            "*.html",
            "*.md",
        ],
    },
    include_package_data=True,
    zip_safe=False,
    entry_points={
        "console_scripts": ["balsamic=BALSAMIC.commands.base:cli"],
    },
)
