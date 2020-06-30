from BALSAMIC import __version__
from setuptools import setup, find_packages

version = __version__

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []

setup(
    name="BALSAMIC",
    version=version,
    url="https://github.com/Clinical-Genomics/BALSAMIC",
    author="Hassan Foroughi Asl",
    author_email='hassan.foroughi@scilifelab.se',
    install_requires=install_requires,
    packages=find_packages("BALSAMIC"),
    package_dir={"": "src"},
    package_data={"assets": ["scripts/*R"], "config": ["*.json"]},  
    include_package_data=True,
    zip_safe=False,
    entry_points={
        'console_scripts': ['balsamic=BALSAMIC.commands.base:cli'],
    },
)
