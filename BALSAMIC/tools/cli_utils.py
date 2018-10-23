import os
import glob
import yaml
from itertools import chain


def add_doc(docstring):
    """
    A decorator for adding docstring. Taken shamelessly from stackexchange. 
    """

    def document(func):
        func.__doc__ = docstring
        return func

    return document


def createDir(path, interm_path=[]):
    '''
    Creates directories by recursively checking if it exists, otherwise increments the number
    '''
    if os.path.isdir(os.path.abspath(path)):
        basepath = os.path.basename(os.path.abspath(path))
        basepath_number = 0
        if "." in basepath:
            basepath_number = int(basepath.split(".")[1])
        basepath_string = basepath.split(".")[0]
        basepath_number += 1
        path = os.path.join(
            os.path.dirname(os.path.abspath(path)),
            ".".join([basepath_string, str(basepath_number)]))
        interm_path.append(path)
        createDir(path, interm_path)
        return interm_path[-1]
    else:
        os.makedirs(os.path.abspath(path), exist_ok=True)
        return interm_path[-1]


def get_packages(yaml_file):
    '''
    return packages found in a conda yaml file

    input: conda yaml file path
    output: list of packages
    '''
    try:
        with open(yaml_file, 'r') as f:
            pkgs = yaml.load(f)['dependencies']
    except OSError:
        print('File not found', yaml_file)

    return (pkgs)


def get_package_split(conda_yamls):
    '''
    Get a list of conda env files, and extract pacakges

    input: conda env files
    output: dict of packages and their version
    '''

    pkgs = [
        "bwa", "bcftools", "cutadapt", "fastqc", "gatk", "manta", "picard",
        "sambamba", "strelka", "samtools", "tabix", "vardic"
    ]

    pkgs = dict(
        [[y.split("=")[0], y.split("=")[1]]
         for y in set(
             chain.from_iterable([get_packages(s) for s in conda_yamls]))
         if y.split("=")[0] in pkgs])

    return (pkgs)
