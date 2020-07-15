"""This file contains consatant variables used by BALSAMIC"""
import sys
from pathlib import Path

import BALSAMIC


#Path to conda folder containing YAML files with verions of software usen un BALSAMIC workflow
CONDA_ENV_PATH = Path(
    Path(sys.modules["BALSAMIC"].__file__).parent.resolve() /
    "conda").as_posix()


#Path to config YAML file to be accessed by Snakemake
CONDA_ENV_YAML = Path(
    Path(sys.modules["BALSAMIC"].__file__).parent.resolve() / "config" /
    "balsamic_env.yaml").as_posix()


#Path to rule files to be accessed by Snakemake
RULE_DIRECTORY = Path(
    sys.modules["BALSAMIC"].__file__).parent.resolve().as_posix() + "/"

#BALSMIC version
BALSAMIC_version = BALSAMIC.__version__

#Configuration of VCF settings
VCF_DICT = {
    "tnsnv": {"mutation": "somatic", "type": "SNV"},
    "manta": {"mutation": "somatic", "type": "SV"},
    "pindel": {"mutation": "somatic", "type": "SV"},
    "cnvkit": {"mutation": "somatic", "type": "CNV"},
    "mutect": {"mutation": "somatic", "type": "SNV"},
    "vardict": {"mutation": "somatic", "type": "SNV"},
    "strelka": {"mutation": "somatic", "type": "SNV"},
    "tnscope": {"mutation": "somatic", "type": "SNV"},
    "vcfmerge": {"mutation": "somatic", "type": "SNV"},
    "dnascope": {"mutation": "germline", "type": "SNV"},
    "tnhaplotyper": {"mutation": "somatic", "type": "SNV"},
    "manta_germline": {"mutation": "germline", "type": "SV"},
    "haplotypecaller": {"mutation": "germline", "type": "SNV"},
    "strelka_germline": {"mutation": "germline", "type": "SNV"},
}

#Configuration of VARDICT settings
VARDICT_SETTINGS = {
    "AD": {
        "tag_value": 5,
        "filter_name" : "balsamic_low_tumor_ad",
        "field" : "INFO"
        },
    "DP": {
        "tag_value" : 100,
        "filter_name" : "balsamic_low_tumor_dp",
        "field" : "INFO",
        },
    "MQ": {
        "tag_value" : 50, 
        "filter_name" : "balsamic_low_mq",
        "field" : "INFO"
        },
    "AF_max": {
        "tag_value" : 1,
        "filter_name" : "balsamic_af_one",
        "field" : "INFO"
        },
    "AF_min": {
        "tag_value": 0.02,
        "filter_name" : "balsamic_low_af",
        "field" : "INFO"
        },
    "varcaller_name" : "VarDict",
    "filter_type" : "general",
    "analysis_type" : "tumor_only",
    "description" : "General purpose filters used for filtering VarDict"
}