# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

""" Script to include variant caller information to the VCF info field and output to new VCF file """

import click
import gzip
from cyvcf2 import VCF, Writer


@click.command()
@click.option(
    "-i",
    "--input_vcf",
    type=click.Path(exists=True),
    required=True,
    help="Input variant called VCF file",
)
@click.option(
    "-o",
    "--output_vcf",
    type=click.Path(exists=False),
    required=True,
    help="New Output VCF with edited INFO",
)
@click.option(
    "-c",
    "--variant_caller",
    type=str,
    required=True,
    help="Variant caller that generated the VCF file",
)
def edit_vcf_info(input_vcf, output_vcf, variant_caller):
    """Add variant-caller text to the INFO field in the VCF file"""
    vcf = VCF(input_vcf)
    vcf.add_info_to_header(
        {
            "ID": "FOUND_IN",
            "Description": "VCF file in which the variant was found",
            "Type": "String",
            "Number": ".",
        }
    )

    new_vcf = Writer(output_vcf, vcf)

    with open(output_vcf, "wb"):
        for variant in vcf:
            variant.INFO["FOUND_IN"] = variant_caller + "|" + output_vcf
            new_vcf.write_record(variant)
    new_vcf.close()
    vcf.close()


if __name__ == "__main__":
    edit_vcf_info()
