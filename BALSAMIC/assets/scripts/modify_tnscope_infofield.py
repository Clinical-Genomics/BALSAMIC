#!/usr/bin/env python
import vcfpy
import click
import sys
import logging
from typing import List, Optional

LOG = logging.getLogger(__name__)


def summarize_ad_to_dp(ad_list):
    """
    Summarizes the AD (allelic depth) field into total DP (read depth).

    Parameters:
    ad_list (list): List of read depths supporting each allele, [ref_depth, alt1_depth, alt2_depth, ...]

    Returns:
    int: Total read depth (DP) across all alleles.
    """
    if ad_list is None:
        return 0  # Return 0 if AD field is not present
    return sum(ad_list)


@click.command()
@click.argument("input_vcf", type=click.Path(exists=True))
@click.argument("output_vcf", type=click.Path())
def process_vcf(input_vcf: str, output_vcf: str):
    """
    Processes the input VCF file and writes the updated information to the output VCF file.

    INPUT_VCF: Path to the input VCF file.
    OUTPUT_VCF: Path to the output VCF file.
    """

    # Open the input VCF file
    reader: vcfpy.Reader = vcfpy.Reader.from_path(input_vcf)

    # Ensure the sample name is 'TUMOR'
    sample_name: str = reader.header.samples.names[0]
    if sample_name != "TUMOR":
        LOG.warning(
            f"Error: The first sample is named '{sample_name}', but 'TUMOR' is expected."
        )
        sys.exit(1)

    # Add AF and DP fields to the header if not already present
    if "AF" not in reader.header.info_ids():
        reader.header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "AF"),
                    ("Number", "A"),
                    ("Type", "Float"),
                    ("Description", "Allele Frequency"),
                ]
            )
        )

    if "DP" not in reader.header.info_ids():
        reader.header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", "DP"),
                    ("Number", "1"),
                    ("Type", "Integer"),
                    ("Description", "Total Depth"),
                ]
            )
        )

    # Open the output VCF file for writing
    with vcfpy.Writer.from_path(output_vcf, reader.header) as writer:
        # Loop through each record (variant)
        for record in reader:
            # Get the TUMOR sample data
            sample_index: int = reader.header.samples.names.index(sample_name)
            tumor_call: vcfpy.Call = record.calls[sample_index]

            # Check and process AD field
            tumor_ad: Optional[List[int]] = tumor_call.data.get(
                "AD", None
            )  # AD is a list [ref_count, alt_count]
            if tumor_ad is None:
                LOG.warning(
                    f"Warning: AD field is missing for record at position {record.POS} on {record.CHROM}"
                )
            else:
                record.INFO["DP"] = summarize_ad_to_dp(tumor_ad)

            # Check and process AF field
            tumor_af: Optional[float] = tumor_call.data.get("AF", None)
            if tumor_af is None:
                LOG.warning(
                    f"Warning: AF field is missing for record at position {record.POS} on {record.CHROM}"
                )
                record.INFO["AF"] = [0.0]  # Default AF to 0.0 if missing
            else:
                record.INFO["AF"] = [tumor_af]  # Wrap AF in a list

            # Write the updated record to the output VCF file
            writer.write_record(record)

    click.echo(f"VCF file processed and saved to {output_vcf}")


if __name__ == "__main__":
    process_vcf()
