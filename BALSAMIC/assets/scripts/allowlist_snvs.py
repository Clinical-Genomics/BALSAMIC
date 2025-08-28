#!/usr/bin/env python
from typing import Dict, List, Tuple
import click
import gzip
import csv

# --- Constants ---
INFO_HEADERS = {
    "AllowlistStatus": "Variant allowlisted based on CLNSIG/ONC or manually curated clinical list",
    "AllowlistedFilters": "Original filters for allowlisted variants that were overridden",
}
VCF_FIELDS = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", "format"]



@click.command()
@click.argument("vcf_path", type=click.Path(exists=True))
@click.option(
    "--allowlist-path",
    type=click.Path(exists=True),
    default=None,
    help="Optional path to allowlist VCF.",
)
@click.argument("output_file", type=click.Path())
def main(vcf_path: str, allowlist_path: str, output_file: str) -> None:
    """VCF processor with allowlist annotation support."""
    variants = read_variants(vcf_path)



if __name__ == "__main__":
    main()
