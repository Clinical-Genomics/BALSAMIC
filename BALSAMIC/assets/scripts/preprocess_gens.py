#!/usr/bin/env python
import click
import io
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Union
from statistics import mean

from BALSAMIC.constants.analysis import SequencingType
from BALSAMIC.constants.tools import GENS_PARAMS
from BALSAMIC.utils.io import read_vcf_file

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.INFO)


@click.group()
@click.option(
    "-o",
    "--output-file",
    required=True,
    type=click.Path(exists=False),
    help="Name of output-file.",
)
@click.option(
    "-s",
    "--sequencing-type",
    required=True,
    type=click.Choice([SequencingType.WGS]),
    help="Sequencing type used.",
)
@click.pass_context
def cli(ctx: click.Context, output_file: str, sequencing_type: SequencingType):
    """GENS pre-processing tool."""
    ctx.ensure_object(dict)
    ctx.obj["output_file"] = output_file
    ctx.obj["sequencing_type"] = sequencing_type


@cli.command()
@click.pass_context
@click.option(
    "-v",
    "--vcf-file-path",
    required=True,
    type=click.Path(exists=True),
    help="Input VCF from germline-caller with SNVs & InDels from DNAscope, called with --given gnomad_af_0.05.vcf ",
)
def calculate_bafs(ctx: click.Context, vcf_file_path: str):
    """
    Processes vcf-file from DNAscope into a bed-file format for GENS, with different number of variants for each zoom-level.

    Args:
        vcf_file: From DNAscope created with --given argument using gnomad_af_min_0.05.vcf.

    Outputs bed-file in file-name specified in output-file.
    """

    LOG.info("Calculating BAFs from VCF.")
    LOG.info("Reading VCF file...")
    vcf_lines: List = read_vcf_file(vcf_file_path)
    LOG.info("Extracting and processing variant info...")
    variants: Dict = get_valid_variants(vcf_lines)
    LOG.info("Writing variant b-allele-frequencies to output...")
    output_file: Path = Path(ctx.obj["output_file"])
    sequencing_type: SequencingType = ctx.obj["sequencing_type"]
    write_b_allele_output(variants, output_file, sequencing_type)


def get_valid_variants(vcf_lines: List[str]) -> Dict:
    """
    Process VCF lines to extract valid variants.

    Args:
        vcf_lines (List[str]): List of VCF lines to be processed.

    Returns:
        dict: A dictionary containing valid variants with variant IDs as keys and variant info as values.

    This function takes a list of VCF lines, processes them, and extracts valid variants.
    Variants are stored in a dictionary with variant IDs as keys and variant information as values.
    If any invalid variants are encountered, appropriate warnings are printed.

    Example usage:
        vcf_lines: List = [...]  # List of VCF lines
        valid_variants: Dict = get_valid_variants(vcf_lines)
    """
    count_invalid_vars: int = 0
    illegal_chromosomes: Dict = {}
    variants = {}
    variant_id: int = 0
    for variant_line in vcf_lines:
        if variant_line.startswith("#"):
            continue

        variant: Dict = extract_variant_info(variant_line)
        if not variant:
            count_invalid_vars += 1
            continue

        v_chrom = variant["chr"]
        if v_chrom not in GENS_PARAMS["ALLOWED_CHR_LIST"]:
            illegal_chromosomes[v_chrom] = illegal_chromosomes.get(v_chrom, 0) + 1
            continue

        variants[variant_id] = variant
        variant_id += 1

    if count_invalid_vars:
        LOG.warning(
            f"Warning: Can't calc AF for a number of variants: {count_invalid_vars}."
        )
    if illegal_chromosomes:
        LOG.warning(
            f"Warning: A number of variants have illegal chromosomes and will be skipped: {illegal_chromosomes}."
        )
    return variants


def extract_variant_info(variant: str) -> Union[Dict, None]:
    """
    Extracts genetic variant information.

    Args:
        variant (str): Tab-separated string representing a genetic variant.

    Returns:
        dict or None: Dictionary with variant details ('chr', 'start', 'ref', 'alt', 'sample', 'af').
        Returns None for uninformative samples or division by zero.

    Raises:
        ValueError: If variant string lacks expected fields.

    Example:
        variant_string = "chr1\t1000\t.\tA\tT\t.\t.\t.\tGT:AD:DP:GQ:PL\t0/1:10,5:15:45:45,0,20"
        extract_variant_info(variant_string)
        {'chr': 'chr1', 'start': '1000', 'ref': 'A', 'alt': 'T', 'sample': '0/1:10,5:15:45:45,0,20', 'af': 0.333333}
    """
    fields: List[str] = variant.split("\t")
    try:
        variant_info = {
            "chr": fields[0],
            "start": fields[1],
            "ref": fields[3],
            "alt": fields[4],
            "sample": fields[9],
        }
    except IndexError:
        error_message = "Invalid variant string format"
        LOG.error(f"{error_message}")
        raise ValueError(f"{error_message}")

    if variant_info["sample"].split(":")[0] == "./.":
        return None

    try:
        allele_depths: List[str] = variant_info["sample"].split(":")[1]
        VD = int(allele_depths.split(",")[1])
        DP = int(variant_info["sample"].split(":")[2])
        variant_info["af"] = round(VD / DP, 6)

    except (ValueError, ZeroDivisionError, IndexError):
        return None

    return variant_info


def write_b_allele_output(
    variants: Dict, output_file: Path, sequencing_type: SequencingType
):
    """
    Writes B-allele frequency (BAF) output to a file for each level of GENS zoom specified in prefix of BAF_SKIP_N.

    Args:
        variants (Dict): A dictionary containing variant information. Each variant should have the following keys:
            - 'start' (int): The start position of the variant.
            - 'chr' (str): The chromosome where the variant is located.
            - 'af' (float): The allele frequency of the variant.

        output_file (Path): The file path where the BAF output will be written.

        sequencing_type (SequencingType): The sequencing type used in analysis.

    Returns:
        None

    Writes BAF information for each variant in the provided `variants` dictionary to the specified `output_file`.
    The output is formatted as follows:
        <prefix>_<chromosome>    <start>    <end>    <allele_frequency>

    Note:
        This function uses a predefined constant BAF_SKIP_N, which determines how many variants to skip before writing.
    """
    SEQ_GENS_PARAMS: Dict = GENS_PARAMS["SEQUENCING_TYPE"][sequencing_type]
    BAF_SKIP_N: Dict = SEQ_GENS_PARAMS["BAF_SKIP_N"]

    with open(output_file.as_posix(), "w") as baf_out:
        for prefix, req_skip_count in BAF_SKIP_N.items():
            skip_count = int(req_skip_count)
            for variant_id in variants:
                variant: Dict = variants[variant_id]
                if skip_count == req_skip_count:
                    v_start: int = int(variant["start"])
                    v_chrom: str = variant["chr"]
                    v_af: float = variant["af"]
                    baf_out.write(
                        f"{prefix}_{v_chrom}\t{v_start}\t{v_start + 1}\t{v_af}\n"
                    )
                    skip_count = 0
                else:
                    skip_count += 1


@cli.command()
@click.pass_context
@click.option(
    "-c",
    "--normalised-coverage-path",
    required=True,
    type=click.Path(exists=True),
    help="Input normalised coverage from GATK DenoiseReadCounts.",
)
def create_coverage_regions(ctx: click.Context, normalised_coverage_path: str) -> None:
    """
    Calculate coverage data.

    Args:
        normalised_coverage_path: Path to normalised coverage file.

    Returns:
        None
    """
    LOG.info("Creating coverage regions for GENS.")
    normalised_coverage_path = Path(normalised_coverage_path)
    output_file: Path = Path(ctx.obj["output_file"])
    sequencing_type: SequencingType = ctx.obj["sequencing_type"]
    SEQ_GENS_PARAMS: Dict = GENS_PARAMS["SEQUENCING_TYPE"][sequencing_type]
    COV_REGION_SIZES: Dict = SEQ_GENS_PARAMS["COV_REGION_SIZES"]
    with open(output_file.as_posix(), "w") as cov_out:
        for prefix, region_size in COV_REGION_SIZES.items():
            LOG.info(
                f"Creating regions for prefix: {prefix}, region_size: {region_size}."
            )
            generate_cov_bed(normalised_coverage_path, region_size, prefix, cov_out)


def write_coverage_region(
    prefix: str,
    region_chr: str,
    region_start: int,
    region_end: int,
    reg_ratios: List[float],
    cov_out: io.TextIOWrapper,
) -> None:
    """
    Write coverage region information to an output file.

    Args:
        prefix (str): Prefix for the output.
        region_chr (str): Chromosome for the region.
        region_start (int): Start position of the region.
        region_end (int): End position of the region.
        reg_ratios (List[float]): List of log2 ratios.
        cov_out (io.TextIOWrapper): Output file.

    Returns:
        None
    """
    mid_point = region_start + (region_end - region_start) // 2
    cov_out.write(
        f"{prefix}_{region_chr}\t{mid_point - 1}\t{mid_point}\t{mean(reg_ratios)}\n"
    )


def extract_coverage_line_values(coverage_line: str) -> Tuple[str, int, int, float]:
    """
    Extract coverage region values from a coverage line.

    Args:
        coverage_line (str): A line containing coverage and genomic position information.

    Returns:
        Tuple[str, int, int, float]: Extracted values (chr, start, end, log2_ratio).
    """
    # Extract coverage region values
    chr_start_stop_ratio: List = coverage_line.strip().split("\t")
    chrom: str = chr_start_stop_ratio[0]
    start, end = int(chr_start_stop_ratio[1]), int(chr_start_stop_ratio[2])
    log2_ratio: float = float(chr_start_stop_ratio[3])
    return chrom, start, end, log2_ratio


def generate_cov_bed(
    normalised_coverage_path: Path,
    region_size_requested: int,
    prefix: str,
    cov_out: io.TextIOWrapper,
) -> None:
    """
    Merge coverage data into coverage regions for GENS.

    Args:
        normalised_coverage_path: Path to normalised coverage file.
        region_size_requested: Size of the coverage region.
        prefix: Prefix for the output.
        cov_out: Output file.

    Returns:
        None
    """

    normalised_coverage: List = normalised_coverage_path.read_text().splitlines()
    minimum_region_size: int = GENS_PARAMS["MINIMUM_REGION_SIZE"]

    region_chrom, chrom, region_start, region_end, end, log2_ratio = None
    first_cov_line: bool = True
    start_new_region: bool = False
    for coverage_line in normalised_coverage:
        if coverage_line.startswith("@") or coverage_line.startswith("CONTIG"):
            continue

        if first_cov_line or start_new_region:
            (
                region_chrom,
                region_start,
                region_end,
                log2_ratio,
            ) = extract_coverage_line_values(coverage_line)
            reg_ratios: List = [log2_ratio]
            first_cov_line: bool = False
            start_new_region: bool = False
            if region_size_requested == minimum_region_size:
                write_coverage_region(
                    prefix, region_chrom, region_start, region_end, reg_ratios, cov_out
                )
                start_new_region: bool = True
            continue
        else:
            chrom, _, end, log2_ratio = extract_coverage_line_values(coverage_line)

        region_size: int = end - region_start + 1
        if region_size == region_size_requested:
            # Region size matches requested region size
            # Step 1: Write region from current line
            # Step 2: Start new region from new line
            reg_ratios.append(log2_ratio)
            write_coverage_region(
                prefix, region_chrom, region_start, end, reg_ratios, cov_out
            )
            start_new_region: bool = True
            continue

        if region_size > region_size_requested:
            # Region size larger due to incomplete genome reference
            # Step 1:  Write region from previous line
            # Step 2: Start new region from current line
            # Conceptual example: window_size = 300:
            # start1 -- 100bases -- end1 + (region_start = start1, region_end = end1)
            # start2 -- 100bases -- end2 + (region_start = start1, region_end = end2)
            # -- gap 500 bases --
            # start3 -- 100bases -- end3
            write_coverage_region(
                prefix, region_chrom, region_start, region_end, reg_ratios, cov_out
            )
            start_new_region: bool = True

        if chrom != region_chrom:
            # New chromosome:
            # Step 1: Write region from previous line
            # Step 2: Start new region from current line
            write_coverage_region(
                prefix, region_chrom, region_start, region_end, reg_ratios, cov_out
            )
            start_new_region: bool = True

        if start_new_region:
            (
                region_chrom,
                region_start,
                region_end,
                log2_ratio,
            ) = extract_coverage_line_values(coverage_line)
            reg_ratios: List = [log2_ratio]
        else:
            region_end: int = end
            reg_ratios.append(log2_ratio)

    # Output last line:
    write_coverage_region(
        prefix, region_chrom, region_start, region_end, reg_ratios, cov_out
    )


if __name__ == "__main__":
    cli()
