#!/usr/bin/env python
import click
import gzip
import io
from pathlib import Path
from typing import Dict, List, Tuple

# Chromosome names allowed in baf files
ALLOWED_CHR_LIST = [
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
    "20",
    "21",
    "22",
    "X",
    "Y",
    "MT",
]

# How many variants to skip per zoom-level
BAF_SKIP_N = {"o": 135, "a": 30, "b": 8, "c": 3, "d": 1}

# The size of each coverage region per zoom-level
COV_WINDOW_SIZES = {"o": 100000, "a": 25000, "b": 5000, "c": 1000, "d": 100}


@click.group()
@click.option(
    "-o",
    "--output-file",
    required=True,
    type=click.Path(exists=False),
    help="Name of output-file.",
)
def cli(output_file):
    """
    GENS pre-processing tool
    """
    pass


@cli.command()
@click.option(
    "-v",
    "--vcf-file",
    required=True,
    type=click.Path(exists=True),
    help="Input VCF from germline-caller with SNVs & InDels from DNAscope, called with --given gnomad_af_0.05.vcf ",
)
def calc_bafs(vcf_file: str):
    """
    Processes vcf-file from DNAscope into a bed-file format for GENS, with different number of variants for each zoom-level.

    Args:
        vcf_file: From DNAscope created with --given argument using gnomad_af_min_0.05.vcf.

    Outputs bed-file in file-name specified in output-file.
    """
    print("Calculating BAFs from VCF...")
    # Step 1: Read file
    vcf_file_path: Path = Path(vcf_file)
    if vcf_file_path.suffix == ".gz":
        with gzip.open(vcf_file_path, "rt") as file:
            vcf_lines = file.read().splitlines()
    else:
        vcf_lines: List = vcf_file.read_text().splitlines()

    # Step 2: Processing each variant
    variants: Dict = process_variants(vcf_lines)

    # Step 3: Writing output to a file
    output_file: str = click.get_current_context().parent.params["output_file"]
    with open(output_file, "w") as baf_out:
        for prefix in BAF_SKIP_N:
            req_skip_count = int(BAF_SKIP_N[prefix])
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


def process_variants(vcf_lines):
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
        valid_variants: Dict = process_variants(vcf_lines)
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
        if v_chrom not in ALLOWED_CHR_LIST:
            illegal_chromosomes[v_chrom] = illegal_chromosomes.get(v_chrom, 0) + 1
            continue

        variants[variant_id] = variant
        variant_id += 1

    if count_invalid_vars:
        print(f"Warning: Can't calc AF for a number of variants: {count_invalid_vars}.")
    if illegal_chromosomes:
        print(
            f"Warning: A number of variants have illegal chromosomes and will be skipped: {illegal_chromosomes}."
        )
    return variants


def extract_variant_info(variant: str):
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
    fields: List = variant.split("\t")
    try:
        variant_info = {
            "chr": fields[0],
            "start": fields[1],
            "ref": fields[3],
            "alt": fields[4],
            "sample": fields[9],
        }
    except IndexError:
        raise ValueError("Invalid variant string format")

    if variant_info["sample"].split(":")[0] == "./.":
        return None

    try:
        allele_depths = variant_info["sample"].split(":")[1]
        VD = int(allele_depths.split(",")[1])
        DP = int(variant_info["sample"].split(":")[2])
        variant_info["af"] = round(VD / DP, 6)
    except (ValueError, ZeroDivisionError, IndexError):
        return None

    return variant_info


@cli.command()
@click.option(
    "-c",
    "--normalised-coverage-path",
    required=True,
    type=click.Path(exists=True),
    help="Input normalised coverage from GATK DenoiseReadCounts.",
)
def calc_cov(normalised_coverage_path: str) -> None:
    """
    Calculate coverage data.

    Args:
        normalised_coverage: Path to normalised coverage file.

    Returns:
        None
    """
    print("Calculating coverage data")
    normalised_coverage_path = Path(normalised_coverage_path)
    output_file = click.get_current_context().parent.params["output_file"]
    with open(output_file, "w") as cov_out:
        for prefix in COV_WINDOW_SIZES:
            window_size = COV_WINDOW_SIZES[prefix]
            generate_cov_bed(normalised_coverage_path, window_size, prefix, cov_out)


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
    chr: str = chr_start_stop_ratio[0]
    start, end = int(chr_start_stop_ratio[1]), int(chr_start_stop_ratio[2])
    log2_ratio: float = float(chr_start_stop_ratio[3])
    return chr, start, end, log2_ratio


def mean(nums: list) -> float:
    """
    Calculate the mean of a list of numbers.

    Args:
        nums: List of numbers.

    Returns:
        Mean value.
    """
    return sum(nums) / len(nums)

def check_if_start_new_region(prefix, window_size, chrom, region_chrom, region_start, region_end, reg_ratios, cov_out):

    start_new_region: bool = False
    region_size: int = region_end - region_start + 1
    if region_size == window_size:
        # Region size matches window size
        write_coverage_region(
            prefix, region_chrom, region_start, region_end, reg_ratios, cov_out
        )
        start_new_region: bool = True

    if chrom != region_chrom:
        # New chromosome
        write_coverage_region(
            prefix, region_chrom, region_start, region_end, reg_ratios, cov_out
        )
        start_new_region: bool = True

    if region_size > window_size:
        # Region size larger due to incomplete genome reference
        # example: window_size = 300:
        # start1 -- 100bases -- end1 + (region_start = start1)
        # start2 -- 100bases -- end2 + (region_end = end2)
        # -- gap 500 bases --
        # start3 -- 100bases -- end3 (region_end = end3)
        # start4 -- 100bases -- end4
        # Step 1: Revert values to before incompleteness and write
        # Step 2: Start current region from previous cov element
        # end_lists = [x, y, z, before_gap (-3), after_gap (-2), current (-1)]
        previous_region_end: int = region_end_list[-3]
        previous_end: int = region_end_list[-2]
        previous_ratio: float = reg_ratios.pop()
        write_coverage_region(
            prefix, region_chrom, region_start, previous_region_end, reg_ratios, cov_out
        )
        region_start: int = region_start_list[-2]
        reg_ratios: List = [previous_ratio]
        region_end_list, region_start_list = [previous_end, end], [region_start, start]
        start_new_region: bool = False

    return start_new_region


def generate_cov_bed(
    normalised_coverage_path: Path,
    window_size: int,
    prefix: str,
    cov_out: io.TextIOWrapper,
) -> None:
    """
    Merge coverage data into coverage regions for GENS.

    Args:
        normalised_coverage: Path to normalised coverage file.
        window_size: Size of the coverage window.
        prefix: Prefix for the output.
        cov_out: Output file.

    Returns:
        None
    """

    normalised_coverage: List = normalised_coverage_path.read_text().splitlines()

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
            if window_size == 100:
                write_coverage_region(
                    prefix, region_chrom, region_start, region_end, reg_ratios, cov_out
                )
                start_new_region: bool = True
            continue
        else:
            chrom, start, end, log2_ratio = extract_coverage_line_values(coverage_line)

        region_size: int = end - region_start + 1
        if region_size == window_size:
            # Region size matches window size
            # Step 1: Write region from current line
            # Step 2: Start new region from new line
            reg_ratios.append(log2_ratio)
            write_coverage_region(
                prefix, region_chrom, region_start, end, reg_ratios, cov_out
            )
            start_new_region: bool = True
            continue

        if region_size > window_size:
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
