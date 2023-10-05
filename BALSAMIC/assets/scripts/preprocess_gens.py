#!/usr/bin/env python
import click
from pathlib import Path
import io


ALLOWED_CHR_LIST = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y', 'MT']
BAF_SKIP_N = {"o": 135,
              "a": 30,
              "b": 8,
              "c": 3,
              "d": 1}
COV_WINDOW_SIZES = {"o": 100000,
                    "a": 25000,
                    "b": 5000,
                    "c": 1000,
                    "d": 100}

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
    help="Input VCF from germline-caller with SNVs & InDels from DNAscope.",
)
def calc_bafs(vcf_file: str):
    """

    Args:
        vcf_file:

    Returns:

    """
    print("Calculating BAFs from VCF...")
    vcf_file = Path(vcf_file)
    vcf_lines = vcf_file.read_text().splitlines()

    # Step 2: Initializing variables
    count_invalid_vars = 0
    illegal_chromosomes = {}
    variants = {}
    variant_id = 0

    # Step 3: Processing each variant
    for variant_line in vcf_lines:
        if variant_line.startswith("#"):
            continue

        variant = extract_variant_info(variant_line)
        if not variant:
            count_invalid_vars += 1
            continue

        v_chrom = variant["chr"]
        if v_chrom not in ALLOWED_CHR_LIST:
            illegal_chromosomes[v_chrom] = illegal_chromosomes.get(v_chrom, 0) + 1
            continue

        variants[variant_id] = variant
        variant_id += 1

    # Step 4: Printing warnings
    if count_invalid_vars:
        print(f"Warning: Can't calc AF for a number of variants: {count_invalid_vars}.")
    if illegal_chromosomes:
        print(f"Warning: A number of variants have illegal chromosomes and will be skipped: {illegal_chromosomes}.")

    # Step 5: Writing output to a file
    output_file = click.get_current_context().parent.params['output_file']
    with open(output_file, 'w') as baf_out:
        for prefix in BAF_SKIP_N:
            req_skip_count = int(BAF_SKIP_N[prefix])
            skip_count = int(req_skip_count)
            for variant_id in variants:
                variant = variants[variant_id]
                if skip_count == req_skip_count:
                    v_start = int(variant["start"])
                    v_chrom = variant["chr"]
                    v_af = variant["af"]
                    baf_out.write(f"{prefix}_{v_chrom}\t{v_start}\t{v_start + 1}\t{v_af}\n")
                    skip_count = 0
                else:
                    skip_count += 1


def extract_variant_info(variant):
    variant_info = {}
    variant = variant.split("\t")
    variant_info["chr"] = variant[0]
    variant_info["start"] = variant[1]
    variant_info["dbsnp"] = variant[2]
    variant_info["ref"] = variant[3]
    variant_info["alt"] = variant[4]
    variant_sample = variant[9]
    variant_info["sample"] = variant_sample
    if variant_sample.split(":")[0] == "./.":
        return None
    allele_depths = variant_sample.split(":")[1]
    #print(variant_sample)
    VD = int(allele_depths.split(",")[1])
    DP = int(variant_sample.split(":")[2])
    try:
        variant_info["af"] = float(round(VD/DP, 6))
    except ZeroDivisionError:
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
def calculate_coverage_data(normalised_coverage_path: str) -> None:
    """
    Calculate coverage data.

    Args:
        normalised_coverage: Path to normalised coverage file.

    Returns:
        None
    """
    print("Calculating coverage data")
    normalised_coverage_path = Path(normalised_coverage_path)
    output_file = click.get_current_context().parent.params['output_file']
    with open(output_file, 'w') as cov_out:
        for prefix in COV_WINDOW_SIZES:
            window_size = COV_WINDOW_SIZES[prefix]
            generate_cov_bed(normalised_coverage_path, window_size, prefix, cov_out)


def mean(nums: list) -> float:
    """
    Calculate the mean of a list of numbers.

    Args:
        nums: List of numbers.

    Returns:
        Mean value.
    """
    return sum(nums) / len(nums)


def generate_cov_bed(
    normalised_coverage_path: Path,
    window_size: int,
    prefix: str,
    cov_out: io.TextIOWrapper
) -> None:
    """
    Generate coverage data.

    Args:
        normalised_coverage: Path to normalised coverage file.
        window_size: Size of the coverage window.
        prefix: Prefix for the output.
        cov_out: Output file.

    Returns:
        None
    """
    region_start, region_end, region_chr, force_end = None, None, None, False
    reg_ratios = []

    normalised_coverage = normalised_coverage_path.read_text().splitlines()
    for line in normalised_coverage:
        if line.startswith('@') or line.startswith('CONTIG'):
            continue
        chr_start_stop_ratio = line.strip().split('\t')
        start, end = int(chr_start_stop_ratio[1]), int(chr_start_stop_ratio[2])
        chr = chr_start_stop_ratio[0]
        log2_ratio = chr_start_stop_ratio[3]

        orig_end = end

        if not region_start:
            region_start, region_end, region_chr = start, end, chr

        if chr == region_chr:
            if start - region_end < window_size:
                reg_ratios.append(float(log2_ratio))
                region_end = end
            else:
                force_end = True
                end = region_end
        else:
            force_end = True
            end = region_end

        if end - region_start + 1 >= window_size or force_end:
            mid_point = region_start + (end - region_start) // 2
            cov_out.write(f"{prefix}_{region_chr}\t{mid_point-1}\t{mid_point}\t{mean(reg_ratios)}\n")
            region_start, region_end, region_chr, reg_ratios = None, None, None, []

        if force_end:
            region_start, region_end, region_chr = start, orig_end, chr
            reg_ratios.append(float(log2_ratio))
            force_end = None


if __name__ == "__main__":
    cli()