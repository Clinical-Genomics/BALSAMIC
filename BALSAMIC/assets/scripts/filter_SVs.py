#!/usr/bin/env python

import vcfpy
import click
import logging
from pathlib import Path


def find_ctg(sample_info):
    """
    Looks for contig in the variant with the highest AF.

    :param sample_info: SAMPLE_PASS_INFO[variant with highest AF]
    :return: bool(if contig exists for max-AF variant or not).
    """
    fields = sample_info.split("|")
    ctg = False
    for field in fields:
        field_name_value = field.split(":")
        field_name = field_name_value[0]
        if field_name == "CTG":
            ctg_value = field_name_value[1]
            if ctg_value != ".":
                ctg = True
    return ctg


def calc_af_get_ctg(pass_sample, pass_info):
    """
    Inputs fields from TIDDIT variant INFO field,
    calculates max AF for merged variants and looks for contig in max-af variant.

    :param pass_sample: TUMOR_PASS_SAMPLE / NORMAL_PASS_SAMPLE
    :param pass_info: TUMOR_PASS_INFO / NORMAL_PASS_INFO

    :return: float(maximum allele frequency), and bool(if contig exists for max-AF variant or not).
    """
    max_af = 0
    final_ctg = False
    any_ctg = False
    for idx, variant in enumerate(pass_sample):
        fields = variant.split("|")

        ctg = find_ctg(pass_info[idx])

        for field in fields:
            field_name_value = field.split(":")
            field_name = field_name_value[0]
            if field_name == "COV":
                B1_cov = field_name_value[1]
                B2_cov = field_name_value[3]
                total_cov = int(B1_cov) + int(B2_cov)
            if field_name == "DV":
                DV = int(field_name_value[1])
            if field_name == "RV":
                RV = int(field_name_value[1])
        try:
            af = float((DV + RV) / total_cov)
        except Exception as e:
            logging.warning(f"Exception: {e} setting AF to 0")
            af = 0
        if ctg:
            any_ctg = True
        if af > max_af:
            max_af = af
            final_ctg = ctg
    if max_af == 0:
        final_ctg = any_ctg
    return max_af, final_ctg


def filter_vcf(vcf_path: Path):
    """
    Adds soft filters based on variant presence in normal.

    :param vcf_path: Path to VCF

    Outputs to standard out.
    """

    vcf = vcfpy.Reader.from_path(vcf_path)

    # Update VCF header
    vcf.header.add_info_line(
        vcfpy.OrderedDict(
            [
                ("ID", "AF_T_MAX"),
                ("Type", "Float"),
                (
                    "Description",
                    "Max AF in tumor, for rows with merged overlapping variants",
                ),
            ]
        )
    )

    vcf.header.add_info_line(
        vcfpy.OrderedDict(
            [
                ("ID", "AF_N_MAX"),
                ("Type", "Float"),
                (
                    "Description",
                    "Max AF in normal, for rows with merged overlapping variants",
                ),
            ]
        )
    )

    vcf.header.add_filter_line(
        vcfpy.OrderedDict(
            [
                ("ID", "normal_variant"),
                ("Description", "AF_T_MAX == 0 and ctg_t == False"),
            ]
        )
    )

    vcf.header.add_filter_line(
        vcfpy.OrderedDict(
            [("ID", "high_normal_af"), ("Description", "AF_N_MAX > 0.25")]
        )
    )

    vcf.header.add_filter_line(
        vcfpy.OrderedDict(
            [
                ("ID", "high_normal_af_fraction"),
                ("Description", "(AF_N_MAX / AF_T_MAX) > 0.25"),
            ]
        )
    )

    vcf.header.add_filter_line(
        vcfpy.OrderedDict(
            [
                ("ID", "in_normal"),
                ("Description", "ctg_n == True and AF_N_MAX == 0 and AF_T_MAX <= 0.25"),
            ]
        )
    )

    writer = vcfpy.Writer.from_path("/dev/stdout", vcf.header)

    # Set soft filters for variants based on presence in the normal sample
    for sv in vcf:
        info = dict([*sv.INFO.items()])

        # Set default values
        ctg_t, ctg_n = False, False
        af_t, af_n = 0, 0

        if "TUMOR_PASS_SAMPLE" in info:
            pass_sample = info["TUMOR_PASS_SAMPLE"]
            pass_info = info["TUMOR_PASS_INFO"]
            af_t, ctg_t = calc_af_get_ctg(pass_sample, pass_info)
        if "NORMAL_PASS_SAMPLE" in info:
            pass_sample = info["NORMAL_PASS_SAMPLE"]
            pass_info = info["NORMAL_PASS_INFO"]
            af_n, ctg_n = calc_af_get_ctg(pass_sample, pass_info)

        # Set filter statuses
        if af_t == 0 and not ctg_t:
            sv.add_filter("normal_variant")
        else:
            # Regardless of CTG, set filter if AF_T / AF_N > 0.25
            if af_t > 0:
                if float(af_n / af_t) > 0.25:
                    normal_af_status = "failed"
                    sv.add_filter("high_normal_af_fraction")
                else:
                    normal_af_status = "pass"

            # Set filter if AF_N > 0.25
            if af_n > 0.25:
                sv.add_filter("high_normal_af")

            # Set filter if CTG_N = True, and AF_N is not pass
            if ctg_n and af_n == 0 and af_t <= 0.25:
                sv.add_filter("in_normal")

        sv.INFO["AF_T_MAX"] = [round(af_t, 4)]
        sv.INFO["AF_N_MAX"] = [round(af_n, 4)]

        writer.write_record(sv)


def get_bnd_id(info):
    """
    Creates and returns unique BND ID based on information in INFO field.

    :param info: Variant INFO field
    :return: unique breakend ID based on SVID and RegionA, and the BND number (1 or 2)
    """
    # example: SV_2414_1
    sv_id = info["TUMOR_PASS_CHROM"][0].split("|")[0]

    # extract sv_id and breakend number
    sv_id_split = sv_id.split("_")
    sv_id_name = "_".join(sv_id_split[0:2])
    sv_id_num = sv_id_split[2]

    # define unique breakend_name
    region_a = str(info["REGIONA"][0]) + "_" + str(info["REGIONA"][1])
    bnd_id = f"{sv_id_name}_{region_a}"
    return bnd_id, sv_id_num


def rescue_t_n_mixed_bnds(vcf_path):
    """
    Removes soft filters and sets PASS to BND variants with mixed soft / PASS filters set.
    Annotates BND variants with unique ID.

    :param vcf_path: Path to VCF

    Outputs to standard out.
    """
    vcf_start = vcfpy.Reader.from_path(vcf_path)

    # First read of VCF-file:
    # define dict with bnd_id - bnd_num - FILTER
    bnd_filter_dict = {}
    for sv in vcf_start:
        info = dict([*sv.INFO.items()])

        svtype = info["SVTYPE"]
        # only relevant for bnd-variants
        if svtype != "BND":
            continue
        # no need to rescue normal variants
        if "TUMOR_PASS_CHROM" not in info:
            continue

        bnd_id, sv_id_num = get_bnd_id(info)

        if bnd_id not in bnd_filter_dict:
            bnd_filter_dict[bnd_id] = {}

        if sv_id_num in bnd_filter_dict[bnd_id]:
            logging.warning(
                f"Conflicting BND-names: {bnd_id}, will not attempt to rescue: {sv}"
            )
            continue

        bnd_filter_dict[bnd_id][sv_id_num] = sv.FILTER

    # Second read fo VCF-file:
    # Annotate and remove BND filters with PASS
    vcf = vcfpy.Reader.from_path(vcf_path)
    vcf.header.add_info_line(
        vcfpy.OrderedDict(
            [
                ("ID", "BND_ID"),
                ("Type", "String"),
                ("Description", "Unique ID for BND-variants, SVID + REGIONA"),
            ]
        )
    )
    vcf.header.add_filter_line(
        vcfpy.OrderedDict(
            [
                ("ID", "PASS"),
                (
                    "Description",
                    "If at least 1 of the 2 BND_ID variants have PASS, set PASS",
                ),
            ]
        )
    )

    writer = vcfpy.Writer.from_path("/dev/stdout", vcf.header)

    for sv in vcf:
        info = dict([*sv.INFO.items()])
        svtype = info["SVTYPE"]
        # only relevant for bnd-variants
        if svtype != "BND":
            writer.write_record(sv)
            continue

        # no need to rescue normal variants
        if "TUMOR_PASS_CHROM" not in info:
            writer.write_record(sv)
            continue

        bnd_id, sv_id_num = get_bnd_id(info)
        sv.INFO["BND_ID"] = [bnd_id]

        filter_bnds = bnd_filter_dict[bnd_id]
        if "1" in filter_bnds and "2" in filter_bnds:
            filter_status1 = filter_bnds["1"]
            filter_status2 = filter_bnds["2"]
            if "PASS" in filter_status1 or "PASS" in filter_status2:
                sv.FILTER = ["PASS"]

        writer.write_record(sv)


@click.group()
def main():
    """SV filtering tool."""
    pass


@main.command("placeholder")
@click.option(
    "-f",
    "--vcf-file",
    required=True,
    type=click.Path(exists=True),
    help="",
)
def do_nothing(vcf_file: Path):
    """For future filtering."""
    pass


@main.command("filter_tiddit_TN")
@click.option(
    "-v",
    "--vcf-file",
    required=True,
    type=click.Path(exists=True),
    help="Input Variant Calling Format(VCF) from merged TIDDIT TUMOR + NORMAL VCF.",
)
@click.option(
    "--method", "-m", required=True, type=click.Choice(["filter", "rescue_bnds"])
)
def filter_tiddit_TN(vcf_file: Path, method: str):
    """Filter TIDDIT variants with Tumor Normal options.
    Outputs filtered VCF in standard-out.

        Args:
            vcf-file: Path; path to vcf-file.

            method: Choose method:

                filter: Set soft-filters for presence of variants in normal.
                rescue_bnds: Rescue BND variants with at least 1 BND with PASS.
    """

    if method == "filter":
        filter_vcf(vcf_file)

    if method == "rescue_bnds":
        rescue_t_n_mixed_bnds(vcf_file)


if __name__ == "__main__":
    main()
