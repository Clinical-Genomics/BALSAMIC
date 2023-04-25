#!/usr/bin/env python

import vcfpy
import click
import logging
from pathlib import Path

LOG = logging.getLogger(__name__)

@click.group()
def cli():
    """
    SV filtering tool.
    """
    pass

@cli.group("tiddit_tn")
@click.option("-v", "--vcf-file", required=True, type=click.Path(exists=True), help="Input Variant Calling Format(VCF) from merged TIDDIT TUMOR + NORMAL VCF.")
@click.pass_context
def tiddit_tn(ctx, vcf_file):
    """
    Manage TIDDIT variants in merged tumor/normal vcf.
    Outputs treated VCF in standard-out.

        Args:
            vcf-file: Path; path to vcf-file.
    """
    ctx.obj = {}
    ctx.obj["vcf_file"] = vcf_file
    pass
@tiddit_tn.command("filter")
@click.pass_context
def filter(ctx: click.Context):
    """
    Add soft-filters based on presence of variants in normal.
    Outputs filtered VCF in standard-out.
    """
    vcf = vcfpy.Reader.from_path(ctx.obj["vcf_file"])

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

        # Extract allele frequency and contig-status for variant in T and N
        if "TUMOR_PASS_SAMPLE" in info:
            pass_sample = info["TUMOR_PASS_SAMPLE"]
            pass_info = info["TUMOR_PASS_INFO"]
            print(type(pass_sample))
            print(type(pass_info))
            af_t, max_af_index = calc_max_af(pass_sample)
            if af_t == 0:
                get_any_contig = True
            else:
                get_any_contig = False
            ctg_t = find_ctg(pass_info, max_af_index, get_any_contig)

        if "NORMAL_PASS_SAMPLE" in info:
            pass_sample = info["NORMAL_PASS_SAMPLE"]
            pass_info = info["NORMAL_PASS_INFO"]
            af_n, max_af_index = calc_max_af(pass_sample)
            if af_n == 0:
                get_any_contig = True
            else:
                get_any_contig = False
            ctg_n = find_ctg(pass_info, max_af_index, get_any_contig)

        # Set filter statuses
        if af_t == 0 and not ctg_t:
            sv.add_filter("normal_variant")
        else:
            # Regardless of CTG, set filter if AF_T / AF_N > 0.25
            if af_t > 0:
                if float(af_n / af_t) > 0.25:
                    sv.add_filter("high_normal_af_fraction")

            # Set filter if AF_N > 0.25
            if af_n > 0.25:
                sv.add_filter("high_normal_af")

            # Set filter if CTG_N = True, AF_N is 0 and AF_T is below 0.25
            if ctg_n and af_n == 0 and af_t <= 0.25:
                sv.add_filter("in_normal")

        sv.INFO["AF_T_MAX"] = [round(af_t, 4)]
        sv.INFO["AF_N_MAX"] = [round(af_n, 4)]

        #writer.write_record(sv)

def is_ctg(variant):
    """

    :return:
    """
    fields = variant.split("|")
    contig = False
    for field in fields:
        field_name_value = field.split(":")
        field_name = field_name_value[0]
        if field_name == "CTG":
            contig_value = field_name_value[1]
            if contig_value != ".":
                contig = True
        return contig
def find_ctg(sample_info, max_af_index, get_any_contig):
    """
    Looks for contig in the variant with the highest AF.

    :param sample_info: SAMPLE_PASS_INFO[variant with highest AF]
    :return: bool(if contig exists for max-AF variant or not).
    """
    contig = False
    if get_any_contig:
        for variant in sample_info:
            contig = is_ctg(variant)
            if contig:
                return contig
    else:
        max_af_variant = sample_info[max_af_index]
        contig = is_ctg(max_af_variant)
    return contig

def calc_max_af(pass_sample):
    """
    Inputs fields from TIDDIT variant INFO field,
    calculates max AF for merged variants and looks for contig in max-af variant.

    :param pass_sample: TUMOR_PASS_SAMPLE / NORMAL_PASS_SAMPLE
    :param pass_info: TUMOR_PASS_INFO / NORMAL_PASS_INFO

    :return: float(maximum allele frequency), and bool(if contig exists for max-AF variant or not).
    """
    max_allele_frequency = 0
    max_af_variant_index = 0
    B1_cov = 0
    B2_cov = 0
    total_cov = 0
    DV = 0
    RV = 0

    for variant_idx, variant in enumerate(pass_sample):
        fields = variant.split("|")
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
            allele_frequency = float((DV + RV) / total_cov)
        except ZeroDivisionError as exc:
            LOG.warning(f"Exception: {exc} setting AF to 0")
            allele_frequency = 0

        if allele_frequency > max_allele_frequency:
            max_af_variant_index = variant_idx
            max_allele_frequency = allele_frequency
    return max_allele_frequency, max_af_variant_index

@tiddit_tn.command("rescue_bnds")
@click.pass_context
def rescue_bnds(ctx: click.Context):
    """
    Rescue BND-variants with at least 1 of 2 BNDs set to PASS.
    Outputs VCF in standard-out.
    """
    vcf_start = vcfpy.Reader.from_path(ctx.obj["vcf_file"])

    # First read of VCF-file:
    # define dict with bnd_id - bnd_num - FILTER
    bnd_filter_dict = {}
    for sv in vcf_start:
        info = dict([*sv.INFO.items()])

        svtype = info["SVTYPE"]
        # only relevant for bnd-variants and tumor-variants
        if svtype == "BND" and "TUMOR_PASS_CHROM" in info:

            bnd_id, sv_id_num = get_bnd_id(info)

            if bnd_id not in bnd_filter_dict:
                bnd_filter_dict[bnd_id] = {}

            if sv_id_num in bnd_filter_dict[bnd_id]:
                logging.warning(
                    f"Conflicting BND-names: {bnd_id}, will not attempt to rescue: {sv}"
                )
                continue

        bnd_filter_dict[bnd_id][sv_id_num] = sv.FILTER

    # Second read of VCF-file:
    # Update VCF header
    vcf = vcfpy.Reader.from_path(vcf_file)
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
    # Annotate and remove BND filters with PASS
    for sv in vcf:
        info = dict([*sv.INFO.items()])
        svtype = info["SVTYPE"]
        # only relevant for bnd and tumor-variants
        if svtype == "BND" and "TUMOR_PASS_CHROM":
            print(type(info))
            bnd_id, sv_id_num = get_bnd_id(info)
            sv.INFO["BND_ID"] = [bnd_id]

            filter_bnds = bnd_filter_dict[bnd_id]
            if "1" in filter_bnds and "2" in filter_bnds and ("PASS" in filter_bnds["1"] or "PASS" in filter_bnds["2"]):
                sv.FILTER = ["PASS"]
        else:
            # simply output non-bnd and normal-variants
            pass
            #writer.write_record(sv)

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
    region_a = f"{info['REGIONA'][0]}_{info['REGIONA'][1]}"
    bnd_id = f"{sv_id_name}_{region_a}"
    return bnd_id, sv_id_num


if __name__ == "__main__":
    cli()
