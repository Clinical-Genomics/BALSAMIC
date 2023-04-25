#!/usr/bin/env python

import vcfpy
from vcfpy import Reader
import click
import logging
from typing import List, Dict, Tuple

LOG = logging.getLogger(__name__)

SV_FILTER_SETTINGS = {
    "tiddit_tumor_normal": {
        "max_tin_fraction": {
            "filter": "high_normal_af_fraction",
            "value": 0.25,
        },
        "max_normal_allele_frequency": {
            "filter": "high_normal_af",
            "value": 0.25,
        }
    }
}
@click.group()
def cli():
    """
    SV filtering tool.
    """
    pass


@cli.group("tiddit_tn")
@click.option(
    "-v",
    "--vcf-file",
    required=True,
    type=click.Path(exists=True),
    help="Input Variant Calling Format(VCF) from merged TIDDIT TUMOR + NORMAL VCF.",
)
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

@tiddit_tn.command("filter")
@click.pass_context
def filter(ctx: click.Context):
    """
    Add soft-filters based on presence of variants in normal.
    Outputs filtered VCF in standard-out.
    """
    vcf: Reader = vcfpy.Reader.from_path(ctx.obj["vcf_file"])
    filter_settings: Dict[str, Dict] = SV_FILTER_SETTINGS["tiddit_tumor_normal"]

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
            [
                ("ID", f"{filter_settings['max_normal_allele_frequency']['filter']}"),
                (
                    "Description",
                    f"AF_N_MAX > {filter_settings['max_normal_allele_frequency']['value']}",
                ),
            ]
        )
    )

    vcf.header.add_filter_line(
        vcfpy.OrderedDict(
            [
                ("ID", f"{filter_settings['max_tin_fraction']['filter']}"),
                (
                    "Description",
                    f"(AF_N_MAX / AF_T_MAX) > {filter_settings['max_tin_fraction']['value']}",
                ),
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
        info: dict = sv.INFO

        # Set default values
        tumor_has_contig, normal_has_contig = False, False
        allele_frequency_tumor, allele_frequency_normal = 0, 0

        # Extract allele frequency and contig-status for variant in T and N
        if "TUMOR_PASS_SAMPLE" in info:
            pass_sample: List[str] = info["TUMOR_PASS_SAMPLE"]
            pass_info: List[str] = info["TUMOR_PASS_INFO"]
            allele_frequency_tumor, max_af_index = get_max_allele_frequency(pass_sample)
            if allele_frequency_tumor == 0:
                get_any_contig = True
            else:
                get_any_contig = False
            tumor_has_contig: bool = look_for_contigs(
                variant_info_pass_field=pass_info,
                max_allele_frequency_variant_index=max_af_index,
                look_in_all_variants=get_any_contig,
            )

        if "NORMAL_PASS_SAMPLE" in info:
            pass_sample: List[str] = info["NORMAL_PASS_SAMPLE"]
            pass_info: List[str] = info["NORMAL_PASS_INFO"]
            allele_frequency_normal, max_af_index = get_max_allele_frequency(pass_sample)
            if allele_frequency_normal == 0:
                get_any_contig: bool = True
            else:
                get_any_contig: bool = False
            normal_has_contig: bool = look_for_contigs(
                variant_info_pass_field=pass_info,
                max_allele_frequency_variant_index=max_af_index,
                look_in_all_variants=get_any_contig,
            )

        # Set filter statuses
        if allele_frequency_tumor == 0 and not tumor_has_contig:
            sv.add_filter("normal_variant")
        else:
            # Regardless of CTG, set filter if AF_T / AF_N > max_tin_fraction
            normal_tumor_af_ratio = (
                float(allele_frequency_normal / allele_frequency_tumor)
                if allele_frequency_tumor > 0
                else 0
            )
            if normal_tumor_af_ratio > filter_settings["max_tin_fraction"]["value"]:
                sv.add_filter("high_normal_af_fraction")

            # Set filter if AF_N > 0.25
            if allele_frequency_normal > filter_settings["max_normal_allele_frequency"]["value"]:
                sv.add_filter("high_normal_af")

            # Set filter if CTG_N = True, AF_N is 0 and AF_T is below 0.25
            if (
                normal_has_contig
                and allele_frequency_normal == 0
                and allele_frequency_tumor <= 0.25
            ):
                sv.add_filter("in_normal")

        sv.INFO["AF_T_MAX"] = [round(allele_frequency_tumor, 4)]
        sv.INFO["AF_N_MAX"] = [round(allele_frequency_normal, 4)]

        writer.write_record(sv)


def has_contig(variant: str) -> bool:
    """
    Parses input string into a nested [key, value] list and looks for presence of a contig.

    Args:
        variant: input string for a variant from the list [SAMPLE]_PASS_INFO

    Returns:
        Bool(True or False) based on if contig was found or not.
    """
    fields: List[str] = variant.split("|")
    for field in fields:
        field_name_value: List[str] = field.split(":")
        field_name: str = field_name_value[0]
        if field_name == "CTG":
            field_value: str = field_name_value[1]
            if field_value != ".":
                return True
    return False


def look_for_contigs(
    variant_info_pass_field: List[str],
    max_allele_frequency_variant_index: int,
    bool_look_in_all_variants: bool,
) -> bool:
    """
    Directs parsing of sample_info to look for contig in either variant with max_af or any variant in sample_info.

    Args:
        variant_info_pass_field: A list of variants from [SAMPLE]_PASS_INFO
        max_allele_frequency_variant_index: Int(index) for variant in sample_info list with the highest allele frequency
        bool_look_in_all_variants: Boolean if it should look for a contig in any of the variants, not only max_af_index

    Returns:
        Bool(True or False) based on if contig was found or not.
    """
    if bool_look_in_all_variants:
        for variant_info in variant_info_pass_field:
            variant_info: List[str] = variant_info
            return has_contig(variant_info)

    max_af_variant_info: List[str] = variant_info_pass_field[max_allele_frequency_variant_index]
    return has_contig(max_af_variant_info)


def get_max_allele_frequency(pass_sample: List[str]) -> Tuple[float, int]:
    """
    Parses [SAMPLE]_PASS_SAMPLE field from TIDDIT info-field, returns the highest allele frequency for
    all merged variants in the sample and the position of this variant in the list.

    Args:
         pass_sample: list of information for variant/s within TUMOR_PASS_SAMPLE / NORMAL_PASS_SAMPLE

    Returns:
         float(maximum allele frequency), and int(position in list of max af variant)
    """
    max_allele_frequency = 0
    max_af_variant_index = 0
    total_cov = 0
    DV = 0
    RV = 0

    for variant_idx, variant in enumerate(pass_sample):
        fields: List[str] = variant.split("|")
        for field in fields:
            field_name_value: List[str] = field.split(":")
            field_name: str = field_name_value[0]
            if field_name == "COV":
                b1_cov = int(field_name_value[1])
                b2_cov = int(field_name_value[3])
                total_cov: int = b1_cov + b2_cov

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
            max_af_variant_index: int = variant_idx
            max_allele_frequency: float = allele_frequency
    return max_allele_frequency, max_af_variant_index


@tiddit_tn.command("rescue_bnds")
@click.pass_context
def rescue_bnds(ctx: click.Context):
    """
    Rescue BND-variants with at least 1 of 2 BNDs set to PASS.
    Outputs VCF in standard-out.
    """
    # First read of VCF-file:
    vcf_start: Reader = vcfpy.Reader.from_path(ctx.obj["vcf_file"])

    # define dict with bnd_id - bnd_num - FILTER
    bnd_filter_dict: dict = {}
    for sv in vcf_start:
        info: dict = sv.INFO

        svtype = info["SVTYPE"]
        # only relevant for bnd-variants and tumor-variants
        if svtype == "BND" and "TUMOR_PASS_CHROM" in info:
            bnd_id, sv_id_num = get_bnd_id(info)

            if bnd_id not in bnd_filter_dict:
                bnd_filter_dict[bnd_id] = {}

            if sv_id_num in bnd_filter_dict[bnd_id]:
                LOG.warning(
                    f"Conflicting BND-names: {bnd_id}, will not attempt to rescue: {sv}"
                )
                continue

            bnd_filter_dict[bnd_id][sv_id_num] = sv.FILTER

    # Second read of VCF-file:
    vcf: Reader = vcfpy.Reader.from_path(ctx.obj["vcf_file"])
    # Update VCF header
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
        info: dict = sv.INFO
        svtype = info["SVTYPE"]
        # only relevant for bnd and tumor-variants
        if svtype == "BND" and "TUMOR_PASS_CHROM" in info:
            bnd_id, sv_id_num = get_bnd_id(variant_info_field=info)
            sv.INFO["BND_ID"] = [bnd_id]

            filter_bnds = bnd_filter_dict[bnd_id]
            # change filter if any of the 2 BNDs are PASS
            if (
                "1" in filter_bnds
                and "2" in filter_bnds
                and ("PASS" in filter_bnds["1"] or "PASS" in filter_bnds["2"])
            ):
                sv.FILTER = ["PASS"]
            writer.write_record(sv)
        else:
            # simply output non-bnd and normal-variants
            writer.write_record(sv)


def get_bnd_id(variant_info_field: dict) -> Tuple[str, str]:
    """
    Creates a unique variant ID for BND variants based on information in the info field.
    This info field may contain multiple BND-variants merged within a sample.
    In this case, the info for the first variant is chosen.

    Args:
         variant_info_field: dict with information for BND variant/s in the VCF.

    Returns:
         str(unique ID for the BND variant) str(the 1st or 2nd BND)
    """
    # example: SV_2414_1, SV_2414 = sv_id_name, 1 = sv_id_num
    sv_id = variant_info_field["TUMOR_PASS_CHROM"][0].split("|")[0]

    # extract sv_id and breakend number
    sv_id_split = sv_id.split("_")
    sv_id_name = "_".join(sv_id_split[0:2])
    sv_id_num = sv_id_split[2]

    # define unique breakend_name
    region_a = f"{variant_info_field['REGIONA'][0]}_{variant_info_field['REGIONA'][1]}"
    bnd_id = f"{sv_id_name}_{region_a}"
    return bnd_id, sv_id_num


if __name__ == "__main__":
    cli()
