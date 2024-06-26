#!/usr/bin/env python

import click
import vcfpy
import logging
from vcfpy import Reader
from pathlib import Path
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
        },
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
                ("Number", "."),
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
                ("Number", "."),
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
    for variant in vcf:
        variant_info: dict = variant.INFO

        # Collect evidence of variant in tumor and normal sample
        evidence_dict: dict = get_tumor_normal_evidence(variant_info)
        allele_frequency_tumor: float = evidence_dict["tumor_max_af"]
        allele_frequency_normal: float = evidence_dict["normal_max_af"]
        tumor_has_contig: bool = evidence_dict["tumor_has_contig"]
        normal_has_contig: bool = evidence_dict["normal_has_contig"]

        # Add AF_MAX to info field
        variant.INFO["AF_T_MAX"] = [round(allele_frequency_tumor, 4)]
        variant.INFO["AF_N_MAX"] = [round(allele_frequency_normal, 4)]

        # Set filter statuses
        if allele_frequency_tumor == 0 and not tumor_has_contig:
            variant.add_filter("normal_variant")
            writer.write_record(variant)
            continue

        # Regardless of CTG, set filter if AF_T / AF_N > max_tin_fraction
        normal_tumor_af_ratio = (
            float(allele_frequency_normal / allele_frequency_tumor)
            if allele_frequency_tumor > 0
            else 0
        )
        if normal_tumor_af_ratio > filter_settings["max_tin_fraction"]["value"]:
            variant.add_filter("high_normal_af_fraction")

        # Set filter if AF_N > 0.25
        if (
            allele_frequency_normal
            > filter_settings["max_normal_allele_frequency"]["value"]
        ):
            variant.add_filter("high_normal_af")

        # Set filter if CTG_N = True, AF_N is 0 and AF_T is below 0.25
        if (
            normal_has_contig
            and allele_frequency_normal == 0
            and allele_frequency_tumor <= 0.25
        ):
            variant.add_filter("in_normal")

        writer.write_record(variant)


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
    look_in_all_variants: bool,
) -> bool:
    """
    Directs parsing of sample_info to look for contig in either variant with max_af or any variant in sample_info.

    Args:
        variant_info_pass_field: A list of variants from [SAMPLE]_PASS_INFO
        max_allele_frequency_variant_index: Int(index) for variant in sample_info list with the highest allele frequency
        look_in_all_variants: Boolean if it should look for a contig in any of the variants, not only max_af_index

    Returns:
        Bool(True or False) based on if contig was found or not.
    """
    if look_in_all_variants:
        for variant_info in variant_info_pass_field:
            return has_contig(variant_info)

    max_af_variant_info: List[str] = variant_info_pass_field[
        max_allele_frequency_variant_index
    ]
    return has_contig(max_af_variant_info)


def calc_allele_frequency(variant_info_sample_field: List[str]) -> float:
    """
    Parses [SAMPLE]_PASS_SAMPLE field from TIDDIT info-field and calculates the allele frequency for the variant.

    Args:
         variant_info_sample_field: list of information for variant from tumor or normal sample

    Returns:
         float(allele frequency)
    """
    total_cov: int = 0
    DV: int = 0
    RV: int = 0
    for field in variant_info_sample_field:
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
    return allele_frequency


def get_max_allele_frequency(variant_info_sample_field: List[str]) -> Tuple[float, int]:
    """
    Finds the highest allele frequency for all merged variants in the sample
    and the position of this variant in the list.

    Args:
         variant_info_sample_field: list of information for variant/s from tumor or normal sample

    Returns:
         float(maximum allele frequency), and int(position in list of max af variant)
    """
    max_allele_frequency: float = 0
    max_af_variant_index: int = 0

    for variant_idx, variant in enumerate(variant_info_sample_field):
        variant_info_fields: List[str] = variant.split("|")
        allele_frequency: float = calc_allele_frequency(
            variant_info_sample_field=variant_info_fields
        )

        if allele_frequency > max_allele_frequency:
            max_af_variant_index: int = variant_idx
            max_allele_frequency: float = allele_frequency
    return max_allele_frequency, max_af_variant_index


def get_tumor_normal_evidence(variant_info_field: dict) -> dict:
    """
    Collects the read evidence for the variant in the tumor and normal sample.
    Gets maximum allele frequency and checks if a contig has been assembled for the variant.

    Args:
         variant_info_field: dict with information for BND variant/s in the VCF.

    Returns:
         Dictionary containing all extracted evidences in the tumor and normal sample.
    """
    # Set default values
    tumor_has_contig, normal_has_contig = False, False
    allele_frequency_tumor, allele_frequency_normal = 0, 0

    # Extract allele frequency and contig-status for variant in T and N
    if "TUMOR_PASS_SAMPLE" in variant_info_field:
        pass_sample: List[str] = variant_info_field["TUMOR_PASS_SAMPLE"]
        pass_info: List[str] = variant_info_field["TUMOR_PASS_INFO"]
        allele_frequency_tumor, max_af_index = get_max_allele_frequency(
            variant_info_sample_field=pass_sample
        )
        if allele_frequency_tumor == 0:
            get_any_contig = True
        else:
            get_any_contig = False
        tumor_has_contig: bool = look_for_contigs(
            variant_info_pass_field=pass_info,
            max_allele_frequency_variant_index=max_af_index,
            look_in_all_variants=get_any_contig,
        )

    if "NORMAL_PASS_SAMPLE" in variant_info_field:
        pass_sample: List[str] = variant_info_field["NORMAL_PASS_SAMPLE"]
        pass_info: List[str] = variant_info_field["NORMAL_PASS_INFO"]
        allele_frequency_normal, max_af_index = get_max_allele_frequency(
            variant_info_sample_field=pass_sample
        )
        if allele_frequency_normal == 0:
            get_any_contig: bool = True
        else:
            get_any_contig: bool = False
        normal_has_contig: bool = look_for_contigs(
            variant_info_pass_field=pass_info,
            max_allele_frequency_variant_index=max_af_index,
            look_in_all_variants=get_any_contig,
        )

    evidence_dict: dict = {
        "tumor_max_af": allele_frequency_tumor,
        "normal_max_af": allele_frequency_normal,
        "tumor_has_contig": tumor_has_contig,
        "normal_has_contig": normal_has_contig,
    }
    return evidence_dict


@tiddit_tn.command("rescue_bnds")
@click.pass_context
def rescue_bnds(ctx: click.Context):
    """
    Rescue BND-variants with at least 1 of 2 BNDs set to PASS.
    Outputs VCF in standard-out.
    """
    vcf_file: Path = Path(ctx.obj["vcf_file"])

    # define dict with bnd_id - bnd_num - FILTER
    bnd_filter_dict: dict = get_bnd_filter_dict(vcf_file)

    vcf: Reader = vcfpy.Reader.from_path(vcf_file)

    # Update VCF header
    vcf.header.add_info_line(
        vcfpy.OrderedDict(
            [
                ("ID", "BND_ID"),
                ("Number", "."),
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
    for variant in vcf:
        variant_info: dict = variant.INFO
        svtype = variant_info["SVTYPE"]
        # only relevant for bnd and tumor-variants
        if svtype == "BND" and "TUMOR_PASS_CHROM" in variant_info:
            bnd_id, sv_id_num = get_bnd_id(variant_info_field=variant_info)
            variant.INFO["BND_ID"] = [bnd_id]

            filter_bnds: Dict[str, str] = bnd_filter_dict[bnd_id]
            # change filter if any of the 2 BNDs are PASS
            if (
                "1" in filter_bnds
                and "2" in filter_bnds
                and ("PASS" in filter_bnds["1"] or "PASS" in filter_bnds["2"])
            ):
                variant.FILTER = ["PASS"]
            writer.write_record(variant)
        else:
            # simply output non-bnd and normal-variants
            writer.write_record(variant)


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
    sv_id: str = variant_info_field["TUMOR_PASS_CHROM"][0].split("|")[0]

    # extract sv_id and breakend number
    sv_id_split: List[str] = sv_id.split("_")
    sv_id_name: str = "_".join(sv_id_split[0:2])
    sv_id_num: str = sv_id_split[2]

    # define unique breakend_name
    region_a: str = (
        f"{variant_info_field['REGIONA'][0]}_{variant_info_field['REGIONA'][1]}"
    )
    bnd_id: str = f"{sv_id_name}_{region_a}"
    return bnd_id, sv_id_num


def get_bnd_filter_dict(vcf_file: Path) -> Dict[str, Dict]:
    """
    Extracts a dictionary of filter statuses for BND (breakend) variants from a VCF file.

    Arguments:
        vcf_file (Path): a pathlib.Path object pointing to the input VCF file containing BND variants.

    Returns:
        bnd_filter_dict (Dict[str, Dict]): a dictionary where each key represents a unique BND variant identifier,
      and each value is another dictionary mapping individual SV IDs to their corresponding filter statuses.
    """
    vcf: Reader = vcfpy.Reader.from_path(vcf_file)

    # define dict with bnd_id - bnd_num - FILTER
    bnd_filter_dict: dict = {}
    for variant in vcf:
        variant_info: dict = variant.INFO

        svtype = variant_info["SVTYPE"]
        # only relevant for bnd-variants and tumor-variants
        if svtype == "BND" and "TUMOR_PASS_CHROM" in variant_info:
            bnd_id, sv_id_num = get_bnd_id(variant_info)

            if bnd_id not in bnd_filter_dict:
                bnd_filter_dict[bnd_id] = {}

            if sv_id_num in bnd_filter_dict[bnd_id]:
                LOG.warning(
                    f"Conflicting BND-names: {bnd_id}, will not attempt to rescue: {variant}"
                )
                continue

            bnd_filter_dict[bnd_id][sv_id_num] = variant.FILTER
    return bnd_filter_dict


if __name__ == "__main__":
    cli()
