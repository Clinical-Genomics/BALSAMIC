#!/usr/bin/env python3

from pathlib import Path
from typing import List, Optional, Set, Tuple

import vcfpy
import click


class RescueReasons:
    CLINVAR_ONC = "ClinvarOnc"
    CLINVAR_PATH = "ClinvarPathogenic"
    CLINVAR_LIKELY_PATH = "ClinvarLikelyPathogenic"
    RESCUE_LIST = "RescueList"


RescueKey = Tuple[str, int, str, str]  # (CHROM, POS, REF, ALT)


def load_rescue_variants(vcf_path: str) -> Set[RescueKey]:
    """Load variants to rescue and return a set of tuples (CHROM, POS, REF, ALT)."""
    reader = vcfpy.Reader.from_path(vcf_path)
    return {
        (record.CHROM, record.POS, record.REF, record.ALT[0].serialize())
        for record in reader
    }


def update_headers(reader) -> None:
    """Update VCF headers."""
    INFO_HEADERS = [
        (
            "RescueFilters",
            "1",
            "String",
            "Original FILTER value moved here because this variant was rescued",
        ),
        (
            "RescueStatus",
            "1",
            "String",
            "Reason(s) for rescuing; pipe-separated (e.g., ClinvarOnc|ClinvarPathogenic)",
        ),
    ]

    for id_, number, type_, desc in INFO_HEADERS:
        reader.header.add_info_line(
            vcfpy.OrderedDict(
                [
                    ("ID", id_),
                    ("Number", number),
                    ("Type", type_),
                    ("Description", desc),
                ]
            )
        )


def get_clinvar_pathogenicity(record: vcfpy.Record) -> str:
    clnsig = record.INFO.get("CLNSIG") or ""
    if isinstance(clnsig, list):
        clnsig = "|".join(clnsig)
    return clnsig.lower()


def determine_clinvar_reasons(record: vcfpy.Record) -> List[str]:
    """Determine rescue reasons based on ClinVar annotations."""
    reasons: List[str] = []
    if "ONC" in record.INFO:
        onc_val = record.INFO["ONC"]
        onc = (
            "|".join(onc_val).lower()
            if isinstance(onc_val, list)
            else str(onc_val).lower()
        )
        # Exclude no_classification_for_the_single_variant, Benign and likely benign
        if "oncogenic" in onc or "tier" in onc or "uncertain" in onc:
            reasons.append(RescueReasons.CLINVAR_ONC)

    clnsig = get_clinvar_pathogenicity(record)
    if "likely_pathogenic" in clnsig:
        reasons.append(RescueReasons.CLINVAR_LIKELY_PATH)
    elif "pathogenic" in clnsig:
        reasons.append(RescueReasons.CLINVAR_PATH)
    return reasons


def record_in_rescue_list(
    record: vcfpy.Record, rescue_keys: Optional[Set[RescueKey]]
) -> bool:
    """Returns True if the record is within the variants to rescue."""
    if not rescue_keys:
        return False
    # Account for multi-ALT (shouldn't happen here but be safe)
    return any(
        (record.CHROM, record.POS, record.REF, alt.serialize()) in rescue_keys
        for alt in record.ALT
    )


def process_record(
    record: vcfpy.Record, rescue_keys: Optional[Set[RescueKey]]
) -> vcfpy.Record:
    """Update record with rescue information and adjusted filter."""
    reasons: List[str] = []
    if record_in_rescue_list(record, rescue_keys):
        reasons.append(RescueReasons.RESCUE_LIST)
    reasons.extend(determine_clinvar_reasons(record))

    if not reasons:
        return record

    if record.FILTER not in ([], ["."], ["PASS"]):
        record.INFO["RescueFilters"] = "|".join(record.FILTER)
        record.FILTER = ["PASS"]
        record.INFO["RescueStatus"] = "|".join(reasons)
    return record


def process_vcf(
    in_vcf_path: str, out_vcf_path: str, rescue_keys: Optional[Set[RescueKey]]
) -> None:
    """Writes a VCF file with rescued variants and additional INFO fields."""
    reader = vcfpy.Reader.from_path(in_vcf_path)
    update_headers(reader)
    with vcfpy.Writer.from_path(out_vcf_path, reader.header) as writer:
        for record in reader:
            writer.write_record(process_record(record, rescue_keys))


@click.command()
@click.option(
    "--rescue-list",
    "rescue_list",
    type=click.Path(exists=True),
    required=False,
    help="Optional: VCF (.vcf or .vcf.gz) of variants to rescue (matched on CHROM, POS, REF, ALT).",
)
@click.option(
    "--vcf",
    "vcf_path",
    type=click.Path(exists=True),
    required=True,
    help="Input VCF (.vcf or .vcf.gz) to process.",
)
@click.option(
    "-o",
    "--out",
    "out_path",
    type=click.Path(dir_okay=False),
    required=True,
    help="Output VCF path (file will be overwritten).",
)
def cli(rescue_list: Optional[Path], vcf_path: Path, out_path: Path) -> None:
    """
    Rescue variants by rescue criteria. If a variant is rescued the filters are set to PASS and
    two new INFO fields are added: RescueFilters (original filters) and RescueStatus (reasons).
    """
    rescue_keys: Optional[Set[RescueKey]] = (
        load_rescue_variants(str(rescue_list)) if rescue_list else None
    )
    process_vcf(str(vcf_path), str(out_path), rescue_keys)


if __name__ == "__main__":
    cli()
