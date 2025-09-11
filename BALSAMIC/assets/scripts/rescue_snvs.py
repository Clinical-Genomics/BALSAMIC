#!/usr/bin/env python3

from pathlib import Path
import vcfpy
import click


class RescueReasons:
    CLINVAR_ONC = "ClinvarOnc"
    CLINVAR_PATH = "ClinvarPathogenic"
    CLINVAR_LIKELY_PATH = "ClinvarLikelyPathogenic"
    RESCUE_LIST = "RescueList"


def load_rescue_variants(vcf_path: str) -> set[tuple]:
    """Load variants to rescue and return an set of tuples (CHROM, POS, REF, ALT)."""
    reader = vcfpy.Reader.from_path(vcf_path)
    return {
        (record.CHROM, record.POS, record.REF, record.ALT[0].serialize())
        for record in reader
    }


def update_headers(reader) -> None:
    """Update vcf headers."""
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

    for id, number, type, desc in INFO_HEADERS:
        reader.header.add_info_line(
            vcfpy.OrderedDict(
                [("ID", id), ("Number", number), ("Type", type), ("Description", desc)]
            )
        )


def determine_clinvar_reasons(record: vcfpy.record) -> list[str]:
    """Determine rescue reasons based on ClinVar annotations."""
    reasons: list[str] = []
    if "ONC" in record.INFO:
        onc = "|".join(record.INFO["ONC"]).lower()
        # Exclude no_classification_for_the_single_variant, Benign and likely benign
        if "oncogenic" in onc or "tier" in onc or "uncertain" in onc:
            reasons.append(RescueReasons.CLINVAR_ONC)
    if "Pathogenic" in record.INFO.get("CLNSIG", ""):
        reasons.append(RescueReasons.CLINVAR_PATH)
    if "Likely_pathogenic" in record.INFO.get("CLNSIG", ""):
        reasons.append(RescueReasons.CLINVAR_LIKELY_PATH)
    return reasons


def record_in_rescue_list(record: vcfpy.Record, rescue_keys: set[tuple] | None) -> bool:
    """Returns True if the record is within the variants to rescue."""
    if not rescue_keys:
        return False
    return any(
        (record.CHROM, record.POS, record.REF, alt.serialize()) in rescue_keys
        for alt in record.ALT
    )  # To account for multiple ALT per row (although this shouldn't be the case)


def process_record(
    record: vcfpy.Record, rescue_keys: set[tuple] | None
) -> vcfpy.Record:
    """Update record with rescue information and adjusted fiter."""
    reasons: list[str] = []
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
    in_vcf_path: str, out_vcf_path: str, rescue_keys: set[tuple] | None
) -> None:
    """Writes a vcf file with rescued variants and additional info fields."""
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
def cli(rescue_list: Path | None, vcf_path: Path, out_path: Path) -> None:
    """
    Rescue variants by rescue criteria. If a variant is rescued the filters are set as PASS and
    two new info fields are added RescueFilters (original filters in input vcf) and RescueStatus
    (list of reasons to rescue).

    Rescue criteria:
      - If --rescue-list is provided: rescue variant match (CHROM, POS, REF, ALT).
      - Clinvar:
        - CLNSIG contains 'Pathogenic' or 'Likely_pathogenic'.
        - ClinVar: ONC field is present and partially matches "oncogenic", "tier", "uncertain"
    """
    rescue_keys = load_rescue_variants(str(rescue_list)) if rescue_list else None
    process_vcf(str(vcf_path), str(out_path), rescue_keys)


if __name__ == "__main__":
    cli()
