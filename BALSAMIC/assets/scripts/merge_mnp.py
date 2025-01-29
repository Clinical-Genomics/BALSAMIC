#!/usr/bin/env python

"""

Copyright (c) Sentieon Inc. All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice, this
    list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Script to combine neighboring PASS variants on the same haplotype,
supporting user-defined merging strategy with the corresponding
arguments.
Original variants that are combined will be marked as "MERGED".
Usage:
  python merge_mnp.py $VCF $FASTA
Usage with codon file (requires merge_by_codon.py strategy file):
  python merge_mnp.py $VCF $FASTA merge_by_codon $CODON_FILE $IGNORE_INDELS

Requirements:
This script requires access to the vcflib library contained in the Sentieon
software package located in $SENTIEON_INSTALL_DIR/lib/python/sentieon.
"""

from __future__ import print_function

import vcflib
import copy
import click
import sys
import pyfaidx
from typing import List, Union, Optional

global vcf, reference


def ifmerge(
    variant1: vcflib.vcf.Variant, variant2: vcflib.vcf.Variant, max_distance: int
) -> bool:
    """
    Determine whether two variants should be merged based on specific conditions.

    Parameters:
        variant1 (vcflib.vcf.Variant): The first variant object.
        variant2 (vcflib.vcf.Variant): The second variant object.
        max_distance (int): The maximum allowable distance between two variants
                            for them to be considered for merging.

    Returns:
        bool: True if the variants should be merged, False otherwise.
    """
    # Check if either variant contains structural variation information
    if "SVTYPE" in variant1.info or "SVTYPE" in variant2.info:
        return False

    # Ensure variants are on the same chromosome and within the maximum distance
    if (
        variant1.chrom != variant2.chrom
        or variant2.pos - variant1.pos > max_distance + len(variant1.ref) - 1
    ):
        return False

    # Extract phasing-related information from samples
    pid1 = variant1.samples[0].get("PID", "")
    pid2 = variant2.samples[0].get("PID", "")
    pgt1 = variant1.samples[0].get("PGT", "")
    pgt2 = variant2.samples[0].get("PGT", "")
    ps1 = variant1.samples[0].get("PS", "")
    ps2 = variant2.samples[0].get("PS", "")

    # Check if variants share the same PID and PGT, or the same PS
    if (pid1 == pid2 and pgt1 == pgt2 and pid1 and pgt1) or (ps1 == ps2 and ps1):
        return True

    return False


def distance(
    v1: Union[List[vcflib.vcf.Variant], vcflib.vcf.Variant], v2: vcflib.vcf.Variant
) -> int:
    """
    Calculate the distance between two VCF variants. If `v1` is a list of variants,
    the function will calculate the distance to the last variant in the list.

    Parameters:
        v1 (Union[List[vcflib.vcf.Variant], vcflib.VCF.Variant]): The first variant(s).
        v2 (vcflib.vcf.Variant): The second variant.

    Returns:
        int: The calculated distance between the two variants. Returns 1000 if variants
             are on different chromosomes or if the distance cannot be determined.
    """
    if isinstance(v1, list):
        # If v1 is a list of variants, return the distance to the last variant in the list
        for v in v1[::-1]:  # Reverse list iteration
            return v2.pos - v.pos - len(v.ref) + 1
        return 1000

    if v1.chrom != v2.chrom:
        # If variants are on different chromosomes, return a large distance
        return 1000

    # Calculate distance between two variants on the same chromosome
    return v2.pos - v1.pos - len(v1.ref) + 1


def merge(
    variant_stack: List[vcflib.vcf.Variant], max_distance: int
) -> List[vcflib.vcf.Variant]:
    """
    Merges a stack of overlapping variants into a single variant.

    This function takes a stack of variants that are close enough (within a specified max_distance)
    and merges them into one variant, combining their reference and alternate alleles, quality scores,
    and sample information.

    Parameters:
        variant_stack (List[vcflib.vcf.Variant]): A list of variants to be merged.
        max_distance (int): The maximum distance between variants for them to be considered for merging.

    Returns:
        List[vcflib.vcf.Variant]: A list of merged variants.
    """

    def _merge(vv: List[vcflib.vcf.Variant]) -> vcflib.vcf.Variant:
        """
        Merges a list of variants into a single variant by combining their properties.

        This function merges reference, alternate alleles, quality, INFO fields, and sample information
        of overlapping variants in the list.

        Parameters:
            vv (List[vcflib.vcf.Variant]): The list of variants to be merged.

        Returns:
            vcflib.vcf.Variant: The merged variant.
        """
        len_vv = len(vv)

        # Check if variants overlap in positions, return None if not
        for i in range(len_vv - 1):
            if vv[i + 1].pos - vv[i].pos < len(vv[i].ref):
                return None

        # Initialize the merged variant with the first variant in the list
        v = copy.deepcopy(vv[0])
        v.id = "."  # Clear the ID for merged variant
        ref = vv[0].ref
        alt = vv[0].alt

        # Combine REF and ALT sequences for all variants
        for i in range(1, len_vv):
            vi = vv[i]
            ref_gap = ""
            if vi.pos - (v.pos + len(ref)) > 0:
                ref_gap = reference[v.chrom][v.pos + len(ref) : vi.pos].seq.upper()
            ref = ref + ref_gap + vi.ref
            alt = [alt[j] + ref_gap + vi.alt[j] for j in range(len(alt))]

        # Calculate the merged quality score
        qual_list = [vi.qual for vi in vv]
        v.qual = sum(qual_list) / len_vv if None not in qual_list else None

        # Merge INFO fields by averaging the values
        for k in v.info.keys():
            vk = [vi.info.get(k) for vi in vv]
            if None in vk:
                v.info[k] = None
            else:
                try:
                    v.info[k] = sum(vk) / len_vv
                except (TypeError, ZeroDivisionError):
                    v.info[k] = None

        # Remove common leading bases between REF and ALT alleles
        i = 0
        while i < len(ref) - 1:
            common = False not in [
                ref[i] == alti[i] if i < len(alti) - 1 else False for alti in alt
            ]
            if not common:
                break
            i += 1
        v.ref = ref[i:]
        v.alt = [alti[i:] for alti in alt]
        v.pos += i

        # Set MNV filter
        all_filters = {flt for vi in vv for flt in vi.filter}
        if len(all_filters) > 1:
            v.filter = ["MERGED_MNV_CONFLICTING_FILTERS"]
        else:
            v.filter = ["MERGED_MNV", f"{all_filters.pop()}"]

        # Mark all constituent variants as "MERGED"
        for vi in vv:
            if "MERGED" not in vi.filter:
                vi.filter.append("MERGED")

        # Merge sample information (AF, AD, AFDP)
        for i in range(len(vcf.samples)):
            t = v.samples[i]
            vv_samples = [vs.samples[i] for vs in vv]
            af_values = [vi.get("AF") for vi in vv_samples]

            # Handle AF (allele frequency)
            if None not in af_values:
                if isinstance(vv_samples[0]["AF"], list):
                    t["AF"] = [
                        sum([vsi["AF"][j] for vsi in vv_samples]) / len_vv
                        for j in range(len(alt))
                    ]
                else:
                    t["AF"] = sum(af_values) / len_vv

            # Handle AD (allele depths)
            t["AD"] = (
                int(sum([vi["AD"][0] for vi in vv_samples]) / len_vv),
                int(sum([vi["AD"][1] for vi in vv_samples]) / len_vv),
            )

            # Handle AFDP (allele frequency depth)
            afdp_values = [vi.get("AFDP") for vi in vv_samples]
            if None not in afdp_values:
                t["AFDP"] = int(sum(afdp_values) / len_vv)

        # Format the final variant
        _ = vcf.format(v)
        for vi in vv:
            _ = vcf.format(vi)

        return v

    vlist = []  # List to store the final merged variants
    to_push = []  # Temporary list to store variants being merged

    # Iterate over the variant stack and merge overlapping variants
    for i in range(len(variant_stack)):
        vi: vcflib.vcf.Variant = variant_stack[i]
        to_merge = [vi]
        for j in range(i + 1, len(variant_stack)):
            vj: vcflib.vcf.Variant = variant_stack[j]
            if ifmerge(to_merge[-1], vj, max_distance):
                to_merge.append(vj)

        # If there are variants to merge, merge them
        if len(to_merge) > 1:
            merged = _merge(to_merge)
            if merged:
                to_push.append(merged)

        # Add merged variants to the final list if they are in order
        while len(to_push):
            if to_push[0].pos <= vi.pos:
                vlist.append(to_push.pop(0))
            else:
                break
        vlist.append(vi)

    # Append any remaining merged variants
    while len(to_push):
        vlist.append(to_push.pop(0))

    return vlist


def process(
    vcf_file: str, ref_file: str, out_file: Optional[str], max_distance: int
) -> None:
    """
    Processes a VCF file, merges neighboring variants into MNVs, and writes the result to an output file or stdout.

    Parameters:
        vcf_file (str): Path to the input VCF file to process.
        ref_file (str): Path to the reference genome file.
        out_file (Optional[str]): Path to the output VCF file. If not specified, the result is printed to stdout.
        max_distance (int): Maximum distance between two variants to be merged.

    This function reads the input VCF file, processes each variant, and merges variants that are close
    enough to each other based on the `max_distance` parameter. The result is written to the specified
    output file or to stdout if no output file is provided.
    """
    global vcf, reference
    vcf = vcflib.VCF(vcf_file)
    reference = pyfaidx.Fasta(ref_file)

    # Open output file (or use stdout if no file is specified)
    out_fh = open(out_file, "w") if out_file else sys.stdout

    # Define and add the MERGED filter to the VCF header if not already present
    new_filters = {
        "MERGED": {
            "Description": '"SNV Merged with neighboring variants"',
            "ID": "MERGED",
        },
        "MERGED_MNV": {
            "Description": '"Created from merged SNVs with same phase-id"',
            "ID": "MERGED_MNV",
        },
        "MERGED_MNV_CONFLICTING_FILTERS": {
            "Description": '"Merged MNV contains SNVs with conflicting filters, such as triallelic_site and in_normal"',
            "ID": "MNV_CONFLICTING_FILTERS",
        }
    }
    for new_filter_id, description_dict in new_filters.items():
        vcf.filters[new_filter_id] = description_dict

    filter_added = False
    for header_line in vcf.headers:
        if header_line.startswith("##FILTER") and not filter_added:
            for new_filter_id, description_dict in new_filters.items():
                id = description_dict["ID"]
                description = description_dict["Description"]
                print(
                    f"##FILTER=<ID={id},Description={description}>",
                    file=out_fh,
                )
            filter_added = True
        print(header_line, file=out_fh)

    # Process variants and merge them if they are within max_distance
    variant_stack = []
    for variant in vcf:
        if not variant_stack:
            # If the stack is empty, add the variant or print if it's an SV
            if "SVTYPE" in variant.info:
                print(variant, file=out_fh)
            else:
                variant_stack.append(variant)
        else:
            # Check if the current variant is close enough to merge
            if distance(variant_stack, variant) <= max_distance:
                variant_stack.append(variant)
            else:
                # Merge variants in the stack and reset it
                merged_variants = merge(variant_stack, max_distance)
                for merged_variant in merged_variants:
                    print(merged_variant, file=out_fh)
                variant_stack = []
                # Print the current variant or add it to the stack
                if "SVTYPE" in variant.info:
                    print(variant, file=out_fh)
                else:
                    variant_stack.append(variant)

    # Process any remaining variants in the stack after the loop
    merged_variants = merge(variant_stack, max_distance)
    for merged_variant in merged_variants:
        print(merged_variant, file=out_fh)

    # Close the output file if it was specified
    if out_file:
        out_fh.close()


@click.command()
@click.argument("vcf_file", type=click.Path(exists=True))
@click.argument("reference", type=click.Path(exists=True))
@click.option(
    "--out_file",
    type=click.Path(),
    default=None,
    help="Output VCF file. If not specified, output will be written to stdout.",
)
@click.option(
    "--max_distance",
    type=int,
    default=5,
    show_default=True,
    help="Maximum distance between two variants to be merged.",
)
def main(
    vcf_file: str, reference: str, out_file: Optional[str], max_distance: int
) -> None:
    """
    Merge phased SNVs into MNVs and output the result to a VCF file or stdout.

    Parameters:
        vcf_file (str): Path to the input VCF file to process.
        reference (str): Path to the reference genome file.
        out_file (Optional[str]): Path to the output VCF file. If not specified, the result will be printed to stdout.
        max_distance (int): Maximum allowed distance between two variants to be merged.

    This command processes the input VCF file and merges phased single nucleotide variants (SNVs)
    into multi-nucleotide variants (MNVs) if they meet the criteria based on `max_distance`.
    """
    # Call the process function to perform the actual merging
    process(vcf_file, reference, out_file, max_distance)


if __name__ == "__main__":
    main()
