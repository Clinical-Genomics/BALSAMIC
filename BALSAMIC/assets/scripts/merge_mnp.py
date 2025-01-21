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

import argparse
import copy
import sys
from pyfaidx import Fasta
from vcflib import VCF

# Global variables
global vcf, reference
MAX_DISTANCE = 20  # Default value

def should_merge(variant1, variant2):
    """Determine whether two variants should be merged."""

    # Do not merge structural variants (SVs)
    if "SVTYPE" in variant1.info or "SVTYPE" in variant2.info:
        return False

    # Variants must be close enough and on the same chromosome
    if variant1.chrom != variant2.chrom or variant2.pos - variant1.pos > MAX_DISTANCE + len(variant1.ref) - 1:
        return False

    # Must share the same PID/PGT or PS (phasing information)
    pid1, pid2 = variant1.samples[0].get("PID", ""), variant2.samples[0].get("PID", "")
    pgt1, pgt2 = variant1.samples[0].get("PGT", ""), variant2.samples[0].get("PGT", "")
    ps1, ps2 = variant1.samples[0].get("PS", ""), variant2.samples[0].get("PS", "")

    if ((pid1 == pid2 and pgt1 == pgt2 and pid1 and pgt1) or (ps1 == ps2 and ps1)):
        return True

    return False

def calculate_distance(variant1, variant2):
    """Calculate the distance between two variants."""
    return variant2.pos - variant1.pos - len(variant1.ref) + 1

def merge_variants(variants):
    """Merge a list of variants into a single variant."""

    def merge_individuals(v_list):
        if len(v_list) < 2:
            return None

        merged_variant = copy.deepcopy(v_list[0])
        merged_variant.id = "."
        ref = v_list[0].ref
        alt = v_list[0].alt

        # Combine the filters from all variants
        combined_filters = set()
        for variant in v_list:
            combined_filters.update(variant.filter)
        if "PASS" in combined_filters and len(combined_filters) > 1:
            combined_filters.remove("PASS")
        merged_variant.filter = list(combined_filters)

        for next_variant in v_list[1:]:
            ref_gap = ""
            gap_length = next_variant.pos - (merged_variant.pos + len(ref))
            if gap_length > 0:
                ref_gap = reference[merged_variant.chrom][merged_variant.pos + len(ref):next_variant.pos].seq.upper()

            ref += ref_gap + next_variant.ref
            alt = [alt[i] + ref_gap + next_variant.alt[i] for i in range(len(alt))]
            next_variant.filter = ["MERGED"]

        merged_variant.qual = (sum(v.qual for v in v_list if v.qual) / len(v_list)) if all(
            v.qual for v in v_list) else None

        merged_variant.ref = ref
        merged_variant.alt = alt

        return merged_variant

    merged_list = []
    to_merge = []

    for variant in variants:
        if not to_merge or should_merge(to_merge[-1], variant):
            to_merge.append(variant)
        else:
            if to_merge:
                merged_variant = merge_individuals(to_merge)
                if merged_variant:
                    merged_list.append(merged_variant)
            to_merge = [variant]

    if to_merge:
        merged_variant = merge_individuals(to_merge)
        if merged_variant:
            merged_list.append(merged_variant)

    return merged_list

def process_vcf(vcf_file, ref_file, out_file):
    """Main processing function for merging VCF variants."""
    global vcf, reference, MAX_DISTANCE

    vcf = VCF(vcf_file)
    reference = Fasta(ref_file)

    if out_file:
        output_handle = open(out_file, "w")
    else:
        output_handle = sys.stdout

    vcf.filters["MERGED"] = {"Description": "Merged with neighboring variants", "ID": "MERGED"}

    # Write headers
    for header in vcf.headers:
        if header.startswith("##FILTER"):
            print('##FILTER=<ID=MERGED,Description="Merged with neighboring variants">', file=output_handle)
        print(header, file=output_handle)

    buffer = []
    for variant in vcf:
        if not buffer or should_merge(buffer[-1], variant):
            buffer.append(variant)
        else:
            merged = merge_variants(buffer)
            for merged_variant in merged:
                print(merged_variant, file=output_handle)
            buffer = [variant]

    # Process remaining variants
    merged = merge_variants(buffer)
    for merged_variant in merged:
        print(merged_variant, file=output_handle)

    if out_file:
        output_handle.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge phased SNVs into MNVs.")
    parser.add_argument("vcf_file", help="Input VCF file")
    parser.add_argument("reference", help="Reference genome file")
    parser.add_argument("--out_file", help="Output VCF file", default=None)
    parser.add_argument("--max_distance", type=int, help="Max distance between variants for merging", default=20)

    args = parser.parse_args()

    if args.max_distance:
        MAX_DISTANCE = args.max_distance

    process_vcf(args.vcf_file, args.reference, args.out_file)