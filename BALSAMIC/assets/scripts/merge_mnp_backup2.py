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

global vcf, reference

def ifmerge(v1, v2, max_distance):

    if "SVTYPE" in v1.info or "SVTYPE" in v2.info:
        return False
    # coordinate off by 1
    if v1.chrom != v2.chrom or v2.pos - v1.pos > max_distance + len(v1.ref) - 1:
        return False
    # share the same PID and PGT
    pid1 = v1.samples[0].get("PID", "")
    pid2 = v2.samples[0].get("PID", "")
    pgt1 = v1.samples[0].get("PGT", "")
    pgt2 = v2.samples[0].get("PGT", "")
    # or the same PS to support use case of VCF coming from whatshap
    ps1 = v1.samples[0].get("PS", "")
    ps2 = v2.samples[0].get("PS", "")
    if (pid1 == pid2 and pgt1 == pgt2 and pid1 != "" and pgt1 != "") or (
            ps1 == ps2 and ps1 != ""
    ):
        return True
    return False


def distance(v1, v2):
    if type(v1) == list:
        for v in v1[::-1]:
            return v2.pos - v.pos - len(v.ref) + 1
        return 1000
    if v1.chrom != v2.chrom:
        return 1000
    return v2.pos - v1.pos - len(v1.ref) + 1


def merge_variants(variants, max_distance):
    def _merge_variant_group(variant_group):
        """
        Merge a group of overlapping variants into a single variant.
        """
        group_length = len(variant_group)

        # Ensure no overlapping variants within the group
        for i in range(group_length - 1):
            if variant_group[i + 1].pos - variant_group[i].pos < len(variant_group[i].ref):
                return None

        # Initialize the merged variant
        merged_variant = copy.deepcopy(variant_group[0])
        merged_variant.id = "."
        ref_seq = variant_group[0].ref
        alt_seqs = variant_group[0].alt

        # Combine REF and ALT sequences
        for i in range(1, group_length):
            current_variant = variant_group[i]
            ref_gap = ""
            if current_variant.pos > merged_variant.pos + len(ref_seq):
                ref_gap = reference[merged_variant.chrom][
                          merged_variant.pos + len(ref_seq):current_variant.pos
                          ].seq.upper()

            ref_seq += ref_gap + current_variant.ref
            alt_seqs = [
                alt_seq + ref_gap + current_variant.alt[j]
                for j, alt_seq in enumerate(alt_seqs)
            ]

        # Update QUAL and INFO fields
        qualities = [v.qual for v in variant_group]
        merged_variant.qual = (
            sum(qualities) / group_length if None not in qualities else None
        )

        for key in merged_variant.info.keys():
            values = [v.info.get(key) for v in variant_group]
            if None in values:
                merged_variant.info[key] = None
            else:
                try:
                    merged_variant.info[key] = sum(values) / group_length
                except:
                    merged_variant.info[key] = None

        # Adjust REF and ALT to remove common prefixes
        prefix_index = 0
        while prefix_index < len(ref_seq) - 1:
            if all(
                    ref_seq[prefix_index] == alt_seq[prefix_index]
                    for alt_seq in alt_seqs if prefix_index < len(alt_seq)
            ):
                prefix_index += 1
            else:
                break
        merged_variant.ref = ref_seq[prefix_index:]
        merged_variant.alt = [alt_seq[prefix_index:] for alt_seq in alt_seqs]
        merged_variant.pos += prefix_index

        # Combine filters and ensure uniqueness
        all_filters = {flt for v in variant_group for flt in v.filter}
        if "PASS" in all_filters and len(all_filters) > 1:
            all_filters.discard("PASS")
        merged_variant.filter = list(all_filters)

        # Explicitly update the filter of the original variants
        for original_variant in variant_group:
            original_variant.filter = ["MERGED"]

        # Update sample fields
        for i, sample in enumerate(vcf.samples):
            sample_data = merged_variant.samples[i]
            group_samples = [v.samples[i] for v in variant_group]

            # Handle allele frequency (AF)
            af_values = [samp.get("AF") for samp in group_samples]
            if None not in af_values:
                if isinstance(af_values[0], list):
                    sample_data["AF"] = [
                        sum(af[j] for af in af_values) / group_length
                        for j in range(len(alt_seqs))
                    ]
                else:
                    sample_data["AF"] = sum(af_values) / group_length

            # Handle allele depth (AD)
            ad_values = [samp["AD"] for samp in group_samples]
            sample_data["AD"] = (
                int(sum(ad[0] for ad in ad_values) / group_length),
                int(sum(ad[1] for ad in ad_values) / group_length),
            )

            # Handle AFDP (allele frequency depth)
            afdp_values = [samp.get("AFDP") for samp in group_samples]
            if None not in afdp_values:
                sample_data["AFDP"] = int(sum(afdp_values) / group_length)

        return merged_variant

    merged_variants = []
    pending_variants = []

    for i, variant in enumerate(variants):
        to_merge = [variant]

        for j in range(i + 1, len(variants)):
            if ifmerge(to_merge[-1], variants[j], max_distance):
                to_merge.append(variants[j])
            else:
                break

        if len(to_merge) > 1:
            merged = _merge_variant_group(to_merge)
            if merged:
                pending_variants.append(merged)

        while pending_variants and pending_variants[0].pos <= variant.pos:
            merged_variants.append(pending_variants.pop(0))

        merged_variants.append(variant)

    while pending_variants:
        merged_variants.append(pending_variants.pop(0))

    return merged_variants

def process(vcf_file, ref_file, out_file, max_distance):
    global vcf, reference
    vcf = vcflib.VCF(vcf_file)
    if out_file:
        out_fh = open(out_file, "w")
    else:
        out_fh = sys.stdout
    reference = pyfaidx.Fasta(ref_file)

    # define header
    vcf.filters["MERGED"] = {
        "Description": '"Merged with neighboring variants"',
        "ID": "MERGED",
    }
    filter = False
    for header_line in vcf.headers:
        if header_line.startswith("##FILTER") and not filter:
            print(
                '##FILTER=<ID=MERGED,Description="Merged with neighboring variants">',
                file=out_fh,
            )
            filter = True
        print(header_line, file=out_fh)
    last = []
    for variant in vcf:
        if len(last) == 0:  # empty stack
            if "SVTYPE" in variant.info:
                # ignore SVs
                print(variant, file=out_fh)
            else:
                last.append(variant)
        else:
            if distance(last, variant) <= max_distance:
                last.append(variant)
            else:
                # perform merge on existing stack and reset it
                last = merge_variants(last, max_distance)
                for v in last:
                    print(v, file=out_fh)
                last = []
                if "SVTYPE" in variant.info:
                    print(variant, file=out_fh)
                else:
                    last.append(variant)
    last = merge_variants(last, max_distance)
    for v in last:
        print(v, file=out_fh)

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
def main(vcf_file, reference, out_file, max_distance):
    """
    Merge phased SNVs into MNVs.

    VCF_FILE: Input VCF file to process.
    REFERENCE: The reference genome file.
    """
    # You can replace 'process' with your actual function implementation
    process(vcf_file, reference, out_file, max_distance)


if __name__ == "__main__":
    main()
