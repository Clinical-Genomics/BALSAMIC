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


def merge(vs):
    def _merge(vv):
        len_vv = len(vv)
        for i in range(len_vv - 1):
            if vv[i + 1].pos - vv[i].pos < len(vv[i].ref):
                return None
        v = copy.deepcopy(vv[0])
        v.id = "."
        ref = vv[0].ref
        alt = vv[0].alt
        for i in range(1, len_vv):
            vi = vv[i]
            ref_gap = ""
            if vi.pos - (v.pos + len(ref)) > 0:
                ref_gap = reference[v.chrom][v.pos + len(ref): vi.pos].seq.upper()
            ref = ref + ref_gap + vi.ref
            alt = [alt[j] + ref_gap + vi.alt[j] for j in range(len(alt))]

        qual_list = [vi.qual for vi in vv]
        v.qual = sum(qual_list) / len_vv if None not in qual_list else None
        for k in v.info.keys():
            vk = [vi.info.get(k) for vi in vv]
            if None in vk:
                v.info[k] = None
            else:
                try:
                    v.info[k] = sum(vk) / len_vv
                except:
                    v.info[k] = None
        # remove common heading letters
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
        #        print(v.alt)
        v.pos += i

        # Combine filters and ensure uniqueness
        all_filters = {flt for vi in vv for flt in vi.filter}
        if "PASS" in all_filters and len(all_filters) > 1:
            all_filters.discard("PASS")
        v.filter = list(all_filters)

        # variants grouped by samples
        for i in range(len(vcf.samples)):
            t = v.samples[i]
            vv_samples = [vs.samples[i] for vs in vv]
            if None in [vi.get("AF") for vi in vv_samples]:
                pass
            elif type(vv_samples[0]["AF"]) == list:
                t["AF"] = [
                    sum([vsi["AF"][j] for vsi in vv_samples]) / len_vv
                    for j in range(len(alt))
                ]
            else:
                t["AF"] = sum([vsi["AF"] for vsi in vv_samples]) / len_vv
            ads = [(vi["AD"][0], vi["AD"][1]) for vi in vv_samples]
            t["AD"] = (
                int(sum([vi["AD"][0] for vi in vv_samples]) / len_vv),
                int(sum([vi["AD"][1] for vi in vv_samples]) / len_vv),
            )
            afdp = [vi.get("AFDP") for vi in vv_samples]
            if None not in afdp:
                t["AFDP"] = int(sum(afdp) / len_vv)
        _ = vcf.format(v)
        for vi in vv:
            _ = vcf.format(vi)
        return v

    vlist = []
    to_push = []
    for i in range(len(vs)):
        vi = vs[i]
        to_merge = [vi]
        for j in range(i + 1, len(vs)):
            vj = vs[j]
            if ifmerge(to_merge[-1], vj):
                to_merge.append(vj)
        if len(to_merge) > 1:
            merged = _merge(to_merge)
            if merged:
                to_push.append(merged)
        while len(to_push):
            if to_push[0].pos <= vi.pos:
                vlist.append(to_push.pop(0))
            else:
                break
        vlist.append(vi)
    while len(to_push):
        vlist.append(to_push.pop(0))
    return vlist

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
                last = merge(last)
                for v in last:
                    print(v, file=out_fh)
                last = []
                if "SVTYPE" in variant.info:
                    print(variant, file=out_fh)
                else:
                    last.append(variant)
    last = merge(last)
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
