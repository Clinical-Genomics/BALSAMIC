#!/usr/bin/env python
# -*- coding: utf-8 -*-

#  Copyright 2022, Khurram Maqbool <khurram.maqbool@scilifelab.se>
#
#  This file is free software: you may copy, redistribute and/or modify it
#  under the terms of the GNU General Public License as published by the
#  Free Software Foundation, either version 2 of the License, or (at your
#  option) any later version.
#
#  This file is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

import vcfpy
import sys
import click
import re
import logging


@click.command()
@click.option(
    "-f",
    "--vcf-file",
    required=True,
    type=click.Path(exists=True),
    help="Input Variant Calling Format(VCF) from ascat or delly for CNVs",
)
@click.option(
    "-c",
    "--cnv-caller",
    default="ascat",
    required=True,
    help="CNV caller ascat or delly",
)
@click.option(
    "-t",
    "--tumor-sample",
    default="TUMOR",
    required=False,
    help="Tumor sample ID from vcf header",
)
@click.option(
    "-n",
    "--normal-sample",
    default="NORMAL",
    required=False,
    help="Normal sample ID from vcf header",
)
def main(vcf_file, cnv_caller, tumor_sample, normal_sample):
    vcf = vcfpy.Reader.from_path(vcf_file)
    if cnv_caller == "ascat" or cnv_caller == "delly":
        try:
            get_header(vcf)
            vcfpy.Writer.from_path("/dev/stdout", vcf.header)
            samples = vcf.header.samples.names
            tumor = tumor_sample
            normal = normal_sample
            for sv in vcf:
                SV = get_sv(sv, cnv_caller, samples, tumor, normal)
                write_sv(SV)
        except Exception:
            sys.exit("Please provide valid VCF file")
    else:
        logging.error("Please provide ascat or delly VCF file for CNVs")
        sys.exit()

# get VCF header
def get_header(vcf):
    vcf_header = vcf.header
    for sample in vcf.header.samples.names:
        if sample == "TUMOUR":
            vcf.header.samples.names[vcf.header.samples.name_to_idx["TUMOUR"]] = "TUMOR"
    dup_line = vcfpy.header.HeaderLine("ALT=<ID", 'DUP,Description="Duplication">')
    vcf_header.add_line(dup_line)
    del_line = vcfpy.header.HeaderLine("ALT=<ID", 'DEL,Description="Deletion">')
    vcf_header.add_line(del_line)

# Write variants
def write_sv(sv):
    print(*sv, sep="\t")

# get INFO from VCF
def get_info(sv, svtype):
    info = dict([*sv.INFO.items()])
    INFO = {}
    for info_key, info_value in info.items():
        if type(info_value) is list:
            INFO[info_key] = ",".join(map(str, info_value))
        else:
            if svtype and info_key == "SVTYPE":
                info_value, n = re.subn("[<,>]", "", svtype)
            INFO[info_key] = str(info_value)
    info = ";".join(
        [
            "{0}={1}".format(info_key, info_value)
            for info_key, info_value in INFO.items()
        ]
    )
    return info

# get genotype calls
def get_calls(sv):
    calls = ""
    CALL = {}
    for call in sv.calls:
        for key, value in call.data.items():
            if type(value) is list:
                CALL[key] = ",".join(map(str, value))
            else:
                CALL[key] = value
        calls += str(":".join(map(str, CALL.values())))
        calls += "\t"
    calls = calls.rstrip()
    return calls

# Process variants
def get_sv(sv, caller, samples, tumor, normal):
    ID = ""
    QUAL = ""
    TYPE = ""
    if sv.ID:
        ID = ",".join(map(str, sv.ID))
    else:
        ID = "."
    if sv.QUAL:
        QUAL = sv.QUAL
    else:
        QUAL = "."
    FORMAT = ":".join(sv.FORMAT)
    calls = get_calls(sv)
    if caller == "ascat":
        if sv.call_for_sample[tumor].data.get("TCN") > sv.call_for_sample[
            normal
        ].data.get("TCN") or sv.call_for_sample[tumor].data.get(
            "MCN"
        ) > sv.call_for_sample[
            normal
        ].data.get(
            "MCN"
        ):
            TYPE = "<DUP>"
            sv.add_filter("PASS")
        elif sv.call_for_sample[tumor].data.get("TCN") < sv.call_for_sample[
            normal
        ].data.get("TCN") or sv.call_for_sample[tumor].data.get(
            "MCN"
        ) < sv.call_for_sample[
            normal
        ].data.get(
            "MCN"
        ):
            TYPE = "<DEL>"
            sv.add_filter("PASS")
        else:
            TYPE = "<CNV>"
            sv.add_filter("GERMLINE")
    elif caller == "delly":
        if len(samples) == 1:
            if (sv.call_for_sample)[tumor].data.get("RDCN") > 2:
                TYPE = "<DUP>"
            elif (sv.call_for_sample)[tumor].data.get("RDCN") < 2:
                TYPE = "<DEL>"
        elif len(samples) == 2:
            if (sv.call_for_sample)[tumor].data.get("RDCN") > (sv.call_for_sample)[
                normal
            ].data.get("RDCN"):
                TYPE = "<DUP>"
            elif (sv.call_for_sample)[tumor].data.get("RDCN") < (sv.call_for_sample)[
                normal
            ].data.get("RDCN"):
                TYPE = "<DEL>"
    INFO = get_info(sv, TYPE)
    return sv.CHROM, sv.POS, ID, sv.REF, TYPE, QUAL, sv.FILTER[0], INFO, FORMAT, calls


if __name__ == "__main__":
    main()
