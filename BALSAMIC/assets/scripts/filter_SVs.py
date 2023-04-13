#!/usr/bin/env python

import vcfpy
import click
import logging

def calc_af_normal(sample_values, sample):
    max_af = 0
    for variant in sample_values:
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
            af = float((DV + RV) / total_cov)
        except Exception as e:
            logging.warning(f"Exception: {e} setting AF to 0")
            af = 0
        if af > max_af:
            max_af = af
    return max_af

def calc_af_tumor(sample_values, sample):
    fields = sample_values[0].split("|")
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
        af = float((DV + RV) / total_cov)
    except Exception as e:
        logging.warning(f"Exception: {e} setting AF to 0")
        af = 0
    return af

def find_ctg(sample_info_values):
    fields = sample_info_values[0].split("|")
    ctg = False
    for field in fields:
        field_name_value = field.split(":")
        field_name = field_name_value[0]
        if field_name == "CTG":
            ctg_value = field_name_value[1]
            if ctg_value != ".":
                ctg = True
    return ctg


def filter_vcf(vcf_path):
    vcf = vcfpy.Reader.from_path(vcf_path)
    vcf_name = os.path.basename(vcf_path)
    write_file = f"filtered_{vcf_name}"
    writer = vcfpy.Writer.from_path(write_file, vcf.header)

    for sv in vcf:
        info = dict([*sv.INFO.items()])

        # Set default values
        ctg_t, ctg_n = False, False
        af_t, af_n = 0, 0

        for info_key, info_value in info.items():
            if info_key == "TUMOR_PASS_SAMPLE":
                af_t = calc_af_tumor(info_value, "tumor")
            if info_key == "NORMAL_PASS_SAMPLE":
                af_n = calc_af_normal(info_value, "normal")
            if info_key == "TUMOR_PASS_INFO":
                ctg_t = find_ctg(info_value)
            if info_key == "NORMAL_PASS_INFO":
                ctg_n = find_ctg(info_value)

        # Set filter statuses
        normal_af_status = ""
        if af_t == 0 and not ctg_t:
            sv.add_filter('normal_variant')
        else:
            # Regardless of CTG, set filter if AF_T / AF_N > 0.25
            if af_t > 0 and af_n > 0:

                if float(af_n / af_t) > 0.25:
                    normal_af_status = "failed"
                    sv.add_filter('high_normal_af_fraction')
                else:
                    normal_af_status = "pass"

                if af_n > 0.25:
                    normal_af_status = "failed"
                    sv.add_filter('high_normal_af')

                if ctg_n and normal_af_status != "pass":
                    sv.add_filter('in_normal')

        writer.write_record(sv)