
import os

def get_vcf_files(config, vcf_type):
  vcf_files = []
  for var_caller in config["vcf"]:
      if config["vcf"][var_caller]["type"] == vcf_type:
          vcf_files.append(config["vcf"][var_caller]["merged"])
  return vcf_files

def get_sample_type(sample, bio_type):
  type_id = []
  for sample_id in sample:
      if sample[sample_id]["type"] == bio_type:
          type_id.append(sample_id)
  return type_id

def get_result_dir(config):
  result_dir = [config["analysis"]["analysis_dir"], config["analysis"]["sample_id"], config["analysis"]["result"]]
  result_dir = os.path.normpath(os.path.join(*result_dir))
  return result_dir
