
def vcf_file_list(config):
  vcf_files = []
  for var_caller in config["vcf"]:
      vcf_files.append(config["vcf"][var_caller]["merged"])
  return vcf_files

def get_sample_type(sample, bio_type):
  type_id = []
  for sample_id in sample:
      if sample[sample_id]["type"] == bio_type:
          type_id.append(sample_id)
  return type_id
