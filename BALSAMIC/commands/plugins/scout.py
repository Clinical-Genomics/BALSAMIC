import os
import logging
import glob
import json
import yaml
import click
import shutil
import sys
import datetime

from BALSAMIC.utils.rule import get_result_dir

LOG = logging.getLogger(__name__)


@click.command("scout", short_help="Creates a scout config yaml file.")
@click.option("--sample-config",
              required=True,
              help="Sample config file. Output of balsamic config sample")
@click.option("--snv-vcf",
              default="vcfmerge",
              help="variant caller to load as vcf_cancer")
@click.option("--tumor", default="TUMOR", help="sample name for tumor sample")
@click.option("--normal",
              default="NORMAL",
              help="sample name for normal sample")
@click.option("--sv-vcf",
              default="manta",
              help="variant caller to load as vcf_cancer_sv")
@click.option("--customer-id",
              required=True,
              help="customer id for scout config")
@click.pass_context
def scout(context, sample_config, snv_vcf, sv_vcf, customer_id, tumor, normal):
    '''
    Create a scout config.yaml file
    '''
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    LOG.info('Adding scout cancer template to delivery directory')

    with open(sample_config, 'r') as fn:
        sample_config = json.load(fn)
    case_name = sample_config['analysis']['case_id']
    capture_kit = os.path.basename(sample_config['panel']['capture_kit'])

    result_dir = get_result_dir(sample_config)
    dst_directory = os.path.join(result_dir, 'scout')
    if not os.path.exists(dst_directory):
        LOG.debug('Creatiing delivery_report directory')
        os.makedirs(dst_directory)

    scout_config_src = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "../../assets/scout_config_template.yaml")
    with open(scout_config_src, 'r') as fn:
        scout_config = yaml.load(fn, Loader=yaml.SafeLoader)

    deliver_wildcards = {
        'bam': {
            'tumor': 'tumor.merged.bam',
            'normal': 'normal.merged.bam'
        },
        'vep': {
            'vcf_cancer_sv': f'SV.somatic.{case_name}.{sv_vcf}.vcf.gz',
            'vcf_cancer': f'SNV.somatic.{case_name}.{snv_vcf}.vcf.gz'
        },
        'qc': 'multiqc_report.html'
    }

    if sample_config['analysis']['analysis_type'] == 'single':
        del deliver_wildcards['bam']['normal']
        scout_config['samples'].pop(1)

    scout_config['owner'] = customer_id
    scout_config['family'] = case_name
    scout_config['family_name'] = case_name

    scout_config['analysis_date'] = str(datetime.datetime.now())

    #scout_config['vcf_cancer_sv'] = os.path.join(result_dir, 'vep', deliver_wildcards['vep']['vcf_cancer_sv'])
    scout_config.pop('vcf_cancer_sv')
    scout_config['vcf_cancer'] = os.path.join(
        result_dir, 'vep', deliver_wildcards['vep']['vcf_cancer'])

    scout_config['multiqc'] = os.path.join(result_dir, 'qc',
                                           deliver_wildcards['qc'])

    # scout sample info for tumor
    scout_config['samples'][0]['bam_path'] = os.path.join(
        result_dir, 'bam', deliver_wildcards['bam']['tumor'])
    scout_config['samples'][0]['capture_kit'] = capture_kit
    scout_config['samples'][0]['sample_id'] = tumor
    scout_config['samples'][0]['sample_name'] = tumor

    if sample_config['analysis']['analysis_type'] == 'paired':
        scout_config['samples'][1]['bam_path'] = os.path.join(
            result_dir, 'bam', deliver_wildcards['bam']['normal'])
        scout_config['samples'][1]['capture_kit'] = capture_kit
        scout_config['samples'][1]['sample_id'] = normal
        scout_config['samples'][1]['sample_name'] = normal

    scout_config_dst = os.path.join(
        dst_directory, sample_config['analysis']['case_id'] + ".scout.yaml")

    LOG.debug('Creating scout config %s', scout_config_dst)
    with open(scout_config_dst, 'w') as f:
        yaml.dump(scout_config, f, default_flow_style=False)
    LOG.info('Scout config template is successfully created: %s',
             scout_config_dst)
