import os
import sys
import logging
import glob
import json
import yaml
import click
import copy
import snakemake
from collections import defaultdict
from yapf.yapflib.yapf_api import FormatFile

from BALSAMIC.utils.cli import get_from_two_key
from BALSAMIC.utils.cli import merge_dict_on_key
from BALSAMIC.utils.cli import get_file_extension
from BALSAMIC.utils.cli import find_file_index
from BALSAMIC.utils.cli import write_json
from BALSAMIC.utils.cli import get_snakefile
from BALSAMIC.utils.cli import CaptureStdout
from BALSAMIC.utils.rule import get_result_dir
from BALSAMIC.utils.exc import BalsamicError

LOG = logging.getLogger(__name__)


@click.command(
    "deliver",
    short_help="Creates a YAML file with output from variant caller and alignment.",
)
@click.option(
    "--sample-config",
    required=True,
    help="Sample config file. Output of balsamic config sample",
)
@click.pass_context
def deliver(context, sample_config):
    """
    cli for deliver sub-command.
    Writes <case_id>.hk in result_directory.
    """
    LOG.info(f"BALSAMIC started with log level {context.obj['loglevel']}.")
    LOG.debug("Reading input sample config")
    with open(sample_config, "r") as fn:
        sample_config_dict = json.load(fn)

    result_dir = get_result_dir(sample_config_dict)
    dst_directory = os.path.join(result_dir, "delivery_report")
    LOG.info("Creatiing delivery_report directory")
    os.makedirs(dst_directory, exist_ok=True)

    yaml_write_directory = os.path.join(result_dir, "delivery_report")
    os.makedirs(yaml_write_directory, exist_ok=True)

    analysis_type = sample_config_dict["analysis"]["analysis_type"]
    sequencing_type = sample_config_dict["analysis"]["sequencing_type"]
    snakefile = get_snakefile(analysis_type, sequencing_type)

    with CaptureStdout() :
        snakemake.snakemake(
            snakefile=snakefile,
            config={"delivery": "True"},
            dryrun=True,
            configfiles=[sample_config],
            quiet=True,
        )

    delivery_file_name = os.path.join(
        yaml_write_directory, sample_config_dict["analysis"]["case_id"] + ".hk"
    )
    delivery_file_raw = os.path.join(
        yaml_write_directory,
        sample_config_dict["analysis"]["case_id"] + "_delivery_raw.hk",
    )
    with open(delivery_file_raw, "r") as fn:
        delivery_file_raw_dict = json.load(fn)

    delivery_file_ready = os.path.join(
        yaml_write_directory,
        sample_config_dict["analysis"]["case_id"] + "_delivery_ready.hk",
    )
    with open(delivery_file_ready, "r") as fn:
        delivery_file_ready_dict = json.load(fn)

    with CaptureStdout() as summary:
        snakemake.snakemake(
            snakefile=snakefile,
            config={"delivery": "True"},
            dryrun=True,
            summary=True,
            configfiles=[sample_config],
            quiet=True,
        )
    summary = [i.split("\t") for i in summary]
    summary_dict = [dict(zip(summary[0], value)) for value in summary[1:]]

    output_files_merged_interm = merge_dict_on_key(
        dict_1=summary_dict, dict_2=delivery_file_raw_dict, by_key="output_file"
    )
    output_files_merged = merge_dict_on_key(
        dict_1=output_files_merged_interm,
        dict_2=delivery_file_ready_dict,
        by_key="output_file",
    )

    delivery_json = dict()
    delivery_json["files"] = list()

    for item in output_files_merged:
        if "date" in item:
            warnings = list()
            interm_dict = copy.deepcopy(item)
            interm_dict["path"] = interm_dict.get("output_file")
            interm_dict["step"] = interm_dict.get("rulename", "unknown")

            file_path_index = find_file_index(interm_dict["path"])
            if len(file_path_index) > 1:
                LOG.warning("More than one index found for %s" % interm_dict["path"])
                LOG.warning("Taking %s index file" % list(file_path_index)[0])
            interm_dict["path_index"] = file_path_index[0] if file_path_index else ""

            interm_dict["format"] = get_file_extension(interm_dict["path"])
            interm_dict["tag"] = ",".join(interm_dict.get("wildcard_name", ["unknown"]))
            interm_dict["id"] = "unknown"

            delivery_id = list()
            delivery_id.append(get_from_two_key(
                interm_dict,
                from_key="wildcard_name",
                by_key="wildcard_value",
                by_value="sample",
                default=None,
            ))

            delivery_id.append(get_from_two_key(
                interm_dict,
                from_key="wildcard_name",
                by_key="wildcard_value",
                by_value="case_name",
                default=None,
            ))

            delivery_id=list(filter(None,delivery_id))
            if len(delivery_id) > 1:
                LOG.error(f"Ambiguous delivery id. Wilcard has both: {delivery_id}")
                raise BalsamicError("Delivery file parsing process failed.")
            
            if delivery_id:
                interm_dict["id"] = delivery_id[0]

            delivery_json["files"].append(interm_dict)

    LOG.debug(f"Writing output file {delivery_file_name}")

    write_json(delivery_json, delivery_file_name)
    with open(delivery_file_name + ".yaml", "w") as fn:
        yaml.dump(delivery_json, fn, default_flow_style=False)


#    for entries in deliveries:
#        delivery_file = entries[2]
#        if os.path.isfile(delivery_file):
#            print(Color(u"[{green}\u2713{/green}]"), entries)
#        else:
#            print(Color(u"[{red}\u2717{/red}]"), entries)
#        break
# print(dag)
