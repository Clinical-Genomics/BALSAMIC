from markdown import Markdown
from jinja2 import Environment, FileSystemLoader
from datetime import datetime
from pathlib import Path

from BALSAMIC import __version__ as balsamic_version
from BALSAMIC.utils.constants import REPORT_MODEL


def report_data_population(collected_qc: dict, meta: dict,
                           lang: str = "sv") -> dict:
    """populates a metadata dictionary that contains qc and case/sample information"""
    meta = {
        **meta,
        **{
            "title": "Kvalitetsrapport",
            "subtitle": "klinisk sekvensering av cancer prover",
            "footnote": "Slut på rapporten",
            "bioinformatic": f"BALSAMIC version {balsamic_version}",
            "qc_table_content": {},
            "coverage_table_content": {}
        }
    }

    meta["qc_table_header"] = [v[lang] for x, v in REPORT_MODEL["qc"].items()]
    meta["coverage_table_header"] = [
        v[lang] for x, v in REPORT_MODEL["coverage"].items()
    ]

    for lims_id, analysis_results in collected_qc.items():
        sample_qc = [meta["sample_map"][lims_id], meta["sample_type"][lims_id]]
        sample_cov = [
            meta["sample_map"][lims_id], meta["sample_type"][lims_id]
        ]

        sample_qc = sample_qc + parse_collected_qc(collected_qc=collected_qc, model_param="qc", lims_id=lims_id)
        sample_cov = sample_cov + parse_collected_qc(collected_qc=collected_qc, model_param="coverage", lims_id=lims_id)

        meta["qc_table_content"][lims_id] = sample_qc
        meta["coverage_table_content"][lims_id] = sample_cov

    return meta

def parse_collected_qc(collected_qc: dict, model_param: str, lims_id: str) -> list:
    """parses collect qc and returns model_param"""
    parsed_qc = list()

    for qc_item, qc_value in REPORT_MODEL[model_param].items():
        decimal_point = qc_value["decimal"]
        qc_to_report = collected_qc[lims_id][qc_item]
        if "as_percent" in qc_value:
            qc_to_report = qc_to_report * 100
        qc_to_report = str(round(qc_to_report, decimal_point))
        parsed_qc.append(qc_to_report)

    return parsed_qc

    

def render_html(meta: dict, html_out: str):
    """renders html report from template"""

    p = Path(__file__).parents[1]
    template_path = Path(p, "assets", "report_template").as_posix()

    report_body = render_body(meta=meta, template_path=template_path)

    md_template = Markdown(
        extensions=['meta', 'tables', 'def_list', 'fenced_code'])

    markdown_text = md_template.convert(source=report_body)

    env = Environment(loader=FileSystemLoader(template_path), autoescape=False)

    template = env.get_template("balsamic_report.html")

    html_report = template.render(body=markdown_text, meta=meta)

    with open(html_out, 'w') as f:
        f.write(html_report)
        return html_out


def render_body(meta: dict,
                template_path: str,
                body_template_md: str = "balsamic_report.md") -> str:
    """renders text body of the report from a markdown template"""
    env = Environment(loader=FileSystemLoader(template_path), autoescape=False)

    template = env.get_template(body_template_md)

    report_body = template.render(meta=meta)

    return report_body
