import os
import subprocess
import json
import click

from pylatex import Document, Section, Subsection, Subsubsection, Tabular, Math, TikZ, Axis, \
    Plot, Figure, Matrix, Alignat, Command, LongTabu, Package
from pylatex.utils import italic, NoEscape, bold

from BALSAMIC.workflows.run_analysis import get_sample_name


@click.command("report", short_help="Report generator for workflow results")
@click.option(
    '-j',
    '--json-report',
    required=True,
    type=click.Path(),
    help='Input JSON file from workflow output')
@click.option(
    '-c',
    '--json-varreport',
    required=True,
    type=click.Path(),
    help='Input JSON file for variant filters')
@click.option(
    '-r',
    '--rulegraph-img',
    #              required=True,
    type=click.Path(),
    help='Input rulegraph from workflow output')
@click.pass_context
def report(context, json_report, json_varreport, rulegraph_img):

    config = json_report
    sample_config = json.load(open(json_report))
    var_config = json.load(open(json_varreport))

    geometry_options = {
        "tmargin": "2.5cm",
        "lmargin": "1cm",
        "margin": "1cm",
        "paperwidth": "210mm",
        "paperheight": "297mm"
    }
    doc = Document(geometry_options=geometry_options)
    doc.packages.append(Package('lscape'))
    doc.packages.append(Package('longtable'))
    doc.packages.append(Package('float'))
    doc.packages.append(Package('caption', options='labelfont=bf'))
    doc.append(
        NoEscape(
            r'\captionsetup[table]{labelsep=space, justification=raggedright, singlelinecheck=off}'
        ))

    doc.preamble.append(
        Command('title', NoEscape(r'BALSAMIC 0.1 \\ \large Developer Report')))
    doc.preamble.append(
        Command('author', 'Patient ID: ' + get_sample_name(config)))
    doc.preamble.append(Command('date', NoEscape(r'\today')))
    doc.append(NoEscape(r'\maketitle'))

    with doc.create(Section(title='Analysis report', numbering=True)):
        with doc.create(
                Subsection(
                    'Summary of variants and variant callers',
                    numbering=True)):
            doc.append(
                "Placeholder for text about BAM alignment metrics and variant callers. Here comes the info on reads, "
                +
                "QC metrics, align metrics, and general sample information. preferabily in table format."
            )
            doc.append("\n")

            summary_tables = ["TMB", "VarClass", "VarCaller", "VarCallerClass"]
            for i in summary_tables:

                shellcmd = [
                    os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), "..",
                        "R_scripts/VariantReport.R")
                ]
                shellcmd.extend([
                    "--infile", sample_config["vcf"]["merged"]["SNV"],
                    "--genomeSize", sample_config["bed"]["genome_size"],
                    "--type", "latex", "--mode", i
                ])

                outTab = subprocess.check_output(shellcmd)
                doc.append(
                    NoEscape(
                        outTab.decode('utf-8').replace("\\centering",
                                                       "\\small")))
            doc.append(NoEscape(r'\normalsize'))

        with doc.create(Subsection("Summary of MVL report", numbering=True)):
            doc.append(
                "Placeholder for general description of MVL settings. A mention to summary "
                +
                "pipeline, summary of MVL settings. Gene coverage for identified genes should go here. Figures!"
            )
            outVar = dict()
            outCov = dict()
            for i in var_config["filters"]:
                shellcmd = [
                    os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), "..",
                        "R_scripts/VariantReport.R")
                ]
                shellcmd.extend([
                    "--infile", sample_config["vcf"]["merged"]["SNV"], "--dp",
                    var_config["filters"][i]["TUMOR"]["DP"], "--tumorad",
                    var_config["filters"][i]["TUMOR"]["AD"], "--inMVL",
                    var_config["filters"][i]["in_mvl"], "--vartype", "SNP",
                    "--varcaller", ",".join(
                        var_config["filters"][i]["variantcaller"]), "--ann",
                    ",".join(var_config["filters"][i]["annotation"]["SNV"]),
                    "--name", var_config["filters"][i]["name"], "--type",
                    "latex"
                ])
                print(" ".join(shellcmd))
                outVar[i] = subprocess.check_output(shellcmd)
                with doc.create(
                        Subsubsection(
                            var_config["filters"][i]["name"], numbering=True)):
                    if outVar[i] != b"FALSE\n":
                        shellcmd.extend(["--exportGene", "T"])
                        print(" ".join(shellcmd))
                        genes = subprocess.check_output(shellcmd).decode(
                            'utf-8')
                        genes = genes.rstrip("\n\r")

                        for s in sample_config["bed"]["exon_cov"]:
                            shellcmd = [
                                os.path.join(
                                    os.path.dirname(os.path.abspath(__file__)),
                                    "..", "R_scripts/CoverageRep.R")
                            ]
                            shellcmd.extend([
                                "--infile",
                                sample_config["bed"]["exon_cov"][s],
                                "--genename", genes, "--name", s.replace("_","\_"), "--type",
                                "latex"
                            ])
                            print(" ".join(shellcmd))
                            outCov[i] = subprocess.check_output(shellcmd)
                            #doc.append(NoEscape(r'\begin{landscape}'))
                            #longtable instead of tabular makes the table span multiple pages, but the header doesn't span. Occasionally
                            #the alignment also is messed up. There must be a hidden package conflict OR general alignment issues.
                            #doc.append(NoEscape(varreport.replace("{tabular}","{longtable}")))
                            doc.append(
                                NoEscape(outVar[i].decode('utf-8').replace(
                                    "\\centering", "\\tiny")))
                            doc.append(
                                NoEscape(outCov[i].decode('utf-8').replace(
                                    "\\centering", "\\tiny")))
                            #doc.append(NoEscape(r'\end{landscape}'))
                    else:
                        doc.append("No variants were found for this filter")

                doc.append(NoEscape(r'\normalsize'))

        with doc.create(Subsection('Analysis pipeline')):
            with doc.create(Figure(position='h!')) as pipeline_img:
                pipeline_img.add_image(rulegraph_img, width='450px')
                pipeline_img.add_caption('Awesome pipeline')

    with doc.create(Section(title="Appendix", numbering=True)):
        with doc.create(Subsection("MVL settings", numbering=True)):
            fmt = "p{2cm}p{3cm}p{3cm}p{3cm}p{3cm}p{3cm}"
            with doc.create(Tabular(fmt)) as data_table:
                header_row1 = [""]
                for i in var_config["filters"]:
                    header_row1.append(var_config["filters"][i]["name"])
                data_table.add_hline()
                data_table.add_row(header_row1, mapper=[bold])
                data_table.add_hline()
                data_table.add_empty_row()
                column = list(var_config["filters"][next(
                    iter(var_config["filters"]))]["TUMOR"].keys())
                for i in column:
                    row = [i]
                    for j in var_config["filters"]:
                        row.append(var_config["filters"][j]["TUMOR"][i])
                    data_table.add_row(row)

                row = ["MVL"]
                for i in var_config["filters"]:
                    row.append(var_config["filters"][i]["in_mvl"])

                row = ["Variantcallers"]
                for i in var_config["filters"]:
                    row.append("\n".join(
                        var_config["filters"][i]["variantcaller"]))
                data_table.add_row(row)
                data_table.add_hline()

    doc.generate_pdf('full', clean_tex=False)
