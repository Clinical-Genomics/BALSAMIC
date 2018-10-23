import os
import subprocess
import json
import click
import glob
import yaml

from datetime import datetime
from itertools import chain
from collections import defaultdict

from pylatex import Document, PageStyle, Section, Subsection, Subsubsection, Tabular, Math, TikZ, Axis, \
    Plot, Figure, Matrix, Alignat, Command, LongTabu, Package, Head, Foot, MiniPage, LargeText, MediumText, \
    LineBreak, SmallText, Tabularx, TextColor, MultiColumn, simple_page_number, NewPage, StandAloneGraphic
from pylatex.utils import italic, NoEscape, bold

from BALSAMIC.workflows.run_analysis import get_sample_name
from BALSAMIC.config.sample import get_config
from BALSAMIC.tools import get_packages, get_package_split


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
    required=False,
    default=get_config('MSK_impact'),
    show_default=True,
    type=click.Choice(['MSK_impact', 'MSK_impact_noStrelka']),
    help='Input JSON file for variant filters')
@click.option(
    '-r',
    '--rulegraph-img',
    type=click.Path(),
    help='Input rulegraph from workflow output.')
@click.pass_context
def report(context, json_report, json_varreport, rulegraph_img):

    config = json_report
    sample_config = json.load(open(json_report))
    var_config = json.load(open(get_config(json_varreport)))

    tex_path = os.path.abspath(
        os.path.join(sample_config["analysis"]["analysis_dir"],
                     "delivery_report"))

    if not rulegraph_img:
        rulegraph_img = sample_config['analysis']['dag']

    os.makedirs(tex_path, exist_ok=True)

    geometry_options = {
        "head": "40pt",
        "headheight": "130pt",
        "headsep": "1cm",
        "margin": "1.5cm",
        "bottom": "1.5cm",
        "includeheadfoot": True
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

    #Add first page style
    first_page = PageStyle("header", header_thickness=1)

    #Add Header
    with first_page.create(Head("C")) as mid_header:

        with mid_header.create(
                MiniPage(width=NoEscape(r"0.2\textwidth"),
                         pos='c')) as logo_wrapper:
            logo_file = os.path.join(
                os.path.dirname(__file__), '..', 'assests/cg.png')
            logo_wrapper.append(
                StandAloneGraphic(
                    image_options="width=50px", filename=logo_file))

        with mid_header.create(
                Tabularx(
                    "p{3cm} p{2cm} X X p{4cm} p{3cm}",
                    width_argument=NoEscape(r"0.8\textwidth"))) as mid_table:
            mid_table.add_row(
                [MultiColumn(6, align='r', data=simple_page_number())])
            mid_table.add_row([
                MultiColumn(
                    6, align='c', data=MediumText("Molecular report on"))
            ])
            mid_table.add_row([
                MultiColumn(
                    6, align='c', data=MediumText(get_sample_name(config)))
            ])
            mid_table.add_empty_row()
            mid_table.add_row([
                'gender', "NA", " ", " ", 'Sample recieved:',
                sample_config['analysis']['date']['sample_received']
            ])
            mid_table.add_row([
                'tumor type', "NA", " ", " ", 'Analysis completion:',
                sample_config['analysis']['date']['analysis_finish']
            ])
            mid_table.add_row([
                'analysis type', "NA", " ", " ", 'PDF Report date:',
                datetime.now().strftime("%Y-%m-%d %H:%M")
            ])
            mid_table.add_row(
                ['sample type', "NA", " ", " ", 'Delivery date', "NA"])
            mid_table.add_row([
                'sample origin', "NA", " ", " ", 'Analysis:',
                r'BALSAMIC v' + sample_config['analysis']['BALSAMIC']
            ])

    doc.preamble.append(first_page)

    #End First page

    #    doc.preamble.append(
    #        Command(
    #            'title',
    #            NoEscape(r'BALSAMIC v' + sample_config["analysis"]["BALSAMIC"] +
    #                     r'\\ \large Developer Report')))
    #    doc.preamble.append(
    #        Command('author', 'Patient ID: ' + get_sample_name(config)))
    #    doc.preamble.append(Command('date', NoEscape(r'\today')))
    #    doc.append(NoEscape(r'\maketitle'))
    doc.change_document_style("header")

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
                    "--type", "latex", "--mode", i, "--outfile",
                    os.path.join(tex_path,
                                 sample_config['analysis']['sample_id'])
                ])
                print(" ".join(shellcmd))

                outTab = subprocess.check_output(shellcmd)
                doc.append(
                    NoEscape(
                        outTab.decode('utf-8').replace("\\centering",
                                                       "\\small")))
            doc.append(NoEscape(r'\normalsize'))
            doc.append(NewPage())

        with doc.create(Subsection("Summary of MVL report", numbering=True)):
            doc.append(
                "Placeholder for general description of MVL settings. A mention to summary "
                +
                "pipeline, summary of MVL settings. Gene coverage for identified genes should go here. Figures!"
            )
            outCov = dict()

            cmd_param = defaultdict(list)
            J = defaultdict(list)
            for i in var_config["filters"]:
                cmd_param["TUMOR_DP"].append(
                    var_config["filters"][i]["TUMOR"]["DP"])
                cmd_param["TUMOR_AD"].append(
                    var_config["filters"][i]["TUMOR"]["AD"])
                cmd_param["TUMOR_AFmax"].append(
                    var_config["filters"][i]["TUMOR"]["AF_max"])
                cmd_param["TUMOR_AFmin"].append(
                    var_config["filters"][i]["TUMOR"]["AF_min"])
                cmd_param["TUMOR_inMVL"].append(
                    var_config["filters"][i]["in_mvl"])
                cmd_param["var_type"].append(",".join(
                    ["SNP", "INDEL", "MNP", "OTHER"]))
                cmd_param["varcaller"].append(",".join(
                    var_config["filters"][i]["variantcaller"]))
                cmd_param["ann"].append(
                    ",".join(var_config["filters"][i]["annotation"]["SNV"]) +
                    "," +
                    ",".join(var_config["filters"][i]["annotation"]["INDEL"]))
                cmd_param["name"].append(i.replace("_", "\_"))
                cmd_param["outfile_tex"].append(tex_path + "/" + i + ".tex")
                cmd_param["outfile_gene"].append(tex_path + "/" + i +
                                                 ".genelist")

            for i in cmd_param:
                J[i] = ";".join(cmd_param[i])

            shellcmd = [
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)), "..",
                    "R_scripts/VariantReport.R")
            ]
            shellcmd.extend([
                "--infile", "'" + sample_config["vcf"]["merged"]["SNV"] + "'",
                "--dp", "'" + J["TUMOR_DP"] + "'", "--tumorad",
                "'" + J["TUMOR_AD"] + "'", "--afmax",
                "'" + J["TUMOR_AFmax"] + "'", "--afmin",
                "'" + J["TUMOR_AFmin"] + "'", "--inMVL",
                "'" + J["TUMOR_inMVL"] + "'", "--exclusiveSets", "TRUE",
                "--vartype", "'" + J["var_type"] + "'", "--varcaller",
                "'" + J["varcaller"] + "'", "--ann", "'" + J["ann"] + "'",
                "--name", "'" + J["name"] + "'", "--type", "latex"
            ])

            subprocess.check_output(
                " ".join(shellcmd +
                         ["--outfile", "'" + J["outfile_tex"] + "'"]),
                shell=True)

            print(" ".join(shellcmd +
                           ["--outfile", "'" + J["outfile_tex"] + "'"]))
            subprocess.check_output(
                " ".join(shellcmd + [
                    "--outfile", "'" + J["outfile_gene"] +
                    "'", "--exportGene", "T"
                ]),
                shell=True)

            for c, i in enumerate(var_config["filters"]):
                with doc.create(
                        Subsubsection(
                            var_config["filters"][i]["name"], numbering=True)):
                    print(cmd_param["outfile_tex"])
                    fname = cmd_param["outfile_tex"][c]
                    if os.stat(fname).st_size > 10:
                        #get gene list
                        with open(cmd_param["outfile_gene"][c]) as myfile:
                            genes = myfile.read().replace('\n', '')

                        with open(fname, 'r') as myfile:
                            data = myfile.read()  #.replace('\n', '')

                        #doc.append(NoEscape(r'\begin{landscape}'))
                        #longtable instead of tabular makes the table span multiple pages, but the header doesn't span. Occasionally
                        #the alignment also is messed up. There must be a hidden package conflict OR general alignment issues.
                        #doc.append(NoEscape(varreport.replace("{tabular}","{longtable}")))
                        doc.append(
                            NoEscape(
                                data.replace("\\centering", "\\scriptsize")))

                        for s in sample_config["bed"]["exon_cov"]:
                            shellcmd = [
                                os.path.join(
                                    os.path.dirname(os.path.abspath(__file__)),
                                    "..", "R_scripts/CoverageRep.R")
                            ]
                            shellcmd.extend([
                                "--infile",
                                sample_config["bed"]["exon_cov"][s],
                                "--genename", genes, "--name",
                                s.replace("_", "\_"), "--type", "latex"
                            ])
                            outCov = subprocess.check_output(shellcmd)

                            doc.append(
                                NoEscape(
                                    outCov.decode('utf-8').replace(
                                        "\\centering", "\\scriptsize")))
                        #doc.append(NoEscape(r'\end{landscape}'))
                    else:
                        doc.append("No variants were found for this filter")


#                doc.append(NoEscape(r'\normalsize'))

            doc.append(NewPage())

        with doc.create(Subsection('Coverage report')):
            for s in sample_config["bed"]["target_cov"]:
                with doc.create(Figure(position='h!')) as cov_img:
                    covplot = ".".join(
                        [os.path.join(tex_path, s), "Coverage.pdf"])
                    shellcmd = [
                        os.path.join(
                            os.path.dirname(os.path.abspath(__file__)), "..",
                            "R_scripts/CoveragePlot.R")
                    ]
                    shellcmd.extend([
                        "--infile", sample_config["bed"]["target_cov"][s],
                        "--outfile", covplot, "--title",
                        s.replace("_", "\_")
                    ])
                    subprocess.check_output(shellcmd)
                    cov_img.add_image(covplot, width='450px')
                    cov_img.add_caption('Coverage report for sample ' +
                                        s.replace("_", "\_"))

            doc.append(NewPage())
        with doc.create(Subsection('Analysis pipeline')):
            with doc.create(Figure(position='h!')) as pipeline_img:
                pipeline_img.add_image(rulegraph_img, width='450px')
                pipeline_img.add_caption('BALSAMIC pipeline')
            doc.append(NewPage())

    with doc.create(Section(title="Appendix", numbering=True)):
        with doc.create(Subsection("MVL settings", numbering=True)):
            fmt = "p{3cm}" * (len(var_config["filters"]) + 1)
            with doc.create(Tabular(fmt)) as data_table:
                header_row1 = [""]
                for i in var_config["filters"]:
                    header_row1.append(var_config["filters"][i]["name"])
                data_table.add_hline()
                data_table.add_row(
                    header_row1, mapper=[bold], color="lightgray")
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

        with doc.create(
                Subsection("Bioinformatic tool in pipeline", numbering=True)):
            doc.append(
                "The following Bioinformatic tools were used in the analysis:\n\n"
            )
            with doc.create(Tabular("p{4cm}p{4cm}")) as data_table:
                data_table.add_hline()
                conda_env = glob.glob(
                    os.path.join(
                        os.path.dirname(os.path.abspath(__file__)), "..",
                        "conda_yaml/*.yaml"))

                pkgs = get_package_split(conda_env)

                data_table.add_row(["Package", "Version"], color="lightgray")
                data_table.add_hline()
                data_table.add_row(
                    ["BALSAMIC", sample_config['analysis']['BALSAMIC']])

                for k, v in pkgs.items():
                    data_table.add_row([k, v])
            doc.append(NewPage())

    print(tex_path)
    doc.generate_tex(os.path.join(tex_path, get_sample_name(config)))
    #    doc.generate_pdf(
    #        os.path.join(tex_path, get_sample_name(config)), clean_tex=False)
    shellcmd = [
        "pdflatex", "-output-directory=" + tex_path,
        os.path.join(tex_path, get_sample_name(config)) + ".tex", "1>",
        "/dev/null"
    ]
    #generate_pdf doesn't run AUX files properly and ends up with incorrect total page numbers. So subprocess for
    #pdflatex is called twice instead.

    print(" ".join(shellcmd))
    subprocess.run(" ".join(shellcmd), shell=True)
    subprocess.run(" ".join(shellcmd), shell=True)
