import os
import subprocess
import json
import click

from pylatex import Document, Section, Subsection, Subsubsection, Tabular, Math, TikZ, Axis, \
    Plot, Figure, Matrix, Alignat, Command, LongTabu, Package
from pylatex.utils import italic, NoEscape, bold


def extractSampleName(config):

    sample_config = json.load(open(config))
    sample_name = sample_config["analysis"]["sample_id"]

    return sample_name


def extractVCFstat():

    variant_stat = ["21", "15479624", "G", "T", "SNV", "300/0.01"]

    return variant_stat


def countReads(config, flag):

    if flag == 'fastq':
        sample_read = config["Metrics"]["fastq"]
    elif flag == 'bam':
        sample_read = config["Metrics"]["BAM"]
    elif flag == 'bam_rmdup':
        sample_read = config["Metrics"]["BAM_rmdup"]
    else:
        sample_read = "NA"

    return sample_read


def countVCF(config, var_caller, flag):

    if flag == 'all':
        vcf_count = config["vcf"][var_caller]["All"]
    elif flag == 'pass':
        vcf_count = config["vcf"][var_caller]["Passed"]
    elif flag == 'annot':
        vcf_count = config["vcf"][var_caller]["Annot"]["Total"]
    elif flag == 'annot_type':
        vcf_count = config["vcf"][var_caller]["Annot"]["misense"] + \
                    "missense and " + \
                    config["vcf"][var_caller]["Annot"]["sense"] + \
                    "sense"
    elif flag == 'pathological':
        vcf_count = config["vcf"][var_caller]["Annot"]["pathological"]
    else:
        vcf_count = str("NA")

    return vcf_count


def varCallerSummary(config, var_caller):

    row = [var_caller]

    row.append(config["vcf"][var_caller]["All"])

    row.append(config["vcf"][var_caller]["Passed"])

    row.append(config["vcf"][var_caller]["Annot"]["misense"])
    row.append(config["vcf"][var_caller]["Annot"]["sense"])
    row.append(config["vcf"][var_caller]["Annot"]["pathological"])

    return row


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
        Command('author', 'Patient ID: ' + extractSampleName(config=config)))
    doc.preamble.append(Command('date', NoEscape(r'\today')))
    doc.append(NoEscape(r'\maketitle'))

    with doc.create(Section(title='Analysis report', numbering=True)):
        with doc.create(
                Subsection('Summary of alignment report', numbering=True)):
            doc.append(
                "Placeholder for text about BAM alignment. Here comes the info on reads, "
                +
                "QC metrics, align metrics, and general sample information. preferabily in table format."
            )
            doc.append("\n")
            fmt = "p{2cm}p{3cm}p{3cm}p{3cm}p{3cm}p{3cm}"
            with doc.create(Tabular(fmt, pos="H")) as data_table:
                header_row1 = [""]
                for i in var_config["filters"]:
                    header_row1.append(var_config["filters"][i]["name"])
                print(header_row1)
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
                    print(row)
                    data_table.add_row(row)
            data_table.add_hline()

        with doc.create(Subsection("Summary of MVL report", numbering=True)):
            doc.append(
                "Placeholder for general description of MVL settings. A mention to summary "
                +
                "pipeline, summary of MVL settings. Gene coverage for identified genes should go here. Figures!"
            )
            outVar = dict()
            outCov = dict()
            for i in ["set_2"]:#var_config["filters"]:
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
                        genes = subprocess.check_output(shellcmd).decode('utf-8')
                        shellcmd = [
                            os.path.join(
                                os.path.dirname(os.path.abspath(__file__)), "..",
                                "R_scripts/CoverageRep.R")
                        ]
                        shellcmd.extend([
                            "--infile", sample_config["bed"]["exon_cov"],
                            "--genename", genes, "--type", "latex"])
                        print(genes)
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
                break

                doc.append(NoEscape(r'\normalsize'))

        with doc.create(
                Subsection(
                    'Summary of variant calling report', numbering=True)):
            doc.append("Variant calling reports and summaries go here")

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
                print(header_row1)
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
                    print(row)
                    data_table.add_row(row)
                
                row = ["MVL"]
                for i in var_config["filters"]:
                    row.append(var_config["filters"][i]["in_mvl"])
                
                row = ["Variantcallers"]
                for i in var_config["filters"]:
                    row.append("\n".join(var_config["filters"][i]["variantcaller"]))
                data_table.add_row(row)
                data_table.add_hline()

    doc.generate_pdf('full', clean_tex=False)


#            for var_caller in sample_config["vcf"]:
#                with doc.create(Subsubsection(var_caller, numbering=False)):
#                    doc.append("Variant calling resulted in " + countVCF(
#                        config=sample_config, var_caller=var_caller, flag="all") +
#                               " variants, of which ")
#                    doc.append(
#                        countVCF(
#                            config=sample_config, var_caller=var_caller, flag="pass")
#                        + " passed our filters (FILTER=PASS). ")
#                    doc.append("From this list " + countVCF(
#                        config=sample_config, var_caller=var_caller, flag="annot") +
#                               " were annotated using VEP and SNPeff. ")
#                    doc.append("In summary, there were " + countVCF(
#                        config=sample_config,
#                        var_caller=var_caller,
#                        flag="annot_type") + " variants, ")
#                    doc.append("where " + countVCF(
#                        config=sample_config,
#                        var_caller=var_caller,
#                        flag="pathological"
#                    ) + " of all predicted to have high functional impact.")
#                    fmt = "X[r] X[r] X[r] X[r] X[r] X[r]"
#                    with doc.create(LongTabu(fmt)) as data_table:
#                        header_row1 = [
#                            "Chrom", "Pos", "RefAllele", "AltAllele", "Type", "DP/AF"
#                        ]
#                        data_table.add_hline()
#                        data_table.add_row(header_row1, mapper=[bold])
#                        data_table.add_hline()
#                        data_table.add_empty_row()
#                        # data_table.end_table_header()
#                        row = extractVCFstat()
#                        for i in range(2):
#                            data_table.add_row(row)
#                        data_table.add_hline()

#        fmt = "X[r] X[r] X[r] X[r] X[r] X[r]"
#        with doc.create(LongTabu(fmt)) as data_table:
#            header_row1 = [
#                "Variant Caller", "Total", "Passed", "Misense", "Sense", "Patholigical"
#            ]
#            data_table.add_hline()
#            data_table.add_row(header_row1, mapper=[bold])
#            data_table.add_hline()
#            data_table.add_empty_row()
#            #data_table.add_caption("tmp")
#            row = extractVCFstat()
#            for var_caller in sample_config["vcf"]:
#                row = varCallerSummary( config=sample_config, var_caller=var_caller)
#                data_table.add_row(row)
#            data_table.add_hline()
