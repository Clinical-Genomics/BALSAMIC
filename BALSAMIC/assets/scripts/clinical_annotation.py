import vcfpy
from pandas import read_csv, DataFrame
import itertools
import click


@click.command()
@click.option(
    "--input",
    "input_vcf",
    required=True,
    type=click.Path(exists=True),
    help="Input VCF file",
)
@click.option(
    "--output-vcf",
    "output_vcf",
    required=True,
    type=click.Path(),
    help="Output annotated VCF file",
)
@click.option(
    "--genes",
    "clinical_genes_file",
    required=True,
    type=click.Path(),
    help="",
)
@click.option(
    "--fusions",
    "known_fusions_file",
    required=True,
    type=click.Path(),
    help="",
)
def main(input_vcf, output_vcf, clinical_genes_file, known_fusions_file):
    """
    Soft-filters a SV vcf file with a list of clinical genes and a fusion list.
    Aberrations outside the provided lists get soft-filtered with 'off_panel'.
    Additional information is added in the info field:
      - CANNO for clinical annotation: fusion_pair,promiscuous_fusion,single_hit_fusion,in_panel
      - CGENES: Clinical genes detected
      - CFUS: Clinical fusions detected
      - GENEA and GENEB, fusion gene partners

    Clinical information:
    - known_fusions_file: The fusion list file contains two columns with two gene symbols. If value in column 2 is missing
    that fusion is considered promiscuous.
    - clinical_genes_file: The gene file is a TSV file. Where the first column is the gene symbol. Only the first column
    is used at the moment.
    """
    annotation = ClinicalAnnotation(clinical_genes_file, known_fusions_file)
    annotate_clinical_aberrations(input_vcf, annotation, output_vcf)


class ClinicalAnnotation:
    """
    Loads and stores clinical gene annotations and known gene fusion data.
    """

    def __init__(self, clinical_genes_file, known_fusions_file):
        self.gene_panel: set[str]
        self.fusion_pairs: set[tuple[str, str]]
        self.promiscuous_fusion_genes: set[str]
        self.genes_in_fusion_pairs: set(str)
        self._load_annotations(clinical_genes_file, known_fusions_file)

    def _load_annotations(self, clinical_genes_file, known_fusions_file) -> None:
        self._load_clinical_genes(clinical_genes_file)
        self._load_fusion_data(known_fusions_file)

    def _load_clinical_genes(self, clinical_genes_file) -> None:
        gene_panel = read_csv(clinical_genes_file, header=0, sep="\t")
        self.gene_panel = set(gene_panel["gene_symbol"])

    def _load_fusion_data(self, known_fusions_file) -> None:
        self.fusion_pairs: set(tuple(str, str)) = set()
        self.promiscuous_fusion_genes: set(str) = set()
        self.genes_in_fusion_pairs: set(str) = set()

        with open(known_fusions_file, "r", encoding="utf-8") as file:
            for line in file:
                genes: list[str] = line.strip().split("\t")
                if len(genes) == 2:
                    self.fusion_pairs.update(
                        [(genes[0], genes[1]), (genes[1], genes[0])]
                    )
                    self.genes_in_fusion_pairs.update([(genes[0], genes[1])])
                elif len(genes) == 1:
                    self.promiscuous_fusion_genes.add(genes[0])


def parse_csq_format(vcf_reader) -> list[str]:
    """Parse csq format and return an list of the fields included."""
    csq_header = vcf_reader.header.get_info_field_info("CSQ")
    return csq_header.description.split("Format: ")[-1].split("|")


def parse_csq_entry(csq_entry, csq_fields) -> dict:
    """Parse a csq entry and return a dictionary."""
    return dict(zip(csq_fields, csq_entry.split("|")))


def include_fusion_partner_info(record, parsed_csq_entries) -> vcfpy.Record:
    "Add GENE A and GENE B in the info field, with information about fusion partners per allele."
    alt = record.ALT[0].serialize()

    if record.INFO["SVTYPE"] != "BND":
        return record

    genes_a = []
    genes_b = []
    for csq_dict in parsed_csq_entries:
        gene_symbol = csq_dict["SYMBOL"]
        # Gene has to have a symbol
        if not gene_symbol:
            continue
        # RefSeq and Ensemble are available and somewhat redundant, we want to keep the latter (also more complete)
        if csq_dict["SOURCE"] != "Ensembl":
            continue
        if csq_dict["Allele"] == alt:
            genes_a.append(gene_symbol)
        else:
            genes_b.append(gene_symbol)
    if genes_a:
        record.INFO["GENEA"] = list(set(genes_a))
    if genes_b:
        record.INFO["GENEB"] = list(set(genes_b))
    return record


def include_vep_info(record, parsed_csq_entries) -> vcfpy.Record:
    genes = set()
    for csq_dict in parsed_csq_entries:
        if csq_dict["SYMBOL"] and csq_dict["SYMBOL"] != "":
            genes.add(csq_dict["SYMBOL"])
    if genes:
        record.INFO["GENES"] = list(genes)
    return record


def include_clinical_aberrations_info(
    record, annotation_set, aberration_type
) -> vcfpy.Record:
    aberration_types = set(record.INFO.get("CANNO") or {})
    clinical_genes = set(record.INFO.get("CGENES") or {})
    genes = set(record.INFO.get("GENES", []))
    if genes and genes.intersection(annotation_set):
        detected_genes = genes.intersection(annotation_set)
        record.INFO["CGENES"] = list(clinical_genes.union(detected_genes))
        record.INFO["CANNO"] = list(aberration_types.union({aberration_type}))
    return record


def include_clinical_fusion_type_info(record, annotation) -> vcfpy.Record:
    # sourcery skip: use-named-expression
    genes_a = set(record.INFO.get("GENEA", []))
    genes_b = set(record.INFO.get("GENEB", []))
    clinical_genes = record.INFO.get("CGENES", [])
    genes_all = set.union(genes_a, genes_b)
    fusions = []
    aberration_type = record.INFO.get("CANNO") or []
    if genes_a and genes_b:
        # Fusions with two hits should match known pairs or promiscuous fusions
        gene_pairs = set(itertools.product(genes_a, genes_b))
        detected_fusion_pairs = list(gene_pairs.intersection(annotation.fusion_pairs))
        if detected_fusion_pairs:
            for ga, gb in detected_fusion_pairs:
                clinical_genes.extend([ga, gb])
                fusions.append(f"{ga}--{gb}")
            aberration_type.append("known_fusion")
    elif genes_a or genes_b:
        # For fusion with only one annotated hit (either 5 or 3), when that hit is within the fusion list
        detected_1hit = list(genes_all.intersection(annotation.genes_in_fusion_pairs))
        if detected_1hit:
            for g in detected_1hit:
                clinical_genes.append(g)
                fusions.append(f"{g}--?")
            aberration_type.append("single_hit_fusion")
    # Promiscuous fusions, regardless of being two hits or one hit
    detected_promiscuous = list(
        genes_all.intersection(annotation.promiscuous_fusion_genes)
    )
    if detected_promiscuous:
        for g in detected_promiscuous:
            clinical_genes.append(g)
            fusions.append(f"{g}--?")
        aberration_type.append("promiscuous_fusion")
    if fusions:
        record.INFO["CGENES"] = clinical_genes
        record.INFO["CFUS"] = fusions
        record.INFO["CANNO"] = aberration_type
    return record


def include_clinical_info(record, annotation) -> vcfpy.Record:
    svtype = record.INFO["SVTYPE"]
    if svtype == "BND":
        include_clinical_fusion_type_info(record, annotation)
    else:
        # Any other type of abberation in the gene panel file
        include_clinical_aberrations_info(record, annotation.gene_panel, "in_panel")
    return record


def soft_filter_off_panel(record) -> None:
    if "CANNO" not in record.INFO:
        record.add_filter("off_panel")


def clean_record_info_fields(record) -> None:
    """Exclude irrelevant information from vcf record after parsing."""
    record.INFO.pop("GENES", None)


def process_record(record, parsed_csq_entries, annotation) -> vcfpy.Record:
    if "CSQ" not in record.INFO:
        return record
    include_vep_info(record, parsed_csq_entries)
    include_fusion_partner_info(record, parsed_csq_entries)
    include_clinical_info(record, annotation)
    soft_filter_off_panel(record)
    clean_record_info_fields(record)
    return record


def update_header(reader) -> None:
    # Add info fields to header
    headers = [
        (
            "CANNO",
            "1",
            "String",
            "Clinical annotation: fusion_pair,promiscuous_fusion,single_hit_fusion,in_panel",
        ),
        ("GENEA", "1", "String", "3' fusion genes"),
        ("GENEB", "1", "String", "5' fusion genes"),
        ("CGENES", "1", "Integer", "Clinical genes detected"),
        ("CFUS", "1", "Integer", "Clinical fusions detected"),
    ]
    for id, number, type, desc in headers:
        reader.header.add_info_line(
            vcfpy.OrderedDict(
                [("ID", id), ("Number", number), ("Type", type), ("Description", desc)]
            )
        )

    # Add filter to header
    reader.header.add_filter_line(
        vcfpy.OrderedDict(
            [
                ("ID", "off_panel"),
                (
                    "Description",
                    "Variant does not involve genes from the gene or fusion lists",
                ),
            ]
        )
    )


def annotate_clinical_aberrations(vcf_path, annotation, output_vcf) -> None:
    reader = vcfpy.Reader.from_path(vcf_path)
    csq_fields = parse_csq_format(reader)
    update_header(reader)
    with vcfpy.Writer.from_path(output_vcf, reader.header) as writer:
        for record in reader:
            if "CSQ" in record.INFO:
                parsed_csq_entries = [
                    parse_csq_entry(entry, csq_fields) for entry in record.INFO["CSQ"]
                ]
                process_record(record, parsed_csq_entries, annotation)
            writer.write_record(record)


if __name__ == "__main__":
    main()
