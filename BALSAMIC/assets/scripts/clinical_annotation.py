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
    Annotate clinical gene aberrations and fusions in VCF and write a TSV file including clinical summary file.
    """
    annotation = ClinicalAnnotation(clinical_genes_file, known_fusions_file)
    annotate_clinical_aberrations(input_vcf, annotation, output_vcf)


class ClinicalAnnotation:
    """
    Loads and stores clinical gene annotations and known gene fusion data.

    Attributes:
        clinical_genes_file (str): Path to the file containing clinical gene panel.
        fusion_file (str): Path to the file containing known fusion gene pairs.
        gene_panel (pd.DataFrame): DataFrame of clinical genes.
        fusion_pairs (Set[tuple[str, str]]): Set of (FiveGene, ThreeGene) fusion gene pairs.
        promiscuous_5 (Set[str]): Genes with only Five' partners.
        promiscuous_3 (Set[str]): Genes with only Three' partners.
    """

    def __init__(self, clinical_genes_file, known_fusions_file):
        self.clinical_genes_file = clinical_genes_file
        self.fusion_file = known_fusions_file

        self.gene_panel: DataFrame
        self.fusion_pairs: set[tuple[str, str]]
        self.promiscuous_5: set[str]
        self.promiscuous_3: set[str]
        self.deletions: set[str]
        self.duplications: set[str]

        self.load_annotations()

    def load_annotations(self) -> None:
        self.load_clinical_genes()
        self.load_fusion_data()

    def load_clinical_genes(self) -> None:
        self.gene_panel = read_csv(self.clinical_genes_file, header=0, sep="\t")
        # self.deletions = set(
        #     self.gene_panel[self.gene_panel["reportDeletion"].astype(bool)]["gene"]
        # )
        # self.duplications = set(
        #     self.gene_panel[self.gene_panel["reportAmplification"].astype(bool)]["gene"]
        # )
        self.all = set(self.gene_panel["gene"])

    def load_fusion_data(self) -> None:
        self.fusion_pairs: set(tuple(str, str)) = set()
        self.promiscuous_fusion_genes: set(str) = set()
        self.genes_in_fusion_pairs: set(str) = set()

        with open(self.fusion_file, "r", encoding="utf-8") as file:
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
    csq_header = vcf_reader.header.get_info_field_info("CSQ")
    csq_fields = csq_header.description.split("Format: ")[-1].split("|")
    return csq_fields


def parse_csq_entry(csq_entry, csq_fields) -> dict:
    return dict(zip(csq_fields, csq_entry.split("|")))


def include_fusion_partner_info(record, parsed_csq_entries):
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


def include_vep_info(record, parsed_csq_entries):
    genes, csq_biotype, csq_impact = set(), set(), set()

    for csq_dict in parsed_csq_entries:
        csq_biotype.add(csq_dict.get("BIOTYPE", ""))
        csq_impact.add(csq_dict.get("IMPACT", ""))
        if csq_dict["SYMBOL"] and csq_dict["SYMBOL"] != "":
            genes.add(csq_dict["SYMBOL"])
    record.INFO["BIOTYPE"] = list(csq_biotype)
    record.INFO["IMPACT"] = list(csq_impact)
    if genes:
        record.INFO["GENE"] = list(genes)
        record.INFO["NGENES"] = str(len(genes))
    return record


def include_clinical_aberrations_info(record, annotation_set, aberration_type) -> bool:
    aberration_types = set(record.INFO.get("CLINICAL_ABERRATION") or {})
    clinical_genes = set(record.INFO.get("CLINICAL_GENES") or {})
    genes = set(record.INFO.get("GENE", []))
    if genes and genes.intersection(annotation_set):
        detected_genes = genes.intersection(annotation_set)
        record.INFO["CLINICAL_GENES"] = list(clinical_genes.union(detected_genes))
        record.INFO["CLINICAL_ABERRATION"] = list(
            aberration_types.union({aberration_type})
        )
    return record


def include_clinical_fusion_type_info(record, annotation) -> set:
    # sourcery skip: use-named-expression
    genes_a = set(record.INFO.get("GENEA", []))
    genes_b = set(record.INFO.get("GENEB", []))
    clinical_genes = record.INFO.get("CLINICAL_GENES", [])
    genes_all = set.union(genes_a, genes_b)
    fusions = []
    aberration_type = record.INFO.get("CLINICAL_ABERRATION") or []
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
        record.INFO["CLINICAL_GENES"] = clinical_genes
        record.INFO["CLINICAL_FUSION"] = fusions
        record.INFO["CLINICAL_ABERRATION"] = aberration_type
    return record


def include_clinical_info(record, annotation):
    svtype = record.INFO["SVTYPE"]
    if svtype == "BND":
        include_clinical_fusion_type_info(record, annotation)
    # elif svtype == "DUP":
    #     include_clinical_aberrations_info(record, annotation.deletions, "duplication")
    # elif svtype == "DEL":
    #     include_clinical_aberrations_info(record, annotation.deletions, "deletion")
    else:
        # Any other in the gene panel file
        include_clinical_aberrations_info(record, annotation.all, "in_gene_list")
    return record


def soft_filter_off_panel(record):
    if not "CLINICAL_ABERRATION" in record.INFO:
        record.add_filter("off_panel")


def process_record(record, parsed_csq_entries, annotation):
    if "CSQ" not in record.INFO:
        return record
    include_vep_info(record, parsed_csq_entries)
    include_fusion_partner_info(record, parsed_csq_entries)
    include_clinical_info(record, annotation)
    soft_filter_off_panel(record)
    return record


def update_header(reader):
    # Add info fields to header
    headers = [
        (
            "CLINICAL_ABERRATION",
            "1",
            "String",
            "fusion_pair,promiscuous_fusion,deletion,duplication,other",
        ),
        ("GENE", "1", "String", "Annotated genes"),
        ("GENEA", "1", "String", "3' fusion genes"),
        ("GENEB", "1", "String", "5' fusion genes"),
        ("NGENES", "1", "Integer", "Number of genes"),
        ("CLINICAL_GENES", "1", "Integer", "Clinical genes"),
        ("BIOTYPE", "1", "String", " "),
        ("IMPACT", "1", "String", ""),
        ("INTRAGENIC", "1", "String", ""),
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


def annotate_clinical_aberrations(vcf_path, annotation, output_vcf):
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
