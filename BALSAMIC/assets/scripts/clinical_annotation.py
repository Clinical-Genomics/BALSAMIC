from typing import List

import vcfpy
from pandas import read_csv, DataFrame
import itertools
import click

SCORE_HIGH = 15
SCORE_MEDIUM = 10
SCORE_LOW = 5


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
    "--output-tsv",
    "output_tsv",
    required=True,
    type=click.Path(),
    help="Output TSV summary file",
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
def main(input_vcf, output_vcf, output_tsv, clinical_genes_file, known_fusions_file):
    """
    Annotate clinical gene aberrations and fusions in VCF and write a TSV file including clinical summary file.
    """
    annotation = ClinicalAnnotation(clinical_genes_file, known_fusions_file)
    annotate_clinical_aberrations(input_vcf, annotation, output_vcf, output_tsv)


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
        self.deletions = set(
            self.gene_panel[self.gene_panel["reportDeletion"].astype(bool)]["gene"]
        )
        self.duplications = set(
            self.gene_panel[self.gene_panel["reportAmplification"].astype(bool)]["gene"]
        )
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


def include_canonical_gene_info(record, parsed_csq_entries):
    gene_info = {}
    for csq_dict in parsed_csq_entries:
        symbol = csq_dict["SYMBOL"]
        # Only use annotation from Ensemble
        if (
            symbol
            and csq_dict["SOURCE"] == "Ensembl"
            and csq_dict["CANONICAL"] == "YES"
        ):
            gene_annot = "|".join(
                [
                    symbol,
                    csq_dict["EXON"],
                    csq_dict["INTRON"],
                    csq_dict["DISTANCE"],
                    csq_dict["Gene"],
                    csq_dict["Feature"],
                    csq_dict["Feature_type"],
                    csq_dict["BIOTYPE"],
                    csq_dict["IMPACT"],
                    csq_dict["Consequence"],
                ]
            )
            if symbol not in gene_info:
                gene_info[symbol] = []
            gene_info[symbol].append(gene_annot)
    if gene_info:
        record.INFO["CANONICAL"] = gene_info.get(symbol)
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


def include_svlen(record):
    svlen = ""
    if "SVLEN" in record.INFO:
        if type(record.INFO["SVLEN"]) == list:
            if len(record.INFO["SVLEN"]) == 1:
                svlen = abs(record.INFO["SVLEN"][0])
        return abs(record.INFO["SVLEN"])
    elif "END" in record.INFO:
        svlen = abs(record.POS - record.INFO["END"])
    record.INFO["SVLEN"] = svlen
    return record


def include_exonic_intronic_flag(record):
    clinical_genes = record.INFO.get("CLINICAL_GENES")
    canonical = record.INFO.get("CANONICAL")
    if clinical_genes and canonical:
        for gene in record.INFO["CLINICAL_GENES"]:
            for canonical in record.INFO["CANONICAL"]:
                fields = canonical.split("|")
                if gene == fields[0]:
                    if fields[1] != "" or fields[2] != "":
                        record.INFO["EXONINTRON"] = "YES"
                        break
    return record


def include_rank(record):
    # POPULATION DATABASES
    # SWEGENAF >= 0.005 or clin_obs >= 0.001 : -10

    # TECHNICAL EVIDENCE
    # MULTIPLE TOOLS: svdb_origin contains 'manta': +3
    # RELIABLE TOOLS: FOUNDBY > 1 : +5

    # VEP ANNOTATION
    # IMPACT == HIGH: +15
    # IMPACT == MODERATE: +10
    # CONSEQUENCE: ??

    # canonical
    # BIOTYPE: protein coding +5
    # EXONINTRON: + 10
    biotype = 0
    impact = 0
    exonintron = 0
    number_tools = 1 if record.INFO.get("FOUNDBY", 0) > 1 else 0
    reliable_tools = 1 if "manta" in record.INFO.get("svdb_origin", "") else 0
    population_dbs = (
        -10
        if (
            record.INFO.get("SWEGENAF", 0) > 0.005
            or record.INFO.get("clin_obs", 0) > 0.001
        )
        else 0
    )
    # Aberration priority:
    clinical_priority_scores = {
        "known_fusion": 15,
        "promiscuous_fusion": 15,
        "single_hit_fusion": 10,
        "duplication": 15,
        "deletion": 15,
        "in_gene_list": 10,
    }
    max_aberration_score = 0
    aberration_score = []
    for aberration_type in record.INFO.get("CLINICAL_ABERRATION") or []:
        aberration_score.append(clinical_priority_scores[aberration_type])
    if aberration_score:
        max_aberration_score = max(aberration_score)

    if "CSQ" in record.INFO:
        csq_biotype = set(record.INFO.get("BIOTYPE") or {})
        csq_impact = set(record.INFO.get("IMPACT") or {})
        biotype = 1 if "protein_coding" in csq_biotype else 0
        exonintron = 10 if record.INFO.get("EXONINTRON", "") == "YES" else 0

        if "HIGH" in csq_impact:
            impact = 3
        elif "MODERATE" in csq_impact:
            impact = 1
        else:
            impact = 0
    rank = [
        max_aberration_score,
        exonintron,
        biotype,
        impact,
        number_tools,
        reliable_tools,
        population_dbs,
    ]
    record.INFO["CLINICAL_SCORE"] = ",".join(map(str, rank))
    record.INFO["RANK"] = sum(rank)
    return record


def include_clinical_info(record, parsed_csq_entries, annotation):
    svtype = record.INFO["SVTYPE"]
    if svtype == "BND":
        include_clinical_fusion_type_info(record, annotation)
    elif svtype == "DUP":
        include_clinical_aberrations_info(record, annotation.deletions, "duplication")
    elif svtype == "DEL":
        include_clinical_aberrations_info(record, annotation.deletions, "deletion")
    # Any other in the gene panel file
    include_clinical_aberrations_info(record, annotation.all, "in_gene_list")
    return record


def process_record(record, parsed_csq_entries, annotation):
    if "CSQ" not in record.INFO:
        return record
    include_vep_info(record, parsed_csq_entries)
    include_fusion_partner_info(record, parsed_csq_entries)
    include_canonical_gene_info(record, parsed_csq_entries)
    include_svlen(record)
    include_clinical_info(record, parsed_csq_entries, annotation)
    include_exonic_intronic_flag(record)
    include_rank(record)
    return record


csq_fields_to_retrieve = [
    "Allele",
    "SYMBOL",
    "HGNC_ID",
    "IMPACT",
    "Consequence",
    "VARIANT_CLASS",
    "Gene",
    "Feature",
    "Feature_type",
    "EXON",
    "Codons",
    "Amino_acids",
    "BIOTYPE",
    "CANONICAL",
    "SOURCE",
    "Protein_position",
    "GENE_PHENO",
]
csq_fields_to_retrieve = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "VARIANT_CLASS",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "CANONICAL",
    "MANE",
    "MANE_SELECT",
    "MANE_PLUS_CLINICAL",
    "TSL",
    "APPRIS",
    "CCDS",
    "ENSP",
    "SWISSPROT",
    "TREMBL",
    "UNIPARC",
    "UNIPROT_ISOFORM",
    "REFSEQ_MATCH",
    "SOURCE",
    "REFSEQ_OFFSET",
    "GIVEN_REF",
    "USED_REF",
    "BAM_EDIT",
    "GENE_PHENO",
    "SIFT",
    "PolyPhen",
    "DOMAINS",
    "miRNA",
    "HGVS_OFFSET",
    "HGVSg",
    "AF",
    "MAX_AF",
    "MAX_AF_POPS",
    "CLIN_SIG",
    "SOMATIC",
    "PHENO",
    "PUBMED",
    "CHECK_REF",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
    "TRANSCRIPTION_FACTORS",
]

info_fields_to_retrieve = [
    "SVTYPE",
    "FOUNDBY",
    "svdb_origin",
    "Cancer_Somatic_Obs",
    "Cancer_Somatic_Frq",
    "swegen_obs",
    "SWEGENAF",
    "clin_obs",
    "Frq",
    "SVLEN",
    "MAPQ",
    "CLINICAL_SCORE",
    "RANK",
    "CLINICAL_FUSION",
    "GENEA",
    "GENEB",
    "NGENES",
    "CLINICAL_ABERRATION",
    "CLINICAL_GENES",
    "CANONICAL",
]
#    "GENE",
# "CLINICAL_DELETION",
# "CLINICAL_DUPLICATION",
# "CLINICAL_GENE_PANEL",
colnames = (
    [
        "CHROM",
        "POS",
        "REF",
        "ALT",
    ]
    + info_fields_to_retrieve
    + csq_fields_to_retrieve
)


def collect_relevant_information(record) -> list[str]:
    general_info = [record.CHROM, record.POS, record.REF, record.ALT[0].serialize()]
    info_fields = [record.INFO.get(field, "") for field in info_fields_to_retrieve]
    info_fields = [",".join(i) if isinstance(i, list) else str(i) for i in info_fields]
    return general_info + info_fields


def update_header(reader):
    headers = [
        (
            "CLINICAL_ABERRATION",
            "1",
            "String",
            "fusion_pair,promiscuous_fusion,deletion,duplication,other",
        ),
        ("CLINICAL_SCORE", "1", "Integer", "Clinical priority score."),
        ("RANK", "1", "Integer", "Ranking"),
        ("CLINICAL_FUSION", "1", "String", "Clinical fusions detected"),
        ("CLINICAL_DELETION", "1", "String", "Reported deletions"),
        ("CLINICAL_DUPLICATION", "1", "String", "Reported duplications"),
        ("CLINICAL_GENE_PANEL", "1", "String", "Detected gene panel match"),
        ("GENE", "1", "String", "Annotated genes"),
        ("GENEA", "1", "String", "3' fusion genes"),
        ("GENEB", "1", "String", "5' fusion genes"),
        ("NGENES", "1", "Integer", "Number of genes"),
        ("CLINICAL_GENES", "1", "Integer", "Clinical genes"),
        ("BIOTYPE", "1", "String", " "),
        ("IMPACT", "1", "String", ""),
        ("CANONICAL", "1", "String", ""),
        ("EXONINTRON", "1", "String", ""),
    ]
    for id, number, type, desc in headers:
        reader.header.add_info_line(
            vcfpy.OrderedDict(
                [("ID", id), ("Number", number), ("Type", type), ("Description", desc)]
            )
        )


def annotate_clinical_aberrations(vcf_path, annotation, output_vcf, output_tsv):
    reader = vcfpy.Reader.from_path(vcf_path)
    csq_fields = parse_csq_format(reader)
    update_header(reader)
    tsv_lines = []
    with vcfpy.Writer.from_path(output_vcf, reader.header) as writer:
        for record in reader:
            if "CSQ" in record.INFO:
                parsed_csq_entries = [
                    parse_csq_entry(entry, csq_fields) for entry in record.INFO["CSQ"]
                ]
                process_record(record, parsed_csq_entries, annotation)
                tsv_lines.append(collect_relevant_information(record))
            writer.write_record(record)
    with open(output_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(colnames) + "\n")
        for line in tsv_lines:
            line = "\t".join(map(str, line))
            f.write(line + "\n")


if __name__ == "__main__":
    main()
