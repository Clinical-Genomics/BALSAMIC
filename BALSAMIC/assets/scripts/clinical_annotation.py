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
        fusion_df = read_csv(self.fusion_file, sep="\t")
        fusion_pairs = fusion_df.dropna(subset=["FiveGene", "ThreeGene"])
        self.fusion_pairs = set(
            zip(fusion_pairs["FiveGene"], fusion_pairs["ThreeGene"])
        )

        self.promiscuous_5 = set(
            fusion_df[fusion_df["FiveGene"].notna() & fusion_df["ThreeGene"].isna()][
                "FiveGene"
            ]
        )
        self.promiscuous_3 = set(
            fusion_df[fusion_df["FiveGene"].isna() & fusion_df["ThreeGene"].notna()][
                "ThreeGene"
            ]
        )


def parse_csq_format(vcf_reader) -> list[str]:
    csq_header = vcf_reader.header.get_info_field_info("CSQ")
    csq_fields = csq_header.description.split("Format: ")[-1].split("|")
    return csq_fields


def get_breakend_side(alt: str) -> int | None:
    """
    Determines if a breakend allele is a 3' or 5' breakend.
    """
    # Check and categorize based on bracket positions
    if alt.startswith("]") and "]" in alt[1:]:
        return 3  # ]chr:pos]X → 3' breakend
    elif alt.startswith("[") and "[" in alt[1:]:
        return 5  # [chr:pos[X → 5' breakend
    elif alt.endswith("[") and "[" in alt[:-1]:
        return 3  # X[chr:pos[ → 3' breakend
    elif alt.endswith("]") and "]" in alt[:-1]:
        return 5  # X]chr:pos] → 5' breakend
    else:
        return None


def switch_breakend(breakend: int | None) -> int | None:
    if breakend == 3:
        return 5
    elif breakend == 5:
        return 3
    else:
        return None


def parse_csq_entry(csq_entry, csq_fields) -> dict:
    return dict(zip(csq_fields, csq_entry.split("|")))


def include_fusion_partner_info(record, parsed_csq_entries):
    alt = record.ALT[0].serialize()
    alt_breakend = get_breakend_side(alt)
    ref_breakend = switch_breakend(alt_breakend)

    if record.INFO["SVTYPE"] != "BND":
        return record

    genes3, genes5 = set(), set()
    for csq_dict in parsed_csq_entries:
        csq_dict["Breakend"] = (
            alt_breakend if csq_dict["Allele"] == alt else ref_breakend
        )
        # print(f"{csq_dict['Allele']}, {csq_dict['Breakend']}, {csq_dict['SYMBOL']}")
        if csq_dict["Breakend"] == 3 and csq_dict["SYMBOL"]:
            genes3.add(csq_dict["SYMBOL"])
        elif csq_dict["Breakend"] == 5 and csq_dict["SYMBOL"]:
            genes5.add(csq_dict["SYMBOL"])
    if genes3:
        record.INFO["GENE3"] = list(genes3)
    if genes5:
        record.INFO["GENE5"] = list(genes5)
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


def include_clinical_aberrations_info(
    record, annotation_set, info_id, clinical_score, aberration_type
) -> bool:
    genes = set(record.INFO.get("GENE", []))
    if genes and genes.intersection(annotation_set):
        record.INFO[info_id] = ",".join(list(genes.intersection(annotation_set)))
        record.INFO["CLINICAL_SCORE"] = clinical_score
        record.INFO["CLINICAL_ABERRATION"] = (
            record.INFO.get("CLINICAL_ABERRATION", "") + "," + aberration_type
        )
    return record


def include_clinical_fusions_info(record, annotation) -> set:
    genes5 = set(record.INFO.get("GENE5", []))
    genes3 = set(record.INFO.get("GENE3", []))
    fusions = ""
    aberration_type = record.INFO.get("CLINICAL_ABERRATION", "")
    if genes5 and genes3:
        # Fusions with two hits should match known pairs or promiscuous fusions
        gene_pairs = set(itertools.product(genes5, genes3))
        detected_fusion_pairs = list(gene_pairs.intersection(annotation.fusion_pairs))
        if detected_fusion_pairs:
            fusions += ",".join([f"{g5}--{g3}" for g5, g3 in detected_fusion_pairs])
            aberration_type += ",known_fusion"
        detected_promiscuous5 = list(genes5.intersection(annotation.promiscuous_5))
        if genes5 and detected_promiscuous5:
            fusions += ",".join([f"{g5}--?" for g5 in detected_promiscuous5])
            aberration_type += ",promiscuous_fusion"
        detected_promiscuous3 = list(genes3.intersection(annotation.promiscuous_3))
        if genes3 and detected_promiscuous3:
            fusions += ",".join([f"?--{g3}" for g3 in detected_promiscuous3])
            aberration_type += ",promiscuous_fusion"
        if fusions:
            clinical_score = 15
            record.INFO["CLINICAL_FUSION"] = fusions
            record.INFO["CLINICAL_SCORE"] = clinical_score
            record.INFO["CLINICAL_ABERRATION"] = aberration_type
            return record
    else:
        # For fusion with only one hit, if that hit is within the fusion list
        if genes5 or genes3:
            annotation_fusion = {
                gene for fusion in annotation.fusion_pairs for gene in fusion
            }
            annotation_fusion = annotation_fusion.union(
                annotation.promiscuous_3, annotation.promiscuous_5
            )
            detected_1hit = list(genes5.union(genes3).intersection(annotation_fusion))
            if detected_1hit:
                aberration_type += ",single_hit_fusion"
                fusions += ",".join(detected_1hit)
                clinical_score = 10
                record.INFO["CLINICAL_FUSION"] = fusions
                record.INFO["CLINICAL_SCORE"] = clinical_score
                record.INFO["CLINICAL_ABERRATION"] = aberration_type
                return record


def include_svlen(record):
    svlen = ""
    if "SVLEN" in record.INFO:
        if type(record.INFO["SVLEN"]) == list:
            if len(record.INFO["SVLEN"]) == 1:
                svlen = abs(record.INFO["SVLEN"][0])
            else:
                print(record.INFO["SVLEN"])
        return abs(record.INFO["SVLEN"])
    elif "END" in record.INFO:
        svlen = abs(record.POS - record.INFO["END"])
    record.INFO["SVLEN"] = svlen
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
    biotype = 0
    impact = 0
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
    clinical_priority = record.INFO.get("CLINICAL_RANK", 0)

    if "CSQ" in record.INFO:
        csq_biotype = set(record.INFO.get("BIOTYPE", []))
        csq_impact = set(record.INFO.get("IMPACT", []))
        biotype = 5 if "protein_coding" in csq_biotype else 0
        if "HIGH" in csq_impact:
            impact = 15
        elif "MODERATE" in csq_impact:
            impact = 10
        else:
            impact = 0
    rank = (
        number_tools
        + reliable_tools
        + population_dbs
        + clinical_priority
        + biotype
        + impact
    )
    record.INFO["RANK"] = rank
    return record


def include_clinical_info(record, parsed_csq_entries, annotation):
    svtype = record.INFO["SVTYPE"]
    if svtype == "BND":
        include_clinical_fusions_info(record, annotation)
    elif svtype == "DUP":
        include_clinical_aberrations_info(
            record, annotation.deletions, "CLINICAL_DUPLICATION", 15, "duplication"
        )
    elif svtype == "DEL":
        include_clinical_aberrations_info(
            record, annotation.deletions, "CLINICAL_DELETION", 15, "deletion"
        )
    if "CLINICAL_SCORE" not in record.INFO:
        # Any other in the gene panel file
        include_clinical_aberrations_info(
            record, annotation.all, "CLINICAL_GENE_PANEL", 5, "other"
        )
    return record


def process_record(record, parsed_csq_entries, annotation):
    if "CSQ" not in record.INFO:
        return record
    include_vep_info(record, parsed_csq_entries)
    include_fusion_partner_info(record, parsed_csq_entries)
    include_svlen(record)
    include_clinical_info(record, parsed_csq_entries, annotation)
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
    "GENE3",
    "GENE5",
    "NGENES",
    "CLINICAL_ABERRATION",
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


def collect_relevant_information(record, parsed_csq_entries, annotation) -> list[str]:
    general_info = [record.CHROM, record.POS, record.REF, record.ALT[0].serialize()]
    info_fields = [record.INFO.get(field, "") for field in info_fields_to_retrieve]
    info_fields = [",".join(i) if isinstance(i, list) else i for i in info_fields]

    relevant_info = []
    for csq_dict in parsed_csq_entries:
        csq_info = [csq_dict.get(field, "") for field in csq_fields_to_retrieve]
        # Only retrieve information for canonical genes in Ensembl
        if (
            csq_dict.get("CANONICAL", "") == "YES"
            and csq_dict.get("SOURCE", "") == "Ensembl"
            and csq_dict.get("SYMBOL") in annotation.all
        ):
            relevant_info.append(general_info + info_fields + csq_info)
    # If there is no canonical ensembl gene, just report one random entry (last one)
    if not relevant_info:
        relevant_info.append(general_info + info_fields + csq_info)
    return relevant_info


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
        ("GENE3", "1", "String", "3' fusion genes"),
        ("GENE5", "1", "String", "5' fusion genes"),
        ("NGENES", "1", "Integer", "Number of genes"),
        ("BIOTYPE", "1", "String", " "),
        ("IMPACT", "1", "String", ""),
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
                tsv_lines.extend(
                    collect_relevant_information(record, parsed_csq_entries, annotation)
                )
            writer.write_record(record)
    with open(output_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(colnames) + "\n")
        for line in tsv_lines:
            line = "\t".join(map(str, line))
            f.write(line + "\n")


if __name__ == "__main__":
    main()
