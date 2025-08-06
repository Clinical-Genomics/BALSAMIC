from BALSAMIC.utils.references import merge_reference_metadata
from pathlib import Path

def test_seeds_from_existing_and_adds_metadata():
    # GIVEN reference dictionary with one reference matching variant observation metadata and one not
    existing = {
        "clinical_snv_observations": "data/clin.vcf.gz",
        "other_resource": "data/other.vcf.gz",
    }

    # WHEN merging reference dict with metadata
    merged = merge_reference_metadata(existing_refs=existing)

    # THEN filepaths will be assigned as Paths in a subdict
    assert merged["clinical_snv_observations"]["file"] == Path("data/clin.vcf.gz")
    assert merged["other_resource"]["file"] == Path("data/other.vcf.gz")

    # THEN metadata applied only where available
    assert merged["clinical_snv_observations"]["category"] == "clinical"
    assert merged["clinical_snv_observations"]["fields"] == ["Frq", "Obs", "Hom"]
    assert "category" not in merged["other_resource"]
    assert "fields" not in merged["other_resource"]


def test_observation_paths_is_optional():
    # GIVEN a reference dictionary but no variant observation dictionary
    existing = {"other_resource": "data/clin.vcf.gz"}

    # WHEN merging with metadata
    merged = merge_reference_metadata(existing_refs=existing)

    # THEN merging should work without issue and only contain new dictionary with filepath
    assert merged["other_resource"]["file"] == Path("data/clin.vcf.gz")
    assert "category" not in merged["other_resource"]
    assert "fields" not in merged["other_resource"]
