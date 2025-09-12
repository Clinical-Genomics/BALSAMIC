import BALSAMIC.utils.references as m
from pathlib import Path
import copy
import pytest


def test_merges_static_metadata(monkeypatch):
    # Arrange
    monkeypatch.setattr(
        m,
        "VARIANT_OBSERVATION_METAVALUES",
        {
            "clinvar": {
                "fields": ["CLNSIG"],
                "ops": ["self"],
                "names": ["CLNSIG"],
                "category": "RESEARCH",
            }
        },
        raising=False,
    )
    refs = {"clinvar": "/ref/clinvar.vcf.gz"}

    # Act
    out = m.add_reference_metadata(refs)

    # Assert
    assert set(out.keys()) == {"clinvar"}
    entry = out["clinvar"]
    assert isinstance(entry["file"], Path)
    assert entry["file"] == Path("/ref/clinvar.vcf.gz")
    # static metadata merged
    assert entry["fields"] == ["CLNSIG"]
    assert entry["ops"] == ["self"]
    assert entry["names"] == ["CLNSIG"]
    assert entry["category"] == "RESEARCH"


def test_handles_key_without_static_metadata(monkeypatch):
    # Arrange: no metadata for 'unknown_key'
    monkeypatch.setattr(m, "VARIANT_OBSERVATION_METAVALUES", {}, raising=False)
    refs = {"unknown_key": "/ref/unknown.vcf.gz"}

    # Act
    out = m.add_reference_metadata(refs)

    # Assert: only 'file' should be present
    assert set(out.keys()) == {"unknown_key"}
    assert out["unknown_key"] == {"file": Path("/ref/unknown.vcf.gz")}


def test_multiple_keys_mix_with_and_without_metadata(monkeypatch):
    monkeypatch.setattr(
        m,
        "VARIANT_OBSERVATION_METAVALUES",
        {"gnomad_variant": {"fields": ["AF"], "ops": ["self"], "names": ["GNOMADAF"]}},
        raising=False,
    )
    refs = {
        "gnomad_variant": "/ref/gnomad.vcf.gz",
        "custom_db": "/ref/custom.vcf.gz",
    }

    out = m.add_reference_metadata(refs)

    assert set(out.keys()) == {"gnomad_variant", "custom_db"}
    assert out["gnomad_variant"]["file"] == Path("/ref/gnomad.vcf.gz")
    # Metadata merged for gnomad_variant
    assert out["gnomad_variant"]["fields"] == ["AF"]
    # No metadata for custom_db
    assert out["custom_db"] == {"file": Path("/ref/custom.vcf.gz")}


def test_input_dict_is_not_mutated(monkeypatch):
    monkeypatch.setattr(m, "VARIANT_OBSERVATION_METAVALUES", {}, raising=False)
    refs = {"a": "/x.vcf.gz"}
    before = copy.deepcopy(refs)

    _ = m.add_reference_metadata(refs)

    assert refs == before  # function should be pure w.r.t input


def test_global_metadata_not_mutated(monkeypatch):
    meta = {"clinvar": {"fields": ["CLNSIG"]}}
    monkeypatch.setattr(m, "VARIANT_OBSERVATION_METAVALUES", meta, raising=False)
    refs = {"clinvar": "/ref/clinvar.vcf.gz"}

    out = m.add_reference_metadata(refs)

    # Ensure function didn't mutate the global meta
    assert meta == {"clinvar": {"fields": ["CLNSIG"]}}
    # Sanity: still merged in the returned structure
    assert out["clinvar"]["fields"] == ["CLNSIG"]
