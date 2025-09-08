from unittest import mock
import logging
import snakemake

from BALSAMIC.constants.analysis import AnalysisWorkflow
from BALSAMIC.utils.cli import get_snakefile

MOCKED_OS_ENVIRON = "os.environ"


def test_workflow_tumor_only_tga_hg19(tumor_only_config, snakemake_runner):
    snakefile = get_snakefile("single", "balsamic")
    rc, out, err, cmd = snakemake_runner(
        snakefile=snakefile,
        configfile=tumor_only_config,
    )
    assert (
        rc == 0
    ), f"Snakemake failed: {' '.join(cmd)}\n--- STDOUT ---\n{out}\n--- STDERR ---\n{err}"
    assert "igh_dux4_detection_tumor_only" not in (out + err)


def test_workflow_tumor_normal_tga_hg19(tumor_normal_config, snakemake_runner):
    snakefile = get_snakefile("paired", "balsamic")
    config_json = tumor_normal_config

    rc, out, err, cmd = snakemake_runner(
        snakefile=snakefile,
        configfile=config_json,
    )

    # workflow ran (dry-run) successfully
    assert (
        rc == 0
    ), f"Snakemake failed: {' '.join(cmd)}\n--- STDOUT ---\n{out}\n--- STDERR ---\n{err}"
    # rule should NOT be present
    assert "igh_dux4_detection_tumor_normal" not in (out + err)


def test_workflow_tumor_only_wgs_hg19(tumor_only_wgs_config, snakemake_runner):
    snakefile = get_snakefile("single", "balsamic")
    config_json = tumor_only_wgs_config

    rc, out, err, cmd = snakemake_runner(
        snakefile=snakefile,
        configfile=config_json,
    )

    assert (
        rc == 0
    ), f"Snakemake failed: {' '.join(cmd)}\n--- STDOUT ---\n{out}\n--- STDERR ---\n{err}"
    # rule SHOULD be present
    assert "igh_dux4_detection_tumor_only" in (out + err)


def test_workflow_tumor_normal_wgs_hg19(tumor_normal_wgs_config, snakemake_runner):
    snakefile = get_snakefile("paired", "balsamic")
    config_json = tumor_normal_wgs_config

    rc, out, err, cmd = snakemake_runner(
        snakefile=snakefile,
        configfile=config_json,
    )

    assert (
        rc == 0
    ), f"Snakemake failed: {' '.join(cmd)}\n--- STDOUT ---\n{out}\n--- STDERR ---\n{err}"
    # rule SHOULD be present
    assert "igh_dux4_detection_tumor_normal" in (out + err)


def test_workflow_qc_tumor_only_hg19(tumor_only_config_qc, snakemake_runner):
    snakefile = get_snakefile("single", AnalysisWorkflow.BALSAMIC_QC)
    config_json = tumor_only_config_qc

    rc, out, err, cmd = snakemake_runner(
        snakefile=snakefile,
        configfile=config_json,
    )

    assert (
        rc == 0
    ), f"Snakemake QC (tumor-only) failed: {' '.join(cmd)}\n--- STDOUT ---\n{out}\n--- STDERR ---\n{err}"


def test_workflow_qc_tumor_normal_hg19(tumor_normal_config_qc, snakemake_runner):
    snakefile = get_snakefile("paired", AnalysisWorkflow.BALSAMIC_QC)
    config_json = tumor_normal_config_qc

    rc, out, err, cmd = snakemake_runner(
        snakefile=snakefile,
        configfile=config_json,
    )

    assert (
        rc == 0
    ), f"Snakemake QC (tumor-normal) failed: {' '.join(cmd)}\n--- STDOUT ---\n{out}\n--- STDERR ---\n{err}"
