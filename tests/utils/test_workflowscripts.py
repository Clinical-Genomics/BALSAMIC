from unittest import mock
from pathlib import Path
import pytest

from BALSAMIC.utils.cli import generate_h5
from BALSAMIC.utils.workflowscripts import get_file_contents
from BALSAMIC.utils.workflowscripts import plot_analysis

def test_get_file_contents():
    # GIVEN a test input file
    test_file = 'tests/test_data/densityplots/dummy_file1.txt'

    # WHEN invoking function
    test_file_built = get_file_contents(test_file, 'umi')
    column_names = ['id', 'AF', 'method']

    # THEN check column names and no. of column matches
    assert all(test_file_built.columns == column_names)
    assert len(test_file_built.columns) == 3


def test_get_wrongfile_contents():
    # GIVEN a test input file
    test_wrongfile = 'tests/test_data/densityplots/dummy_wrongfile.txt'

    # WHEN invoking function
    with pytest.raises(ValueError):
        test_wrongfile_built = get_file_contents(test_wrongfile, 'umi')
        assert len(test_wrongfile_built.columns) != 3


def test_plot_analysis(tmp_path_factory):
    # GIVEN a dummy log file
    dummy_log_file = Path(
        "tests/test_data/dummy_run_logs/BALSAMIC.T_panel.bwa_mem.123.sh_31415926535.err"
    )
    dummy_h5 = "tests/test_data/dummy_run_logs/BALSAMIC.T_panel.bwa_mem.123.h5"
    dummy_path = tmp_path_factory.mktemp("dummy_pdf_path")
    dummy_pdf_name = dummy_path / "BALSAMIC.T_panel.bwa_mem.123.pdf"
    dummy_pdf_name.touch()

    # WHEN calling plot_analysis
    actual_pdf_file = plot_analysis(dummy_log_file, dummy_h5, dummy_pdf_name)

    assert Path(actual_pdf_file).exists()


def test_plot_analysis_bad_h5(tmp_path_factory):
    # GIVEN a dummy log file
    dummy_log_file = Path(
        "tests/test_data/dummy_run_logs/BALSAMIC.T_panel.bwa_mem.123.sh_31415926535.err"
    )
    dummy_h5 = "tests/test_data/dummy_run_logs/bad_format.h5"
    dummy_path = tmp_path_factory.mktemp("dummy_pdf_path")
    dummy_pdf_name = dummy_path / "plot_file.pdf"
    dummy_pdf_name.touch()

    # WHEN calling plot_analysis
    actual_pdf_file = plot_analysis(dummy_log_file, dummy_h5, dummy_pdf_name)

    assert actual_pdf_file is None
