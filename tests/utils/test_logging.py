"""Test Balsamic logging utility methods."""
import logging
import coloredlogs

from pathlib import Path

import pytest

from BALSAMIC.utils.logging import add_file_logging
from BALSAMIC.constants.constants import LogLevel


def test_writes_messages_to_file(tmp_path: Path, caplog: pytest.LogCaptureFixture):
    """The helper should attach a FileHandler that receives log records."""

    # GIVEN logging instance
    LOG = logging.getLogger(__name__)

    coloredlogs.DEFAULT_FIELD_STYLES = {
        "asctime": {"color": "green"},
        "hostname": {"color": "magenta"},
        "levelname": {"color": "yellow", "bold": True},
        "programname": {"color": "cyan"},
        "name": {"color": "blue"},
    }
    coloredlogs.install(
        level=LogLevel.INFO,
        fmt="%(programname)s %(hostname)s %(asctime)s %(name)s pid:%(process)d [%(levelname)s] %(message)s",
    )

    log_file = tmp_path / "run.log"

    # GIVEN logfile path to store log messages
    add_file_logging(str(log_file), logger_name=__name__)

    # WHEN printing log information
    LOG.info("hello world")

    # THEN logfile should be created
    assert log_file.exists()

    # THEN message should exist inside the logfile
    assert "hello world" in log_file.read_text()

    # WHEN printing log debug info
    LOG.debug("invisible")

    # THEN message should not be stored in logfile
    assert "invisible" not in log_file.read_text()


def test_add_file_logging_does_not_duplicate_handlers(tmp_path: Path):
    """Test that logger function correctly only adds one filehandler in case its accidentally added twice."""

    # GIVEN logging instance
    log_file = tmp_path / "run.log"

    logger = logging.getLogger("test_logger")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    # WHEN adding filehandler
    add_file_logging(log_file, logger_name="test_logger")

    # THEN 1 file handler should exist
    assert sum(isinstance(h, logging.FileHandler) for h in logger.handlers) == 1

    # WHEN adding a second filderhandler
    add_file_logging(log_file, logger_name="test_logger")

    # THEN 1 file handler should exist
    assert sum(isinstance(h, logging.FileHandler) for h in logger.handlers) == 1
