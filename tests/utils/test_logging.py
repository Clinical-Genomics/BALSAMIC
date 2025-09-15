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


def test_add_file_logging_uses_existing_formatter(tmp_path: Path):
    """Test that add_file_logging copies the formatter from an existing handler."""

    # GIVEN logging instance
    log_file = tmp_path / "run.log"
    logger = logging.getLogger("formatted_logger")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    # GIVEN a StreamHandler with a custom formatter
    stream_handler = logging.StreamHandler()
    custom_formatter = logging.Formatter("CUSTOM FORMAT: %(message)s")
    stream_handler.setFormatter(custom_formatter)
    logger.addHandler(stream_handler)

    # WHEN adding file logging
    add_file_logging(str(log_file), logger_name="formatted_logger")

    # THEN the FileHandler that was added should be found
    file_handlers = [h for h in logger.handlers if isinstance(h, logging.FileHandler)]
    assert len(file_handlers) == 1

    file_handler = file_handlers[0]

    # THEN the formatter should be correctly copied to the file logger
    assert isinstance(file_handler.formatter, logging.Formatter)
    assert file_handler.formatter._fmt == "CUSTOM FORMAT: %(message)s"

    # THEN the messages written to logger ends up in file with the right format
    logger.info("formatter test")
    contents = log_file.read_text()
    assert "CUSTOM FORMAT: formatter test" in contents
