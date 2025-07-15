import logging
from BALSAMIC.constants.analysis import LogFile
import os
import re


def set_log_filename(analysis_dir: str, run_start: bool = False):
    """
    Determine the appropriate logfile name in a given analysis directory based on existing log versions.

    This function inspects the given directory for files following the naming pattern:
    - `balsamic.log`
    - `balsamic.log.1`, `balsamic.log.2`, etc.

    Args:
        analysis_dir (str): Path to the directory where log files are stored.
        run_start (bool): If False (default), return the latest existing logfile name.
                          If True, return the next incremented logfile name suitable for a new log.

    Returns:
        str: The name of the logfile (either the latest existing or the next version, based on `run_start`).

    Notes:
        - `LogFile.LOGNAME` is used as the base name (e.g., 'balsamic.log').
        - The version-less logfile (e.g., `balsamic.log`) is treated as version 0.
    """
    latest_version = -1

    # Pattern: match balsamic.log or balsamic.log.1, .2, etc.
    pattern = re.compile(rf"^{re.escape(LogFile.LOGNAME)}(?:\.(\d+))?$")

    for filename in os.listdir(analysis_dir):
        match = pattern.match(filename)
        if match:
            version = int(match.group(1)) if match.group(1) else 0
            if version > latest_version:
                latest_version = version

    if run_start:
        # Suggest the next version name
        if latest_version == 0 and not os.path.exists(
            os.path.join(analysis_dir, LogFile.LOGNAME)
        ):
            return LogFile.LOGNAME
        else:
            return f"{LogFile.LOGNAME}.{latest_version + 1}"
    else:
        # Return the latest existing logfile name
        return (
            LogFile.LOGNAME
            if latest_version == 0
            else f"{LogFile.LOGNAME}.{latest_version}"
        )


def add_file_logging(log_file: str, logger_name: str = None):
    """Adds a file handler to the specified logger without modifying its format.

    Args:
        log_file (str): Path to the log file.
        logger_name (str, optional): Name of the logger to configure. If None, uses the root logger.
    """
    logger = logging.getLogger(logger_name) if logger_name else logging.getLogger()

    # Avoid adding duplicate file handlers
    if any(isinstance(h, logging.FileHandler) for h in logger.handlers):
        return

    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.INFO)

    # Try to copy the formatter from an existing handler
    for handler in logger.handlers:
        if handler.formatter:
            file_handler.setFormatter(handler.formatter)
            break
    else:
        # Fallback to a default formatter if no other handler has a formatter
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
        )

    logger.addHandler(file_handler)
