import logging
from BALSAMIC.constants.analysis import LogFile
import os
import re
from pathlib import Path


def set_log_filename(analysis_dir: str, run_start: bool = False, config_case=False):
    """
    Determine the appropriate logfile name in a given analysis directory based on existing log versions.

    This function inspects the given directory for files following the naming pattern:
    - `balsamic.run.log`, `balsamic.run.log.1`, `balsamic.run.log.2`, ...
    - or, if `config_case=True`: `balsamic.config.log`, `balsamic.config.log.1`, ...

    Args:
        analysis_dir (str): Path to the directory where log files are stored.
        run_start (bool): If False (default), return the latest existing logfile name.
                          If True, return the next incremented logfile name suitable for a new log.
        config_case (bool): If True, use the config logfile base name (`balsamic.config.log`)
                            instead of the run logfile name.

    Returns:
        str: The name of the logfile (either the latest existing or the next version, based on `run_start`).

    Notes:
        - `LogFile.RUN_LOGNAME` and `LogFile.CONFIG_LOGNAME` are used as base names.
        - The version-less logfile (e.g., `balsamic.run.log`) is treated as version 0.
    """
    if config_case:
        logname = LogFile.CONFIG_LOGNAME
    else:
        logname = LogFile.RUN_LOGNAME

    if not Path(analysis_dir, logname).is_file():
        return Path(analysis_dir, logname).as_posix()

    # Pattern: match balsamic.log or balsamic.log.1, .2, etc.
    pattern = re.compile(rf"^{re.escape(logname)}(?:\.(\d+))?$")

    latest_version = -1

    for filename in os.listdir(analysis_dir):
        match = pattern.match(filename)
        if match:
            version = int(match.group(1)) if match.group(1) else 0
            if version > latest_version:
                latest_version = version

    if run_start:
        return Path(f"{analysis_dir}/{logname}.{latest_version + 1}").as_posix()
    else:
        # Return the latest existing logfile name
        logfile = logname if latest_version == 0 else f"{logname}.{latest_version}"
        return Path(f"{analysis_dir}/{logfile}").as_posix()


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
