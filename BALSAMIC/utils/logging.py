import logging


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
