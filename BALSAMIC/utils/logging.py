import logging

def add_file_logging(log_file: str, logger_name: str = None):
    """Adds a file handler to the specified logger without modifying its format.

    Args:
        log_file (str): Path to the log file.
        logger_name (str, optional): Name of the logger to configure. If None, uses the root logger.
    """
    logger = logging.getLogger(logger_name) if logger_name else logging.getLogger()

    # Avoid adding duplicate handlers
    if not any(isinstance(h, logging.FileHandler) for h in logger.handlers):
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)

        # Copy formatter from an existing handler if available
        if logger.handlers:
            file_handler.setFormatter(logger.handlers[0].formatter)

        logger.addHandler(file_handler)