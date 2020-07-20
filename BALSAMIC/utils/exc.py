class BalsamicError(Exception):
    """Base exception for the BALSAMIC."""

    def __init__(self, message):
        super().__init__(message)
