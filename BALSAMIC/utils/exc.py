class BalsamicError(Exception):
    """Base exception for the BALSAMIC."""

    def __init__(self, message):
        super(BalsamicError, self).__init__()
        self.message = message
