class BalsamicError(Exception):
    """Base exception for the BALSAMIC."""

    def __init__(self, message):
        super().__init__(message)


class WorkflowRunError(BalsamicError):
    """
    Exception for handling workflow errors.
    Raise this exception when workflow or rules fails to run or execute
    """
