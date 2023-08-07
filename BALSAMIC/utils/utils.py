"""Helper functions."""


def remove_unnecessary_spaces(string: str) -> str:
    """Return a string removing unnecessary empty spaces."""
    return " ".join(string.split())
