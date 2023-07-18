"""Custom class types."""
from enum import Enum


class StrEnum(str, Enum):
    """StrEnum data class."""

    def __str__(self) -> str:
        return str.__str__(self)
