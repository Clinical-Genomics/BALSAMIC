"""QC validation metrics model."""
import logging
from typing import Optional, Any, List, Annotated

from pydantic import BaseModel, AfterValidator

from BALSAMIC.constants.metrics import VALID_OPS

LOG = logging.getLogger(__name__)


class MetricCondition(BaseModel):
    """Defines the metric condition model.

    Attributes:
        norm (str, optional)       : Validation condition (e.g., "eq", "lt").
        threshold (Any, optional)  : Validation cutoff or expected value.
    """

    norm: Optional[str] = None
    threshold: Optional[Any] = None


class Metric(BaseModel):
    """Defines the metric attributes model.

    Attributes:
        header (str, optional)                : Data.
        id (str, required)                    : Unique sample identifier (sample_id, case_id or project_id).
        input (str, required)                 : Input file.
        name (str, required)                  : Metric name.
        step (str, required)                  : Step that generated the metric.
        value (Any, required)                 : Metric value (can be float or str).
        condition (MetricCondition, required) : Metric validation condition.
    """

    header: Optional[str] = None
    id: str
    input: str
    name: str
    step: str
    value: Any
    condition: Optional[MetricCondition]


def validate_metric(metric: Metric):
    """Checks if a metric meets its filtering condition."""
    if metric.condition:
        norm: Optional[str] = metric.condition.norm
        threshold: Optional[Any] = metric.condition.threshold
        value: Any = metric.value

        # Validate the norm operator
        if norm not in VALID_OPS:
            raise ValueError(f"Unsupported operation: {norm}")

        # Attempt validation using the operator
        try:
            if not VALID_OPS[norm](value, threshold):
                raise ValueError(
                    f"QC metric {metric.name}: {value} validation has failed. "
                    f"(Condition: {norm} {threshold}, ID: {metric.id})."
                )
        except TypeError:
            raise ValueError(
                f"Type mismatch in QC metric {metric.name}: {value} and {threshold} "
                f"are not compatible with operator {norm}. (ID: {metric.id})."
            )

    LOG.info(f"QC metric {metric.name}: {metric.value} meets its condition.")
    return metric


class MetricValidation(BaseModel):
    """Defines the metric validation model.

    Attributes:
        metrics (List[Metric], required) : Metric model to validate.
    """

    metrics: List[Annotated[Metric, AfterValidator(validate_metric)]]
