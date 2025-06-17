"""QC validation metrics model."""
import logging
from typing import Optional, Any, List, Annotated, Callable

from pydantic import BaseModel, AfterValidator

from BALSAMIC.constants.metrics import VALID_OPS, METRIC_WARNINGS

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

def validate_metric(metric: Metric) -> Metric:
    """
    Checks if a metric meets its filtering condition.

    Raises:
    ValueError
        If the operator is unsupported, the operands are incompatible,
        or the metric fails validation and is **not** listed in
        ``METRIC_WARNINGS``
    """

    cond = metric.condition
    if cond is None:
        LOG.info("QC metric %s: %s (no condition).", metric.name, metric.value)
        return metric

    op: Optional[Callable[[Any, Any], bool]] = VALID_OPS.get(cond.norm)
    if op is None:
        raise ValueError(f"Unsupported operation: {cond.norm!r}")

    try:
        passed = op(metric.value, cond.threshold)
    except TypeError as exc:
        raise ValueError(
            f"Type mismatch for QC metric {metric.name}: "
            f"{metric.value!r} {cond.norm} {cond.threshold!r} (ID: {metric.id})"
        ) from exc

    if passed:
        LOG.info(
            "QC metric %s: %s meets its condition (%s %s, ID: %s).",
            metric.name, metric.value, cond.norm, cond.threshold, metric.id
        )
        return metric

    # ── Validation failed ──────────────────────────────────────────────
    msg = (
        f"QC metric {metric.name}: {metric.value} failed "
        f"(condition: {cond.norm} {cond.threshold}, ID: {metric.id})."
    )
    if metric.name in METRIC_WARNINGS:
        LOG.info(msg)
    else:
        raise ValueError(msg)

    return metric

class MetricValidation(BaseModel):
    """Defines the metric validation model.

    Attributes:
        metrics (List[Metric], required) : Metric model to validate.
    """

    metrics: List[Annotated[Metric, AfterValidator(validate_metric)]]
