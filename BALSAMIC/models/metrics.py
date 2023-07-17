"""QC validation metrics model."""
import logging
from typing import Optional, Any, List

from pydantic import BaseModel, validator

from BALSAMIC.constants.metrics import VALID_OPS

LOG = logging.getLogger(__name__)


class MetricCondition(BaseModel):
    """Defines the metric condition model.

    Attributes:
        norm (string, optional)     : Validation condition.
        threshold (float, optional) : Validation cut off.
    """

    norm: Optional[str] = None
    threshold: Optional[float] = None


class Metric(BaseModel):
    """Defines the metric attributes model.

    Attributes:
        header (str, optional)                : Data.
        id (str, required)                    : Unique sample identifier (sample_id, case_id or project_id).
        input (str, required)                 : Input file.
        name (str, required)                  : Metric name.
        step (str, required)                  : Step that generated the metric.
        value (Any, required)                 : Metric value.
        condition (MetricCondition, required) : Metric validation condition.
    """

    header: Optional[str]
    id: str
    input: str
    name: str
    step: str
    value: Any = ...
    condition: Optional[MetricCondition] = ...

    @validator("name")
    def validate_name(cls, name, values):
        """Updates the name if the source is FastQC."""
        if "fastqc-percent_duplicates" in name:
            return "PERCENT_DUPLICATION_R" + values["input"].split("_")[-2]
        return name


class MetricValidation(BaseModel):
    """Defines the metric validation model.

    Attributes:
        metrics (List[Metric], required) : Metric model to validate.
    """

    metrics: List[Metric]

    @validator("metrics", each_item=True)
    def validate_metrics(cls, metric):
        """Checks if a metric meets its filtering condition."""
        if metric.condition and not VALID_OPS[metric.condition.norm](
            metric.value, metric.condition.threshold
        ):
            raise ValueError(
                f"QC metric {metric.name}: {metric.value} validation has failed. "
                f"(Condition: {metric.condition.norm} {metric.condition.threshold}, ID: {metric.id})."
            )
        LOG.info(f"QC metric {metric.name}: {metric.value} meets its condition.")
        return metric
