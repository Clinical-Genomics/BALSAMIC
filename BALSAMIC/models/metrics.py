"""QC validation metrics model."""
import logging
from typing import Optional, Any, List

from pydantic import BaseModel, validator

from BALSAMIC.constants.metrics import VALID_OPS

LOG = logging.getLogger(__name__)


class MetricConditionModel(BaseModel):
    """Defines the metric condition model

    Attributes:
        norm: string (optional); validation condition
        threshold: float (optional); validation cut off
    """

    norm: Optional[str] = None
    threshold: Optional[float] = None


class MetricModel(BaseModel):
    """Defines the metric attributes model

    Attributes:
        header: str (optional); data
        id: str (required); unique sample identifier (sample_id, case_id or project_id)
        input: str (required); input file
        name: str (required); metric name
        step: str (required); step that generated the metric
        value: Any (required and can take None as a value); metric value
        condition: MetricConditionModel (required and can take None as a value); metric validation condition
    """

    header: Optional[str]
    id: str
    input: str
    name: str
    step: str
    value: Any = ...
    condition: Optional[MetricConditionModel] = ...

    @validator("name")
    def validate_name(cls, name, values):
        """Updates the name if the source is FastQC"""

        if "fastqc-percent_duplicates" in name:
            return "PERCENT_DUPLICATION_R" + values["input"].split("_")[-2]

        return name


class MetricValidationModel(BaseModel):
    """Defines the metric validation model

    Attributes:
        metrics: List[MetricModel] (required); metric model to validate

    Raises:
        ValueError: when a metric does not meet its validation requirements
    """

    metrics: List[MetricModel]

    @validator("metrics", each_item=True)
    def validate_metrics(cls, metric):
        """Checks if a metric meets its filtering condition"""

        if metric.condition and not VALID_OPS[metric.condition.norm](
            metric.value, metric.condition.threshold
        ):
            raise ValueError(
                f"QC metric {metric.name}: {metric.value} validation has failed. "
                f"(Condition: {metric.condition.norm} {metric.condition.threshold}, ID: {metric.id})."
            )

        LOG.info(f"QC metric {metric.name}: {metric.value} meets its condition.")

        return metric
