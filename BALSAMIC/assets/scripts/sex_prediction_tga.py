import click
import pandas as pd
import os
import warnings
from typing import List, Optional, Dict, Any, Tuple

from BALSAMIC.utils.io import write_json


def process_coverage_data(filepath: str) -> Optional[pd.DataFrame]:
    """
    Process coverage data to keep only information about sex chromosomes and returns a DataFrame.

    Args:
        filepath (str): Path to tsv coverage file

    Returns:
        Optional[pd.DataFrame]: Processed DataFrame or None if insufficient data.
    """

    colnames: List[str] = ["chromosome", "start", "end", "gene", "depth", "log2"]
    data: pd.DataFrame = pd.read_csv(filepath, sep="\t", header=None, names=colnames)
    filtered_data = data[data["chromosome"].isin(["X", "Y"])].copy()

    if filtered_data.shape[0] < 10:
        warnings.warn(
            "Cannot compute sex due to insufficient data, ignoring", UserWarning
        )
        return None

    filtered_data["key"] = (
        filtered_data["chromosome"]
        + "_"
        + filtered_data["start"].astype(str)
        + "_"
        + filtered_data["end"].astype(str)
    )

    filtered_data = filtered_data[["key", "chromosome", "depth"]].set_index("key")
    filtered_data["depth"] = pd.to_numeric(filtered_data["depth"], errors="coerce")
    return filtered_data


def get_stats(df: pd.DataFrame) -> dict[str, float]:
    """Compute statistical summaries of depth grouped by chromosome.

    Args:
        df (pd.DataFrame): DataFrame containing coverage data.

    Returns:
        Dict[str, float]: Dictionary containing extracted statistics such as mean, median, min, max, etc.
    """
    stats = (
        df.groupby("chromosome")["depth"]
        .agg(
            mean="mean", median="median", min="min", max="max", std="std", count="count"
        )
        .reset_index()
    )
    # Create the stats dictionary
    stats_dict = {}
    for index, row in stats.iterrows():
        chromosome = row["chromosome"]
        for metric in ["mean", "median", "min", "max", "std", "count"]:
            stats_dict[f"{chromosome}_{metric}"] = row[metric]
    return stats_dict


def retrieve_file_info(filepath: str) -> (str, str):
    """Extract sample name and CNN type from file name."""
    filename: str = os.path.basename(filepath)
    sample_name: str = filename.split(".")[0]
    cnn_type: str = "antitarget" if "antitarget" in filename else "target"
    return sample_name, cnn_type


def get_predicted_sex(y_x_frac: float) -> Dict[str, str]:
    """Predict sex and confidence based on the Y/X ratio.

    Args:
        y_x_frac (float): Ratio of Y to X.

    Returns:
        Dict[str, str]: Predicted sex and confidence level.
    """
    if y_x_frac > 0.2:
        if y_x_frac >= 0.9:
            return {"predicted_sex": "male", "confidence": "high"}
        elif y_x_frac >= 0.7:
            return {"predicted_sex": "male", "confidence": "medium"}
        else:
            return {"predicted_sex": "male", "confidence": "low"}
    elif y_x_frac <= 0.15:
        if y_x_frac >= 0.05:
            return {"predicted_sex": "female", "confidence": "medium"}
        elif y_x_frac >= 0.00001:
            return {"predicted_sex": "female", "confidence": "high"}
        else:
            return {"predicted_sex": "female", "confidence": "low"}
    else:
        return {"predicted_sex": "unknown", "confidence": "low"}


def predict_sex(cnn_file: str) -> Dict[str, Any]:
    """Predict the sex of a sample based on CNN file data.

    Args:
        cnn_file (str): Path to the CNN file containing coverage data.

    Returns:
        Dict[str, Any]: A dictionary containing sample information, data statistics,
                        and predicted sex based on mean and median Y/X coverage ratios.
    """
    sample_name, cnn_type = retrieve_file_info(cnn_file)

    predicted_sex = {
        "sample_name": sample_name,
        "cnn_type": cnn_type,
    }

    data_df: pd.DataFrame = process_coverage_data(cnn_file)

    if not isinstance(data_df, pd.DataFrame):
        predicted_sex["sex_prediction"] = "failed"
        predicted_sex["data_dict"] = "NA"
        return predicted_sex

    data_dict: dict = get_stats(data_df)

    # Categorize data amount based on X and Y chromosome counts
    if data_dict["X_count"] < 10 or data_dict["Y_count"] < 10:
        data_dict["data_amount"] = "low"
    elif data_dict["X_count"] < 50 or data_dict["Y_count"] < 50:
        data_dict["data_amount"] = "medium"
    else:
        data_dict["data_amount"] = "high"

    # Compute mean and median ratios
    data_dict["Y_mean/X_mean"] = round(data_dict["Y_mean"] / data_dict["X_mean"], 5)
    data_dict["Y_median/X_median"] = round(
        data_dict["Y_median"] / data_dict["X_median"], 5
    )

    predicted_sex["sex_prediction"] = {
        "by_mean": get_predicted_sex(data_dict["Y_mean/X_mean"]),
        "by_median": get_predicted_sex(data_dict["Y_median/X_median"]),
    }

    predicted_sex["data_dict"] = data_dict

    return predicted_sex


def calculate_prediction_score(sex_prediction: dict) -> int:
    """Convert sex prediction confidence strings into scores.

    Args:
        sex_prediction (dict): Sex prediction dictionary

    Returns:
       int: Arbitrary score for sex-prediction confidence.
    """
    frac_conf_score = {
        "high": 4,
        "medium": 3,
        "low": 1,
        "NA": 0,
    }
    data_conf_score = {
        "high": 6,
        "medium": 5,
        "low": 1,
        "NA": 0,
    }
    data_type_score = {
        "target": 1,
        "antitarget": 0,
    }
    frac_conf: str = sex_prediction["frac_conf"]
    data_conf: str = sex_prediction["data_conf"]
    data_type: str = sex_prediction["data_type"]
    score = (
        frac_conf_score[frac_conf]
        + data_conf_score[data_conf]
        + data_type_score[data_type]
    )
    return score


def get_prediction(prediction: Dict[str, Any]) -> Tuple[str, str, str, str]:
    """Extract sex predictions and confidence levels from a prediction dictionary.

    Args:
        prediction (Dict[str, Any]): A dictionary containing sex prediction details.

    Returns:
        Tuple[str, str, str, str]: A tuple containing:
            - Predicted sex by mean ratio (or "unknown" if prediction failed).
            - Confidence level for the mean ratio prediction (or "NA" if prediction failed).
            - Predicted sex by median ratio (or "unknown" if prediction failed).
            - Confidence level for the median ratio prediction (or "NA" if prediction failed).
    """
    if "failed" in prediction["sex_prediction"]:
        return "unknown", "NA", "unknown", "NA"

    return (
        prediction["sex_prediction"]["by_mean"]["predicted_sex"],
        prediction["sex_prediction"]["by_mean"]["confidence"],
        prediction["sex_prediction"]["by_median"]["predicted_sex"],
        prediction["sex_prediction"]["by_median"]["confidence"],
    )


def consolidate_sex_predictions(
    sex_predictions: List[Dict[str, Any]]
) -> Dict[str, Any]:
    """Consolidate multiple sex predictions into a final prediction.

    Args:
        sex_predictions (List[Dict[str, Any]]): A list of dictionaries containing sex predictions.
            Each dictionary should have:
                - "sex" (str): Predicted sex.
                - "score" (int): Score for the prediction.
                - "data_type" (str): Type of data used for the prediction (e.g., "target").

    Returns:
        Dict[str, Any]: A dictionary containing:
            - "sex" (str): The consolidated predicted sex.
            - "sex_score" (int): The highest score among the predictions.
            - "prediction_confidence" (str): Confidence level based on the score ("low", "medium", "high", or "uncertain").
    """
    final_sex = "unknown"
    final_score = 0

    for sex_prediction in sex_predictions:
        sex = sex_prediction["sex"]
        score = sex_prediction["score"]
        data_type = sex_prediction["data_type"]

        # Set sex based on highest score
        if score > final_score:
            final_sex = sex
            final_score = score
        elif score == final_score:
            # If score is the same, prioritize target data
            if data_type == "target":
                final_sex = sex

    score_confidence_levels = {
        "high": 9,
        "medium": 8,
        "low": 4,
    }

    final_confidence = "uncertain"
    if final_score > score_confidence_levels["low"]:
        final_confidence = "low"
    if final_score > score_confidence_levels["medium"]:
        final_confidence = "medium"
    if final_score > score_confidence_levels["high"]:
        final_confidence = "high"

    return {
        "predicted_sex": final_sex,
        "predicted_sex_score": final_score,
        "prediction_confidence": final_confidence,
    }


def summarise_sample_sex_prediction(
    target_predicted_sex: Dict[str, Any],
    antitarget_predicted_sex: Dict[str, Any],
    sample_type: str,
) -> Dict[str, Any]:
    """Summarize sex predictions for a sample based on target and antitarget data.

    Args:
        target_predicted_sex (Dict[str, Any]): Predicted sex information from target data.
        antitarget_predicted_sex (Dict[str, Any]): Predicted sex information from antitarget data.
        sample_type (str): Type of the sample (e.g., "tumor" or "normal").

    Returns:
        Dict[str, Any]: A dictionary summarizing the sex predictions, including individual
                        predictions, confidence scores, and final consolidated prediction.
    """
    sample_predicted_sex = {sample_type: {}}

    if "failed" in target_predicted_sex["sex_prediction"]:
        target_data_amount = "NA"
    else:
        target_data_amount = target_predicted_sex["data_dict"]["data_amount"]

    (
        mean_target_predicted_sex,
        mean_target_predicted_sex_conf,
        median_target_predicted_sex,
        median_target_predicted_sex_conf,
    ) = get_prediction(target_predicted_sex)

    if "failed" in antitarget_predicted_sex["sex_prediction"]:
        antitarget_data_amount = "NA"
    else:
        antitarget_data_amount = antitarget_predicted_sex["data_dict"]["data_amount"]

    (
        mean_antitarget_predicted_sex,
        mean_antitarget_predicted_sex_conf,
        median_antitarget_predicted_sex,
        median_antitarget_predicted_sex_conf,
    ) = get_prediction(antitarget_predicted_sex)

    sample_predicted_sex[sample_type]["predicted_sex"] = [
        {
            "sex": mean_target_predicted_sex,
            "frac_conf": mean_target_predicted_sex_conf,
            "data_conf": target_data_amount,
            "data_type": "target",
        },
        {
            "sex": median_target_predicted_sex,
            "frac_conf": median_target_predicted_sex_conf,
            "data_conf": target_data_amount,
            "data_type": "target",
        },
        {
            "sex": mean_antitarget_predicted_sex,
            "frac_conf": mean_antitarget_predicted_sex_conf,
            "data_conf": antitarget_data_amount,
            "data_type": "antitarget",
        },
        {
            "sex": median_antitarget_predicted_sex,
            "frac_conf": median_antitarget_predicted_sex_conf,
            "data_conf": antitarget_data_amount,
            "data_type": "antitarget",
        },
    ]

    for idx, sex_prediction in enumerate(
        sample_predicted_sex[sample_type]["predicted_sex"]
    ):
        # Convert to score from sex prediction confidence strings
        score = calculate_prediction_score(sex_prediction)
        sample_predicted_sex[sample_type]["predicted_sex"][idx].update({"score": score})

    sample_predicted_sex[sample_type] = consolidate_sex_predictions(
        sample_predicted_sex[sample_type]["predicted_sex"]
    )

    # Add raw data
    sample_predicted_sex[sample_type]["target_predicted_sex"] = target_predicted_sex
    sample_predicted_sex[sample_type][
        "antitarget_predicted_sex"
    ] = antitarget_predicted_sex

    return sample_predicted_sex


def case_sex_prediction(predicted_sex: Dict) -> Dict:
    """Add case level sex prediction, setting sex to conflicting if tumor and normal doesn't match."""
    predicted_sex["case_sex"] = {}
    tumor_sex: str = predicted_sex["tumor"]["predicted_sex"]
    if "normal" in predicted_sex:
        normal_sex: str = predicted_sex["normal"]["predicted_sex"]
        if tumor_sex != normal_sex:
            case_sex = "conflicting"
            warnings.warn(
                f"Conflicting sex prediction found, tumor sex: {tumor_sex}, is not the same as normal sex: {normal_sex}",
                UserWarning,
            )
        else:
            case_sex: str = tumor_sex
    else:
        case_sex = tumor_sex

    predicted_sex["case_sex"] = case_sex
    return predicted_sex


@click.command()
@click.option(
    "--target-cnn-tumor",
    type=click.Path(exists=True),
    required=True,
    help="Path to the target CNN tumor file.",
)
@click.option(
    "--antitarget-cnn-tumor",
    type=click.Path(exists=True),
    required=True,
    help="Path to the antitarget CNN tumor file.",
)
@click.option(
    "--output",
    type=click.Path(writable=True),
    required=True,
    help="Path to the output file to be created.",
)
@click.option(
    "--target-cnn-normal",
    type=click.Path(exists=True),
    default=None,
    help="Optional path to the target CNN normal file.",
)
@click.option(
    "--antitarget-cnn-normal",
    type=click.Path(exists=True),
    default=None,
    help="Optional path to the antitarget CNN normal file.",
)
def predict_sex_main(
    target_cnn_tumor: str,
    antitarget_cnn_tumor: str,
    output: str,
    target_cnn_normal: Optional[str],
    antitarget_cnn_normal: Optional[str],
) -> None:
    """Process tumor and optional normal CNN files to predict sex and output JSON results.

    Args:
        target_cnn_tumor (str): Path to the target CNN tumor file.
        antitarget_cnn_tumor (str): Path to the antitarget CNN tumor file.
        output (str): Path to the output JSON file.
        target_cnn_normal (Optional[str]): Path to the target CNN normal file.
        antitarget_cnn_normal (Optional[str]): Path to the antitarget CNN normal file.
    """
    # Predict tumor sex based on CNVkit output
    tumor_target_predicted_sex: dict = predict_sex(target_cnn_tumor)
    tumor_antitarget_predicted_sex: dict = predict_sex(antitarget_cnn_tumor)
    predicted_sex: dict = summarise_sample_sex_prediction(
        tumor_target_predicted_sex, tumor_antitarget_predicted_sex, "tumor"
    )

    if target_cnn_normal:
        # Predict normal sex based on CNVkit output
        normal_target_predicted_sex: dict = predict_sex(target_cnn_normal)
        normal_antitarget_predicted_sex: dict = predict_sex(antitarget_cnn_normal)
        normal_predicted_sex: dict = summarise_sample_sex_prediction(
            normal_target_predicted_sex, normal_antitarget_predicted_sex, "normal"
        )
        predicted_sex.update(normal_predicted_sex)

    # Create case-level prediction (compare tumor and normal sex)
    predicted_sex: dict = case_sex_prediction(predicted_sex)

    # Write json report
    write_json(predicted_sex, output)


if __name__ == "__main__":
    predict_sex_main()
