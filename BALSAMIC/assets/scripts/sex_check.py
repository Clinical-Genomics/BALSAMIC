import click
import pandas as pd
import json
import os




def read_cov(filepath):
    with open(filepath, "r") as rf:
        rows = rf.readlines()
        data_rows = [r.strip("\n").split("\t") for r in rows]
    return data_rows

def process_data(data_rows):
    filtered_data = [sublist for sublist in data_rows if sublist[0] in {"Y", "X"}]
    if len(filtered_data) < 10:
        print("Cannot compute sex, ignoring")
        return None
    dicty = {}
    for l in filtered_data:
        c = l[0]
        p = l[1]
        s = l[2]
        d = l[4]
        n = f"{c}_{p}_{s}"
        dicty[n] = {}
        dicty[n]["chrom"] = c
        dicty[n]["depth"] = d
    df = pd.DataFrame.from_dict(dicty).T
    df["depth"] = pd.to_numeric(df["depth"], errors="coerce")
    return df

def get_stats(df):
    stats = df.groupby("chrom")["depth"].agg(
        mean="mean",
        median="median",
        min="min",
        max="max",
        std="std",
        count="count"
    ).reset_index()
    return stats

def extract_stats(stats):
    stats = stats.T.to_dict()
    sample_dict = {}
    sample_dict["X_mean"] = stats[0]["mean"]
    sample_dict["X_median"] = stats[0]["median"]
    sample_dict["X_min"] = stats[0]["min"]
    sample_dict["X_max"] = stats[0]["min"]
    sample_dict["X_std"] = stats[0]["std"]
    sample_dict["X_count"] = stats[0]["count"]
    sample_dict["Y_mean"] = stats[1]["mean"]
    sample_dict["Y_median"] = stats[1]["median"]
    sample_dict["Y_min"] = stats[1]["min"]
    sample_dict["Y_max"] = stats[1]["max"]
    sample_dict["Y_std"] = stats[1]["std"]
    sample_dict["Y_count"] = stats[1]["count"]
    return sample_dict

def retrieve_file_info(cnn_file):
    filename = os.path.basename(cnn_file)
    sample_name = filename.split(".")[0]
    if "antitarget" in filename:
        cnn_type = "antitarget"
    else:
        cnn_type = "target"
    return sample_name, cnn_type

def get_predicted_sex(y_x_frac):
    sex_prediction = {}

    if y_x_frac > 0.2:
        # Above 0.2
        sex_prediction["predicted_sex"] = "male"
        if y_x_frac >= 0.5:
            # above 0.5
            sex_prediction["predicted_sex_confidence"] = "high"
        elif y_x_frac >= 0.3 and y_x_frac < 0.5:
            # between 0.3 and 0.5
            sex_prediction["predicted_sex_confidence"] = "medium"
        else:
            # between 0.2 and 0.3
            sex_prediction["predicted_sex_confidence"] = "low"
    elif y_x_frac > 0.1 and y_x_frac < 0.2:
        # Between 0.1 and 0.2
        sex_prediction["predicted_sex"] = "unknown"
        sex_prediction["predicted_sex_confidence"] = "low"
    else:
        # Below 0.1
        sex_prediction["predicted_sex"] = "female"
        if y_x_frac >= 0.08:
            # Between 0.08 and 0.1
            sex_prediction["predicted_sex_confidence"] = "low"
        elif y_x_frac >= 0.05:
            # Between 0.05 and 0.08
            sex_prediction["predicted_sex_confidence"] = "medium"
        else:
            # Between 0 and 0.05
            sex_prediction["predicted_sex_confidence"] = "high"

    return sex_prediction

def predict_sex(cnn_file):
    sample_name, cnn_type = retrieve_file_info(cnn_file)

    predicted_sex = {}
    predicted_sex["sample_name"] = sample_name
    predicted_sex["cnn_type"] = cnn_type

    data_rows = read_cov(cnn_file)
    data_df = process_data(data_rows)
    if not isinstance(data_df, pd.DataFrame):
        predicted_sex["sex_prediction"] = "failed"
        predicted_sex["data_dict"] = "NA"
        return predicted_sex

    stats = get_stats(data_df)
    data_dict = extract_stats(stats)


    if data_dict["X_count"] < 10 or data_dict["Y_count"] < 10:
        data_dict["data_amount"] = "low"
    elif data_dict["X_count"] < 50 or data_dict["Y_count"] < 50:
        data_dict["data_amount"] = "medium"
    else:
        data_dict["data_amount"] = "high"

    data_dict["Y_mean/X_mean"] = round(data_dict["Y_mean"] / data_dict["X_mean"], 5)
    data_dict["Y_median/X_median"] = round(data_dict["Y_median"] / data_dict["X_median"], 5)

    predicted_sex["sex_prediction"] = {}
    predicted_sex["sex_prediction"]["by_mean"] = get_predicted_sex(data_dict["Y_mean/X_mean"])
    predicted_sex["sex_prediction"]["by_median"] = get_predicted_sex(data_dict["Y_median/X_median"])

    predicted_sex["data_dict"] = data_dict

    return predicted_sex

def get_prediction(prediction):
    if "failed" in prediction["sex_prediction"]:
        return "NA", "NA", "NA", "NA"
    else:
        return prediction["sex_prediction"]["by_mean"]["predicted_sex"], prediction["sex_prediction"]["by_mean"]["predicted_sex_confidence"], prediction["sex_prediction"]["by_median"]["predicted_sex"], prediction["sex_prediction"]["by_median"]["predicted_sex_confidence"]


def calculate_prediction_score(sex_prediction):
    frac_conf_score = {
        "high": 3,
        "medium": 2,
        "low": 1,
        "NA": 0,
    }
    data_conf_score = {
        "high": 5,
        "medium": 4,
        "low": 1,
        "NA": 0,
    }
    data_type_score = {
        "target": 1,
        "antitarget": 0,
    }
    frac_conf = sex_prediction["frac_conf"]
    data_conf = sex_prediction["data_conf"]
    data_type = sex_prediction["data_type"]
    score = frac_conf_score[frac_conf] + data_conf_score[data_conf] + data_type_score[data_type]
    return score

def consolidate_sex_predictions(sex_predictions):

    final_sex = "NA"
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
            # If score is the same, prioritise target data
            if data_type == "target":
                final_sex = sex

    score_confidence_levels = {
        "high": 6,
        "medium": 5,
        "low": 3,
    }
    final_confidence = "uncertain"
    if final_score > score_confidence_levels["low"]:
        final_confidence = "low"
    if final_score > score_confidence_levels["medium"]:
        final_confidence = "medium"
    if final_score > score_confidence_levels["high"]:
        final_confidence = "high"

    return {"sex": final_sex, "sex_score": score, "prediction_confidence": final_confidence}

def summarise_sample_sex_prediction(target_predicted_sex, antitarget_predicted_sex, sample_type):
    sample_predicted_sex = {sample_type: {}}

    if "failed" in target_predicted_sex["sex_prediction"]:
        target_data_amount = "NA"
    else:
        target_data_amount = target_predicted_sex["data_dict"]["data_amount"]

    mean_target_predicted_sex, mean_target_predicted_sex_conf, median_target_predicted_sex,  median_target_predicted_sex_conf = get_prediction(target_predicted_sex)

    if "failed" in antitarget_predicted_sex["sex_prediction"]:
        antitarget_data_amount = "NA"
    else:
        antitarget_data_amount = antitarget_predicted_sex["data_dict"]["data_amount"]

    mean_antitarget_predicted_sex, mean_antitarget_predicted_sex_conf, median_antitarget_predicted_sex, median_antitarget_predicted_sex_conf = get_prediction(antitarget_predicted_sex)

    sample_predicted_sex[sample_type]["predicted_sex"] = [
        {"sex": mean_target_predicted_sex, "frac_conf": mean_target_predicted_sex_conf, "data_conf": target_data_amount, "data_type": "target"},
        {"sex": median_target_predicted_sex, "frac_conf": median_target_predicted_sex_conf, "data_conf": target_data_amount, "data_type": "target"},
        {"sex": mean_antitarget_predicted_sex, "frac_conf": mean_antitarget_predicted_sex_conf, "data_conf": antitarget_data_amount, "data_type": "antitarget"},
        {"sex": median_antitarget_predicted_sex, "frac_conf": median_antitarget_predicted_sex_conf, "data_conf": antitarget_data_amount, "data_type": "antitarget"},
    ]

    for idx, sex_prediction in enumerate(sample_predicted_sex[sample_type]["predicted_sex"]):
        score = calculate_prediction_score(sex_prediction)
        sample_predicted_sex[sample_type]["predicted_sex"][idx].update({"score": score})

    sample_predicted_sex[sample_type]["final_sex"] = consolidate_sex_predictions(sample_predicted_sex[sample_type]["predicted_sex"])

    # Add raw data
    sample_predicted_sex[sample_type]["target_predicted_sex"] = target_predicted_sex
    sample_predicted_sex[sample_type]["antitarget_predicted_sex"] = antitarget_predicted_sex

    return sample_predicted_sex

def case_sex_prediction(predicted_sex):
    predicted_sex["case_sex"] = {}
    tumor_sex = predicted_sex["tumor"]["final_sex"]["sex"]
    if "normal" in predicted_sex:
        normal_sex = predicted_sex["normal"]["final_sex"]["sex"]
        if tumor_sex != normal_sex:
            case_sex = "conflicting"
        else:
            case_sex = tumor_sex
    else:
        case_sex = tumor_sex

    predicted_sex["case_sex"] = case_sex
    return predicted_sex

def write_json(json_obj: dict, path: str) -> None:
    """Write JSON format data to an output file."""
    try:
        with open(path, "w") as fn:
            json.dump(json_obj, fn, indent=4)
    except OSError as error:
        raise OSError(f"Error while writing JSON file: {path}, error: {error}")

@click.command()
@click.option('--target-cnn-tumor', type=click.Path(exists=True), required=True,
              help="Path to the target CNN tumor file.")
@click.option('--antitarget-cnn-tumor', type=click.Path(exists=True), required=True,
              help="Path to the antitarget CNN tumor file.")
@click.option('--output', type=click.Path(writable=True), required=True, help="Path to the output file to be created.")
@click.option('--target-cnn-normal', type=click.Path(exists=True), default=None,
              help="Optional path to the target CNN normal file.")
@click.option('--antitarget-cnn-normal', type=click.Path(exists=True), default=None,
              help="Optional path to the antitarget CNN normal file.")
def process_files(target_cnn_tumor, antitarget_cnn_tumor, output, target_cnn_normal, antitarget_cnn_normal):

    tumor_target_predicted_sex = predict_sex(target_cnn_tumor)
    tumor_antitarget_predicted_sex = predict_sex(antitarget_cnn_tumor)
    predicted_sex = summarise_sample_sex_prediction(tumor_target_predicted_sex, tumor_antitarget_predicted_sex, "tumor")

    if target_cnn_normal:
        normal_target_predicted_sex = predict_sex(target_cnn_normal)
        normal_antitarget_predicted_sex = predict_sex(antitarget_cnn_normal)
        normal_predicted_sex = summarise_sample_sex_prediction(normal_target_predicted_sex, normal_antitarget_predicted_sex, "normal")
        predicted_sex.update(normal_predicted_sex)

    # Create case-level prediction (compare tumor and normal sex)
    predicted_sex = case_sex_prediction(predicted_sex)

    # Write json report
    write_json(predicted_sex, output)





if __name__ == '__main__':
    process_files()