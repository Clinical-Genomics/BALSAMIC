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

def predict_sex(cnn_file):
    data_rows = read_cov(cnn_file)
    data_df = process_data(data_rows)
    stats = get_stats(data_df)
    sample_dict = extract_stats(stats)
    sample_name, cnn_type = retrieve_file_info(cnn_file)

    predicted_sex = {}
    predicted_sex["sample_name"] = sample_name
    predicted_sex["cnn_type"] = cnn_type

    if sample_dict["X_count"] < 10 or sample_dict["Y_count"] < 10:
        predicted_sex["data_amount"] = "low"
    elif sample_dict["X_count"] < 50 or sample_dict["Y_count"] < 50:
        predicted_sex["data_amount"] = "medium"
    else:
        predicted_sex["data_amount"] = "high"

    predicted_sex["Y_mean/X_mean"] = round(sample_dict["Y_mean"] / sample_dict["X_mean"], 5)



    sex_prediction = {}

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

    antitarget_cnn_tumor





if __name__ == '__main__':
    process_files()