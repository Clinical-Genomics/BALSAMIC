import json
import sys

multiple_json = [{
    'file_name': 'multiqc_picard_insertSize.json',
    'required_metrics': ['MEAN_INSERT_SIZE']
}, {
    'file_name': 'multiqc_picard_dups.json',
    'required_metrics': ['PERCENT_DUPLICATION']
}, {
    'file_name':
    'multiqc_picard_HsMetrics.json',
    'required_metrics': [
        'MEAN_TARGET_COVERAGE', 'MEDIAN_TARGET_COVERAGE',
        'PCT_TARGET_BASES_50X', 'PCT_TARGET_BASES_100X',
        'PCT_TARGET_BASES_200X', 'PCT_TARGET_BASES_500X',
        'PCT_TARGET_BASES_1000X', 'FOLD_80_BASE_PENALTY'
    ]
}]


def get_qc_metrics(case_name, output):
    """Reads picard metrics for a particular case and returns a new dict

    Args:
        case_name: given casename

    Returns:
        metrics_dict: dictionary

    """
    qc_data = {}
    file_path = '/home/proj/production/cancer/cases/' + case_name + '/analysis/qc/multiqc_data/'

    # Loop through json files
    for single_json in multiple_json:
        file_name = file_path + single_json['file_name']
        with open(file_name, 'r') as f:
            json_file = json.load(f)
        for k in single_json['required_metrics']:
            for j in json_file.keys():
                sampleid = j.split('_')[1]
                if sampleid not in qc_data:
                    qc_data[sampleid] = {}
                qc_data[sampleid][k] = json_file[j][k]

    with open(output, 'w') as output_file:
        json.dump(qc_data, output_file, indent=2)


if __name__ == '__main__':
    get_qc_metrics(*sys.argv[1:])
