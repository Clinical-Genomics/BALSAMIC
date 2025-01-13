Custom QC tools
======================

This is by no means an exhaustive list of all quality control steps done in balsamic, but contains only a subset of the custom developed scripts for quality control.

Sex check
======================

Since 'Balsamic 17.0.0', the given sample sex is being compared to the predicted sample sex. This was added to detect the situations where a sample mix-up has occurred with a sample of the opposite sex.

The sample sex prediction is done in 2 different ways, depending on if the case is TGA or WGS, and the final result is output into the file: `analysis/qc/sex_prediction.json`

In the final quality control step in balsamic the sexes of the samples included in the case are compared to the given sex in the case config. If the predicted sex of any sample is conflicting with the case in the config, the quality control fails.

Read more about the development and testing of this feature in this `Github issue <https://github.com/Clinical-Genomics/BALSAMIC/issues/1517>`_.


Sex prediction of WGS cases
----------------------------

For WGS cases the sex of the samples are predicted using the custom python script in balsamic: `sex_prediction_wgs.py` based on the files:

-   tumor.[sample-id]\_X_cov_per_base.txt
-   tumor.[sample-id]\_Y_cov_per_base.txt
-   normal.[sample-id]\_X_cov_per_base.txt (optional for TN cases)
-   normal.[sample-id]\_Y_cov_per_base.txt (optional for TN cases)

Which contains the coverages of each base in chromosome X and Y respectively. The sex of the sample is then determined simply by calculating the fraction of the median coverage of the Y-bases and X-bases, and if the fraction of [median per base coverage on Y] / [median per base coverage on X] > 0.1 the sample is predicted to be male, and if it's between 0.1 and 0.08, it's set as unknown and otherwise as female.

Most often the fraction for females is ~ 0.0, and for males ~ 1.0.

All data used in the prediction are stored in the final `sex_prediction.json`

This JSON file is formatted to contain one sub-dict per sample, "tumor, and "normal". Each sample level dict contains:
- predicted_sex: str (the predicted sample sex)
- median_x_coverage: int (the median coverage over x chromosome)
- median_y_coverage: int (the median coverage over y chromosome)
- y_x_median_frac: float (the fraction of median per base y coverage divided by the median per base x coverage)

Sex prediction of TGA cases
-----------------------------

For all TGA cases regardless of panel, the intermediate CNN files from CNVkit is used:

-   [sample-id].targetcoverage.cnn
-   [sample-id].antitargetcoverage.cnn

The first file contains the regions of the panel bedfile, and the second contains the regions of the genome outside the targets defined the panel bedfile.

Both of these files contain the average coverage within these genomic target regions and anti-target regions, and these are used to calculate mean and median coverages across these regions in the X and Y chromosomes. Based on this, similarly to the WGS TO method, fractions of relative coverage across the Y and X chromosomes are calculated and depending on where the fractions fall on certain thresholds, the sex is determined to be male or female.

For each sample the both the mean and median fractions from both files, a total of 4 fractions, are consolidated to one final sample prediction. For each Y / X fraction the sex is predicted based on the thresholds:

+---------------------------+-----------------+--------------------------+
| Y / X fraction thresholds | Predicted Sex  | Predicted Sex Confidence |
+===========================+=================+==========================+
| X > 0.9                   | male           | high                     |
+---------------------------+-----------------+--------------------------+
| 0.9 < X <= 0.7            | male           | medium                   |
+---------------------------+-----------------+--------------------------+
| 0.7 < X <= 0.2            | male           | low                      |
+---------------------------+-----------------+--------------------------+
| 0.2 < X <= 0.15           | unknown        | low                      |
+---------------------------+-----------------+--------------------------+
| 0.15 < X <= 0.05          | female         | medium                   |
+---------------------------+-----------------+--------------------------+
| 0.05 < X <= 0.00001       | female         | high                     |
+---------------------------+-----------------+--------------------------+
| X < 0.00001               | female         | low                      |
+---------------------------+-----------------+--------------------------+

Then the amount of data for each file in X and Y is taken into account to create for each sex prediction a final score, where high confidence and larger amounts of data in the prediction are given more points, and the target file which was shown to be less noisy is prioritised, and the prediction with the highest score is set as the final sample sex prediction.

All data used in the prediction are stored in the final `sex_prediction.json`

This JSON file is formatted to contain one sub-dict per sample, "tumor, and "normal". Each sample level dict contains:
- predicted_sex: str (the final predicted sample sex)
- predicted_sex_score: int (the final score of the prediction)
- prediction_confidence: str (the confidence of the prediction; low | medium | high)
- target_predicted_sex: dict (sex prediction based only on the target coverage cnn file)
- antitarget_predicted_sex: dict (sex prediction based only on the antitarget coverage cnn file)

Each target and antitarget prediction dictionary also contains the raw data used in the prediction, such as:
- Y_mean/X_mean: float (the fraction of the mean coverage over the regions in Y divided by the mean coverage over the regions in X)
- Y_median/X_median: float (the fraction of the median coverage over the regions in Y divided by the median coverage over the regions in X)
