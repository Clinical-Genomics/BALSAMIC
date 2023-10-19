
from BALSAMIC.constants.analysis import SequencingType

# GENS parameters
GENS_PARAMS = {
    "MINIMUM_WINDOW_SIZE": 100,
    "ALLOWED_CHR_LIST": [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
        "X",
        "Y",
        "MT",
    ],
    "SEQUENCING_TYPE": {
        SequencingType.WGS: {
            "BAF_SKIP_N": {"o": 135, "a": 30, "b": 8, "c": 3, "d": 1},
            "COV_WINDOW_SIZES": {"o": 100000, "a": 25000, "b": 5000, "c": 1000, "d": 100},
        }
    }
}
