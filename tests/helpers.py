"""Some helper functions and classes to fixtures"""

import json
from pathlib import Path


class ConfigHelper:
    def __init__(self):
        self.case_id = str()

    def read_config(self, balsamic_config):
        with open(balsamic_config, "r") as f:
            sample_config = json.load(f)

        self.case_id = sample_config["analysis"]["case_id"]
        self.analysis_dir = sample_config["analysis"]["analysis_dir"]
        self.result_dir = sample_config["analysis"]["result"]
        self.delivery_dir = Path(self.result_dir, "delivery_report").as_posix()


class Map(dict):
    """Mock class to use dot notation to access values of a dictionary"""

    def __init__(self, *args, **kwargs):
        super(Map, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.items():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.items():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(Map, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(Map, self).__delitem__(key)
        del self.__dict__[key]
