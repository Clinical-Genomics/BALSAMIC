"""Some helper functions and classes to fixtures"""

import json
from pathlib import Path

class ConfigHelper: 
    def __init__(self):
        self.case_id=str()

    def read_config(self, balsamic_config): 
        with open(balsamic_config, 'r') as f: 
            sample_config = json.load(f) 

        self.case_id = sample_config["analysis"]["case_id"] 
        self.analysis_dir = sample_config['analysis']['analysis_dir'] 
        self.result_dir = sample_config['analysis']['result'] 
        self.delivery_dir = Path(self.result_dir, "delivery_report").as_posix() 
