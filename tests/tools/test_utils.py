import unittest
import json

from BALSAMIC.tools import get_ref_path, iterdict


class TestUtils(unittest.TestCase):
    """docstring for TestUtils"""

    def setUp(self):
        self.test_ref_json_path = "tests/test_data/references/reference.json"

    def test_get_ref_path(self):
        # GIVEN a sample json file path
        test_ref = self.test_ref_json_path

        # WHEN giving a path for json file,
        test_ref_json = get_ref_path(test_ref)

        # THEN It will read the file and return a dict with updated absolute path
        self.assertTrue(isinstance(test_ref_json, dict))
        self.assertTrue(test_ref_json['path']['genomefa'].startswith('/'))

    def test_iterdict(self):
        """ GIVEN a dict for iteration """
        test_dict = json.load(open(self.test_ref_json_path, 'r'))

        #WHEN passing dict to this function
        dict_gen = iterdict(test_dict)

        #THEN it will create dict generator, we can iterate it, get the key, values as string
        for key, value in dict_gen:
            assert isinstance(key, str) == True
            assert isinstance(value, str) == True
