import unittest

from BALSAMIC.tools import get_ref_path


class TestUtils(unittest.TestCase):
    """docstring for TestUtils"""

    def test_get_ref_path(self):
        test_ref = "test_data/references/reference.json"

        test_ref_json = get_ref_path(test_ref)

        self.assertTrue(isinstance(test_ref_json, dict))
        self.assertTrue(test_ref_json['path']['genomefa'].startswith('/'))
