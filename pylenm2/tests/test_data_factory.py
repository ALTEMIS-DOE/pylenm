import unittest
import pandas as pd
import numpy as np
from pylenm2.data.data_factory import PylenmDataFactory
from pylenm2.utils import constants as c

class TestPylenmDataFactory(unittest.TestCase):

    # @classmethod
    # def setUpClass(cls):
    def setUp(self):
        # print(f"Running test: {self.__class__.__name__}.{self._testMethodName}")
        self.pdf = PylenmDataFactory()
        
        self.dummy_data = pd.DataFrame({
            'STATION_ID': ['1', '2', '3'],
            'ANALYTE_NAME': ['A', 'B', 'C'],
            'RESULT': ["result1", "result2", "result3"],
            'RESULT_UNITS': ['unit1', 'unit2', 'unit3'],
            'COLLECTION_DATE': ['coll_date1', 'coll_date2', 'coll_date3'],
        })

        self.dummy_construction_data = pd.DataFrame({
            'STATION_ID': ['1', '2', '3'],
            'GROUND_ELEVATION': ['1', '2', '3'],
            'LONGITUDE': ['1', '2', '3'],
            'LATITUDE': ['1', '2', '3'],
            'WELL_USE': ['active1', 'active2', 'active3'],
            'AQUIFER': ['aq1', 'aq2', 'aq3'],
            'TOTAL_DEPTH': [1.0, 2.0, 3.0],
        })

    # def setUp(self):
    #     print(f"Running test: {self.__class__.__name__}.{self._testMethodName}")

    def test_init(self):
        self.assertIsNone(self.pdf.data)
        self.assertIsNone(self.pdf.construction_data)
        # self.assertEqual(self.pdf.__jointData, [None, 0])

    # def test_has_columns(self):
    #     data = self.dummy_data
    #     required_cols = ['STATION_ID', 'ANALYTE_NAME']
    #     print(f"[SAPLOGS] {self.pdf._has_columns(data, required_cols) = }")
    #     self.assertTrue(self.pdf._has_columns(data, required_cols))

    # def test_is_valid(self):
    #     data = self.dummy_data
    #     required_cols = ['STATION_ID', 'ANALYTE_NAME']
    #     self.assertTrue(self.pdf._is_valid(data, required_cols))

    def test_is_valid_data(self):
        data = self.dummy_data
        self.assertTrue(self.pdf.is_valid_data(data))

    def test_is_valid_construction_data(self):
        construction_data = self.dummy_construction_data
        self.assertTrue(self.pdf.is_valid_construction_data(construction_data))

    def test_setData(self):
        data = self.dummy_data
        self.pdf.setData(data)
        self.assertIsNotNone(self.pdf.data)
        self.assertIsInstance(self.pdf.unit_dictionary, dict)

    def test_setConstructionData(self):
        construction_data = self.dummy_construction_data
        self.pdf.setConstructionData(construction_data)
        self.assertIsNotNone(self.pdf.construction_data)

    def test_get_unit(self):
        data = self.dummy_data
        self.pdf.setData(data)
        self.assertEqual(self.pdf.get_unit('A'), 'unit1')

    def test_is_set_jointData(self):
        self.assertFalse(self.pdf.is_set_jointData(10))


if __name__ == '__main__':
    unittest.main()
