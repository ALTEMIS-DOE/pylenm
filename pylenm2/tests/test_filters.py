import os
import unittest
import numpy as np
import pandas as pd
from pylenm2 import PylenmDataModule
from pylenm2.data import filters
from pylenm2.utils import utilities
from pylenm2.utils import constants as c

class TestDataFunctions(unittest.TestCase):

    def setUp(self):

        url_1 = '../../notebooks2/data/FASB_Data_thru_3Q2015_Reduced_Demo.csv'
        url_2 = '../../notebooks2/data/FASB Well Construction Info.xlsx'

        # Set low_memory=True only if the file size is greater than 50MB
        url_1_file_size_mb = utilities.get_file_size(url_1)
        if url_1_file_size_mb > 50:
            concentration_data = pd.read_csv(url_1, low_memory=True)
        else:
            concentration_data = pd.read_csv(url_1, low_memory=False)

        construction_data = pd.read_excel(url_2, engine='openpyxl')    # for .xlsx files

        self.pdm = PylenmDataModule(
            data=concentration_data, 
            construction_data=construction_data,
            verbose=False,
            logger_level=99,    # setting high value to turn off logging
        )

    def test_simplify_data(self):
        simplified_data = filters.simplify_data(self.pdm.data)
        self.assertIsInstance(simplified_data, pd.DataFrame)
        self.assertEqual(simplified_data.shape[1], len(c.REQUIRED_DATA_COLUMNS) + 1)  # +1 for COLLECTION_TIME

        simplified_data_with_extra_column = filters.simplify_data(self.pdm.data, columns=['STATION_TYPE'])
        self.assertIsInstance(simplified_data_with_extra_column, pd.DataFrame)
        self.assertEqual(simplified_data_with_extra_column.shape[1], len(c.REQUIRED_DATA_COLUMNS) + 2)  # +2 for COLLECTION_TIME and STATION_TYPE

    def test_filter_by_column(self):
        filtered_data = filters.filter_by_column(self.pdm.data, 'STATION_ID', ['FSB133D', 'FSB122D', 'FOB 15D'])
        self.assertIsInstance(filtered_data, pd.DataFrame)
        self.assertEqual(filtered_data.shape[0], 1648)

        filtered_data_with_no_matches = filters.filter_by_column(self.pdm.data, 'STATION_ID', ['4'])
        self.assertIsNone(filtered_data_with_no_matches)

    def test_filter_stations(self):
        filtered_stations = filters.filter_stations(self.pdm, units=['A'])
        self.assertIsInstance(filtered_stations, list)
        self.assertEqual(len(filtered_stations), 31)

        filtered_stations_with_no_matches = filters.filter_stations(self.pdm, ['X'])
        self.assertIsInstance(filtered_stations_with_no_matches, list)
        self.assertEqual(len(filtered_stations_with_no_matches), 0)

    def test_query_data(self):
        queried_data = filters.query_data(
            self.pdm, 
            station_name='FM-A7U', 
            analyte_name='DEPTH_TO_WATER',
        )
        self.assertIsInstance(queried_data, pd.DataFrame)
        self.assertEqual(queried_data.shape[0], 2)

        queried_data_with_no_matches = filters.query_data(self.pdm, '4', 'A')
        self.assertIsNone(queried_data_with_no_matches)

if __name__ == '__main__':
    unittest.main()