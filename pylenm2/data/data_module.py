import numpy as np
import pandas as pd

from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

df_logger = logger_config.setup_logging(
    module_name=__name__,
    level=logging.INFO,
    # level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)


# class PylenmDataFactory(object):
class PylenmDataModule(object):
    """Class object that initilaizes Pylenm given data.
    """
    
    def __init__(
        self, 
        data: pd.DataFrame=None,
        construction_data: pd.DataFrame=None,
        verbose=True,
        logger_level:int=logging.WARNING,
    ) -> None:
        """Initializes pylenm with a Pandas DataFrame

        Args:
            data (pd.DataFrame): Data to be imported.
        """

        self.logger = logger_config.setup_logging(module_name=__name__, level=logger_level)

        if data is not None: 
            self.set_data(data, verbose=verbose)
        else:
            self.data = None
        
        if construction_data is not None:
            self.set_construction_data(construction_data, verbose=verbose)
        else:
            self.construction_data = None
        
        self.__jointData = [None, 0]    # TODO: Check if we can make this a single value


    def __has_columns(self, data: pd.DataFrame, required_cols=list()):
        data_cols = set([x.upper() for x in data.columns])
        return set(required_cols).issubset(data_cols)


    def __is_valid(self, data: pd.DataFrame, required_cols=list()):
        """Validates data for Pylenm usage.

        Args:
            data (pd.DataFrame): Data to be validated.
            required_cols (list, optional): List of required columns. Defaults to REQUIRED_DATA_COLUMNS.

        Returns:
            bool: True if data is valid, False otherwise.
        """

        if data is None:
            self.logger.error("`data` is None!")
            return False
        
        if not isinstance(data, pd.DataFrame):
            self.logger.error("`data` must be a Pandas DataFrame!")
            return False
        
        if not self.__has_columns(data=data, required_cols=required_cols):
            missing_cols = list(
                set(required_cols).difference(set(data.columns))
            )
            self.logger.error(f"Missing these columns from the data: {missing_cols}.")
            return False

        return True


    def is_valid_data(self, data: pd.DataFrame):
        return self.__is_valid(data, required_cols=c.REQUIRED_DATA_COLUMNS)
    

    def is_valid_construction_data(self, data: pd.DataFrame):
        return self.__is_valid(data, required_cols=c.REQUIRED_CONSTRUCTION_DATA_COLUMNS)

    
    def __set_units(self):
        # analytes = list(np.unique(self.data[['ANALYTE_NAME']]))
        # mask1 = ~self.data[['ANALYTE_NAME','RESULT_UNITS']].duplicated()
        # res = self.data[['ANALYTE_NAME','RESULT_UNITS']][mask1]
        # mask2 = ~self.data[['ANALYTE_NAME']].duplicated()
        # res = res[mask2]        # NOTE: Generates warning => UserWarning: Boolean Series key will be reindexed to match DataFrame index. # TODO: Check later.
        
        res = self.data[['ANALYTE_NAME','RESULT_UNITS']].drop_duplicates(
            subset="ANALYTE_NAME",
        )
        unit_dictionary = pd.Series(res.RESULT_UNITS.values,index=res.ANALYTE_NAME).to_dict()
        self.unit_dictionary = unit_dictionary
    

    def get_unit(self, analyte_name):
        """Returns the unit of the analyte you specify. Example: 'DEPTH_TO_WATER' may return 'ft'

        Args:
            analyte_name (str): ame of the analyte to be processed

        Returns:
            str: unit of analyte
        """
        return self.unit_dictionary[analyte_name]


    # SETTING DATA
    def set_data(self, data: pd.DataFrame, verbose: bool=False) -> None:
        """Saves the dataset into pylenm.

        Args:
            data (pd.DataFrame): Dataset to be imported.
            verbose (bool, optional): Prints success message. Defaults to True.

        Returns:
            None
        """
        if self.is_valid_data(data):
            cols_upper = [x.upper() for x in data.columns]
            data.columns = cols_upper
            self.data = data
            
            if "COLLECTION_DATE" in self.data.columns:
                self.update_collection_date_and_time()
            
            self.logger.info("Successfully imported the data!")
            
            if(verbose):
                print('Successfully imported the data!\n')
            
            self.__set_units()
        
        else:
            self.logger.info("Failed to import data!")
            self.logger.warning(c.PYLENM_DATA_REQUIREMENTS)


    def set_construction_data(
        self, 
        construction_data: pd.DataFrame, 
        verbose: bool=False,
    ) -> None:
        """Imports the additional station information as a separate DataFrame.

        Args:
            construction_data (pd.DataFrame): Data with additonal details.
            verbose (bool, optional): Prints success message. Defaults to True.

        Returns:
            None
        """
        if self.is_valid_construction_data(construction_data):
            cols_upper = [x.upper() for x in list(construction_data.columns)]
            construction_data.columns = cols_upper
            self.construction_data = construction_data.set_index(['STATION_ID'])
            self.logger.info("Successfully imported the construction data!")
            
            if(verbose):
                print('Successfully imported the construction data!\n')
        
        else:
            self.logger.info("Failed to import the construction data!")
            self.logger.warning(c.PYLENM_CONSTRUCTION_DATA_REQUIREMENTS)


    def set_jointData(self, data, lag):
        self.__jointData[0] = data
        self.__jointData[1] = lag

    
    @property
    def jointData(self):
        return self.__jointData
    

    # def jointData_is_set(self, lag):
    def is_set_jointData(self, lag):
        """Checks to see if getJointData function was already called and saved for given lag.

        Args:
            lag (int): number of days to look ahead and behind the specified date (+/-)

        Returns:
            bool: True if JointData was already calculated, False, otherwise.
        """
        # if(str(type(self.__jointData[0])).lower().find('dataframe') == -1):
        if not isinstance(self.__jointData[0], pd.DataFrame):
            return False
        
        if self.__jointData[1] != lag:
            return False
        
        return True


    def get_data(self):
        return self.data
    

    def get_construction_data(self):
        return self.construction_data


    def update_collection_date_and_time(self):
        # Formating the date to Datetime
        is_dayfirst = c.COLLECTION_DATE_FORMAT.startswith("%d")
        self.data.COLLECTION_DATE = pd.to_datetime(self.data.COLLECTION_DATE,  format="mixed", dayfirst=is_dayfirst)
        
        # Creating a new column for COLLECTION_TIME
        self.data["COLLECTION_TIME"] = self.data.COLLECTION_DATE.dt.strftime(c.COLLECTION_TIME_FORMAT)