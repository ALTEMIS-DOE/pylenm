import os
import re
import numpy as np
import pandas as pd

import pylenm2
import pylenm2.data
import pylenm2.data.data_module
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

filters_logger = logger_config.setup_logging(
    module_name=__name__,
    # level=logging.INFO,
    level=logging.ERROR,
    logfile_dir=c.LOGFILE_DIR,
)


def simplify_data(
        data, 
        inplace=False, 
        columns=None, 
        save_csv=False, 
        file_name= 'data_simplified', 
        save_dir='data/',
    ):
    """Removes all columns except 'COLLECTION_DATE', 'STATION_ID', 'ANALYTE_NAME', 'RESULT', and 'RESULT_UNITS', i.e. the `REQUIRED_DATA_COLUMNS`!
        
        If the user specifies additional columns in addition to the ones listed above, those columns will be kept.
        The function returns a dataframe and has an optional parameter to be able to save the dataframe to a csv file.

    Args:
        data (pd.DataFrame, pylenm2.PylenmDataModule, optional): data to simplify.
        inplace (bool, optional): save data to current working dataset. Defaults to False.
        columns (list, optional): list of any additional columns on top of  ['COLLECTION_DATE', 'STATION_ID', 'ANALYTE_NAME', 'RESULT', and 'RESULT_UNITS'] to be kept in the dataframe. Defaults to None.
        save_csv (bool, optional): flag to determine whether or not to save the dataframe to a csv file. Defaults to False.
        file_name (str, optional): name of the csv file you want to save. Defaults to 'data_simplified'.
        save_dir (str, optional): name of the directory you want to save the csv file to. Defaults to 'data/'.

    Returns:
        pd.DataFrame
    """
    # if(str(type(data)).lower().find('dataframe') == -1):
    #     data = self.data
    # else:
    #     data = data
    if isinstance(data, pd.DataFrame):
        data_df = data.copy(deep=True)
    elif isinstance(data, (pylenm2.PylenmDataModule, pylenm2.data.data_module.PylenmDataModule)):
        data_df = data.data.copy(deep=True)
    else:
        filters_logger.error("`data` must be either a pandas DataFrame or PylenmDataModule!")
        # raise ValueError("`data` must be either a pandas DataFrame or PylenmDataModule!")
        return None
    
    required_data_cols = c.REQUIRED_DATA_COLUMNS + ["COLLECTION_TIME"]
    if columns==None:
        # sel_cols = ['COLLECTION_DATE','STATION_ID','ANALYTE_NAME','RESULT','RESULT_UNITS']
        # sel_cols = c.REQUIRED_DATA_COLUMNS
        sel_cols = required_data_cols
    else:
        # hasColumns = all(item in list(data.columns) for item in columns)
        # if(hasColumns):

        if set(columns).issubset(data_df.columns):
            # sel_cols = ['COLLECTION_DATE','STATION_ID','ANALYTE_NAME','RESULT','RESULT_UNITS'] + columns
            # sel_cols = c.REQUIRED_DATA_COLUMNS + columns
            sel_cols = required_data_cols + columns
        else:
            extra_columns = set(columns).difference(data_df.columns)
            filters_logger.error(f'Following specified column(s) do not exist in the data: {extra_columns}')
            raise ValueError("Specified column(s) do not exist in the data!")

    data_df = data_df[sel_cols]

    # # Formating the date to Datetime
    # # data_df.COLLECTION_DATE = pd.to_datetime(data_df.COLLECTION_DATE, format=c.COLLECTION_DATE_FORMAT)
    # try:
    #     # For when time is specified in the HH:MM format
    #     data_df.COLLECTION_DATE = pd.to_datetime(
    #         data_df.COLLECTION_DATE, 
    #         format=f"{c.COLLECTION_DATE_FORMAT} {c.COLLECTION_TIME_FORMAT}",
    #     )
    # except ValueError as ve:
    #     filters_logger.warning(f"COLLECTION_DATE data does not match the specified '{c.COLLECTION_DATE_FORMAT} {c.COLLECTION_TIME_FORMAT}' format. Checking '{c.COLLECTION_DATE_FORMAT}' format now.")
        
    #     data_df.COLLECTION_DATE = pd.to_datetime(
    #         data_df.COLLECTION_DATE, 
    #         format=c.COLLECTION_DATE_FORMAT,
    #     )
    
    # # Creating a new column for COLLECTION_TIME
    # data_df["COLLECTION_TIME"] = data_df.COLLECTION_DATE.dt.strftime(c.COLLECTION_TIME_FORMAT)

    data_df = data_df.sort_values(by="COLLECTION_DATE")
    dup = data_df[data_df.duplicated(['COLLECTION_DATE', 'STATION_ID','ANALYTE_NAME', 'RESULT'])]
    data_df = data_df.drop(dup.index)
    data_df = data_df.reset_index().drop('index', axis=1)
    
    if(save_csv):
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        data_df.to_csv(save_dir + file_name + '.csv')
        print('Successfully saved "' + file_name +'.csv" in ' + save_dir)
    
    if(inplace):
        if isinstance(data, pd.DataFrame):
            data.drop(data.index, inplace=True)
            data.drop(columns=list(set(data.columns).difference(sel_cols)), inplace=True)
            data[sel_cols] = data_df[sel_cols]
    
        elif isinstance(data, (pylenm2.PylenmDataModule, pylenm2.data.data_module.PylenmDataModule)):
            data.set_data(data_df, verbose=False)
        
        else:
            raise pylenm2.UnreachableCodeError("Code execution should never reach here!")
    
    return data_df


def filter_by_column(data, col, equals=[]):
    """Filters construction data based on one column. You only specify ONE column to filter by, but can selected MANY values for the entry.

    TODO: Handle Error returns better!

    Args:
        data (pd.DataFrame): dataframe to filter.
        col (str, optional): column to filter. Example: col='STATION_ID'. Defaults to None.
        equals (list, optional): values to filter col by. Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].

    Returns:
        pd.DataFrame: returns filtered dataframe
    """
    if (data is None):
        filters_logger.error('ERROR: DataFrame was not provided to this function.')
        # return 'ERROR: DataFrame was not provided to this function.'
        return None
    else:
        # if(str(type(data)).lower().find('dataframe') == -1):
        if not isinstance(data, pd.DataFrame):
            filters_logger.error('ERROR: Data provided is not a pandas DataFrame.')
            # return 'ERROR: Data provided is not a pandas DataFrame.'
            return None
        else:
            data = data
    
    # DATA VALIDATION
    if (col==None):
        filters_logger.error('ERROR: Specify a column name to filter by.')
        # return 'ERROR: Specify a column name to filter by.'
        return None

    # data_cols = list(data.columns)
    # if((col in data_cols)==False): # Make sure column name exists 
    if col not in data.columns: # Make sure column name exists 
        filters_logger.error(f'Error: Column name {col} does not exist')
        # return 'Error: Column name "{}" does not exist'.format(col)
        return None
    
    if (equals==[]):
        filters_logger.error(f"ERROR: Specify a value that {col} should equal to.")
        # return 'ERROR: Specify a value that "{}" should equal to'.format(col)
        return None
    
    # data_val = list(data[col])
    # for value in equals:
    #     if((value in data_val)==False):
    #         return 'ERROR: No value equal to "{}" in "{}".'.format(value, col)
    values_not_in_data = set(equals).difference(set(data[col]))
    if len(values_not_in_data) > 0:
        filters_logger.warning(f'WARNING: {values_not_in_data} do not exist in {col}.')
        # return f'ERROR: {values_not_in_data} do not exist in {col}.'
        return None

    # QUERY
    # final_data = pd.DataFrame()
    # for value in equals:
    #     current_data = data[data[col]==value]
    #     final_data = pd.concat([final_data, current_data])
    final_data = data[[item in equals for item in data[col]]]
    return final_data


def filter_stations(data_pylenm_dm, units):
    """Returns a list of the station names filtered by the unit(s) specified.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        units (list): Letter of the station to be filtered (e.g. [‘A’] or [‘A’, ‘D’])

    Returns:
        list: station names filtered by the unit(s) specified
    """
    
    data = data_pylenm_dm.data
    
    if units==None:
        units = ['A', 'B', 'C', 'D']
    elif not isinstance(units, (list, tuple)):
        units = [units]
    
    def getUnits():
        stations = list(np.unique(data.STATION_ID))
        stations = pd.DataFrame(stations, columns=['STATION_ID'])
        
        for index, row in stations.iterrows():
            mo = re.match('.+([0-9])[^0-9]*$', row.STATION_ID)
            last_index = mo.start(1)
            stations.at[index, 'unit'] = row.STATION_ID[last_index+1:]
            u = stations.unit.iloc[index]
            
            if(len(u)==0): # if has no letter, use D
                stations.at[index, 'unit'] = 'D'
            if(len(u)>1): # if has more than 1 letter, remove the extra letter
                if(u.find('R')>0):
                    stations.at[index, 'unit'] = u[:-1]
                else:
                    stations.at[index, 'unit'] = u[1:]
            
            u = stations.unit.iloc[index]
            
            if(u=='A' or u=='B' or u=='C' or u=='D'):
                pass
            else:
                stations.at[index, 'unit'] = 'D'
        return stations
    
    df = getUnits()
    
    res = df.loc[df.unit.isin(units)]
    
    return list(res.STATION_ID)


def query_data(data_pylenm_dm, station_name, analyte_name):
    """Filters data by passing the data and specifying the station_name and analyte_name

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        station_name (str): name of the station to be processed
        analyte_name (str): name of the analyte to be processed

    Returns:
        pd.DataFrame: filtered data based on query conditons
    """
    data = data_pylenm_dm.data
    
    query = data[data.STATION_ID == station_name]
    query = query[query.ANALYTE_NAME == analyte_name]
    
    if query.shape[0]==0:
        return None
        # return 0        # TODO: Handle this better!
        # return pd.DataFrame(columns=data.columns)   # TODO: Use this once you make sure that the above return value is not being used anywhere else.
    else:
        return query
