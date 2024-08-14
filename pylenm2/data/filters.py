import os
import pandas as pd

import pylenm2
import pylenm2.data
import pylenm2.data.data_factory
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

filters_logger = logger_config.setup_logging(
    module_name=__name__,
    # level=logging.INFO,
    level=logging.DEBUG,
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
        data (pd.DataFrame, optional): data to simplify. Defaults to None.
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
        data_df = data
    elif isinstance(data, pylenm2.PylenmDataFactory):
        data_df = data.data
    else:
        filters_logger.error("`data` must be either a pandas DataFrame or PylenmDataFactory!")
        raise ValueError("`data` must be either a pandas DataFrame or PylenmDataFactory!")

        
    if columns==None:
        # sel_cols = ['COLLECTION_DATE','STATION_ID','ANALYTE_NAME','RESULT','RESULT_UNITS']
        sel_cols = c.REQUIRED_DATA_COLUMNS
    else:
        # hasColumns = all(item in list(data.columns) for item in columns)
        # if(hasColumns):

        if set(columns).issubset(data_df.columns):
            # sel_cols = ['COLLECTION_DATE','STATION_ID','ANALYTE_NAME','RESULT','RESULT_UNITS'] + columns
            sel_cols = c.REQUIRED_DATA_COLUMNS + columns
        else:
            extra_columns = set(columns).difference(data_df.columns)
            filters_logger.error(f'Following specified column(s) do not exist in the data: {extra_columns}')
            raise ValueError("Specified column(s) do not exist in the data!")

    data_df = data_df[sel_cols]
    data_df.COLLECTION_DATE = pd.to_datetime(data_df.COLLECTION_DATE)
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
    
        elif isinstance(data, pylenm2.PylenmDataFactory):
            data.setData(data_df, verbose=False)
        
        else:
            raise pylenm2.UnreachableCodeError("Code execution should never reach here!")
    
    return data_df


def filter_by_column(self, data=None, col=None, equals=[]):
    """Filters construction data based on one column. You only specify ONE column to filter by, but can selected MANY values for the entry.

    Args:
        data (pd.DataFrame, optional): dataframe to filter. Defaults to None.
        col (str, optional): column to filter. Example: col='STATION_ID'. Defaults to None.
        equals (list, optional): values to filter col by. Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].

    Returns:
        pd.DataFrame: returns filtered dataframe
    """
    if(data is None):
        return 'ERROR: DataFrame was not provided to this function.'
    else:
        if(str(type(data)).lower().find('dataframe') == -1):
            return 'ERROR: Data provided is not a pandas DataFrame.'
        else:
            data = data
    # DATA VALIDATION
    if(col==None):
        return 'ERROR: Specify a column name to filter by.'
    data_cols = list(data.columns)
    if((col in data_cols)==False): # Make sure column name exists 
        return 'Error: Column name "{}" does not exist'.format(col)
    if(equals==[]):
        return 'ERROR: Specify a value that "{}" should equal to'.format(col)
    data_val = list(data[col])
    for value in equals:
        if((value in data_val)==False):
            return 'ERROR: No value equal to "{}" in "{}".'.format(value, col)

    # QUERY
    final_data = pd.DataFrame()
    for value in equals:
        current_data = data[data[col]==value]
        final_data = pd.concat([final_data, current_data])
    return final_data


def filter_wells(self, units):
    """Returns a list of the well names filtered by the unit(s) specified.

    Args:
        units (list): Letter of the well to be filtered (e.g. [‘A’] or [‘A’, ‘D’])

    Returns:
        list: well names filtered by the unit(s) specified
    """
    data = self.data
    if(units==None):
        units= ['A', 'B', 'C', 'D']
    def getUnits():
        wells = list(np.unique(data.STATION_ID))
        wells = pd.DataFrame(wells, columns=['STATION_ID'])
        for index, row in wells.iterrows():
            mo = re.match('.+([0-9])[^0-9]*$', row.STATION_ID)
            last_index = mo.start(1)
            wells.at[index, 'unit'] = row.STATION_ID[last_index+1:]
            u = wells.unit.iloc[index]
            if(len(u)==0): # if has no letter, use D
                wells.at[index, 'unit'] = 'D'
            if(len(u)>1): # if has more than 1 letter, remove the extra letter
                if(u.find('R')>0):
                    wells.at[index, 'unit'] = u[:-1]
                else:
                    wells.at[index, 'unit'] = u[1:]
            u = wells.unit.iloc[index]
            if(u=='A' or u=='B' or u=='C' or u=='D'):
                pass
            else:
                wells.at[index, 'unit'] = 'D'
        return wells
    df = getUnits()
    res = df.loc[df.unit.isin(units)]
    return list(res.STATION_ID)


def query_data(self, well_name, analyte_name):
    """Filters data by passing the data and specifying the well_name and analyte_name

    Args:
        well_name (str): name of the well to be processed
        analyte_name (str): name of the analyte to be processed

    Returns:
        pd.DataFrame: filtered data based on query conditons
    """
    data = self.data
    query = data[data.STATION_ID == well_name]
    query = query[query.ANALYTE_NAME == analyte_name]
    if(query.shape[0] == 0):
        return 0
    else:
        return query