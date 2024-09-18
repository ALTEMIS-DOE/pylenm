import os
import numpy as np
import pandas as pd
from tqdm import tqdm

from pylenm2.data import filters
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

fetchers_logger = logger_config.setup_logging(
    module_name=__name__,
    level=logging.INFO,
    # level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)


def _get_individual_analyte_df(data, dates, analyte):
    """<Function docstring> TODO: write function docstring

    Args:
        data (pd.DataaFrame): input dataframe.
        dates (pd.Series): input dates.
            TODO: Confirm the data type.
        analyte (str): name of the analyte.
    
    Returns:
        pd.DataFrame: Sample analytes based on the input arguments.
    """

    sample = data[analyte]
    sample_analyte = pd.DataFrame(sample, index=dates, columns=sample.index)
    
    for station in sample.index:
        sample_analyte[station] = sample[station]
    
    return sample_analyte


# Helper function to return start and end date for a date and a lag (+/- days)
# def __getLagDate(self, date, lagDays=7):
def _getLagDate(date, lagDays=7):
    date = pd.to_datetime(date)
    dateStart = date - pd.DateOffset(days=lagDays)
    dateEnd = date + pd.DateOffset(days=lagDays)
    return dateStart, dateEnd


def getCleanData(data_pylenm_dm, analytes):
    """Creates a table filling the data from the concentration dataset for a given analyte list where the columns are multi-indexed as follows [analytes, station names] and the index is all of the dates in the dataset. Many NaN should be expected.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        analytes (list): list of analyte names to use

    Returns:
        pd.DataFrame
    """
    curr = data_pylenm_dm.data[['STATION_ID', 'COLLECTION_DATE', 'ANALYTE_NAME', 'RESULT']]
    
    # Get data from only the specified analytes (TODO: improve code)
    main = pd.DataFrame()
    for ana in analytes:
        main = pd.concat([main, curr[curr.ANALYTE_NAME==ana]])

    # Get the pivot data results
    piv = main.pivot_table(index=['COLLECTION_DATE'],columns=['ANALYTE_NAME', 'STATION_ID'], values='RESULT', aggfunc=np.mean)
    
    piv.index = pd.to_datetime(piv.index)
    piv.sort_index(inplace=True)
    
    return piv


def getCommonDates(
        # self, 
        data_pylenm_dm,
        analytes, 
        lag=[3,7,10],
    ) -> pd.DataFrame:
    """Creates a table which counts the number of stations within a range specified by a list of lag days.
    TODO: Very slow. Check if efficiency can be improved.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        analytes (list): list of analyte names to use
        lag (list, optional): list of days to look ahead and behind the specified date (+/-). Defaults to [3,7,10].

    Returns:
        pd.DataFrame
    """
    # piv = self.getCleanData(analytes)
    piv = getCleanData(
        data_pylenm_dm=data_pylenm_dm,
        analytes=analytes,
    )
    dates = piv.index
    names=['Dates', 'Lag']
    tuples = [dates, lag]
    finalData = pd.DataFrame(index=pd.MultiIndex.from_product(tuples, names=names), columns=['Date Ranges', 'Number of stations'])
    
    for date in dates:
        for i in lag:
            # dateStart, dateEnd = self.__getLagDate(date, lagDays=i)
            dateStart, dateEnd = _getLagDate(date, lagDays=i)
            mask = (piv.index > dateStart) & (piv.index <= dateEnd)
            result = piv[mask].dropna(axis=1, how='all')
            numStations = len(list(result.columns.get_level_values(1).unique()))
            dateRange = str(dateStart.date()) + " - " + str(dateEnd.date())
            finalData.loc[date, i]['Date Ranges'] = dateRange
            finalData.loc[date, i]['Number of stations'] = numStations
    
    return finalData


def getJointData(
        data_pylenm_dm, 
        analytes, 
        lag=3,
    ) -> pd.DataFrame:
    """Creates a table filling the data from the concentration dataset for a given analyte list where the columns are multi-indexed as follows [analytes, station names] and the index is the date ranges secified by the lag.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        analytes (list): list of analyte names to use
        lag (int, optional): number of days to look ahead and behind the specified date (+/-). Defaults to 3.

    Returns:
        pd.DataFrame
    """
    # if(self.jointData_is_set(lag=lag)==True):
    if(data_pylenm_dm.is_set_jointData(lag=lag)):
        finalData = data_pylenm_dm.__jointData[0]
        return finalData
    
    piv = getCleanData(
        data_pylenm_dm=data_pylenm_dm,
        analytes=analytes,
    )
    
    dates = piv.index
    dateRanges = []
    for date in dates:
        dateStart, dateEnd = _getLagDate(date, lagDays=lag)
        dateRange = str(dateStart.date()) + " - " + str(dateEnd.date())
        dateRanges.append(dateRange)
    
    finalData = pd.DataFrame(columns=piv.columns, index=dateRanges)
    # numLoops = len(dates)
    # everySomePercent = []
    print("Generating data with a lag of {} days.".format(lag).upper())
    # print("Progress:")
    # for x in list(np.arange(1, 100, 1)):
    #     everySomePercent.append(round((x/100)*numLoops))
    
    # for date, iteration in zip(dates, range(numLoops)):
    for date in tqdm(dates, desc="Progress"):
        # if(iteration in everySomePercent):
        #     print(str(round(iteration/numLoops*100)) + "%", end=', ')
        
        dateStart, dateEnd = _getLagDate(date, lagDays=lag)
        dateRange = str(dateStart.date()) + " - " + str(dateEnd.date())
        mask = (piv.index > dateStart) & (piv.index <= dateEnd)
        result = piv[mask].dropna(axis=1, how='all')
        resultCollapse = pd.concat([result[col].dropna().reset_index(drop=True) for col in result], axis=1)
        
        # HANDLE MULTIPLE VALUES
        if(resultCollapse.shape[0]>1):
            resultCollapse = pd.DataFrame(resultCollapse.mean()).T
        resultCollapse = resultCollapse.rename(index={0: dateRange})
        for ana_station in resultCollapse.columns:
            finalData.loc[dateRange, ana_station] =  resultCollapse.loc[dateRange, ana_station]
    
        # Save data to the pylenm global variable
        # data_pylenm_dm.__set_jointData(data=finalData, lag=lag)
        data_pylenm_dm.set_jointData(data=finalData, lag=lag)
    
    for col in finalData.columns:
        finalData[col] = finalData[col].astype('float64')
    
    print("Completed")
    return finalData


def get_analyte_details(
        data_pylenm_dm, 
        analyte_name, 
        filter=False, 
        col=None, 
        equals=[], 
        save_to_file = False, 
        save_dir='analyte_details',
    ):
    """Returns a csv file saved to save_dir with details pertaining to the specified analyte. Details include the station names, the date ranges and the number of unique samples.

    TODO: Handle error returns better!

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        analyte_name (str): name of the analyte to be processed
        filter (bool, optional): whether to filter the data. Defaults to False.
        col (str, optional): column to filter. Example: col='STATION_ID'. Defaults to None.
        equals (list, optional): values to filter col by. Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].
        save_to_file (bool, optional): whether to save data to file. Defaults to False.
        save_dir (str, optional): name of the directory you want to save the csv file to. Defaults to 'analyte_details'.

    Returns:
        pd.DataFrame: Table with station information
    """
    data = data_pylenm_dm.data
    data = data[data.ANALYTE_NAME == analyte_name].reset_index().drop('index', axis=1)
    data = data[~data.RESULT.isna()]
    data = data.drop(['ANALYTE_NAME', 'RESULT', 'RESULT_UNITS'], axis=1)
    data.COLLECTION_DATE = pd.to_datetime(data.COLLECTION_DATE)
    if(filter):
        filter_res = filters.filter_by_column(
            data=data_pylenm_dm.construction_data, 
            col=col, 
            equals=equals,
        )
        if('ERROR:' in str(filter_res)):    # TODO: Handle Error returns better
            fetchers_logger.error("Ran into ERROR when calling filter_by_column()!")
            return filter_res
        
        query_stations = list(data.STATION_ID.unique())
        filter_stations = list(filter_res.index.unique())
        intersect_stations = list(set(query_stations) & set(filter_stations))
        if(len(intersect_stations)<=0):
            fetchers_logger.error('ERROR: No results for this query with the specifed filter parameters.')
            return 'ERROR: No results for this query with the specifed filter parameters.'
        data = data[data['STATION_ID'].isin(intersect_stations)]        

    info = []
    stations = np.unique(data.STATION_ID.values)
    for station in stations:
        current = data[data.STATION_ID == station]
        startDate = current.COLLECTION_DATE.min().date()
        endDate = current.COLLECTION_DATE.max().date()
        numSamples = current.duplicated().value_counts()[0]
        info.append({
            'Station Name': station, 
            'Start Date': startDate, 
            'End Date': endDate,
            'Date Range (days)': endDate-startDate ,
            'Unique samples': numSamples,
        })
    details = pd.DataFrame(info)
    details.index = details['Station Name']
    details = details.drop('Station Name', axis=1)
    details = details.sort_values(by=['Start Date', 'End Date'])
    details['Date Range (days)'] = (details['Date Range (days)']/ np.timedelta64(1, 'D')).astype(int)
    if(save_to_file):
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        details.to_csv(save_dir + '/' + analyte_name + '_details.csv')
    return details


def get_data_summary(
        data_pylenm_dm, 
        analytes=None, 
        sort_by='date', 
        ascending=False, 
        filter=False, 
        col=None, 
        equals=[],
    ):
    """Returns a dataframe with a summary of the data for certain analytes. Summary includes the date ranges and the number of unique samples and other statistics for the analyte results.

    TODO: Handle error returns better!

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        analytes (list, optional): list of analyte names to be processed. If left empty, a list of all the analytes in the data will be used. Defaults to None.
        sort_by (str, optional): {‘date’, ‘samples’, ‘stations’} sorts the data by either the dates by entering: ‘date’, the samples by entering: ‘samples’, or by unique station locations by entering ‘stations’. Defaults to 'date'.
        ascending (bool, optional): flag to sort in ascending order.. Defaults to False.
        filter (bool, optional): flag to indicate filtering. Defaults to False.
        col (str, optional): column to filter. Example: col='STATION_ID'. Defaults to None.
        equals (list, optional): values to filter col by. Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].

    Returns:
        pd.DataFrame: Table with station information
    """
    data = data_pylenm_dm.data
    if(analytes == None):
        analytes = data.ANALYTE_NAME.unique()
    data = data.loc[data.ANALYTE_NAME.isin(analytes)].drop(['RESULT_UNITS'], axis=1)
    data = data[~data.duplicated()] # remove duplicates
    data.COLLECTION_DATE = pd.to_datetime(data.COLLECTION_DATE)
    data = data[~data.RESULT.isna()]
    if(filter):
        filter_res = filters.filter_by_column(
            data=data_pylenm_dm.construction_data, 
            col=col, 
            equals=equals,
        )
        if('ERROR:' in str(filter_res)):    # TODO: Handle Error returns better
            fetchers_logger.error("Ran into ERROR when calling filter_by_column()!")
            return filter_res
        
        query_stations = list(data.STATION_ID.unique())
        filter_stations = list(filter_res.index.unique())
        intersect_stations = list(set(query_stations) & set(filter_stations))
        if(len(intersect_stations)<=0):
            fetchers_logger.error('ERROR: No results for this query with the specifed filter parameters.')
            return 'ERROR: No results for this query with the specifed filter parameters.'
        data = data[data['STATION_ID'].isin(intersect_stations)]

    info = []
    for analyte_name in analytes:
        query = data[data.ANALYTE_NAME == analyte_name]
        startDate = min(query.COLLECTION_DATE)
        endDate = max(query.COLLECTION_DATE)
        numSamples = query.shape[0]
        stationCount = len(query.STATION_ID.unique())
        stats = query.RESULT.describe().drop('count', axis=0)
        stats = pd.DataFrame(stats).T
        stats_col = [x for x in stats.columns]

        result = {
            'Analyte Name': analyte_name, 
            'Start Date': startDate, 
            'End Date': endDate, 
            'Date Range (days)':endDate-startDate, 
            '# unique stations': stationCount, 
            '# samples': numSamples, 
            'Unit': data_pylenm_dm.get_unit(analyte_name), 
        }
        for num in range(len(stats_col)):
            result[stats_col[num]] = stats.iloc[0][num] 

        info.append(result)

    details = pd.DataFrame(info)
    details.index = details['Analyte Name']
    details = details.drop('Analyte Name', axis=1)
    if(sort_by.lower() == 'date'):
        details = details.sort_values(by=['Start Date', 'End Date', 'Date Range (days)'], ascending=ascending)
    elif(sort_by.lower() == 'samples'):
        details = details.sort_values(by=['# samples'], ascending=ascending)
    elif(sort_by.lower() == 'stations'):
        details = details.sort_values(by=['# unique stations'], ascending=ascending)

    return details


def get_station_analytes(
        data_pylenm_dm, 
        station_names=None, 
        filter=False, 
        col=None, 
        equals=[],
    ):
    """Displays the analyte names available at given station locations.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        station_names (str, List[str], optional): names of the station. If left empty, all stations are returned. The input can either be a str, List[str], or None. Defaults to None.
        filter (bool, optional): flag to indicate filtering. Defaults to False.
        col (str, optional): column to filter. Example: col='STATION_ID'. Defaults to None.
        equals (list, optional): values to filter col by. Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].

    Returns:
        None
    """
    data = data_pylenm_dm.data
    bb = "\033[1m"
    be = "\033[0m"
    
    if(filter):
        filter_res = filters.filter_by_column(
            data=data_pylenm_dm.construction_data, 
            col=col, 
            equals=equals,
        )
        
        if('ERROR:' in str(filter_res)):
            fetchers_logger.error("Ran into ERROR when calling filter_by_column()!")
            return filter_res
        
        query_stations = list(data.STATION_ID.unique())
        filter_stations = list(filter_res.index.unique())
        intersect_stations = list(set(query_stations) & set(filter_stations))
        if(len(intersect_stations)<=0):
            fetchers_logger.error('ERROR: No results for this query with the specifed filter parameters.')
            return 'ERROR: No results for this query with the specifed filter parameters.'
        data = data[data['STATION_ID'].isin(intersect_stations)]
    
    if(station_names==None):
        stations = list(data.STATION_ID.unique())
    elif isinstance(station_names, str):
        stations = [station_names]
    elif isinstance(station_names, list):
        if len(station_names) > 0:
            stations = station_names
        else:
            stations = list(data.STATION_ID.unique())
    else:
        fetchers_logger.error("Invalid input for station_names! The input should either be a str, List[str], or None.")
        raise ValueError("Invalid input for station_names! The input should either be a str, List[str], or None.")
    
    station_analytes = dict()
    for station in stations:
        analytes = sorted(list(data[data.STATION_ID==station].ANALYTE_NAME.unique()))
        
        fetchers_logger.debug(f"{bb}{str(station)}{be}: {str(analytes)}")

        station_analytes[station] = analytes

    return station_analytes