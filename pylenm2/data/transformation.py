import numpy as np
import pandas as pd

from typing import Tuple

from pylenm2.data import fetchers
from pylenm2.stats import metrics
from pylenm2.stats import preprocess
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

transformation_logger = logger_config.setup_logging(
    module_name=__name__,
    # level=logging.INFO,
    level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)


def interpolate_well_data(
        # self, 
        data_pylenm_dm, 
        well_name, 
        analytes, 
        frequency='2W',
    ) -> pd.DataFrame:
    """Resamples the data based on the frequency specified and interpolates the values of the analytes.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        well_name (str): name of the well to be processed.
        analytes (list): list of analyte names to use
        frequency (str, optional): {‘D’, ‘W’, ‘M’, ‘Y’} frequency to interpolate. See https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html for valid frequency inputs. (e.g. ‘W’ = every week, ‘D ’= every day, ‘2W’ = every 2 weeks). Defaults to '2W'.

    Returns:
        pd.DataFrame
    """
    # data = self.data
    data = data_pylenm_dm.data
    inter_series = {}
    query = data[data.STATION_ID == well_name]
    
    for analyte in analytes:
        series = query[query.ANALYTE_NAME == analyte]
        series = (series[['COLLECTION_DATE', 'RESULT']])
        series.COLLECTION_DATE = pd.to_datetime(series.COLLECTION_DATE)
        series.index = series.COLLECTION_DATE
        original_dates = series.index
        breakpoint()
        series = series.drop('COLLECTION_DATE', axis=1)
        series = series.rename({'RESULT': analyte}, axis=1)
        upsampled = series.resample(frequency).mean()
        interpolated = upsampled.interpolate(method='linear', order=2)
        inter_series[analyte] = interpolated
    
    join = inter_series[analytes[0]]
    join = join.drop(analytes[0], axis=1)
    
    for analyte in analytes:
        join = join.join(inter_series[analyte])
    
    join = join.dropna()
    
    return join


def interpolate_wells_by_analyte(
        data_pylenm_dm, 
        analyte, 
        frequency='2W', 
        rm_outliers=True, 
        z_threshold=3,
    ) -> pd.DataFrame:
    """Resamples analyte data based on the frequency specified and interpolates the values in between. NaN values are replaced with the average value per well.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        analyte (_type_): analyte name for interpolation of all present wells.
        frequency (str, optional): {‘D’, ‘W’, ‘M’, ‘Y’} frequency to interpolate. See https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html for valid frequency inputs. (e.g. ‘W’ = every week, ‘D ’= every day, ‘2W’ = every 2 weeks). Defaults to '2W'.
        rm_outliers (bool, optional): flag to remove outliers in the data. Defaults to True.
        z_threshold (int, optional): z_score threshold to eliminate outliers. Defaults to 3.

    Returns:
        pd.DataFrame: interpolated data
    """
    # data = data_pylenm_dm.data

    df_t, dates = _transform_time_series( 
        data_pylenm_dm=data_pylenm_dm,
        analytes=[analyte], 
        resample=frequency, 
        rm_outliers=True, 
        z_threshold=z_threshold,
    )

    res_interp = fetchers._get_individual_analyte_df(
        data=df_t, 
        dates=dates, 
        analyte=analyte,
    )
    res_interp = res_interp.dropna(axis=1, how='all')
    
    return res_interp


# IN THE WORKS
def _transform_time_series(
        data_pylenm_dm, 
        analytes=[], 
        resample='2W', 
        rm_outliers=False, 
        z_threshold=4,
    ) -> Tuple[pd.DataFrame, pd.DatetimeIndex]:
    """<Function docstring> TODO: write function docstring.
    TODO: The function can be optimized a lot for a faster performance. Come back to this once everything is done.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        ... TODO: Write the args docstring.

    Returns:
        Tuple[pd.DataFrame, pd.DatetimeIndex]: Returns a tuple of the dataframe and the dates. TODO: Write more details about the returned data.
    """

    data = data_pylenm_dm.data
    
    def transform_time_series_by_analyte(data, analyte_name):
        """<Nested-Function docstring> TODO: write function docstring.

        Args:
            data (pd.DataFrame): input dataframe.
            analyte_name (str): name of the analyte.
        
        Returns:
            pd.DataFrame: return sample dataframe for the analyte.
        """
        # wells_analyte = np.unique(data[data.ANALYTE_NAME == analyte_name].STATION_ID)

        # all_dates = np.unique(data.COLLECTION_DATE)
        # Create array of equally spaced dates
        start_date = pd.Timestamp(data.COLLECTION_DATE.min())
        end_date = pd.Timestamp(data.COLLECTION_DATE.max())
        date_delta = (end_date - start_date) + pd.Timedelta(days=1)    # to include the end date as well
        t = np.linspace(start_date.value, end_date.value, date_delta.days)
        t = pd.to_datetime(t).date
        # t = pd.Series(t)
        # t = t.apply(lambda x: x.replace(minute=0, hour=0, second=0, microsecond=0, nanosecond=0))

        # condensed = data[data.ANALYTE_NAME == analyte_name].groupby(['STATION_ID','COLLECTION_DATE']).mean()    # NOTE: Breaks the code
        # condensed = data[data.ANALYTE_NAME == analyte_name].groupby(['STATION_ID','COLLECTION_DATE'])['RESULT'].mean().to_frame('RESULT')     # NOTE: Works. Result must have (well, date) as index
        condensed = data[data.ANALYTE_NAME == analyte_name].groupby(['STATION_ID','COLLECTION_DATE'], as_index=False)['RESULT'].mean()  # NOTE: Much better approach.

        analyte_df_resample = condensed.pivot(columns="COLLECTION_DATE", index="STATION_ID", values="RESULT")

        # analyte_df_resample = pd.DataFrame(index=wells_analyte, columns=t)
        analyte_df_resample.sort_index(inplace=True)
        
        # for well in wells_analyte:    # NOTE: Performs pivot. Implemented more efficiently above.
        #     for date in condensed.loc[well].index:
        #         analyte_df_resample.at[well, pd.to_datetime(date).date()] = condensed.loc[well,date].RESULT
        
        analyte_df_resample = analyte_df_resample.astype('float').T
        analyte_df_resample = analyte_df_resample.interpolate(method='linear')
        return analyte_df_resample


    data_analyte_groups = data.groupby("ANALYTE_NAME", as_index=False)
    transformed_ts_data_analyte = data_analyte_groups.apply(
        lambda subdf: transform_time_series_by_analyte(data=subdf, analyte_name=subdf.ANALYTE_NAME.iloc[0])
    )


    # Save each analyte data
    cutoff_dates = []
    analyte_data = []
    for analyte in analytes:
        ana_data = transform_time_series_by_analyte(data, analyte)

        if(rm_outliers):
            col_num = ana_data.shape[1]

            for col in range(col_num):

                try:
                    ana_data.iloc[:,col] = preprocess.remove_outliers(
                        ana_data.iloc[:,col], 
                        z_threshold=z_threshold,
                        nan_policy='omit',  # Omit nan values while computing.
                    )
                except Exception as e:
                    transformation_logger.error(e)

            ana_data = ana_data.interpolate(method='linear')
        
        ana_data.index = pd.to_datetime(ana_data.index)
        
        # Resample
        ana_data_resample = ana_data.resample(resample).mean()
        
        # Save data
        analyte_data.append(ana_data_resample)
        
        # Determine cuttoff point for number of NaNs in dataset
        passes_limit = []
        for date in ana_data_resample.index:
            limit = 0.7 * ana_data_resample.shape[1]
            curr = ana_data_resample.isna().loc[date,:].value_counts()
            if('False' in str(curr)):
                curr_total = ana_data_resample.isna().loc[date,:].value_counts()[0]
                if curr_total > limit:
                    passes_limit.append(date)
        passes_limit = pd.to_datetime(passes_limit)
        
        cutoff_dates.append(passes_limit.min())
    
    start_index = pd.Series(cutoff_dates).max()

    # Get list of shared wells amongst all the listed analytes
    combined_well_list = []
    for x in range(len(analytes)):
        combined_well_list = combined_well_list + list(analyte_data[x].columns)
    
    combined_count = pd.Series(combined_well_list).value_counts()
    shared_wells = list(
        combined_count[
            list(pd.Series(combined_well_list).value_counts()==len(analytes))
        ].index
    )

    # Vectorize data
    vectorized_df = pd.DataFrame(columns=analytes, index=shared_wells)

    for analyte, num in zip(analytes, range(len(analytes))):
        for well in shared_wells:
            analyte_data_full = analyte_data[num][well].fillna(analyte_data[num][well].mean())
            vectorized_df.at[well, analyte] = analyte_data_full[start_index:].values

    dates = ana_data_resample[start_index:].index
    
    return vectorized_df, dates


def add_dist_to_source(
        XX, 
        source_coordinate=c.DEFAULT_SOURCE_COORDINATES, 
        col_name='dist_to_source',
    ) -> pd.DataFrame:
    """adds column to data with the distance of a record to the source coordinate

    Args:
        XX (pd.DataFrame): data with coordinate information
        source_coordinate (list, optional): source coordinate. Defaults to [436642.70,3681927.09].
        col_name (str, optional): name to assign new column. Defaults to 'dist_to_source'.

    Returns:
        pd.DataFrame: returns original data with additional column with the distance.
    """
    x1,y1 = source_coordinate
    distances = []
    for i in range(XX.shape[0]):
        x2,y2 = XX.iloc[i][0], XX.iloc[i][1]
        # distances.append(self.dist([x1,y1],[x2,y2]))
        distances.append(metrics.dist([x1,y1],[x2,y2]))
    XX[col_name] = distances
    return XX