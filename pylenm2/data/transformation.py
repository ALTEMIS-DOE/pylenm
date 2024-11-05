import numpy as np
import pandas as pd
from tqdm import tqdm
from typing import Tuple

from pylenm2.data import fetchers
from pylenm2.stats import metrics
from pylenm2.stats import preprocess
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

transformation_logger = logger_config.setup_logging(
    module_name=__name__,
    level=logging.INFO,
    # level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)


def interpolate_station_data(
        # self, 
        data_pylenm_dm, 
        station_name, 
        analytes, 
        frequency='2W',
    ) -> pd.DataFrame:
    """Resamples the data based on the frequency specified and interpolates the values of the analytes.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        station_name (str): name of the station to be processed.
        analytes (list): list of analyte names to use
        frequency (str, optional): {‘D’, ‘W’, ‘M’, ‘Y’} frequency to interpolate. See https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html for valid frequency inputs. (e.g. ‘W’ = every week, ‘D ’= every day, ‘2W’ = every 2 weeks). Defaults to '2W'.

    Returns:
        pd.DataFrame
    """
    # data = self.data
    data = data_pylenm_dm.data
    inter_series = {}
    query = data[data.STATION_ID == station_name]
    
    for analyte in analytes:
        series = query[query.ANALYTE_NAME == analyte]
        series = (series[['COLLECTION_DATE', 'RESULT']])
        series.COLLECTION_DATE = pd.to_datetime(series.COLLECTION_DATE)
        series.index = series.COLLECTION_DATE
        original_dates = series.index
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


def interpolate_stations_by_analyte(
        data_pylenm_dm, 
        analyte, 
        frequency='2W', 
        rm_outliers=True, 
        z_threshold=3,
    ) -> pd.DataFrame:
    """Resamples analyte data based on the frequency specified and interpolates the values in between. NaN values are replaced with the average value per station.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        analyte (_type_): analyte name for interpolation of all present stations.
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
    TODO: Get the logic confirmed by Zexuan.

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
        # stations_analyte = np.unique(data[data.ANALYTE_NAME == analyte_name].STATION_ID)

        # # Create array of equally spaced dates
        # start_date = pd.Timestamp(data.COLLECTION_DATE.min())
        # end_date = pd.Timestamp(data.COLLECTION_DATE.max())
        # date_delta = (end_date - start_date) + pd.Timedelta(days=1)    # to include the end date as station
        # t = np.linspace(start_date.value, end_date.value, date_delta.days)
        # t = pd.to_datetime(t).date
        # # t = pd.Series(t)
        # # t = t.apply(lambda x: x.replace(minute=0, hour=0, second=0, microsecond=0, nanosecond=0))

        # condensed = data[data.ANALYTE_NAME == analyte_name].groupby(['STATION_ID','COLLECTION_DATE']).mean()    # NOTE: Breaks the code
        # condensed = data[data.ANALYTE_NAME == analyte_name].groupby(['STATION_ID','COLLECTION_DATE'])['RESULT'].mean().to_frame('RESULT')     # NOTE: Works. Result must have (station, date) as index
        condensed = data[data.ANALYTE_NAME == analyte_name].groupby(['STATION_ID','COLLECTION_DATE'], as_index=False)['RESULT'].mean()  # NOTE: Much better approach.

        analyte_df_resample = condensed.pivot(columns="COLLECTION_DATE", index="STATION_ID", values="RESULT")

        # analyte_df_resample = pd.DataFrame(index=stations_analyte, columns=t)
        analyte_df_resample.sort_index(inplace=True)
        
        # for station in stations_analyte:    # NOTE: Performs pivot. Implemented more efficiently above.
        #     for date in condensed.loc[station].index:
        #         analyte_df_resample.at[station, pd.to_datetime(date).date()] = condensed.loc[station,date].RESULT
        
        analyte_df_resample = analyte_df_resample.astype('float').T
        analyte_df_resample = analyte_df_resample.interpolate(method='linear')
        return analyte_df_resample


    # data_analyte_groups = data.groupby("ANALYTE_NAME", as_index=False)
    # transformed_ts_data_analyte = data_analyte_groups.apply(
    #     lambda subdf: transform_time_series_by_analyte(
    #         data=subdf, analyte_name=subdf.ANALYTE_NAME.iloc[0]
    #     )
    # )


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

    # Get list of shared stations amongst all the listed analytes
    combined_station_list = []
    for x in range(len(analytes)):
        combined_station_list = combined_station_list + list(analyte_data[x].columns)
    
    combined_count = pd.Series(combined_station_list).value_counts()
    shared_stations = list(
        combined_count[
            list(pd.Series(combined_station_list).value_counts()==len(analytes))
        ].index
    )

    # Vectorize data
    vectorized_df = pd.DataFrame(columns=analytes, index=shared_stations)

    for analyte, num in zip(analytes, range(len(analytes))):
        for station in shared_stations:
            analyte_data_full = analyte_data[num][station].fillna(analyte_data[num][station].mean())
            vectorized_df.at[station, analyte] = analyte_data_full[start_index:].values

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


## by K. Whiteaker, kwhit@alum.mit.edu
def time_average_all_stations(
        # self, 
        data_pylenm_dm, 
        analyte, 
        period='1W', 
        rm_outliers=True, 
        std_thresh=2.2, 
        lowess_frac=0.1,
    ):
    """Transforms all analyte data from exact measurements, at potentially 
    different dates for each station, into periodic averaged data (ex. weekly 
    averages) at the same dates for all stations.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        analyte (_type_): analyte name
        period (str, optional): {‘D’, ‘W’, ‘M’, ‘Y’} time period over which to average. See https://pandas.pydata.org/pandas-docs/stable/user_guide/timeseries.html for valid frequency inputs. (e.g. ‘2W’ = every 2 weeks). Defaults to '1W'
        rm_outliers (bool, optional): flag to remove outliers in the data via LOWESS fit. Defaults to True
        std_thresh (int, optional): number of standard deviations in (observation - fit) outside of which is considered an outlier. Defaults to 2.2
        lowess_frac (float, optional): fraction of total data points considered in local fits of LOWESS smoother. A smaller value captures more local behaviour, and may be required for large datasets. Defaults to 0.1

    Returns:
        pd.DataFrame: dataframe with columns = stations, rows = dates, and data representing the average for that station between each date
    """
    # import data and filter for the chosen analyte
    data_concentration = data_pylenm_dm.data
    data_concentration = data_concentration[data_concentration.ANALYTE_NAME == analyte]
    data_construction = data_pylenm_dm.construction_data

    # create reshaped_df with columns=all available stations and rows=all measurement dates, to iterate over when conducting periodic averaging
    stations = np.unique(data_concentration[data_concentration.ANALYTE_NAME == analyte].STATION_ID)
    times = pd.to_datetime(np.unique(data_concentration.COLLECTION_DATE))
    reshaped_df = pd.DataFrame(index=times, columns=stations)
    
    # remove outliers from each station's data using remove_outlier_lowess(), and save the remaining data in reshaped_df
    for station in stations:
        data_in = data_concentration[data_concentration.STATION_ID == station].set_index('COLLECTION_DATE').RESULT
        data_in.index = pd.to_datetime(data_in.index)
        if rm_outliers:
            data_xOutliers = preprocess.remove_outliers_lowess(
                data_in, lowess_frac, std_thresh, 
            )
        else:
            data_xOutliers = data_in
        reshaped_df[station].loc[data_xOutliers.index] = data_xOutliers
    
    # find first and last measurement dates
    date_start = pd.Timestamp(times.min())
    date_end = pd.Timestamp(times.max())
    # make an array of datetimes, 1 period (ex. 1 week) apart, that spanning from before the first measurement date to after the last measurement date
    period_timesteps = pd.date_range(
        min(times) - pd.Timedelta(period), 
        max(times) + pd.Timedelta(period), 
        freq=period
    ).to_numpy()

    # conduct periodic averaging: take periodic (ex. 1 week) averages of each station and save the result to averaged_df
    averaged_df = pd.DataFrame(
        index=period_timesteps, 
        columns=stations, 
        dtype=float,
    )
    
    for period_step in tqdm(range(period_timesteps.size-1), desc="Timesteps", total=period_timesteps.size-1):  # for each period...
        # if period_step%100==0:
        #     print(period_step)  # print progress
        period_step_bin = []
        relevant_times = times[(times > period_timesteps[period_step]) & (times < period_timesteps[period_step+1])]
        
        for relevant_time in relevant_times:  # for each existing measurement that occurred within this period...
            timestep_data = reshaped_df.loc[relevant_time]  # identify all data points at this timestep
            period_step_bin.append(timestep_data)
        
        if len(period_step_bin)>0:  # if there are any measurements in this period
            period_data = pd.concat(period_step_bin, axis=1)
            period_station_avg = period_data.mean(axis=1, skipna=True)  # for each station, take the average over this period
            averaged_df.loc[period_timesteps[period_step]] = period_station_avg
    
    return averaged_df