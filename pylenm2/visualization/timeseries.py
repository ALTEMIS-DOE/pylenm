import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from dateutil.relativedelta import relativedelta

from pyproj import Proj, Transformer

import pylenm2
from pylenm2.data import filters, fetchers, transformation
from pylenm2.utils import constants as c

import logging
from pylenm2 import logger_config

timeseries_logger = logger_config.setup_logging(
    module_name=__name__,
    level=logging.INFO,
    # level=logging.DEBUG,
    logfile_dir=c.LOGFILE_DIR,
)


def plot_all_time_series_simple(
        # self, 
        data_pylenm_dm,
        analyte_name=None, 
        start_date=None, 
        end_date=None, 
        title='Dataset: Time ranges', 
        x_label='Station', 
        y_label='Year',
        min_days=10, 
        x_min_lim=-5, 
        x_max_lim = 170, 
        y_min_date='1988-01-01', 
        y_max_date='2020-01-01', 
        return_data=False, 
        filter=False, 
        col=None, 
        equals=[],
    ):
    """Plots the start and end date of analyte readings for differnt locations/sensors/stations.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        analyte_name (str, optional): analyte to examine. Defaults to None.
        start_date (str, optional): start date of horizontal time to show alignment. Defaults to None.
        end_date (str, optional): end date of horizontal time to show alignment.. Defaults to None.
        title (str, optional): plot title. Defaults to 'Dataset: Time ranges'.
        x_label (str, optional): x axis label. Defaults to 'Station'.
        y_label (str, optional): y axis label. Defaults to 'Year'.
        min_days (int, optional): minimum number of days required to plot the time series . Defaults to 10.
        x_min_lim (int, optional): x axis starting point. Defaults to -5.
        x_max_lim (int, optional): x axis ending point. Defaults to 170.
        y_min_date (str, optional): y axis starting date. Defaults to '1988-01-01'.
        y_max_date (str, optional): y axis ending date. Defaults to '2020-01-01'.
        return_data (bool, optional): flag to return data. Defaults to False.
        filter (bool, optional): flag to indicate filtering. Defaults to False.
        col (str, optional): column to filter. Example: col='STATION_ID'. Defaults to None.
        equals (list, optional): values to filter col by. Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].
    """
    # data = self.simplify_data()
    data = filters.simplify_data(data=data_pylenm_dm)
    
    # Filtering
    # -----------------
    if(filter):
        filter_res = filters.filter_by_column(
            data=data_pylenm_dm.construction_data, 
            col=col, 
            equals=equals,
        )
        # if('ERROR:' in str(filter_res)):
        if filter_res is None:
            timeseries_logger.debug("Ran into ERROR when calling filter_by_column()!")
            return filter_res
        
        query_stations = list(data.STATION_ID.unique())
        filter_stations = list(filter_res.index.unique())
        intersect_stations = list(set(query_stations) & set(filter_stations))
        if(len(intersect_stations)<=0):
            timeseries_logger.warning('ERROR: No results for this query with the specified filter parameters.')
            # return 'ERROR: No results for this query with the specifed filter parameters.'
            return None
        
        data = data[data['STATION_ID'].isin(intersect_stations)]

    # Preparing data
    # -----------------
    if(analyte_name!=None):
        data = data[data.ANALYTE_NAME == analyte_name]
    
    stations = data.STATION_ID.unique()
    stations_dateRange=pd.DataFrame(columns=['STATION_ID','START_DATE','END_DATE'])
    for i in range(len(stations)):
        stationName=stations[i]
        stationNamedData=data[data['STATION_ID']==stations[i]]
        minDate=min(stationNamedData['COLLECTION_DATE'])
        maxDate=max(stationNamedData['COLLECTION_DATE'])
        stations_dateRange.loc[stations_dateRange.shape[0]]=[stationName,minDate,maxDate]

    stations_dateRange["RANGE"] = stations_dateRange.END_DATE - stations_dateRange.START_DATE
    
    # stations_dateRange.RANGE = stations_dateRange.RANGE.astype('timedelta64[D]').astype('int')
    # stations_dateRange = stations_dateRange[stations_dateRange.RANGE>min_days]
    stations_dateRange = stations_dateRange[stations_dateRange.RANGE.dt.days>min_days]
    
    stations_dateRange.sort_values(by=["RANGE","END_DATE","START_DATE"], ascending = (False, False, True), inplace=True)
    stations_dateRange.reset_index(inplace=True)
    stations_dateRange.drop('index', axis=1, inplace=True)
    stations = np.array(stations_dateRange.STATION_ID)

    # Plotting
    # -----------------
    fig, ax = plt.subplots(1, 1, sharex=False,figsize=(20,6),dpi=300)

    ax.set_xticks(range(len(stations)))
    ax.set_xticklabels(stations, rotation='vertical', fontsize=6)

    ax.plot(stations_dateRange['START_DATE'], c='blue', marker='o',lw=0, label='Start date')
    ax.plot(stations_dateRange['END_DATE'], c='red', marker='o',lw=0, label='End date')

    ax.hlines([max(stations_dateRange['END_DATE'])], x_min_lim, x_max_lim, colors='purple', label='Selected end date')
    if(start_date==None):
        ax.hlines([min(stations_dateRange['START_DATE'])], x_min_lim, x_max_lim, colors='green', label='Selected start date')
    else:
        ax.hlines([pd.to_datetime(start_date)], x_min_lim, x_max_lim, colors='green', label='Selected start date')

    x_label = x_label + ' (count: ' + str(stations_dateRange.shape[0])+ ')'
    ax.set_xlabel(x_label, fontsize=20)
    ax.set_ylabel(y_label, fontsize=20)   
    ax.set_xlim([x_min_lim, x_max_lim])
    ax.set_ylim([pd.to_datetime(y_min_date), pd.to_datetime(y_max_date)]) 
    ax.plot([], [], ' ', label="Time series with at least {} days".format(min_days))
    ax.legend()
    
    if(analyte_name!=None):
        title = title + ' (' + analyte_name + ')'
    fig.suptitle(title, fontsize=20)
    for i in range(stations_dateRange.shape[0]):
        ax.vlines(i,stations_dateRange.loc[i,'START_DATE'],stations_dateRange.loc[i,'END_DATE'],colors='k')
    
    if(return_data):
        return stations_dateRange


def plot_all_time_series(
        # self, 
        data_pylenm_dm,
        analyte_name=None, 
        title='Dataset: Time ranges', 
        x_label='Station', 
        y_label='Year', 
        x_label_size=8, 
        marker_size=30,
        min_days=10, 
        x_min_lim=None, 
        x_max_lim=None, 
        y_min_date=None, 
        y_max_date=None, 
        sort_by_distance=True, 
        source_coordinate=[436642.70,3681927.09], 
        log_transform=False, 
        cmap=mpl.cm.rainbow, 
        drop_cols=[], 
        return_data=False, 
        filter=False, 
        col=None, 
        equals=[], 
        cbar_min=None, 
        cbar_max=None, 
        reverse_y_axis=False, 
        fontsize = 20, 
        figsize=(20,6), 
        dpi=300, 
        y_2nd_label=None,
    ):      # TODO: Getting error in the notebook!
    """Plots the start and end date of analyte readings for differnt locations/sensors/stations with colored concentration reading.

    Args:
        data_pylenm_dm (pylenm2.PylenmDataModule): PylenmDataModule object containing the concentration and construction data.
        analyte_name (str, optional): analyte to examine. Defaults to None.
        title (str, optional): plot title. Defaults to 'Dataset: Time ranges'.
        x_label (str, optional): x axis label. Defaults to 'Station'.
        y_label (str, optional): y axis label. Defaults to 'Year'.
        x_label_size (int, optional): x axis label font size. Defaults to 8.
        marker_size (int, optional): point size for time series. Defaults to 30.
        min_days (int, optional): minimum number of days required to plot the time series . Defaults to 10.
        x_min_lim (int, optional): x axis starting point. Defaults to None.
        x_max_lim (int, optional): x axis ending point. Defaults to None.
        y_min_date (str, optional): y axis starting date. Defaults to None.
        y_max_date (str, optional): y axis ending date. Defaults to None.
        sort_by_distance (bool, optional): flag to sort by distance from source center. Defaults to True.
        source_coordinate (list, optional): Easting, Northing coordinate of source center. Defaults to [436642.70,3681927.09].
        log_transform (bool, optional): flag to toggle log base 10 transformation. Defaults to False.
        cmap (cmap, optional): color map for plotting. Defaults to mpl.cm.rainbow.
        drop_cols (list, optional): columns, usually stations, to exclude. Defaults to [].
        return_data (bool, optional): flag to return data. Defaults to False.
        filter (bool, optional): flag to indicate filtering. Defaults to False.
        col (str, optional): column to filter. Example: col='STATION_ID'. Defaults to None.
        equals (list, optional): values to filter col by. Examples: equals=['FAI001A', 'FAI001B']. Defaults to [].
        cbar_min (float, optional): color bar lower boundary. Defaults to None.
        cbar_max (float, optional): color bar upper boundary. Defaults to None.
        reverse_y_axis (bool, optional): flag that reverses y axis. Defaults to False.
        fontsize (int, optional): plot font size. Defaults to 20.
        figsize (tuple, optional): matplotlib style figure size. Defaults to (20,6).
        dpi (int, optional): DPI of figure. Defaults to 300.
        y_2nd_label (str, optional): color bar label manual override. Defaults to None.
    """
    
    # Get clean data
    # dt = self.getCleanData([analyte_name])
    dt = fetchers.getCleanData(
        data_pylenm_dm, 
        analytes=[] if analyte_name is None else [analyte_name],
    )
    dt = dt[analyte_name]
    
    # try:
    # Filter data (if necessary)
    if(filter):
        filter_res = filters.filter_by_column(
            data=data_pylenm_dm.get_construction_data(), 
            col=col, 
            equals=equals,
        )
        # if('ERROR:' in str(filter_res)):
        if filter_res is None:
            timeseries_logger.debug("Ran into ERROR when calling filter_by_column()!")
            return filter_res
        
        query_stations = list(dt.columns.unique())
        filter_stations = list(filter_res.index.unique())
        intersect_stations = list(set(query_stations) & set(filter_stations) & set(dt.columns))
        if(len(intersect_stations)<=0):
            timeseries_logger.warning('ERROR: No results for this query with the specifed filter parameters.')
            # return 'ERROR: No results for this query with the specifed filter parameters.'
            return None

        dt = dt[intersect_stations]

    # Get station information
    station_info = data_pylenm_dm.get_construction_data()
    shared_stations = list(set(station_info.index) & set(dt.columns))
    dt = dt[shared_stations]
    station_info = station_info.T[shared_stations]
    dt = dt.reindex(sorted(dt.columns), axis=1)
    station_info = station_info.reindex(sorted(station_info.columns), axis=1)
    station_info = station_info.T
    
    # Add distance to source
    transformer = Transformer.from_crs("epsg:4326", "epsg:26917") # Latitude/Longitude to UTM
    UTM_x, UTM_y = transformer.transform(station_info.LATITUDE, station_info.LONGITUDE)
    X = np.vstack((UTM_x,UTM_y)).T
    station_info = pd.DataFrame(X, index=list(station_info.index),columns=['Easting', 'Northing'])
    station_info = transformation.add_dist_to_source(
        XX=station_info, 
        source_coordinate=source_coordinate,
    )

    # Sort by distance (if necessary)
    if(sort_by_distance):
        station_info.sort_values(
            by=['dist_to_source'], 
            ascending=True, 
            inplace=True,
        )
    
    dt = dt[station_info.index]

    # except Exception as e:
    #     timeseries_logger.error(f"Error occurred but ignoring the error - {e}")

    # Iterpolate/drop data
    dt = dt.interpolate()
    dt = dt.drop(drop_cols, axis=1) # DROP BAD ONES 
    
    # Take log transform of the data (if necessary)
    if(log_transform):
        dt[dt <= 0] = 0.00000001
        dt = np.log10(dt)
    
    # Plot data
    stations = dt.columns
    # if(cbar_min==None):
    #     cbar_min = dt.min().min()
    # if(cbar_max==None):
    #     cbar_max = dt.max().max()
    cbar_min = dt.min().min() if cbar_min is None else cbar_min 
    cbar_max = dt.max().max() if cbar_max is None else cbar_max
    norm = mpl.colors.Normalize(vmin=cbar_min, vmax=cbar_max)

    fig, ax = plt.subplots(
        1, 2, sharex=False, figsize=figsize, dpi=dpi, 
        gridspec_kw={'width_ratios': [40, 1]},
    )
    ax[0].set_xticks(range(len(stations)))
    ax[0].set_xticklabels(stations, rotation='vertical', fontsize=x_label_size)

    for col in stations:
        curr_start = dt[col].first_valid_index()
        curr_end =  dt[col].last_valid_index()
        length = len(list(dt[col].loc[curr_start:curr_end].index))
        color_vals = list(dt[col].loc[curr_start:curr_end])
        color_vals = [cmap(norm((x))) for x in color_vals]
        
        x = col
        ys = list(dt[col].loc[curr_start:curr_end].index)
        ax[0].scatter([str(x)]*length, ys, c=color_vals, marker='o',lw=0,s=marker_size, alpha=0.75)

    # Set plot labels
    x_label = x_label + ' (count: ' + str(dt.shape[1])+ ')'
    ax[0].set_xlabel(x_label, fontsize=fontsize)
    ax[0].set_ylabel(y_label, fontsize=fontsize)
    
    # Set plot limits
    # if(x_min_lim==None):
    #     x_min_lim = -5
    # if(x_max_lim==None):
    #     x_max_lim = len(stations)+5
    x_min_lim = -5 if x_min_lim is None else x_min_lim
    x_max_lim = len(stations)+5 if x_max_lim is None else x_max_lim
    ax[0].set_xlim([x_min_lim, x_max_lim])
    
    if (y_min_date==None):
        y_min_date = dt.index.min() + relativedelta(years=-1)
    if (y_max_date==None):
        y_max_date = dt.index.max() + relativedelta(years=1)
    ax[0].set_ylim([pd.to_datetime(y_min_date), pd.to_datetime(y_max_date)]) 
    ax[0].plot([], [], ' ', label="Time series with at least {} days".format(min_days))
    ax[0].set_facecolor((0, 0, 0, 0.1))
            
    # Add colorbar
    label_cb = "Concentration ({})".format(data_pylenm_dm.get_unit(analyte_name))
    if(log_transform):
        label_cb = "Log " + label_cb
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=ax[1], 
        orientation='vertical',
    )
    if(y_2nd_label!=None):
        label_cb = y_2nd_label
    cbar.set_label(label=label_cb, size=fontsize)

    ax[0].tick_params(axis='y', labelsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize)
    plt.tight_layout()
    if(reverse_y_axis):
        ax[0].invert_yaxis()
    
    if(analyte_name!=None):
        title = title + ' (' + analyte_name + ')'
    fig.suptitle(title, fontsize=fontsize, y=1.05)
    
    # Return data if necessary
    if(return_data):
        return dt
